
library(numDeriv)

##-----------------------------------------------
generate_animals <- function(
    width,
    height,
    density,
    distribution = c("uniform", "clustered"),
    cluster_radius = 5,
    mean_cluster_size = 10)
     {

  distribution <- match.arg(distribution)
  area <- width * height

  # Total number of animals implied by density
  N <- rpois(1, lambda = density * area)

  if (distribution == "uniform") {
    pts <- data.frame(
      x = runif(N, 0, width),
      y = runif(N, 0, height),
      cluster_id = 1
    )
  } else {
    # Number of cluster centres needed
    n_clusters <- max(1, ceiling(N / mean_cluster_size))

    centres <- data.frame(
      x = runif(n_clusters, 0, width),
      y = runif(n_clusters, 0, height)
    )

    # Allocate animals among clusters
    cluster_id <- sample(seq_len(n_clusters), N, replace = TRUE)
    x = rnorm(N, mean = centres$x[cluster_id], sd = cluster_radius)
    y = rnorm(N, mean = centres$y[cluster_id], sd = cluster_radius)

    # rejection sampling to handle points outside the simulated area
    outside <- x < 0 | x > width | y < 0 | y > height
    while (any(outside)) {
      idx <- which(outside)
      x[idx] <- rnorm(length(idx), centres$x[cluster_id[idx]], cluster_radius)
      y[idx] <- rnorm(length(idx), centres$y[cluster_id[idx]], cluster_radius)
      outside <- x < 0 | x > width | y < 0 | y > height
    }
    pts <- data.frame(
      x = x, y = y, cluster = cluster_id
    )
  }
  pts
}
##---------------------------------------
## Fit multiple detection functions and return the best by AIC
##---------------------------------------
select_best_ds <- function(dist_ds, w, binned = FALSE, breaks = NULL) {
  # Fit detection functions to distance data using ds() in the Distance package
  # cycle through candidate models and select best using AIC
  # Candidate models: key + adjustment + formula combinations
  candidates <- list(
    list(key = "hn",   adjustment = "cos", formula = ~1),
    list(key = "hn",   adjustment = NULL, formula = ~gs),
    list(key = "hr",   adjustment = "cos", formula = ~1),
    list(key = "unif", adjustment = "cos", formula = ~1)
  )

  fits <- list()
  aics <- c()

  for (i in seq_along(candidates)) {
    cand <- candidates[[i]]
    fit_i <- tryCatch({
      if (binned && !is.null(breaks)) {
        dist_ds<- bin_distances_cut(dist_ds, cutpoints = breaks)
        ds(dist_ds, key = cand$key, adjustment = cand$adjustment,
           transect = "point", formula = cand$formula)
      } else {
        ds(dist_ds, key = cand$key, adjustment = cand$adjustment,
           transect = "point", truncation = w, formula = cand$formula)
      }
    }, error = function(e) NULL)

    if (!is.null(fit_i)) {
      fits[[length(fits) + 1]] <- fit_i
      aics <- c(aics, AIC(fit_i)$AIC)
    }
  }

  if (length(fits) == 0) stop("All detection function models failed to converge.")
  fits[[which.min(aics)]]
}

##---------------------------------------
## standard CTDS - one replicate (multiple cameras)
##---------------------------------------
sim_ctds <- function(D_true, sigma_closest, sigma_true=NULL, w, fov,
                     n_cam, width, height,
                     distribution = "uniform", perfect_fov = TRUE,
                     cluster_radius = 30, mean_cluster_size = 5,
                     min_total = 10, binned = FALSE, breaks = NULL) {
  require(Distance)
  animals<- generate_animals(width = width,
                             height = height,
                             density = D_true,
                             distribution = distribution,
                             cluster_radius = cluster_radius,
                             mean_cluster_size = mean_cluster_size)

  # Camera locations, kept >= w from every edge so each sector stays
  # fully inside the region for any bearing (requires region wider than 2w).
  cx <- runif(n_cam, w, width - w)
  cy <- runif(n_cam, w, height - w)


  #  bearing <- runif(n_cam, 0, 360)
  bearing <- rep(270, n_cam) ##  All cams facing "south"

  counts <- integer(n_cam)
  dist_list <- vector("list", n_cam)
  for (k in seq_len(n_cam)) {

    dk <- sample_sector(animals, origin = c(cx[k], cy[k]), bearing = bearing[k],
                        radius = w, angle = fov)
    if(nrow(dk) == 0L) {
      dist_list[[k]]<- data.frame(Sample.Label = paste0("C",k),
                                  distance = NA,
                                  gs=NA,
                                  size=NA)


    } else {
      dk<- dk[order(dk$r),]
      nk<- nrow(dk)
      r_closest <- dk$r[1]
      if (rbinom(1, 1, hn_func(r_closest, sigma_true)) == 1L) {
        # Camera triggered by closest individual
        if(!perfect_fov){
          # Further animals may be obscured. detection is HN
          p <- c(1.0, hn_func(dk$r[2:nk], sigma_true)) # r[1] already detected
          dk<- dk[rbinom(nk, 1, p) == 1, , drop = FALSE]
          dist_list[[k]]<- data.frame(Sample.Label = paste0("C",k),
                                      distance = dk$r,
                                      gs = nrow(dk),
                                      size=1)
          counts[k]<- nrow(dk)
        } else {
          dist_list[[k]] <- data.frame(Sample.Label = paste0("C",k),
                                       distance = dk$r,
                                       gs = nrow(dk),
                                       size=1)
          counts[k]<- nrow(dk)
        }
      } else {
        dist_list[[k]]<- data.frame(Sample.Label = paste0("C",k),
                                    distance = NA,
                                    gs=NA,
                                    size=NA)
      }
    }
  }


  dist <- list_rbind(dist_list)
  dist<- dist |> mutate(object = row_number(), Region.Label = "CTDS",
                        Effort = 1, Area = 1)
  dist_ds<- dist |> filter(!is.na(distance))


  if (nrow(dist_ds) < min_total) return(NULL)

  tryCatch({
    fit <- select_best_ds(dist_ds, w = w, binned = binned, breaks = breaks)


    est <- dht2(fit, flatfile = dist, strat_formula = ~1, er_est="P2",
                sample_fraction = fov/360)

    # CV of the detection function
    fit_summ <- summary(fit)
    cv_det <- fit_summ$ds$average.p.se / fit_summ$ds$average.p

    data.frame(
      n_total      = est$n,
      mean_n       = mean(counts),
      D            = est$Abundance,
      se           = est$Abundance_se,
      cv_encounter = est$ER_CV,
      cv_detection = cv_det,
      lcl          = est$LCI,
      ucl          = est$UCI,
      cover        = est$LCI <= D_true & D_true <= est$UCI
    )
  }, error = function(e) {
    data.frame(
      n_total      = sum(counts),
      mean_n       = mean(counts),
      D            = NA_real_,
      se           = NA_real_,
      cv_encounter = NA_real_,
      cv_detection = NA_real_,
      lcl          = NA_real_,
      ucl          = NA_real_,
      cover        = NA
    )
  })
}

##---------------------------------------
## standard CTDS - (time lapse cameras)
##---------------------------------------
sim_lapse <- function(D_true, sigma_true, w, fov,
                     n_cam, width, height,
                     distribution = "uniform", perfect_fov = TRUE,
                     cluster_radius = 30, mean_cluster_size = 5,
                     min_total = 10, binned = FALSE, breaks = NULL) {
  require(Distance)
  animals<- generate_animals(width = width,
                             height = height,
                             density = D_true,
                             distribution = distribution,
                             cluster_radius = cluster_radius,
                             mean_cluster_size = mean_cluster_size)

  # Camera locations, kept >= w from every edge so each sector stays
  # fully inside the region for any bearing (requires region wider than 2w).
  cx <- runif(n_cam, w, width - w)
  cy <- runif(n_cam, w, height - w)

  #  bearing <- runif(n_cam, 0, 360)
  bearing <- rep(270, n_cam) ##  All cams facing "south"

  counts <- integer(n_cam)
  dist_list <- vector("list", n_cam)
  for (k in seq_len(n_cam)) {

    dk <- sample_sector(animals, origin = c(cx[k], cy[k]), bearing = bearing[k],
                        radius = w, angle = fov)
    if(nrow(dk) == 0L) {
      ## no individuals in sector
      dist_list[[k]]<- data.frame(Sample.Label = paste0("C",k),
                                  distance = NA,
                                  gs=NA,
                                  size=NA)

    } else if(!perfect_fov){
          # animals may be obscured. detection is HN
          nk<- nrow(dk)
          p <- hn_func(dk$r, sigma_true)
          dk<- dk[rbinom(nk, 1, p) == 1, , drop = FALSE]
          if(nrow(dk) >= 1) {
            dist_list[[k]]<- data.frame(Sample.Label = paste0("C",k),
                                        distance = dk$r,
                                        gs = nrow(dk),
                                        size=1)
            counts[k]<- nrow(dk)
          } else {
              dist_list[[k]]<- data.frame(Sample.Label = paste0("C",k),
                                          distance = NA,
                                          gs=NA,
                                          size=NA)
          }
        } else {
          # No detection error
          dist_list[[k]] <- data.frame(Sample.Label = paste0("C",k),
                                       distance = dk$r,
                                       gs = nrow(dk),
                                       size=1)
          counts[k]<- nrow(dk)
      }
  }

  dist <- list_rbind(dist_list)
  dist<- dist |> mutate(object = row_number(), Region.Label = "CTDS",
                        Effort = 1, Area = 1)
  dist_ds<- dist |> filter(!is.na(distance))

  if (nrow(dist_ds) < min_total) return(NULL)

  # Fit multiple detection functions and keep the one with lowest AIC
  # estimate density using dht2 catching errors

  tryCatch({
    fit <- select_best_ds(dist_ds, w = w, binned = binned, breaks = breaks)

    est <- dht2(fit, flatfile = dist, strat_formula = ~1, er_est="P2",
                sample_fraction = fov/360)

    # CV of the detection function
    fit_summ <- summary(fit)
    cv_det <- fit_summ$ds$average.p.se / fit_summ$ds$average.p

    data.frame(
      n_total      = est$n,
      mean_n       = mean(counts),
      D            = est$Abundance,
      se           = est$Abundance_se,
      cv_encounter = est$ER_CV,
      cv_detection = cv_det,
      lcl          = est$LCI,
      ucl          = est$UCI,
      cover        = est$LCI <= D_true & D_true <= est$UCI
    )
  }, error = function(e) {
    data.frame(
      n_total      = sum(counts),
      mean_n       = mean(counts),
      D            = NA_real_,
      se           = NA_real_,
      cv_encounter = NA_real_,
      cv_detection = NA_real_,
      lcl          = NA_real_,
      ucl          = NA_real_,
      cover        = NA
    )
  })
}
##---------------------------------------
## Closest distance
##---------------------------------------
sim_closest <- function(D_true, sigma_closest, sigma_true, w, fov,
                        n_cam, width, height,
                        distribution = "uniform", perfect_fov = TRUE,
                        cluster_radius = 30, mean_cluster_size = 5,
                        min_total = 10, binned = FALSE, breaks = NULL) {

  animals<- generate_animals(width = width,
                             height = height,
                             density = D_true,
                             distribution = distribution,
                             cluster_radius = cluster_radius,
                             mean_cluster_size = mean_cluster_size)
  # Camera locations, kept >= w from every edge so each sector stays
  # fully inside the region for any bearing (requires region wider than 2w).
  cx <- runif(n_cam, w, width - w)
  cy <- runif(n_cam, w, height - w)

  bearing <- rep(270, n_cam)  ##  All cams facing "south"

  counts <- integer(n_cam)
  dist_vec <- rep(NA_real_, n_cam)
  for (k in seq_len(n_cam)) {

    dk <- sample_sector(animals, origin = c(cx[k], cy[k]), bearing = bearing[k],
                        radius = w, angle = fov)
    if(nrow(dk) == 0L) next
    dk<- dk[order(dk$r),]
    r_closest<- dk$r[1]
    nk<- nrow(dk)
    if (rbinom(1, 1, hn_func(r_closest, sigma_closest)) == 1L) {
      # Camera triggered by closest individual
      if(!perfect_fov){
        # Further animals may be obscured. detection is HN
        p <- c(1.0, hn_func(dk$r[2:nk], sigma_true)) # r[1] already detected
        dk<- dk[rbinom(nk, 1, p) == 1, , drop = FALSE]
        counts[k]<- nrow(dk) # retain detected for group size
        dist_vec[k]<- r_closest # for detection function
      } else {
        counts[k]<- nrow(dk)
        dist_vec[k]<- r_closest
      }
    }
  }
  if (length(dist_vec[!is.na(dist_vec)]) < min_total) return(NULL)

  if(binned & !is.null(breaks)) {
    tabs<- table_counts(counts, dist_vec, breaks = breaks)
    bin_counts<- tabs$bin_list
    grps<- tabs$grps
    fit <- fit_detection_hn(bin_counts, w=w, gs=grps, binned=TRUE, breaks=breaks)
  } else {
    grps<- counts[!is.na(dist_vec)]
    dist <- dist_vec[!is.na(dist_vec)]
    fit <- fit_detection_hn(dist, w=w, gs=grps)
  }

  est <- estimate_density_closest(counts, w = w, angle = fov, fit = fit)

  data.frame(
    n_total      = est$n,
    mean_n       = mean(counts),
    D            = est$D,
    se           = est$se,
    cv_encounter = est$cv_encounter,
    cv_detection = est$cv_detection,
    lcl          = est$lcl,
    ucl          = est$ucl,
    cover        = est$lcl <= D_true & D_true <= est$ucl
  )
}
##---------------------------------------------------------------------
sim_once <- function(D_true, sigma_true, w, fov, n_cam, width, height,
                     distribution = "uniform", perfect_fov = TRUE,
                     cluster_radius = 30, mean_cluster_size = 5, min_total = 10,
                     binned = FALSE, breaks = NULL) {

  animals<- generate_animals(width = width,
                             height = height,
                             density = D_true,
                             distribution = distribution,
                             cluster_radius = cluster_radius,
                             mean_cluster_size = mean_cluster_size)

  # Camera locations, kept >= w from every edge so each sector stays
  # fully inside the region for any bearing (requires region wider than 2w).
  cx <- runif(n_cam, w, width - w)
  cy <- runif(n_cam, w, height - w)
  #  bearing <- runif(n_cam, 0, 360)
  bearing <- rep(270, n_cam) ##  All cams facing "south"

  counts <- integer(n_cam)
  dist_list <- vector("list", n_cam)
  for (k in seq_len(n_cam)) {

    dk <- sample_sector(animals, origin = c(cx[k], cy[k]), bearing = bearing[k],
                        radius = w, angle = fov, gr = hn_func, sigma = sigma_true)
    if(nrow(dk) == 0L) next
    counts[k] <- nrow(dk)  # this is now group size
    dist_list[[k]] <- dk$r

  }

  dist <- unlist(dist_list)
  if (length(dist) < min_total) return(NULL)

  if(binned & !is.null(breaks)) {
    bin_counts<- make_bins(dist, breaks = breaks)
    fit <- fit_detection_hn(bin_counts, w=w, gs=1, binned=TRUE, breaks=breaks)
  } else {
    fit <- fit_detection_hn(dist, w=w, gs=1)
  }

  est <- estimate_density_ctds(counts, w = w, angle = fov, fit = fit)

  data.frame(
    n_total      = est$n,
    mean_n       = mean(counts),
    D            = est$D,
    se           = est$se,
    cv_encounter = est$cv_encounter,
    cv_detection = est$cv_detection,
    lcl          = est$lcl,
    ucl          = est$ucl,
    cover        = est$lcl <= D_true & D_true <= est$ucl
  )
}


##--------------------------------------------------------------

sample_sector <- function(points,
                          origin = c(0, 0),
                          bearing = 0,
                          radius,
                          angle,
                          gr = NULL,
                          ...) {
  # Sample (detect) animals falling within a camera FOV sector.
  #   points  : data.frame with columns x, y (e.g. from generate_animals)
  #   origin  : c(x, y) camera location
  #   bearing : FOV direction in degrees
  #   radius  : truncation distance of the sector (w)
  #   angle   : full angular width of the FOV in degrees (split +/- about bearing)
  #   gr      : optional detection function gr(r, ...)
  # Returns detected points with polar coordinates r (distance from camera)
  # and theta (angle relative to bearing, radians), plus the detection
  # probability p when gr is supplied.
  dx <- points$x - origin[1]
  dy <- points$y - origin[2]
  r <- sqrt(dx^2 + dy^2)


  # Angle of each point relative to the bearing, wrapped to (-pi, pi]
  theta <- atan2(dy, dx) - bearing * pi / 180
  theta <- atan2(sin(theta), cos(theta))


  half <- (angle * pi / 180) / 2
  inside <- r <= radius & abs(theta) <= half


  out <- points[inside, , drop = FALSE]
  out$r <- r[inside]
  out$theta <- theta[inside]

  # detection function:
  if (!is.null(gr)) {
    p <- gr(out$r, ...)
    out$p <- p
    out <- out[rbinom(nrow(out), 1, p) == 1, , drop = FALSE]
  }

  out
}
##-----------------------------------------------------
table_counts<- function(counts, dists, breaks) {
  # function to bin distances
  # counts : number of animals occuring in FOV (group size)
  # dists  : distances
  # breaks : bin endpoints
  #
  ii<- which(counts > 0)
  c_nonzero<- counts[ii]
  dd<- dists[ii]

  fc<- table(c_nonzero)
  n<- length(fc)
  blist<- vector("list", n)
  grps<- as.integer(names(fc))

  for(i in seq_len(n)) {
    ii<- which(c_nonzero == grps[i])
    d2<- dd[ii]
    blist[[i]]<- make_bins(d2, breaks=bin_bp)
  }
  list(bin_list=blist, grps=grps)
}

##-------------------------------------------------------
make_bins <- function(dist, breaks, right = FALSE, include_lowest = TRUE) {
  ## bin data keeping site hierarchy
  bin <- cut(dist, breaks = breaks, right = right, include.lowest = include_lowest)
  bin_levels <- levels(bin)
  # observed counts
  counts <- as.data.frame(table(bin = bin, useNA = "no"),
                          stringsAsFactors = FALSE)

  full <- data.frame(bin = bin_levels, stringsAsFactors = FALSE)
  merged <- merge(full, counts, by = "bin", all.x = TRUE, sort=FALSE)
  merged$Freq[is.na(merged$Freq)] <- 0L
  merged$bin<- factor(merged$bin, levels=bin_levels)
  # cast to matrix
  mat <- xtabs(Freq ~ bin, data = merged)
  return(mat)
}

##----------------------------------------
fit_detection_hn <- function(dist, w, gs = 1, binned = FALSE, breaks = NULL) {
  # Fit a half-normal detection function to observed sector distances by
  # conditional MLE
  # Returns sigma with an SE on the log scale (from the Hessian).
  if(binned & ! is.null(breaks)) {
    mle <- optim(
      par = log(w / 2),
      fn = nll.cond.binned.hn,
      counts = dist, gs = gs, breaks = breaks,
      method = "Brent", lower = -5, upper = 5,
      hessian = TRUE)
  } else {
    mle <- optim(
      par = log(w / 2),
      fn = nll.cond.point.hn,
      x = dist, w = w, gs = gs,
      method = "Brent", lower = -5, upper = 5,
      hessian = TRUE)
  }
  list(
    sigma        = exp(mle$par),
    log_sigma    = mle$par,
    se_log_sigma = sqrt(1 / mle$hessian[1, 1]),
    mle          = mle
  )
}

##----------------------------------------
estimate_density_ctds <- function(counts, w, angle, fit, level = 0.95) {
  # Density from several independent camera sectors of equal area
  # using design-based appraoch
  # counts : count of total detections per camera (length K, may include zeros)
  # fit    : output of fit_detection_hn() on the pooled distances
  # Encounter-rate variance is estimated empirically from the cameras
  # counts using Fewster et al: Biometrics (2009)
  # The detection-function variance is added via the delta method.

  K <- length(counts)
  if (K < 2) stop("need at least 2 cameras for an empirical variance")

  n_total <- sum(counts)
  theta<- angle/360 * pi * w^2
  m <- function(sigma, w, gs) {
    # product of availability, given group size (gs) and detection
    integrate(function(r) availability_cont(r, w, gs) * hn_func(r, sigma), 0 , w)$value
  }

  a <- theta * m(fit$sigma, w, gs=1)        # per-camera effective area
  D <- n_total / (K * a)

  # Encounter-rate variance (P2)
  R <- mean(counts)
  var_er <- sum((counts - R)^2) / (K * (K-1))
  cv2_er<- var_er/R^2

  # Detection-function variance (delta method)
  pbar_fun <- function(lsigma, w){m(exp(lsigma), w, 1)}
  gradient <- numDeriv::grad(pbar_fun, x = fit$log_sigma, w = w)
  se_det <- gradient * fit$se_log_sigma
  cv2_det <- (se_det/m(fit$sigma,w,gs=1))^2

  cv2 <- cv2_er + cv2_det
  se_D <- D * sqrt(cv2)

  z <- qnorm(1 - (1 - level) / 2)
  c_mult <- exp(z * sqrt(log(1 + cv2)))

  list(
    D = D, n = n_total, K = K, eff_area = a, sigma = fit$sigma,
    se = se_D, cv = sqrt(cv2),
    cv_encounter = sqrt(cv2_er), cv_detection = sqrt(cv2_det),
    lcl = D / c_mult, ucl = D * c_mult
  )
}

##----------------------------------------
estimate_density_closest <- function(counts, w, angle, fit, level = 0.95) {
  # Density from several independent camera sectors of equal area.
  # counts : count of total indivdiuals in FOV (=group size: length K, may include zeros)
  # fit    : output of fit_detection_hn() on the pooled closest distances
  # Encounter-rate variance is estimated empirically from the cameras
  # counts using Fewster et al: Biometrics (2009)
  # The detection-function variance is added via the delta method.

  K <- length(counts)
  if (K < 2) stop("need at least 2 cameras for an empirical variance")

  n_total <- sum(counts)
  theta<- angle/360 * pi * w^2  # sector area
  m <- function(sigma, w, gs) {
    # product of availability, given group size (gs) and detection
    integrate(function(r) availability_cont(r, w, gs) * hn_func(r, sigma), 0 , w)$value
  }

  ak<- rep(NA_real_, K)
  nk<- rep(NA_real_, K)
  c_nonzero<- counts
  c_nonzero[c_nonzero < 1]<- 1 # Zero counts get group size 1 for areas

  # effective detection probability now depends on group size
  for(k in seq_len(K)) {
    ak[k] <- m(fit$sigma, w, gs=c_nonzero[k])
    nk[k] <- counts[k]/ak[k]
  }

  D<- sum(nk)/(K * theta)

  # Encounter-rate variance (R2)
  R <- mean(counts)
  var_er <- sum((counts - R)^2) / (K * (K-1))
  cv2_er<- var_er/R^2

  # Detection-function variance (delta method)
  pbar_fun <- function(lsigma, w){m(exp(lsigma), w, 1)}
  gradient <- numDeriv::grad(pbar_fun, x = fit$log_sigma, w = w)
  se_det <- gradient * fit$se_log_sigma
  cv2_det <- (se_det/m(fit$sigma,w,gs=1))^2

  cv2 <- cv2_er + cv2_det
  se_D <- D * sqrt(cv2)

  z <- qnorm(1 - (1 - level) / 2)
  c_mult <- exp(z * sqrt(log(1 + cv2)))

  list(
    D = D, n = n_total, K = K, eff_area = mean(ak), sigma = fit$sigma,
    se = se_D, cv = sqrt(cv2),
    cv_encounter = sqrt(cv2_er), cv_detection = sqrt(cv2_det),
    lcl = D / c_mult, ucl = D * c_mult
  )
}



##---------------------------------------
## standard CTDS with closest distance trigger
##---------------------------------------
run_density_sim <- function(n_rep = 5,
                            D_true = 0.01,
                            sigma_closest = 4,
                            sigma_true= 12,
                            w = 12,
                            fov = 40,
                            perfect_fov = FALSE,
                            n_cam = 60,
                            width = 500,
                            height = 500,
                            distribution = c("uniform", "clustered"),
                            cluster_radius = 30,
                            mean_cluster_size = 5,
                            binned = FALSE,
                            breaks = NULL,
                            progress = TRUE) {
  # Each replicate places n_cam random cameras; the encounter-rate variance is
  # estimated across those cameras. Under clustering the empirical variance
  # captures over-dispersion, so CI coverage should recover toward nominal.
  distribution <- match.arg(distribution)
  res <- vector("list", n_rep)
  pb <- if (progress) utils::txtProgressBar(max = n_rep, style = 3) else NULL
  for (i in seq_len(n_rep)) {
    res[[i]] <- sim_ctds(D_true, sigma_closest, sigma_true, w, fov,
                         n_cam, width, height,
                         distribution = distribution,
                         perfect_fov = perfect_fov,
                         cluster_radius = cluster_radius,
                         mean_cluster_size = mean_cluster_size,
                         binned = binned,
                         breaks = breaks)
    if (progress) utils::setTxtProgressBar(pb, i)
  }
  if (progress) close(pb)
  out <- bind_rows(res)
  attr(out, "truth") <- list(D_true = D_true, sigma_true = sigma_true,
                             w = w, fov = fov, n_cam = n_cam,
                             distribution = distribution,
                             model = "CTDS")
  out
}

##---------------------------------------
## closest distance only with group size adjustment
##---------------------------------------
run_density_closest <- function(n_rep = 5,
                            D_true = 0.01,
                            sigma_closest = 4,
                            sigma_true = 12,
                            w = 15,
                            fov = 40,
                            perfect_fov = FALSE,
                            n_cam = 60,
                            width = 500,
                            height = 500,
                            distribution = c("uniform", "clustered"),
                            cluster_radius = 30,
                            mean_cluster_size = 5,
                            binned = FALSE,
                            breaks = NULL,
                            progress = TRUE) {
  # Each replicate places n_cam random cameras; the encounter-rate variance is
  # estimated across those cameras. Under clustering the empirical variance
  # captures over-dispersion, so CI coverage should recover toward nominal.
  distribution <- match.arg(distribution)
  res <- vector("list", n_rep)
  pb <- if (progress) utils::txtProgressBar(max = n_rep, style = 3) else NULL
  for (i in seq_len(n_rep)) {
    res[[i]] <- sim_closest(D_true, sigma_closest, sigma_true, w, fov,
                            n_cam, width, height,
                            distribution = distribution,
                            perfect_fov = perfect_fov,
                            cluster_radius = cluster_radius,
                            mean_cluster_size = mean_cluster_size,
                            binned = binned,
                            breaks = breaks)
    if (progress) utils::setTxtProgressBar(pb, i)
  }
  if (progress) close(pb)
  out <- bind_rows(res)
  attr(out, "truth") <- list(D_true = D_true, sigma_true = sigma_true,
                             w = w, fov = fov, n_cam = n_cam,
                             distribution = distribution,
                             model = "Closest")
  out
}
##---------------------------------------
## standard CTDS with time lapse (no sensor)
##---------------------------------------
run_density_lapse <- function(n_rep = 5,
                                D_true = 0.01,
                                sigma_true = 12,
                                w = 15,
                                fov = 40,
                                perfect_fov = FALSE,
                                n_cam = 60,
                                width = 500,
                                height = 500,
                                distribution = c("uniform", "clustered"),
                                cluster_radius = 30,
                                mean_cluster_size = 5,
                                binned = FALSE,
                                breaks = NULL,
                                progress = TRUE) {
  # Each replicate places n_cam random cameras; the encounter-rate variance is
  # estimated across those cameras. Under clustering the empirical variance
  # captures over-dispersion, so CI coverage should recover toward nominal.
  distribution <- match.arg(distribution)
  res <- vector("list", n_rep)
  pb <- if (progress) utils::txtProgressBar(max = n_rep, style = 3) else NULL
  for (i in seq_len(n_rep)) {
    res[[i]] <- sim_lapse(D_true, sigma_true, w, fov,
                            n_cam, width, height,
                            distribution = distribution,
                            perfect_fov = perfect_fov,
                            cluster_radius = cluster_radius,
                            mean_cluster_size = mean_cluster_size,
                            binned = binned,
                            breaks = breaks)
    if (progress) utils::setTxtProgressBar(pb, i)
  }
  if (progress) close(pb)
  out <- bind_rows(res)
  attr(out, "truth") <- list(D_true = D_true, sigma_true = sigma_true,
                             w = w, fov = fov, n_cam = n_cam,
                             distribution = distribution,
                             model = "lapse")
  out
}

##---------------------------------------
## Summaries and plots
##---------------------------------------
summarise_density_sim <- function(res) {
  truth <- attr(res, "truth")
  tibble(
    Model          = truth$model,
    distribution   = truth$distribution,
    n_cam          = truth$n_cam,
    replicates     = length(res$D[!is.na(res$D)]),
    D_true         = truth$D_true,
    D_hat          = mean(res$D, na.rm=TRUE),
    bias           = 100 *(mean(res$D, na.rm=TRUE) - truth$D_true) / truth$D_true,
    mean_SE        = mean(res$se, na.rm=TRUE),
    mean_cv_enc    = mean(res$cv_encounter, na.rm=TRUE),
    mean_cv_det    = mean(res$cv_detection, na.rm=TRUE),
    mean_CV        = mean(res$se/res$D, na.rm=TRUE),
    coverage_95    = mean(res$cover,na.rm=TRUE)
  )
}

##----------------------------------------------
plot_density_sim <- function(res, n_show = 100) {
  truth <- attr(res, "truth")

  p1 <- ggplot(res, aes(D)) +
    geom_histogram(bins = 30, fill = "grey70", color = "white") +
    geom_vline(xintercept = truth$D_true, color = "firebrick", linewidth = 1) +
    geom_vline(xintercept = mean(res$D), color = "steelblue",
               linetype = 2, linewidth = 1) +
    labs(x = expression(hat(D)), y = "Count",
         title = "Sampling distribution of density estimates",
         subtitle = "red = truth, blue dashed = mean estimate") +
    theme_bw()

  samp <- res |>
    slice_sample(n = min(n_show, nrow(res))) |>
    arrange(D) |>
    mutate(ord = row_number())

  p2 <- ggplot(samp, aes(ord, D, color = cover)) +
    geom_hline(yintercept = truth$D_true, color = "firebrick") +
    geom_pointrange(aes(ymin = lcl, ymax = ucl), fatten = 1) +
    scale_color_manual(values = c(`TRUE` = "grey40", `FALSE` = "orange")) +
    labs(x = "Replicate (ordered)", y = expression(hat(D)),
         title = paste0("95% CIs across ", nrow(samp), " replicates"),
         color = "covers truth") +
    theme_bw()

  p1 / p2
}

##------------------------------------
generate_cam_locs<- function(n_cam, w, width, height, bearing, camera_layout = c("random","grid")) {
  if (camera_layout == "grid") {
    # Build a near-square grid of >= n_cam points within the interior, spaced
    # evenly and inset by w from each edge, then take the first n_cam of them.
    n_col <- ceiling(sqrt(n_cam * (width - 2 * w) / (height - 2 * w)))
    n_row <- ceiling(n_cam / n_col)
    gx <- seq(w, width - w, length.out = n_col)
    gy <- seq(w, height - w, length.out = n_row)
    grid <- expand.grid(x = gx, y = gy)[seq_len(n_cam), ]
    cx <- grid$x
    cy <- grid$y
  } else {
    cx <- runif(n_cam, w, width - w)
    cy <- runif(n_cam, w, height - w)
  }
  bearing = rep(bearing, n_cam)
  data.frame(cx=cx, cy=cy, bearing = bearing)
}
##----------------------------------------
sector_polygon <- function(id, origin, bearing, radius, angle, n = 30) {
  # Build a closed polygon (data.frame of vertices) for one camera FOV sector,
  # matching the geometry in sample_sector(): bearing in degrees (0 = +x axis,
  # CCW), angle is the full FOV width split symmetrically about the bearing.
  # Returns data.frame(id, x, y) with the apex followed by the arc.
  half <- (angle * pi / 180) / 2
  b <- bearing * pi / 180
  a <- seq(b - half, b + half, length.out = n)
  data.frame(
    id = id,
    x = c(origin[1], origin[1] + radius * cos(a)),
    y = c(origin[2], origin[2] + radius * sin(a))
  )
}

##----------------------------------------
plot_animals_sectors <- function(animals,
                                 cams,
                                 radius,
                                 angle,
                                 bearing,
                                 width = NULL,
                                 height = NULL,
                                 show_detected = FALSE,
                                 detected_color = "firebrick",
                                 gr = NULL,
                                 ...) {
  # Plot a 2D animal point pattern with a set of camera FOV sectors overlaid.
  #   animals        : data.frame with columns x, y (e.g. from generate_animals)
  #   cams           : data.frame with columns cx, cy, bearing (one row per camera);
  #                    an id column is added if absent
  #   radius         : sector detection radius (w); angle: full FOV width in degrees
  #   width, height  : optional plot limits (defaults to the animal extent)
  #   show_detected  : when TRUE, animals falling inside any sector are coloured.
  #                    Detection is assessed with sample_sector() using the same
  #                    geometry; pass gr (and its params via ...) to thin by a
  #                    detection function, otherwise detection is perfect within view.
  #   detected_color : colour used for detected animals, sector fill, and camera points.
  if (is.null(cams$id)) cams$id <- seq_len(nrow(cams))
  if (is.null(width))  width  <- max(animals$x)
  if (is.null(height)) height <- max(animals$y)

  sectors <- cams |>
    purrr::pmap(\(id, cx, cy, bearing, ...)
                sector_polygon(id, c(cx, cy), bearing,
                               radius = radius, angle = angle)) |>
    dplyr::bind_rows()

  if (show_detected) {
    # Flag animals occurring within a sector
    animals$.row <- seq_len(nrow(animals))
    detected_rows <- cams |>
      purrr::pmap(\(id, cx, cy, bearing, ...)
                  sample_sector(animals, origin = c(cx, cy), bearing = bearing,
                                radius = radius, angle = angle, gr = gr, ...)$.row) |>
      unlist() |>
      unique()
    animals$detected <- animals$.row %in% detected_rows
  }

  p <- ggplot(animals, aes(x, y))

  if (show_detected) {
    p <- p +
      geom_polygon(data = sectors, aes(group = id),
                   fill = detected_color, alpha = 0.15, color = detected_color) +
      geom_point(aes(color = detected, alpha = detected), size = 1) +
      scale_color_manual(values = c(`FALSE` = "grey40", `TRUE` = detected_color),
                         labels = c(`FALSE` = "no", `TRUE` = "yes")) +
      scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.9), guide = "none") +
      labs(color = "Within FOV")
  } else {
    p <- p +
      geom_point(alpha = 0.4, size = 1) +
      geom_polygon(data = sectors, aes(group = id),
                   fill = detected_color, alpha = 0.3, color = detected_color)
  }

  p +
    geom_point(data = cams, aes(cx, cy), color = detected_color, size = 1.5) +
    coord_equal(xlim = c(0, width), ylim = c(0, height)) +
    labs(title = "Animal distribution with camera sectors",
         x = "x (m)", y = "y (m)") +
    theme_bw()
}

