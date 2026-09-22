
library(numDeriv)

## ----------  TOP ----------------------------------------

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
    # Number of cluster centres
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


## ----- Model selection by AIC ----------------------------
select_best_ds <- function(dist_ds, w, binned = FALSE, breaks = NULL) {
  # Fit detection functions to distance data using ds() in the Distance package
  # cycle through candidate models and select best using AIC
  # Candidate models: key + adjustment + formula combinations
  candidates <- list(
      list(key = "hn",   adjustment = "cos"),
      list(key = "hr",   adjustment = "cos"),
      list(key = "unif", adjustment = "cos")
    )

  fits <- list()
  aics <- c()

  for (i in seq_along(candidates)) {
    cand <- candidates[[i]]
    fit_i <- tryCatch({
      if (binned && !is.null(breaks)) {
        dist_ds<- bin_distances_cut(dist_ds, cutpoints = breaks)
        suppressWarnings(
          suppressMessages(
            ds(dist_ds, key = cand$key, adjustment = cand$adjustment, transect = "point")
            )
          )
      } else {
        suppressWarnings(
          suppressMessages(
            ds(dist_ds, key = cand$key, adjustment = cand$adjustment,
               transect = "point", truncation = w)
            )
          )
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


##---- Simulation functions -----------------------

sim_ctds <- function(D_true, sigma_closest, sigma_true=NULL,
                     w, fov, n_cam, width, height,
                     distribution = "uniform", perfect_fov = TRUE,
                     cluster_radius = 30, mean_cluster_size = 5,
                     min_total = 10, binned = FALSE, breaks = NULL) {
  ## standard CTDS for sensor based cameras assuming camera is triggered
  ## by the closest individual

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
      if (rbinom(1, 1, hn_func(r_closest, sigma_closest)) == 1L) {
        # Camera triggered by closest individual
        if(!perfect_fov){
          # Further animals may be obscured. detection is sigma_true
          p <- c(1.0, hn_func(dk$r[2:nk], sigma_true)) # r[1] already detected
          dk<- dk[rbinom(nk, 1, p) == 1, , drop = FALSE]
          dist_list[[k]]<- data.frame(Sample.Label = paste0("C",k),
                                      distance = dk$r,
                                      gs = nrow(dk),
                                      size=1)
          counts[k]<- nrow(dk)
        } else {
          # No detection error in FOV
          dist_list[[k]] <- data.frame(Sample.Label = paste0("C",k),
                                       distance = dk$r,
                                       gs = nrow(dk),
                                       size=1)
          counts[k]<- nrow(dk)
        }
      } else {
        # camera not triggered (animals present but not detected)
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
      mean_n       = mean(counts[counts > 0]),
      D            = est$Abundance,
      se           = est$Abundance_se,
      cv_encounter = est$ER_CV,
      cv_detection = cv_det,
      cv           = est$Abundance_se/est$Abundance,
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
      cv           = NA_real_,
      lcl          = NA_real_,
      ucl          = NA_real_,
      cover        = NA
    )
  })
}

##----- Closest distance -----------------------

sim_closest <- function(D_true, sigma_closest, sigma_true,
                        w, fov, n_cam, width, height,
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
        # No detection error. all individuals in FOV observed
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
    mean_n       = mean(counts[counts>0]),
    D            = est$D,
    se           = est$se,
    cv_encounter = est$cv_encounter,
    cv_detection = est$cv_detection,
    cv           = est$se/est$D,
    lcl          = est$lcl,
    ucl          = est$ucl,
    cover        = est$lcl <= D_true & D_true <= est$ucl
  )
}

## ----- Time lapse CTDS -----------------

sim_lapse <- function(D_true, sigma_true, w, fov, n_cam, width, height,
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
                        radius = w, angle = fov)
    if(nrow(dk) == 0L) next
      else if(!perfect_fov){
      # animals may be obscured. detection is HN
      nk<- nrow(dk)
      p <- hn_func(dk$r, sigma_true)
      dk<- dk[rbinom(nk, 1, p) == 1, , drop = FALSE]
      if(nrow(dk) >= 1) {
        dist_list[[k]]<- dk$r
        counts[k]<- nrow(dk)
        }
      } else {
      # No detection error. All individual in sector recorded
      dist_list[[k]] <- dk$r
      counts[k]<- nrow(dk)
    }
  }

  dist <- unlist(dist_list)
  if (length(dist) < min_total) return(NULL)

  if(binned & !is.null(breaks)) {
    bin_counts<- make_bins(dist, breaks = breaks)
    fit <- fit_detection_hn(bin_counts, w=w, gs=1, binned=TRUE, breaks=breaks)
  } else {
    fit <- fit_detection_hn(dist, w=w, gs=1, binned = FALSE)
  }

  est <- estimate_density_ctds(counts, w = w, angle = fov, fit = fit)

  data.frame(
    n_total      = est$n,
    mean_n       = mean(counts[counts>0]),
    D            = est$D,
    se           = est$se,
    cv_encounter = est$cv_encounter,
    cv_detection = est$cv_detection,
    cv           = est$se/est$D,
    lcl          = est$lcl,
    ucl          = est$ucl,
    cover        = est$lcl <= D_true & D_true <= est$ucl
  )
}

## ----- Sample sector ----------------------------------------

sample_sector <- function(points,
                          origin = c(0, 0),
                          bearing = 0,
                          radius,
                          angle) {
  # Sample (detect) animals falling within a camera FOV sector.
  #   points  : data.frame with columns x, y (e.g. from generate_animals)
  #   origin  : c(x, y) camera location
  #   bearing : FOV direction in degrees
  #   radius  : truncation distance of the sector (w)
  #   angle   : full angular width of the FOV in degrees (split +/- about bearing)
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

  out
}

##---- Binning -------------------------------------
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


bin_distances_bw <- function(data, bin_width, truncation) {
  breaks <- seq(0, truncation, by = bin_width)

  data |>
    dplyr::filter(distance <= truncation) |>
    dplyr::mutate(
      distbegin = breaks[findInterval(distance, breaks)],
      distend   = distbegin + bin_width
    )
}


bin_distances_cut <- function(data, cutpoints) {
  truncation <- max(cutpoints)

  data |>
    filter(distance <= truncation) |>
    mutate(distbegin = cutpoints[findInterval(distance, cutpoints)],
           distend   = cutpoints[findInterval(distance, cutpoints) + 1]
    )
}


##---------------------------------------
## Availability and LL functions
##---------------------------------------

hn_func<- function(x, sigma) {exp(-x^2/(2*sigma^2))}


availability_cont <- function(x, w, n=1) {
  # equation  5 for continuous distances
  return((2*x*n)/w^2 * (1 - (x/w)^2)^(n-1))
}

availability_bins <- function(bin_start, bin_end, w, n=1) {
  # equation 6 for binned distances
  return((1 - bin_start^2/w^2)^n - (1 - bin_end^2/w^2)^n)
}

##---------------------------------------
bin_probs_hn <- function(breaks, sigma, gs) {
  # probabilities for each bin given half-normal detection
  # Need to integrate Eq. 5 over each bin interval to remove
  # bias for larger bins.
  integrand <- function(r, sigma, w, gs) {
    # product of availability, given group size (gs) and detection
    availability_cont(r, w, gs) * hn_func(r, sigma)
  }
  K<- length(breaks) - 1
  bin_probs <- numeric(K)
  w<- max(breaks)
  for (j in 1:K) {
    # integrate det func over each bin interval (lower, upper)
    bin_probs[j] <- integrate(integrand, sigma=sigma, w=w, gs=gs,
                              lower = breaks[j], upper = breaks[j+1])$value
  }
  return(bin_probs)
}

##---- conditional likelihood for continuous data ----

nll.cond.point.hn <- function(parm, x, w, gs){
  # HN detection function
  sigma <- exp(parm)
  p <- hn_func(x, sigma)
  A<- availability_cont(x, w, gs)
  intergrand<- function(r, sigma, w, gs) {
    availability_cont(r, w, gs) * hn_func(r, sigma)
  }
  gs_unique <- unique(gs)
  pbar_unique <- vapply(gs_unique, function(g) {
    integrate(intergrand, 0, w, sigma = sigma, w = w, gs = g)$value
  }, 1)
  pbar <- pbar_unique[match(gs, gs_unique)]
  if (any(!is.finite(pbar)) || any(pbar <= 0)) return(1e10)
  nll <- (-1)*sum(log(p*A/pbar))
  if (!is.finite(nll)) return(1e10)
  return(nll)
}

##---- conditional likelihood for binned data -----
nll.cond.binned.hn <- function(parm, counts, gs, breaks){
  #gs is now a vector
  sigma <- exp(parm)
  n_group_sizes<- length(gs)
  if(!is.list(counts)) counts<- list(counts)
  if(n_group_sizes != length(counts)) stop("error")
  nll<- rep(NA, n_group_sizes)
  for(i in 1:n_group_sizes) {
    grp_counts<- counts[[i]]
    pd<- bin_probs_hn(breaks, sigma, gs[i])
    cp<- pd/sum(pd) # cp must sum to 1
    nll[i] <- sum(grp_counts * log(cp))
  }
  if(any(!is.finite(nll))) return(1e10)
  return(-sum(nll))
}

##---- Detection function-----------------------

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

##----- Estimate density ctds--------------------------

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
  pbar <- function(sigma, w, gs) {
    # product of availability, given group size (gs) and detection
    integrate(function(r) availability_cont(r, w, gs) * hn_func(r, sigma), 0 , w)$value
  }

  a <- theta * pbar(fit$sigma, w = w, gs = 1) # per-camera effective area
  D <- n_total / (K * a)

  # Encounter-rate variance (P2)
  R <- mean(counts)
  var_er <- sum((counts - R)^2) / (K * (K-1))
  cv2_er<- var_er/R^2

  # Detection-function variance (delta method)
  pbar_fun <- function(lsigma, w){pbar(exp(lsigma), w = w, gs = 1)}
  gradient <- numDeriv::grad(pbar_fun, x = fit$log_sigma, w = w)
  se_det <- gradient * fit$se_log_sigma
  cv2_det <- (se_det/pbar(fit$sigma, w, gs = 1))^2

  cv2 <- cv2_er + cv2_det
  se_D <- D * sqrt(cv2)

  z <- qnorm(1 - (1 - level) / 2)
  c_mult <- exp(z * sqrt(log(1 + cv2)))

  list(
    D = D,
    n = n_total,
    K = K,
    eff_area = a,
    sigma = fit$sigma,
    se = se_D,
    cv = sqrt(cv2),
    cv_encounter = sqrt(cv2_er),
    cv_detection = sqrt(cv2_det),
    lcl = D / c_mult,
    ucl = D * c_mult
  )
}

##---- Estimate Density Closest-----------------------------
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
  pbar <- function(sigma, w, gs) {
    # product of availability, given group size (gs) and detection
    integrate(function(r) availability_cont(r, w, gs) * hn_func(r, sigma), 0 , w)$value
  }

  ak<- rep(NA_real_, K)
  nk<- rep(NA_real_, K)
  c_nonzero<- counts
  c_nonzero[c_nonzero < 1]<- 1 # Zero counts get group size 1 for areas
  med_gs<- median(counts[counts > 0])  # median group size (for variance calcs)

    # effective detection probability now depends on group size
  for(k in seq_len(K)) {
    ak[k] <- pbar(fit$sigma, w, gs=c_nonzero[k])
    nk[k] <- counts[k]/ak[k]
  }

  D<- sum(nk)/(K * theta)

  # Encounter-rate variance (R2)
  R <- mean(counts)
  var_er <- sum((counts - R)^2) / (K * (K-1))
  cv2_er<- var_er/R^2

  # Detection-function variance (delta method)
  pbar_fun <- function(lsigma, w, gs){pbar(exp(lsigma), w, gs)}
  gradient <- numDeriv::grad(pbar_fun, x = fit$log_sigma, w = w, gs = med_gs)
  se_det <- gradient * fit$se_log_sigma
  cv2_det <- (se_det/pbar(fit$sigma,w = w,gs = med_gs))^2

  cv2 <- cv2_er + cv2_det
  se_D <- D * sqrt(cv2)

  z <- qnorm(1 - (1 - level) / 2)
  c_mult <- exp(z * sqrt(log(1 + cv2)))

  list(
    D = D,
    n = n_total,
    K = K,
    eff_area = mean(ak),
    sigma = fit$sigma,
    se = se_D,
    cv = sqrt(cv2),
    cv_encounter = sqrt(cv2_er),
    cv_detection = sqrt(cv2_det),
    lcl = D / c_mult,
    ucl = D * c_mult
  )
}

##---- standard CTDS with closest distance trigger ------------------

run_density_ctds<- function(n_rep = 5,
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
    res[[i]] <- sim_ctds(D_true, sigma_closest, sigma_true,
                         w, fov, n_cam, width, height,
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


##---- closest distance only with group size adjustment----------------

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

##---- standard CTDS with time lapse (no sensor) -----------------------

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


##---- Summaries and plots -------------------------------------

summarise_density_sim <- function(res) {
  truth <- attr(res, "truth")
  tibble(
    Model          = truth$model,
    distribution   = truth$distribution,
    n_cam          = truth$n_cam,
    replicates     = length(res$D[!is.na(res$D)]),
    n_total        = mean(res$n_total, na.rm=TRUE),
    mean_n         = mean(res$mean_n, na.rm=TRUE),
    D_true         = truth$D_true,
    D_hat          = median(res$D, na.rm=TRUE),
    bias           = 100 *(median(res$D, na.rm=TRUE) - truth$D_true) / truth$D_true,
    mean_cv_enc    = mean(res$cv_encounter, na.rm=TRUE),
    mean_cv_det    = mean(res$cv_detection, na.rm=TRUE),
    mean_CV        = mean(res$se/res$D, na.rm=TRUE),
    coverage_95    = mean(res$cover,na.rm=TRUE)
  )
}

