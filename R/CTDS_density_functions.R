
library(numDeriv)

## ----------  TOP ----------------------------------------

#' Generate a simulated animal population
#'
#' Simulates animal locations within a rectangular study region, either
#' uniformly at random or aggregated in Gaussian clusters. The number of
#' animals is Poisson with mean \code{density * width * height}.
#'
#' @param width,height Numeric. Dimensions of the rectangular study region.
#' @param density Numeric. Expected animal density (animals per unit area).
#' @param distribution Character. Spatial distribution of animals: \code{"uniform"}
#'   (complete spatial randomness) or \code{"clustered"} (Gaussian clusters).
#' @param cluster_radius Numeric. Standard deviation of animal displacements
#'   around cluster centres (clustered distribution only).
#' @param mean_cluster_size Numeric. Rough mean number of animals per cluster;
#'   used to set the number of cluster centres.
#'
#' @return A \code{data.frame} with columns \code{x}, \code{y} and a cluster
#'   identifier (\code{cluster_id} for uniform, \code{cluster} for clustered).
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
#' Select the best point-transect detection function by AIC
#'
#' Fits a set of candidate detection functions (half-normal, hazard-rate and
#' uniform keys, each with cosine adjustments) to radial distance data using
#' \code{\link[Distance]{ds}}, and returns the model with the lowest AIC.
#' Models that fail to converge are dropped.
#'
#' @param dist_ds A \code{data.frame} of detected objects in \code{Distance}
#'   flatfile format, containing at least a \code{distance} column (or
#'   \code{distbegin}/\code{distend} when binned).
#' @param w Numeric. Truncation distance (sector radius).
#' @param binned Logical. If \code{TRUE}, distances are first binned using
#'   \code{breaks} and fitted as binned point-transect data.
#' @param breaks Numeric vector of bin cutpoints, required when
#'   \code{binned = TRUE}.
#'
#' @return A fitted \code{ds} model object with the lowest AIC among the
#'   candidates. Errors if all candidate models fail.
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

#' Simulate one camera trap distance sampling (CTDS) replicate
#'
#' Simulates sensor-triggered CTDS: an animal population is generated,
#' \code{n_cam} cameras are placed at random locations (all facing south),
#' and a camera is triggered by the closest individual with probability given
#' by a half-normal detection function. Density is estimated by fitting a
#' point-transect detection function with \code{\link[Distance]{ds}} (chosen
#' by AIC) and estimating abundance with \code{\link[Distance]{dht2}}.
#'
#' @param D_true Numeric. True animal density (animals per unit area).
#' @param sigma_closest Numeric. Scale of the half-normal function governing
#'   detection (triggering) of the closest individual.
#' @param sigma_true Numeric or \code{NULL}. Scale of the half-normal
#'   detection function for other individuals within the field of view; only
#'   used when \code{perfect_fov = FALSE}.
#' @param w Numeric. Truncation distance (sector radius).
#' @param fov Numeric. Camera field of view in degrees (full angular width).
#' @param n_cam Integer. Number of camera locations.
#' @param width,height Numeric. Dimensions of the rectangular study region.
#' @param distribution Character. \code{"uniform"} or \code{"clustered"};
#'   passed to \code{\link{generate_animals}}.
#' @param perfect_fov Logical. If \code{TRUE}, all animals in the sector are
#'   recorded once the camera triggers; if \code{FALSE}, further animals are
#'   detected with half-normal probability using \code{sigma_true}.
#' @param cluster_radius,mean_cluster_size Passed to
#'   \code{\link{generate_animals}} for clustered distributions.
#' @param min_total Integer. Minimum total detections required; replicates
#'   with fewer return \code{NULL}.
#' @param binned Logical. If \code{TRUE}, distances are binned before fitting.
#' @param breaks Numeric vector of bin cutpoints, used when
#'   \code{binned = TRUE}.
#'
#' @return A one-row \code{data.frame} with the number of detections
#'   (\code{n_total}), mean count per triggered camera (\code{mean_n}),
#'   density estimate (\code{D}) with standard error, CV components
#'   (\code{cv_encounter}, \code{cv_detection}, \code{cv}), confidence limits
#'   (\code{lcl}, \code{ucl}) and whether the interval covers \code{D_true}
#'   (\code{cover}); \code{NA}s if model fitting fails. Returns \code{NULL} if
#'   fewer than \code{min_total} detections are obtained.
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

#' Simulate one closest-distance CTDS replicate
#'
#' Simulates sensor-triggered camera surveys where only the distance to the
#' closest (triggering) individual at each camera is used to fit the detection
#' function, with availability adjusted for the number of individuals (group
#' size) in the field of view. Density is estimated with
#' \code{\link{estimate_density_closest}}.
#'
#' @inheritParams sim_ctds
#'
#' @return A one-row \code{data.frame} with the same columns as
#'   \code{\link{sim_ctds}}; returns \code{NULL} if fewer than
#'   \code{min_total} cameras record a triggering individual.
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

#' Simulate one time-lapse CTDS replicate
#'
#' Simulates time-lapse (non-sensor) camera surveys: at a snapshot moment all
#' animals within each camera sector are recorded (subject to detection error
#' when \code{perfect_fov = FALSE}), and distances to all detected animals
#' are used to fit a half-normal detection function. Density is estimated with
#' \code{\link{estimate_density_ctds}}.
#'
#' @inheritParams sim_ctds
#'
#' @return A one-row \code{data.frame} with the same columns as
#'   \code{\link{sim_ctds}}; returns \code{NULL} if fewer than
#'   \code{min_total} distances are recorded.
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

#' Subset points falling within a camera sector
#'
#' Detects which points (e.g. animals) fall within a circular sector defined
#' by a camera location, viewing direction, radius and angular width.
#'
#' @param points A \code{data.frame} with columns \code{x} and \code{y}, e.g.
#'   from \code{\link{generate_animals}}.
#' @param origin Numeric vector of length 2. Camera location \code{c(x, y)}.
#' @param bearing Numeric. FOV direction in degrees (0 = east, increasing
#'   counter-clockwise in mathematical convention as used by \code{atan2}).
#' @param radius Numeric. Sector radius (truncation distance \code{w}).
#' @param angle Numeric. Full angular width of the FOV in degrees, split
#'   \eqn{\pm} about \code{bearing}.
#'
#' @return The subset of \code{points} inside the sector, with added columns
#'   \code{r} (radial distance from the camera) and \code{theta} (angle
#'   relative to the bearing, radians).
sample_sector <- function(points,
                          origin = c(0, 0),
                          bearing = 0,
                          radius,
                          angle) {

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
#' Tabulate binned distances by group size
#'
#' Bins the triggering (closest) distances separately for each observed group
#' size, for use with \code{\link{nll.cond.binned.hn}}.
#'
#' @param counts Integer vector of group sizes (number of animals in the FOV)
#'   per camera; zeros are ignored.
#' @param dists Numeric vector of closest distances per camera, parallel to
#'   \code{counts} (may contain \code{NA} for non-triggered cameras).
#' @param breaks Numeric vector of bin cutpoints passed to
#'   \code{\link{make_bins}}.
#'
#' @return A list with elements \code{bin_list} (a list of binned count
#'   vectors, one per distinct group size) and \code{grps} (the corresponding
#'   group sizes).
table_counts<- function(counts, dists, breaks) {

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

#' Bin distances into counts
#'
#' Cuts distances into bins defined by \code{breaks} and returns a count per
#' bin, including zero counts for empty bins.
#'
#' @param dist Numeric vector of distances.
#' @param breaks Numeric vector of bin cutpoints passed to \code{\link{cut}}.
#' @param right Logical. Should intervals be closed on the right? Default
#'   \code{FALSE} (closed on the left).
#' @param include_lowest Logical. Should the lowest break be included in the
#'   first bin? Default \code{TRUE}.
#'
#' @return A one-dimensional \code{xtabs} object (named numeric vector) of
#'   counts per bin, with bins in the order given by \code{breaks}.
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


#' Bin exact distances into equal-width intervals
#'
#' Converts exact distances into binned distances with columns
#' \code{distbegin} and \code{distend} (the format expected by
#' \code{\link[Distance]{ds}} for binned data), truncating at
#' \code{truncation}.
#'
#' @param data A \code{data.frame} with a \code{distance} column.
#' @param bin_width Numeric. Width of each distance bin.
#' @param truncation Numeric. Maximum distance retained.
#'
#' @return \code{data} restricted to distances up to \code{truncation}, with
#'   added columns \code{distbegin} and \code{distend} giving the bin
#'   endpoints.
bin_distances_bw <- function(data, bin_width, truncation) {
  breaks <- seq(0, truncation, by = bin_width)

  data |>
    dplyr::filter(distance <= truncation) |>
    dplyr::mutate(
      distbegin = breaks[findInterval(distance, breaks)],
      distend   = distbegin + bin_width
    )
}


#' Bin exact distances into intervals given by cutpoints
#'
#' Like \code{\link{bin_distances_bw}} but with arbitrary bin endpoints,
#' truncating at the largest cutpoint.
#'
#' @param data A \code{data.frame} with a \code{distance} column.
#' @param cutpoints Numeric vector of bin endpoints (must include 0 and the
#'   truncation distance).
#'
#' @return \code{data} restricted to distances up to \code{max(cutpoints)},
#'   with added columns \code{distbegin} and \code{distend} giving the bin
#'   endpoints.
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

#' Half-normal detection function
#'
#' @param x Numeric vector of radial distances.
#' @param sigma Numeric. Scale parameter of the half-normal function.
#'
#' @return Detection probability \eqn{\exp(-x^2 / (2\sigma^2))} at each
#'   distance.
hn_func<- function(x, sigma) {exp(-x^2/(2*sigma^2))}


#' Availability density of the closest of n individuals (continuous)
#'
#' Density of the distance to the closest of \code{n} uniformly distributed
#' individuals within a circular sector of radius \code{w} (equation 5 of the
#' accompanying manuscript). With \code{n = 1} this is the triangular
#' availability density \eqn{2x/w^2}.
#'
#' @param x Numeric vector of radial distances.
#' @param w Numeric. Sector radius (truncation distance).
#' @param n Numeric. Group size (number of individuals in the sector).
#'
#' @return The availability density evaluated at \code{x}.
availability_cont <- function(x, w, n=1) {
  return((2*x*n)/w^2 * (1 - (x/w)^2)^(n-1))
}

#' Availability probability of the closest of n individuals (binned)
#'
#' Probability that the closest of \code{n} uniformly distributed individuals
#' within a sector of radius \code{w} falls in the distance bin
#' \eqn{[bin\_start, bin\_end)} (equation 6 of the accompanying manuscript).
#'
#' @param bin_start,bin_end Numeric. Lower and upper bin endpoints.
#' @param w Numeric. Sector radius (truncation distance).
#' @param n Numeric. Group size (number of individuals in the sector).
#'
#' @return The probability of the closest individual falling in the bin.
availability_bins <- function(bin_start, bin_end, w, n=1) {
  return((1 - bin_start^2/w^2)^n - (1 - bin_end^2/w^2)^n)
}
##---------------------------------------
#' Expected bin probabilities under a half-normal detection function
#'
#' Integrates the product of availability and the half-normal detection
#' function over each distance bin, giving the (unnormalised) probability of
#' a detection falling in each bin.
#'
#' @param breaks Numeric vector of bin cutpoints; the largest value is taken
#'   as the truncation distance \code{w}.
#' @param sigma Numeric. Scale of the half-normal detection function.
#' @param gs Numeric. Group size used in the availability function.
#'
#' @return Numeric vector of length \code{length(breaks) - 1} with the
#'   integrated probability for each bin.
bin_probs_hn <- function(breaks, sigma, gs) {

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

#' Conditional negative log-likelihood for exact distances (half-normal)
#'
#' Negative log-likelihood of observed radial distances conditional on
#' detection, under a half-normal detection function. Group-size-specific
#' availability is allowed via \code{gs}; the normalising constant
#' (mean detection probability) is computed once per distinct group size.
#'
#' @param parm Numeric. Log of the half-normal scale parameter \code{sigma}.
#' @param x Numeric vector of observed radial distances.
#' @param w Numeric. Truncation distance (sector radius).
#' @param gs Numeric vector of group sizes, parallel to \code{x}; use a
#'   scalar or constant vector for no group-size adjustment.
#'
#' @return The conditional negative log-likelihood; returns \code{1e10} if
#'   the likelihood is not finite (for use with \code{\link{optim}}).
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

#' Conditional negative log-likelihood for binned distances (half-normal)
#'
#' Negative multinomial log-likelihood of binned detection counts conditional
#' on detection, under a half-normal detection function. Expected bin
#' probabilities are obtained by integrating availability times detection
#' over each bin (see \code{\link{bin_probs_hn}}) and normalising.
#'
#' @param parm Numeric. Log of the half-normal scale parameter \code{sigma}.
#' @param counts A list of binned count vectors (e.g. from
#'   \code{\link{make_bins}}), one per group size; a single vector is
#'   coerced to a list of length one.
#' @param gs Numeric vector of group sizes, one per element of \code{counts}.
#' @param breaks Numeric vector of bin cutpoints; the largest value is taken
#'   as the truncation distance.
#'
#' @return The conditional negative log-likelihood; returns \code{1e10} if
#'   not finite (for use with \code{\link{optim}}).
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

#' Fit a half-normal detection function by conditional MLE
#'
#' Fits a half-normal detection function to observed radial distances by
#' maximising the conditional likelihood (given detection), using
#' \code{\link{optim}} with the Brent method on the log scale. Works with
#' either exact distances or binned counts.
#'
#' @param dist For exact data, a numeric vector of radial distances; for
#'   binned data, a list of binned count vectors (see
#'   \code{\link{nll.cond.binned.hn}}).
#' @param w Numeric. Truncation distance (sector radius).
#' @param gs Numeric. Group size(s) used for availability: a vector parallel
#'   to \code{dist} for exact data, or one value per element of \code{dist}
#'   for binned data. Default \code{1} (no group-size adjustment).
#' @param binned Logical. Is \code{dist} binned? Requires \code{breaks}.
#' @param breaks Numeric vector of bin cutpoints, required when
#'   \code{binned = TRUE}.
#'
#' @return A list with elements \code{sigma} (estimated half-normal scale),
#'   \code{log_sigma}, \code{se_log_sigma} (standard error on the log scale,
#'   from the Hessian) and \code{mle} (the full \code{\link{optim}} result).
fit_detection_hn <- function(dist, w, gs = 1, binned = FALSE, breaks = NULL) {

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

#' Estimate density from time-lapse camera sectors
#'
#' Design-based density estimate from independent camera sectors of equal
#' area. Encounter-rate variance is estimated empirically across cameras
#' following Fewster et al. (2009, Biometrics, estimator P2); detection-
#' function variance is added via the delta method, and log-normal confidence
#' intervals are constructed.
#'
#' @param counts Integer vector of detection counts per camera (length
#'   \code{K >= 2}, may include zeros).
#' @param w Numeric. Truncation distance (sector radius).
#' @param angle Numeric. Camera field of view in degrees (full angular width).
#' @param fit Output of \code{\link{fit_detection_hn}} fitted to the pooled
#'   distances.
#' @param level Numeric. Confidence level for the interval; default 0.95.
#'
#' @return A list with the density estimate \code{D}, total detections
#'   \code{n}, number of cameras \code{K}, effective per-camera area
#'   \code{eff_area}, estimated \code{sigma}, standard error \code{se},
#'   overall CV and its components (\code{cv}, \code{cv_encounter},
#'   \code{cv_detection}), and confidence limits \code{lcl}, \code{ucl}.
estimate_density_ctds <- function(counts, w, angle, fit, level = 0.95) {

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
#' Estimate density from closest-distance camera sectors
#'
#' Design-based density estimate where each camera contributes the group size
#' (number of individuals in the FOV) of the triggered cluster, and the
#' effective detection probability of the closest individual depends on group
#' size through the availability function. Encounter-rate variance is
#' estimated empirically across cameras (Fewster et al. 2009, Biometrics,
#' estimator R2); detection-function variance is added via the delta method
#' evaluated at the median group size, and log-normal confidence intervals
#' are constructed.
#'
#' @param counts Integer vector of group sizes per camera (length
#'   \code{K >= 2}, may include zeros).
#' @param w Numeric. Truncation distance (sector radius).
#' @param angle Numeric. Camera field of view in degrees (full angular width).
#' @param fit Output of \code{\link{fit_detection_hn}} fitted to the pooled
#'   closest distances (with group-size adjustment).
#' @param level Numeric. Confidence level for the interval; default 0.95.
#'
#' @return A list with the same elements as \code{\link{estimate_density_ctds}};
#'   \code{eff_area} is the mean per-camera effective area.
estimate_density_closest <- function(counts, w, angle, fit, level = 0.95) {

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

#' Run a sensor-triggered CTDS simulation study
#'
#' Repeats \code{\link{sim_ctds}} \code{n_rep} times. Each replicate places
#' \code{n_cam} random cameras; the encounter-rate variance is estimated
#' across those cameras, so under clustering the empirical variance captures
#' over-dispersion and CI coverage should recover toward nominal.
#'
#' @param n_rep Integer. Number of simulation replicates.
#' @inheritParams sim_ctds
#' @param progress Logical. Show a text progress bar? Default \code{TRUE}.
#'
#' @return A \code{data.frame} with one row per successful replicate (columns
#'   as in \code{\link{sim_ctds}}), with an attribute \code{"truth"} holding
#'   the true parameter values and model label (\code{"CTDS"}).
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

#' Run a closest-distance CTDS simulation study
#'
#' Repeats \code{\link{sim_closest}} \code{n_rep} times. Variance is estimated
#' empirically across cameras as in \code{\link{run_density_ctds}}.
#'
#' @inheritParams run_density_ctds
#'
#' @return A \code{data.frame} with one row per successful replicate (columns
#'   as in \code{\link{sim_closest}}), with an attribute \code{"truth"}
#'   holding the true parameter values and model label (\code{"Closest"}).
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

#' Run a time-lapse CTDS simulation study
#'
#' Repeats \code{\link{sim_lapse}} \code{n_rep} times. Variance is estimated
#' empirically across cameras as in \code{\link{run_density_ctds}}.
#'
#' @inheritParams run_density_ctds
#'
#' @return A \code{data.frame} with one row per successful replicate (columns
#'   as in \code{\link{sim_lapse}}), with an attribute \code{"truth"} holding
#'   the true parameter values and model label (\code{"lapse"}).
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

#' Summarise a density simulation study
#'
#' Computes performance summaries (bias, CV components and confidence-interval
#' coverage) from the output of \code{\link{run_density_ctds}},
#' \code{\link{run_density_closest}} or \code{\link{run_density_lapse}}.
#'
#' @param res A \code{data.frame} from one of the \code{run_density_*}
#'   functions, carrying a \code{"truth"} attribute.
#'
#' @return A one-row \code{tibble} with the model label, distribution, number
#'   of cameras, number of successful replicates, mean detections
#'   (\code{n_total}, \code{mean_n}), true and median estimated density,
#'   percent relative \code{bias}, mean encounter/detection/overall CVs and
#'   95\% interval \code{coverage_95}.
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

