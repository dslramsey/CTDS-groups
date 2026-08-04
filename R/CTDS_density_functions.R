
##----------------------------------------
simulate_animals <- function(n,
                             xlim = c(0, 1),
                             ylim = c(0, 1),
                             distribution = c("uniform", "clustered"),
                             n_clusters = 5,
                             cluster_sd = NULL) {
  # Simulate a 2D point pattern of n animals within the rectangle
  # defined by xlim/ylim.
  #   distribution = "uniform"   : complete spatial randomness
  #   distribution = "clustered" : Thomas-style process; parents (cluster
  #                                centres) are uniform, offspring are Gaussian
  #                                scattered about a parent with sd cluster_sd
  # Returns a data.frame(x, y, cluster); cluster is NA for the uniform case.
  distribution <- match.arg(distribution)

  if (distribution == "uniform") {
    return(data.frame(
      x = runif(n, xlim[1], xlim[2]),
      y = runif(n, ylim[1], ylim[2]),
      cluster = NA_integer_
    ))
  }

  # Clustered: default scatter is 5% of the smaller region dimension
  if (is.null(cluster_sd)) {
    cluster_sd <- 0.05 * min(diff(xlim), diff(ylim))
  }

  # Cluster centres (parents) placed uniformly
  cx <- runif(n_clusters, xlim[1], xlim[2])
  cy <- runif(n_clusters, ylim[1], ylim[2])

  # Assign each animal to a cluster, then scatter around its parent.
  cluster_id <- sample.int(n_clusters, n, replace = TRUE)
  x <- rnorm(n, cx[cluster_id], cluster_sd)
  y <- rnorm(n, cy[cluster_id], cluster_sd)

  # Vectorized rejection: resample any points that fall outside the rectangle
  outside <- x < xlim[1] | x > xlim[2] | y < ylim[1] | y > ylim[2]
  while (any(outside)) {
    idx <- which(outside)
    x[idx] <- rnorm(length(idx), cx[cluster_id[idx]], cluster_sd)
    y[idx] <- rnorm(length(idx), cy[cluster_id[idx]], cluster_sd)
    outside <- x < xlim[1] | x > xlim[2] | y < ylim[1] | y > ylim[2]
  }

  data.frame(x = x, y = y, cluster = cluster_id)
}
##-----------------------------------------------
generate_animals <- function(
    width,
    height,
    density,
    distribution = c("uniform", "clustered"),
    cluster_radius = 50,
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


##----------------------------------------
sample_sector <- function(points,
                          origin = c(0, 0),
                          bearing = 0,
                          radius,
                          angle,
                          gr = NULL,
                          ...) {
  # Sample (detect) animals falling within a camera FOV sector.
  #   points  : data.frame with columns x, y (e.g. from simulate_animals)
  #   origin  : c(x, y) camera location
  #   bearing : central viewing direction in degrees (0 = +x axis, CCW)
  #   radius  : detection radius of the sector (w)
  #   angle   : full angular width of the FOV in degrees (split +/- about bearing)
  #   gr      : optional detection function gr(r, ...) returning P(detect) in
  #             [0, 1]. When supplied, animals inside the sector are thinned by
  #             a Bernoulli(gr(r)) draw. When NULL, detection is perfect within
  #             the sector.
  #   ...     : additional parameters passed to gr() (e.g. sigma).
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

  # Thin by the detection function: keep each in-view animal with prob gr(r)
  if (!is.null(gr)) {
    p <- gr(out$r, ...)
    out$p <- p
    out <- out[rbinom(nrow(out), 1, p) == 1, , drop = FALSE]
  }

  out
}

##----------------------------------------
fit_detection_hn <- function(dist, w) {
  # Fit a half-normal detection function to observed sector distances by
  # conditional MLE (point/sector transect; group size 1). Reuses
  # nll.cond.point.hn, whose conditional likelihood does not depend on the
  # sector angle. Returns sigma with an SE on the log scale (from the Hessian).
  mle <- optim(
    par = log(w / 2),
    fn = nll.cond.point.hn,
    x = dist, w = w, gs = 1,
    method = "Brent", lower = -5, upper = 5,
    hessian = TRUE
  )
  list(
    sigma        = exp(mle$par),
    log_sigma    = mle$par,
    se_log_sigma = sqrt(1 / mle$hessian[1, 1]),
    mle          = mle
  )
}

##----------------------------------------
estimate_density <- function(dist, w, angle, fit, level = 0.95) {
  # Estimate density of individuals from sector distance-sampling detections.
  #   dist  : detected distances (from sample_sector)
  #   w     : sector radius; angle: full FOV width in degrees
  #   fit   : output of fit_detection_hn()
  # Effective sampled area a = theta * integral_0^w r * gr(r) dr, i.e. the
  # sector area A = 0.5*theta*w^2 times the mean detection probability. Density
  # D = n / a. SE combines a Poisson encounter component with detection-function
  # uncertainty via the delta method (assumes a single homogeneous sector).
  n <- length(dist)
  theta <- angle * pi / 180
  m <- function(sigma) integrate(function(r) r * hn_func(r, sigma), 0, w)$value

  a <- theta * m(fit$sigma)
  D <- n / a

  # d log(a) / d log(sigma) by central difference
  h <- 1e-4
  dloga_dlogsig <- (log(theta * m(exp(fit$log_sigma + h))) -
                      log(theta * m(exp(fit$log_sigma - h)))) / (2 * h)

  cv2 <- 1 / n + (dloga_dlogsig * fit$se_log_sigma)^2
  se_D <- D * sqrt(cv2)

  # Lognormal CI
  z <- qnorm(1 - (1 - level) / 2)
  c_mult <- exp(z * sqrt(log(1 + cv2)))

  list(
    D = D, n = n, eff_area = a, sigma = fit$sigma,
    se = se_D, cv = sqrt(cv2),
    lcl = D / c_mult, ucl = D * c_mult
  )
}

##----------------------------------------
estimate_density_multi <- function(counts, dist, w, angle, fit, level = 0.95) {
  # Density from several independent camera sectors of equal geometry.
  #   counts : per-camera detection counts (length K, may include zeros)
  #   dist   : pooled detected distances across all cameras (for the fit)
  #   fit    : output of fit_detection_hn() on the pooled distances
  # Per-camera effective area a = theta * integral_0^w r*gr(r) dr; total
  # effort is K*a, so D = sum(counts) / (K*a). The encounter-rate variance is
  # estimated EMPIRICALLY from the variation in counts among cameras
  # (var(mean count) = s^2 / K), so it captures over-dispersion under
  # clustering rather than assuming Poisson. The detection-function
  # contribution is added via the delta method.
  K <- length(counts)
  if (K < 2) stop("need at least 2 cameras for an empirical variance")

  n_total <- sum(counts)
  theta <- angle * pi / 180
  m <- function(sigma) integrate(function(r) r * hn_func(r, sigma), 0, w)$value

  a <- theta * m(fit$sigma)        # per-camera effective area
  D <- n_total / (K * a)

  # Encounter-rate component: CV^2 of the mean count across cameras
  R <- mean(counts)
  cv2_er <- (K/(n_total^2 * (K - 1))) * sum((counts - R)^2)


  # Detection-function component via delta method on log(sigma)
  h <- 1e-4
  dloga <- (log(theta * m(exp(fit$log_sigma + h))) -
              log(theta * m(exp(fit$log_sigma - h)))) / (2 * h)
  cv2_det <- (dloga * fit$se_log_sigma)^2

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

##---------------------------------------
## One replicate (multiple cameras)
##---------------------------------------
sim_once <- function(D_true, sigma_true, w, fov, n_cam, width, height,
                     distribution = "uniform", cluster_radius = 30,
                     mean_cluster_size = 5, min_total = 10,
                     camera_layout = c("random", "grid")) {
  camera_layout <- match.arg(camera_layout)
  # Simulate one field, survey it with n_cam random cameras, fit and estimate.
  # Returns a one-row data.frame, or NULL if too few pooled detections to fit.
  #D_true<- D_km2/1e6 # D_true is density per m2
  n_animals <- round(D_true * width * height)

  # animals <- simulate_animals(n_animals, xlim, ylim, distribution,
  #                             n_clusters = n_clusters, cluster_sd = cluster_sd)
  animals<- generate_animals(width = width,
                             height = height,
                             density = D_true,
                             distribution = distribution,
                             cluster_radius = cluster_radius,
                             mean_cluster_size = mean_cluster_size)

  # Camera locations, kept >= w from every edge so each sector stays
  # fully inside the region for any bearing (requires region wider than 2w).
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
#  bearing <- runif(n_cam, 0, 360)
  bearing <- rep(270, n_cam)

  counts <- integer(n_cam)
  dist_list <- vector("list", n_cam)
  for (k in seq_len(n_cam)) {
    dk <- sample_sector(animals, origin = c(cx[k], cy[k]), bearing = bearing[k],
                        radius = w, angle = fov, gr = hn_func, sigma = sigma_true)
    counts[k] <- nrow(dk)
    dist_list[[k]] <- dk$r
  }
  dist <- unlist(dist_list)
  if (length(dist) < min_total) return(NULL)

  fit <- fit_detection_hn(dist, w = w)
  est <- estimate_density_multi(counts, dist, w = w, angle = fov, fit = fit)

  data.frame(
    n_total      = est$n,
    mean_n       = mean(counts),
    sigma        = fit$sigma,
    D            = est$D,
    se           = est$se,
    cv_encounter = est$cv_encounter,
    cv_detection = est$cv_detection,
    lcl          = est$lcl,
    ucl          = est$ucl,
    cover        = est$lcl <= D_true & D_true <= est$ucl
  )
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
##---------------------------------------
## Many replicates
##---------------------------------------
run_density_sim <- function(n_rep = 5,
                            D_true = 0.01,
                            sigma_true = 4,
                            w = 12,
                            fov = 40,
                            n_cam = 60,
                            width = 500,
                            height = 500,
                            distribution = c("uniform", "clustered"),
                            cluster_radius = 30,
                            mean_cluster_size = 5,
                            camera_layout = c("random", "grid"),
                            progress = TRUE) {
  # Each replicate places n_cam random cameras; the encounter-rate variance is
  # estimated across those cameras. Under clustering the empirical variance
  # captures over-dispersion, so CI coverage should recover toward nominal.
  distribution <- match.arg(distribution)
  camera_layout <- match.arg(camera_layout)
  res <- vector("list", n_rep)
  pb <- if (progress) utils::txtProgressBar(max = n_rep, style = 3) else NULL
  for (i in seq_len(n_rep)) {
    res[[i]] <- sim_once(D_true, sigma_true, w, fov, n_cam, width, height,
                         distribution = distribution,
                         cluster_radius = cluster_radius,
                         mean_cluster_size = mean_cluster_size,
                         camera_layout = camera_layout)
    if (progress) utils::setTxtProgressBar(pb, i)
  }
  if (progress) close(pb)
  out <- bind_rows(res)
  attr(out, "truth") <- list(D_true = D_true, sigma_true = sigma_true,
                             w = w, fov = fov, n_cam = n_cam,
                             distribution = distribution)
  out
}


##---------------------------------------
## Summaries and plots
##---------------------------------------
summarise_density_sim <- function(res) {
  truth <- attr(res, "truth")
  tibble(
    distribution   = truth$distribution,
    n_cam          = truth$n_cam,
    replicates     = nrow(res),
    D_true         = truth$D_true,
    mean_Dhat      = mean(res$D),
    pct_bias       = 100 * (mean(res$D) - truth$D_true) / truth$D_true,
    emp_SD         = sd(res$D),
    mean_est_SE    = mean(res$se),
    coverage_95    = mean(res$cover)
  )
}

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
                                 gr = NULL,
                                 ...) {
  # Plot a 2D animal point pattern with a set of camera FOV sectors overlaid.
  #   animals : data.frame with columns x, y (e.g. from generate_animals)
  #   cams    : data.frame with columns cx, cy, bearing (one row per camera);
  #             an id column is added if absent
  #   radius  : sector detection radius (w); angle: full FOV width in degrees
  #   width, height : optional plot limits (defaults to the animal extent)
  #   show_detected : when TRUE, animals falling inside any sector are coloured.
  #             Detection is assessed with sample_sector() using the same
  #             geometry; pass gr (and its params via ...) to thin by a
  #             detection function, otherwise detection is perfect within view.
  if (is.null(cams$id)) cams$id <- seq_len(nrow(cams))
  if (is.null(width))  width  <- max(animals$x)
  if (is.null(height)) height <- max(animals$y)

  sectors <- cams |>
    purrr::pmap(\(id, cx, cy, bearing, ...)
                sector_polygon(id, c(cx, cy), bearing,
                               radius = radius, angle = angle)) |>
    dplyr::bind_rows()

  if (show_detected) {
    # Flag animals detected by any sector. Track original row ids so the same
    # animal seen by multiple cameras is not double-counted.
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
                   fill = "firebrick", alpha = 0.15, color = "firebrick") +
      geom_point(aes(color = detected, alpha = detected), size = 1) +
      scale_color_manual(values = c(`FALSE` = "grey40", `TRUE` = "firebrick"),
                         labels = c(`FALSE` = "no", `TRUE` = "yes")) +
      scale_alpha_manual(values = c(`FALSE` = 0.3, `TRUE` = 0.9), guide = "none") +
      labs(color = "detected")
  } else {
    p <- p +
      geom_point(alpha = 0.4, size = 1) +
      geom_polygon(data = sectors, aes(group = id),
                   fill = "firebrick", alpha = 0.3, color = "firebrick")
  }

  p +
    geom_point(data = cams, aes(cx, cy), color = "firebrick", size = 1.5) +
    coord_equal(xlim = c(0, width), ylim = c(0, height)) +
    labs(title = "Animal distribution with camera sectors",
         x = "x (m)", y = "y (m)") +
    theme_bw()
}
