
##--------------------------------------
## Data generation
##---------------------------------------
get_closest_dist<- function(n, w, theta_max=40) {
  # Could just simulate d directly but for completeness...
  Pi <- base::pi
  theta_max<- theta_max * Pi/180
  rho <- sqrt(runif(n)) * w
  theta <- runif(n, 0, theta_max)
  x <- rho * cos(theta)
  y <- rho * sin(theta)
  d<- sqrt(x^2 + y^2)
  return(min(d))
}
##--------------------------------------------------------------
generate_closest_dist<- function(entries, grp_size, w) {
  # for generating closest distances with no detection component
  counts <- rep(grp_size, entries)
  d<- sapply(counts, function(x) get_closest_dist(n=x, w = w))
  return(d)
}
##--------------------------------------
get_detections<- function(n, gr, w, closest=TRUE, theta_max=40, any=FALSE,
                          clustered=FALSE, ...) {
  # Could just simulate d directly but for completeness...
  Pi <- base::pi
  theta_max<- theta_max * Pi/180 # Camera FOV in radians
  if(n > 1 && clustered)
    pts<- dcluster_sector(n, w, theta_max = theta_max)
  else
    pts<- dunif_sector(n, w, theta_max)
  r<- sort(pts$r)
  if(any) {
    p<- gr(r, ...)
    detected<- rbinom(n, 1, p)
    if(any(detected == 1)) {
      if(closest) return(r[1])
      else return(r)
    } else return(0)
  } else {
    p<- gr(r[1], ...)
    detected<- rbinom(1, 1, p)
    if(detected == 1) {
      if(closest) return(r[1])
      else return(r)
    } else return(0)
  }
}
##----------------------------------------
dunif_sector <- function(n, w, theta_max) {
   # theta_max should be in radians
  rho <- sqrt(runif(n)) * w
  theta <- runif(n, 0, theta_max)
  x <- rho * cos(theta)
  y <- rho * sin(theta)
  r<- sqrt(x^2 + y^2)
  data.frame(x=x,y=y,r=r,theta=theta)
}
##----------------------------------------
inside_sector <- function(x, y, w, theta_max) {
  r <- sqrt(x^2 + y^2)
  theta <- atan2(y, x)
  theta <- ifelse(theta < 0, theta + theta_max, theta)
  r <= w && theta >= 0 && theta <= theta_max
}
##----------------------------------------
dcluster_sector <- function(n, w, theta_max, n_clusters = 1, cluster_sd = 2,
                            max_iter = 10000) {
  # Generate cluster centres uniformly within the circle
  # theta_max should be in radians
  clusters <- dunif_sector(n_clusters, w, theta_max)
  cluster_id <- sample(seq_len(n_clusters), n, replace = TRUE)

  out <- data.frame(
    x = numeric(n),
    y = numeric(n),
    cluster = cluster_id
  )
  for (i in seq_len(n)) {
    ok <- FALSE
    iter <- 0
    while (!ok && iter < max_iter) {
      iter <- iter + 1
      cx <- clusters$x[cluster_id[i]]
      cy <- clusters$y[cluster_id[i]]
      x <- rnorm(1, cx, cluster_sd)
      y <- rnorm(1, cy, cluster_sd)
      if (inside_sector(x, y, w, theta_max)) {
        ok <- TRUE
        out$x[i] <- x
        out$y[i] <- y
      }
    }
    if (!ok) {
      stop("Failed to generate point inside circle. Try reducing cluster_sd.")
    }
  }
  out$r <- sqrt(out$x^2 + out$y^2)
  out$theta <- atan2(out$y, out$x)
  out
}

#-----------------------------------------------------------
bin_matrix <- function(data, dist_col, breaks,
                            right = FALSE, include_lowest = TRUE) {
  ## bin data keeping site hierarchy
  dist <- data[[dist_col]]
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
##------------------------------------------------------------------------
generate_detections_r <- function(entries_per_group, gr, w, closest = TRUE, any=FALSE,
                                clustered=FALSE, ...) {
  # used when clustered=TRUE (C++ only handles uniform case)
  # Generation of detection data given a detection function gr(),  returns either the
  # closest distance of a group (closest=TRUE) or all distances. Detection can be conditional on the
  # closest distance only (any=FALSE) or any of the distances in the group (any=TRUE)
  # Handles continuous and binned distances
  G <- length(entries_per_group)
  out <- vector("list", G)
  for (g in seq_len(G)) {
    # per-entry chunks (accumulate per entry)
    e_chunks <- vector("list", entries_per_group[g])
    d_chunks <- vector("list", entries_per_group[g])
    for (j in seq_len(entries_per_group[g])) {
      v <- get_detections(g, gr, w, ..., closest = closest, any=any, clustered = clustered)
      v <- v[v > 0]
      if (length(v) == 0L) next
      e_chunks[[j]] <- rep.int(j, length(v))
      d_chunks[[j]] <- as.numeric(v)
    }
    entry <- unlist(e_chunks, use.names = FALSE)
    distance <- unlist(d_chunks, use.names = FALSE)
    if(is.null(entry) | is.null(distance)) next
    df <- data.frame(entry = entry, group = g, distance = distance)
    out[[g]] <- df
  }
  out
}

##------------------------------------------------------------------------
generate_detections <- function(entries_per_group, gr, w, closest = TRUE,
                                binned = FALSE, breaks = NULL,
                                clustered = FALSE, any = FALSE, ...) {
  # Generation of detection data given a detection function gr().
  #
  # Detection modes:
  #   any=FALSE (default): group detected if the closest individual detected through gr();
  #                        return closest distance (closest=TRUE) or all distances.
  #   any=TRUE           : each individual detected independently using gr(r_k);
  #                        return closest detected distance (closest=TRUE) or all
  #                        distances.
  #
  # Uses a fast C++ backend for the non-clustered case; falls back to pure R
  # when clustered=TRUE.

  if (binned && is.null(breaks)) stop("binning requested but no breaks supplied")

  if (clustered) {
    out<- generate_detections_r(entries_per_group, gr, w, closest,
                                  any, clustered, ...)
  } else {
  # single-argument closure so C++ only sees gr(r)
    gr_fn <- function(x) gr(x, ...)

    out <- generate_detections_cpp(
      entries_per_group = as.integer(entries_per_group),
      gr      = gr_fn,
      w       = w,
      closest = closest,
      any     = any
    )
}
  if (binned) {
    for (g in seq_along(out)) {
      if (!is.null(out[[g]])) {
        out[[g]] <- bin_matrix(out[[g]], "distance", breaks = breaks)
      }
    }
  }
  out
}

##---------------------------------------
## Availability and LL functions
##---------------------------------------

hn_func<- function(x, sigma) {exp(-x^2/(2*sigma^2))}


availability_bin_mid <- function(delta, midpoint, w) {
  # Equation 10 - not used directly
  return((2 * delta * midpoint)/w^2)
}

availability_cont <- function(x, w, n=1) {
  # equation  5
  return((2*x*n)/w^2 * (1 - (x/w)^2)^(n-1))
}

availability_bins <- function(bin_start, bin_end, w, n=1) {
  # equation 6
 return((1 - bin_start^2/w^2)^n - (1 - bin_end^2/w^2)^n)
}


##---------------------------------------
bin_probs_hn <- function(breaks, sigma, gs) {
  # probabilities for each bin given half-normal detection
  # Need to integrate over each bin interval
  integrand <- function(r, sigma, w, gs) {
    # product of availability, given group size and detection
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

##------------------------------------------------------------
# (3) Define conditional likelihood for continuous data
nll.cond.point.hn <- function(parm, x, w, gs){
  # HN detection function
  sigma <- exp(parm)
  p <- hn_func(x, sigma)
  A<- availability_cont(x, w, gs)
  intergrand<- function(r, sigma, w, gs) {
    availability_cont(r, w, gs) * hn_func(r, sigma)
  }
  pbar <- integrate(intergrand, 0, w, sigma=sigma, w=w, gs=gs)$value
  LL <- sum(log(p*A/pbar))
  return(-LL)
}

##------------------------------
nll.cond.binned.hn <- function(parm, counts, gs, breaks){
  #gs is now a vector
  sigma <- exp(parm)
  n_group_sizes<- length(gs)
  if(!is.list(counts)) counts<- list(counts)
  if(n_group_sizes != length(counts)) stop("error")
  LL<- rep(NA, n_group_sizes)
  for(i in 1:n_group_sizes) {
    grp_counts<- counts[[i]]
    pd<- bin_probs_hn(breaks, sigma, gs[i])
    cp<- pd/sum(pd) # cp must sum to 1
    LL[i] <- sum(grp_counts * log(cp))
  }
  return(-sum(LL))
}

