
##--------------------------------------
## Data generation
##---------------------------------------
get_closest_dist<- function(n, w) {
  # just simulating r directly using polar identity...
  r <- sqrt(runif(n)) * w
  return(min(r))
}
##----------------------------------------
dunif_circle <- function(n, w) {
  # polar identity
  r <- sqrt(runif(n)) * w
  return(r)
}
##--------------------------------------------------------------
generate_closest_dist<- function(entries, grp_size, w) {
  # for generating closest distances with no detection component
  counts <- rep(grp_size, entries)
  d<- sapply(counts, function(x) get_closest_dist(n = x, w = w))
  return(d)
}
##--------------------------------------
get_detections<- function(n, gr, w, closest=TRUE, ...) {
  # Could just simulate d directly but for completeness...
  r<- sort(dunif_circle(n, w))
  p<- gr(r[1], ...)
    detected<- rbinom(1, 1, p)
    if(detected == 1) {
      if(closest) return(r[1])
      else return(r)
    } else return(0)
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
generate_detections_r <- function(entries_per_group, gr, w, closest = TRUE,
                                  binned = FALSE, breaks = NULL, ...) {

  # Generation of detection data given a detection function gr(),
  # returns either the closest distance of a group (closest=TRUE) or all distances.
  # Handles continuous and binned distances
  if (binned && is.null(breaks)) stop("binning requested but no breaks supplied")

  G <- length(entries_per_group)
  out <- vector("list", G)
  for (g in seq_len(G)) {
    # per-entry chunks (accumulate per entry)
    e_chunks <- vector("list", entries_per_group[g])
    d_chunks <- vector("list", entries_per_group[g])
    for (j in seq_len(entries_per_group[g])) {
      v <- get_detections(g, gr, w, closest = closest, ...)
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
  if (binned) {
    for (g in seq_along(out)) {
      if (!is.null(out[[g]])) {
        out[[g]] <- bin_matrix(out[[g]], "distance", breaks = breaks)
      }
    }
  }
  out
}

##------------------------------------------------------------------------
generate_detections <- function(entries_per_group, gr, w, closest = TRUE,
                                binned = FALSE, breaks = NULL, ...) {
  # Generation of detection data given a detection function gr().
  # Uses a fast C++ backend .

  if (binned && is.null(breaks)) stop("binning requested but no breaks supplied")

  # single-argument closure so C++ only sees gr(r)
    gr_fn <- function(x) gr(x, ...)

    out <- generate_detections_cpp(
      entries_per_group = as.integer(entries_per_group),
      gr      = gr_fn,
      w       = w,
      closest = closest
    )

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

