# ---------------------------------------------------------------------------
# Sub-population bias under a pooled detection function
# ---------------------------------------------------------------------------
# Two sub-populations with identical true density that differ only in how animals
# are distributed in space: one uniform (homogeneous Poisson), one clustered.
# Clustering raises the number of individuals recorded per detected snapshot moment
# without changing abundance, so the two sub-populations differ in group-size
# composition while being identical in the quantity being estimated.
#
# Both su-bpopulations are surveyed together and analysed with a single pooled
# detection function, as would be common in real-world studies.
#
# Observation process is as in the main simulations (see results 3.1): the closest individual
# triggers the camera following a half-normal with scale sigma_closest (here set to 4m); any
# further individuals in the FOV are then detected following a half-normal with
# scale sigma_fov (here set to 12m), so group size is subject to detection error.
#
# Models (all share one pooled detection function, sub-population = stratum):
#   CTDS         distances to all detected individuals; key selected by AIC
#                from half-normal, hazard-rate and uniform keys with cosine
#                adjustments (select_best_ds)
#   MCDS         as CTDS but with sub-population ID as a covariate on the scale
#                parameter; key selected by AIC from half-normal and
#                hazard-rate. No adjustment terms: ds() cannot enforce
#                monotonicity alongside covariates so adjusted covariate models are
#                not guaranteed to be valid detection functions. A uniform key
#                has no scale parameter for a covariate to act on and is
#                likewise excluded.
#   Closest      distance to the closest detected individual only, with
#                availability adjusted for the number of individuals detected
#   Time lapse   camera records at every snapshot moment irrespective of
#                triggering, with the standard availability calculation
#
#
# Outputs: outputs/res_subpopulation.rds
#          outputs/subpop_field.png, subpop_abundance.png, subpop_rel_bias.png
# ---------------------------------------------------------------------------

library(tidyverse)
library(Distance)

source("R/CTDS_density_functions.R")

##---------------------------------------
## Settings
##---------------------------------------
n_rep         <- 100          # replicates; ~2 h at these settings, reduce to test
D_true        <- 0.01         # true density of BOTH subpopulations (m^-2)
n_cam         <- 1000         # cameras (= snapshot moments) per subpopulation
w             <- 15           # truncation distance (m)
fov           <- 40           # camera viewing angle (degrees)
sigma_closest <- 4            # scale of the trigger function (m)
sigma_fov     <- 12           # scale for detection of further individuals (m)
clus_radius   <- 10           # SD for cluster radius (m), as in the main simulations
clus_size     <- 30           # mean cluster size, as used for D = 0.05
width         <- 500
height        <- 500

subpops <- c("Uniform", "Clustered")     # Region.Label values
area   <- width * height                 # area of each subpopulation (m^2)
sector <- (fov / 360) * pi * w^2         # camera sector area (m^2)
N_true <- D_true * area                  # true abundance of each subpopulation

##---------------------------------------
## Helpers local to this script
##---------------------------------------
# run a Distance fit quietly, returning NULL on failure
quiet_ds <- function(expr) {
  tryCatch(suppressWarnings(suppressMessages({
    utils::capture.output(m <- expr); m
  })), error = function(e) NULL)
}

# abundance by sub-population from a fitted ds model
abund_ds <- function(fit, flatfile) {
  e <- as.data.frame(dht2(fit, flatfile = flatfile, strat_formula = ~Region.Label,
                          er_est = "P2", sample_fraction = fov / 360))
  e <- e[e$Region.Label %in% subpops, ]
  data.frame(subpop = as.character(e$Region.Label),
             N = e$Abundance, lcl = e$LCI, ucl = e$UCI)
}

##---------------------------------------
## One replicate: survey both sub-populations, analyse jointly
##---------------------------------------
generate_survey <- function() {

  flat <- laps <- clos <- vector("list", length(subpops))

  for (i in seq_along(subpops)) {
    s <- subpops[i]

    animals <- if (s == "Uniform")
      generate_animals(width = width, height = height, density = D_true,
                       distribution = "uniform")
    else
      generate_animals(width = width, height = height, density = D_true,
                       distribution = "clustered",
                       cluster_radius = clus_radius,
                       mean_cluster_size = clus_size)

    # cameras at random locations, kept >= w from every edge, all facing south
    # note that some camera trap wedges can overlap. this does not invalidate or bias
    # results but non-independence does arise
    cx <- runif(n_cam, w, width  - w)
    cy <- runif(n_cam, w, height - w)

    rows <- lrows <- vector("list", n_cam)
    counts <- integer(n_cam)
    r_min  <- rep(NA_real_, n_cam)

    for (k in seq_len(n_cam)) {
      lab   <- paste0(s, "_C", k)
      blank <- data.frame(Sample.Label = lab, distance = NA, gs = NA,
                          size = NA, Region.Label = s)

      dk <- sample_sector(animals, origin = c(cx[k], cy[k]), bearing = 270,
                          radius = w, angle = fov)

      ## --- time lapse: camera records regardless of triggering, and every
      ##     individual present is detected conditional on sigma_fov
      if (nrow(dk)) {
        lk <- dk[rbinom(nrow(dk), 1, hn_func(dk$r, sigma_fov)) == 1, , drop = FALSE]
        lrows[[k]] <- if (nrow(lk))
          data.frame(Sample.Label = lab, distance = lk$r, gs = nrow(lk),
                     size = 1, Region.Label = s) else blank
      } else {
        lrows[[k]] <- blank
      }

      ## --- sensor-triggered protocols
      if (!nrow(dk)) { rows[[k]] <- blank; next }
      dk <- dk[order(dk$r), ]
      nk <- nrow(dk)
      r_closest <- dk$r[1]

      # camera triggered by the closest individual
      if (rbinom(1, 1, hn_func(r_closest, sigma_closest)) != 1L) {
        rows[[k]] <- blank; next
      }
      # further individuals may be obscured; the closest is already detected p = 1.0
      if (nk > 1L) {
        p  <- c(1.0, hn_func(dk$r[2:nk], sigma_fov))
        dk <- dk[rbinom(nk, 1, p) == 1, , drop = FALSE]
      }

      rows[[k]]   <- data.frame(Sample.Label = lab, distance = dk$r,
                                gs = nrow(dk), size = 1, Region.Label = s)
      counts[k]   <- nrow(dk)       # group size = number of individuals detected
      r_min[k]    <- r_closest      # distance to the closest individual
    }

    flat[[i]] <- list_rbind(rows)
    laps[[i]] <- list_rbind(lrows)
    clos[[i]] <- data.frame(subpop = s, count = counts, rmin = r_min)
  }

  list(flat = list_rbind(flat) |> mutate(object = row_number(), Effort = 1,
                                        Area = area),
       laps = list_rbind(laps) |> mutate(object = row_number(), Effort = 1,
                                        Area = area),
       clos = list_rbind(clos))
}

##---------------------------------------
## One replicate: generate, then analyse four ways
##---------------------------------------
sim_once <- function() {

  g    <- generate_survey()
  flat <- g$flat; laps <- g$laps; clos <- g$clos

  fd <- filter(flat, !is.na(distance))     # triggered images, all distances
  ld <- filter(laps, !is.na(distance))     # time-lapse images, all distances
  if (nrow(fd) < 20 || nrow(ld) < 20) return(NULL)

  out <- list()

  ## ---- CTDS ------------------------------
  fit <- quiet_ds(select_best_ds(fd, w = w))
  if (!is.null(fit)) {
    out$ctds <- abund_ds(fit, flat) |>
      mutate(model = "CTDS", key = fit$ddf$name.message)
  }

  ## ---- MCDS -----------------------------

  cand <- list(); cand_nm <- character(0)
  for (ky in c("hn", "hr")) {
    m <- quiet_ds(ds(fd, transect = "point", key = ky, formula = ~Region.Label,
                     truncation = w, adjustment = NULL))
    if (!is.null(m)) {
      cand[[length(cand) + 1L]] <- m
      cand_nm <- c(cand_nm, m$ddf$name.message)
    }
  }
  if (length(cand)) {
    j   <- which.min(vapply(cand, function(m) AIC(m)$AIC, numeric(1)))
    fit <- cand[[j]]
    out$mcds <- abund_ds(fit, flat) |>
      mutate(model = "MCDS", key = cand_nm[j])
  }

  ## ---- Closest --------------------------
  try({
    ok  <- !is.na(clos$rmin)
    fit <- fit_detection_hn(clos$rmin[ok], w = w, gs = clos$count[ok])
    out$closest <- list_rbind(lapply(subpops, function(s) {
      e <- estimate_density_closest(clos$count[clos$subpop == s],
                                    w = w, angle = fov, fit = fit)
      data.frame(subpop = s, N = e$D * area, lcl = e$lcl * area,
                 ucl = e$ucl * area)
    })) |> mutate(model = "Closest", key = "half-normal", sigma = fit$sigma)
  }, silent = TRUE)

  ## ---- Time lapse -----------------------

  try({
    fit <- fit_detection_hn(ld$distance, w = w, gs = 1)
    out$lapse <- list_rbind(lapply(subpops, function(s) {
      cams_s <- unique(laps$Sample.Label[laps$Region.Label == s])
      ck <- as.numeric(table(factor(ld$Sample.Label[ld$Region.Label == s],
                                    levels = cams_s)))
      e  <- estimate_density_ctds(ck, w = w, angle = fov, fit = fit)
      data.frame(subpop = s, N = e$D * area,
                 lcl = e$lcl * area, ucl = e$ucl * area)
    })) |> mutate(model = "Time lapse", key = "half-normal", sigma = fit$sigma)
  }, silent = TRUE)

  if (!length(out)) return(NULL)
  res <- bind_rows(out)

  ## observed group-size composition, for context
  mg <- tapply(clos$count[clos$count > 0], clos$subpop[clos$count > 0], mean)
  nd <- table(factor(fd$Region.Label, levels = subpops))

  res |>
    mutate(N_true  = N_true,
           mean_gs = as.numeric(mg[subpop]),
           n_det   = as.numeric(nd[subpop]),
           cover   = lcl <= N_true & N_true <= ucl)
}

##---------------------------------------
## Replicate sims
##---------------------------------------

set.seed(2026)

res <- vector("list", n_rep)
pb  <- txtProgressBar(max = n_rep, style = 3)
for (i in seq_len(n_rep)) {
  r <- sim_once()
  if (!is.null(r)) { r$rep <- i; res[[i]] <- r }
  setTxtProgressBar(pb, i)
}
close(pb)

lev <- c("CTDS", "MCDS", "Closest", "Time lapse")
res <- list_rbind(res) |>
  mutate(model    = factor(model, levels = lev),
         subpop   = factor(subpop, levels = subpops),
         rel_bias = (N - N_true) / N_true)

dir.create("outputs", showWarnings = FALSE)
write_rds(res, "outputs/res_subpopulation.rds")

# you can write in res here and investigate with summaries and plots below
res <- readRDS("outputs/res_subpopulation.rds")

##---------------------------------------
## Summaries
##---------------------------------------
cat("\n true abundance of each subpopulation:", N_true, "\n")

res |>
  summarise(reps        = n(),
            N_mean      = round(mean(N)),
            mean_bias   = round(mean(rel_bias), 3),
            cv          = round(sd(N) / mean(N), 3),
            coverage    = round(mean(cover, na.rm = TRUE), 2),
            .by = c(model, subpop)) |>
  arrange(model, subpop) |> as.data.frame() |> print()


cat("\n Fitted scale parameters \n")
res |> filter(!is.na(sigma)) |>
  summarise(sigma = round(median(sigma), 2), .by = model) |>
  as.data.frame() |> print()

cat("\n key function selected \n")
res |> filter(model %in% c("CTDS", "MCDS")) |>
  count(model, key) |> arrange(model, desc(n)) |> as.data.frame() |> print()

##---------------------------------------
## Figures
##---------------------------------------
theme_ms <- theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size = 13),
        axis.text  = element_text(size = 10),
        strip.text = element_text(face = "bold", size = 11))

## relative bias

bias_dat <- bind_rows(
  res |> select(model, rep, subpop, rel_bias),
  res |> summarise(N = sum(N), k = n(), .by = c(model, rep)) |>
    filter(k == length(subpops)) |>
    mutate(subpop = "Total (pooled)",
           rel_bias = (N - 2 * N_true) / (2 * N_true)) |>
    select(model, rep, subpop, rel_bias)) |>
  mutate(subpop = factor(subpop, levels = c(subpops, "Total (pooled)"),
                         labels = c("A (uniform)","B (clustered)","Total (pooled)")))

bias_dat |> summarise(bias = mean(rel_bias), .by = c(model,subpop))


p_bias <- ggplot(bias_dat, aes(subpop, rel_bias, fill = subpop)) +
  geom_boxplot(outlier.size = 0.5, outliers = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dotted") +
  facet_wrap(~model, nrow = 1) +
  coord_cartesian(ylim = c(-0.6, 0.6)) +
  scale_fill_manual(values = c("grey75", "grey45", "steelblue")) +
  labs(x = NULL, y = expression(bold("Relative bias of " * hat(N)))) +
  theme_ms + theme(axis.text.x = element_text(angle = 25, hjust = 1))

p_bias

