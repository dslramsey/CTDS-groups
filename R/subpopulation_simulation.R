# ---------------------------------------------------------------------------
# Subpopulation bias under a pooled detection function
# ---------------------------------------------------------------------------
# Two subpopulations with identical true density that differ only in how animals
# are distributed in space: one uniform (homogeneous Poisson), one clustered.
# Clustering raises the number of individuals recorded per detected snapshot moment
# without changing abundance, so the two subpopulations differ in group-size
# composition while being identical in the quantity being estimated.
#
# Both subpopulations are surveyed together and analysed with a single pooled
# detection function, as would be common in real-world studies,
# as small per-subpopulation sample sizes usually require.
# Any difference between the two estimates is therefore bias: the truth is the
# same abundance in both.
#
# Observation process is as in the main simulations (see results 3.1): the closest individual
# triggers the camera following a half-normal with scale sigma_closest (here set to 4m); any
# further individuals in the FOV are then detected following a half-normal with
# scale sigma_fov (here set to 12m), so group size is subject to detection error.
#
# Models (all share one pooled detection function, subpopulation = stratum):
#   CTDS         distances to all detected individuals; key selected by AIC
#                from half-normal, hazard-rate and uniform keys with cosine
#                adjustments (select_best_ds)
#   MCDS         as CTDS but with subpopulation as a covariate on the scale
#                parameter; key selected by AIC from half-normal and
#                hazard-rate. No adjustment terms: ds() cannot enforce
#                monotonicity alongside covariates ("Monotonicity cannot be
#                enforced with covariates"), so adjusted covariate models are
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
n_rep         <- 2000          # replicates; ~2 h at these settings, reduce to test
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

# abundance by subpopulation from a fitted ds model
abund_ds <- function(fit, flatfile) {
  e <- as.data.frame(dht2(fit, flatfile = flatfile, strat_formula = ~Region.Label,
                          er_est = "P2", sample_fraction = fov / 360))
  e <- e[e$Region.Label %in% subpops, ]
  data.frame(subpop = as.character(e$Region.Label),
             N = e$Abundance, lcl = e$LCI, ucl = e$UCI)
}

##---------------------------------------
## One replicate: survey both subpopulations, analyse jointly
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

  ## ---- CTDS: one pooled key selected by AIC ------------------------------
  fit <- quiet_ds(select_best_ds(fd, w = w))
  if (!is.null(fit)) {
    out$ctds <- abund_ds(fit, flat) |>
      mutate(model = "CTDS", key = fit$ddf$name.message)
  }

  ## ---- MCDS: subpopulation as a covariate on the scale -------------------
  ## Key function selected by AIC, as for model CTDS, but without adjustment
  ## terms. ds() turns the monotonicity constraint off whenever the detection
  ## function has covariates (monotonicity = ifelse(formula == ~1, "strict",
  ## "none")) and refuses to turn it back on, so a covariate model carrying
  ## cosine adjustments can be non-monotonic or exceed one and cannot be
  ## constrained. The CTDS candidate set has no covariates and so keeps its
  ## adjustment terms under a strict monotonicity constraint.
  ## A uniform key has no scale parameter for a covariate to act on and so
  ## cannot be included either.
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

  ## ---- Closest: pooled sigma, availability adjusted by group size --------
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

  ## ---- Time lapse: pooled sigma, standard availability -------------------
  ## estimate_density_ctds() is the manuscript's own design-based estimator and
  ## adds the detection-function variance (delta method) to the empirical
  ## encounter-rate variance, so the interval is built the same way as for the
  ## other models. Counts are per camera and include cameras that saw nothing.
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
## Replicate 2,000 times
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
# res <- readRDS("outputs/res_subpopulation.rds")

##---------------------------------------
## Summaries
##---------------------------------------
cat("\n===== true abundance of each subpopulation:", N_true, "=====\n")

cat("\n===== group-size composition and sample size =====\n")
res |>
  summarise(mean_group_size = round(mean(mean_gs), 2),
            detections      = round(mean(n_det)), .by = subpop) |>
  as.data.frame() |> print()

cat("\n===== estimated abundance by subpopulation =====\n")
cat("mean_bias is the primary summary; median_bias is shown for comparison\n")
res |>
  summarise(reps        = n(),
            N_mean      = round(mean(N)),
            N_median    = round(median(N)),
            p10         = round(quantile(N, 0.1)),   # spread of N-hat across
            p90         = round(quantile(N, 0.9)),   # reps, NOT a confidence interval
            mean_bias   = round(mean(rel_bias), 3),
            median_bias = round(median(rel_bias), 3),
            cv          = round(sd(N) / mean(N), 3),
            skew        = round(mean((N - mean(N))^3) / sd(N)^3, 2),
            coverage    = round(mean(cover, na.rm = TRUE), 2),
            .by = c(model, subpop)) |>
  arrange(model, subpop) |> as.data.frame() |> print()

cat("\n===== ratio of the two subpopulations (truth = 1) =====\n")
res |>
  select(model, rep, subpop, N) |>
  pivot_wider(names_from = subpop, values_from = N) |>
  mutate(ratio = Clustered / Uniform) |>
  summarise(reps         = sum(!is.na(ratio)),
            ratio_mean   = round(mean(ratio, na.rm = TRUE), 3),
            ratio_median = round(median(ratio, na.rm = TRUE), 3),
            lcl          = round(quantile(ratio, 0.1, na.rm = TRUE), 3),
            ucl          = round(quantile(ratio, 0.9, na.rm = TRUE), 3),
            .by = model) |>
  arrange(model) |> as.data.frame() |> print()

cat("\n===== fitted scale parameter (true sensor sigma =", sigma_closest, ") =====\n")
res |> filter(!is.na(sigma)) |>
  summarise(sigma = round(median(sigma), 2), .by = model) |>
  as.data.frame() |> print()

cat("\n===== key function selected =====\n")
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

## estimated abundance against the truth
p_abund <- ggplot(res, aes(subpop, N, fill = model)) +
  geom_boxplot(outlier.size = 0.5, outliers = FALSE) +
  geom_hline(yintercept = N_true, linetype = "dashed", linewidth = 0.7) +
  facet_wrap(~model, nrow = 1) +
  labs(x = NULL, y = expression(bold("Estimated abundance " * hat(N)))) +
  theme_ms
# ggsave("outputs/subpop_abundance.png", p_abund, width = 10, height = 5, dpi = 300)

## relative bias, in the style of existing Figure 3
## relative bias for each subpopulation and for the two pooled. Pooling
## robustness protects the total, not the parts, so the third box is the
## one that should sit on zero even when the first two do not.
bias_dat <- bind_rows(
  res |> select(model, rep, subpop, rel_bias),
  res |> summarise(N = sum(N), k = n(), .by = c(model, rep)) |>
    filter(k == length(subpops)) |>
    mutate(subpop = "Total (pooled)",
           rel_bias = (N - 2 * N_true) / (2 * N_true)) |>
    select(model, rep, subpop, rel_bias)) |>
  mutate(subpop = factor(subpop, levels = c(subpops, "Total (pooled)")))

p_bias <- ggplot(bias_dat, aes(subpop, rel_bias, fill = subpop)) +
  geom_boxplot(outlier.size = 0.5, outliers = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dotted") +
  facet_wrap(~model, nrow = 1) +
  coord_cartesian(ylim = c(-0.6, 0.6)) +
  scale_fill_manual(values = c("grey75", "grey45", "steelblue")) +
  labs(x = NULL, y = expression(bold("Relative bias of " * hat(N)))) +
  theme_ms + theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave("outputs/subpop_rel_bias.png", p_bias, width = 10, height = 5.2, dpi = 300)

cat("\n===== relative bias: subpopulations and the pooled total =====\n")
print(as.data.frame(bias_dat |>
  summarise(mean_bias = round(mean(rel_bias), 3),
            median_bias = round(median(rel_bias), 3), .by = c(model, subpop)) |>
  pivot_wider(names_from = subpop, values_from = c(mean_bias, median_bias)) |>
  arrange(model)))

#### OPTIONAL Plot for supp material, example populations ####
## example realisations of the two subpopulations, with camera sectors
set.seed(7)
sector_poly <- function(id, ox, oy, bearing, radius, angle, n = 40) {
  half <- (angle * pi / 180) / 2
  b    <- bearing * pi / 180
  a    <- seq(b - half, b + half, length.out = n)
  data.frame(id = id, x = c(ox, ox + radius * cos(a)),
             y = c(oy, oy + radius * sin(a)))
}
cam_demo <- data.frame(id = 1:25,
                       x = runif(25, w, width - w), y = runif(25, w, height - w))
sectors  <- pmap(cam_demo, function(id, x, y)
  sector_poly(id, x, y, 270, w, fov)) |> list_rbind()

mgs   <- res |> summarise(m = mean(mean_gs), .by = subpop) |> deframe()
field <- bind_rows(
  generate_animals(width, height, D_true, "uniform") |> select(x, y) |>
    mutate(panel = sprintf("(A) Uniform\nmean group size in view = %.1f", mgs["Uniform"])),
  generate_animals(width, height, D_true, "clustered",
                   cluster_radius = clus_radius, mean_cluster_size = clus_size) |>
    select(x, y) |>
    mutate(panel = sprintf("(B) Clustered\nmean group size in view = %.1f", mgs["Clustered"])))

p_field <- ggplot(field, aes(x, y)) +
  geom_point(alpha = 0.25, size = 0.3) +
  geom_polygon(data = expand_grid(panel = unique(field$panel), sectors),
               aes(group = id), fill = "firebrick", alpha = 0.35,
               colour = "firebrick", linewidth = 0.25) +
  facet_wrap(~panel) +
  coord_equal(xlim = c(0, width), ylim = c(0, height)) +
  labs(x = "x (m)", y = "y (m)") +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 12),
        strip.text = element_text(face = "bold", size = 11))
ggsave("outputs/subpop_field.png", p_field, width = 9, height = 5, dpi = 300)

print(p_field); print(p_abund); print(p_bias)
