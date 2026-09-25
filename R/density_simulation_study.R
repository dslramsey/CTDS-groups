# Simulation study: bias and coverage of standard CTDS distance-sampling
# compared with closest individual detetcion, using multiple independent
# camera placements per replicate.
# Density estimated using conditional likelihood and an empirical
# encounter-rate variance estimator.
#
# For each replicate we simulate one animal field with known global density,
# then place n_cam cameras at random or systematic locations.

library(tidyverse)
library(Distance)


source("R/group_size_funcs.R")
source("R/CTDS_density_functions.R")


##---------------------------------------
## Uniform density field
##---------------------------------------
n_rep<- 100                 # increase for better inference
D<- 0.1                     # Density (m2)
ncams<- 500                 # Camera sectors (and snapshot moments)
fov<- 40                    # actually camera angle
w<- 15                      # truncation distance (m)
sigma_closest<- 4           # True sigma for HN detection of the camera sensor
sigma_true<- 12             # True sigma for HN detection of animals in camera FOV
perfect_fov<- FALSE         # Whether FOV detection is turned on (FALSE) or not (TRUE)
bearing <- 270              # Camera orientation
width<- 500                 # Sampling frame area width and height (m)
height<- 500
clus_rad <- 10              # Average radius of clusters (m)
clus_size<- 60              # average size of clustersd set to 30 for D = 0.005 and 0.01,
                            # or 60 for D=0.05 or 0.1

delta<- 2.5
bin_bp <- seq(0, w, by = delta)  # bin breakpoints

  res_unif <- run_density_ctds(n_rep = n_rep,
                              D_true = D,
                              sigma_closest = sigma_closest,
                              sigma_true = sigma_true,
                              w = w,
                              fov = fov,
                              perfect_fov = perfect_fov,
                              n_cam = ncams,
                              width = width,
                              height = height,
                              distribution = "uniform",
                              cluster_radius = clus_rad,
                              mean_cluster_size =  clus_size,
                              binned = FALSE,
                              breaks = NULL)

  #write_rds(res_unif, "outputs/res_unif_ctds_high.rds")
  summarise_density_sim(res_unif)

  res_unif_cl <- run_density_closest(n_rep = n_rep,
                                  D_true = D,
                                  sigma_closest = sigma_closest,
                                  sigma_true = sigma_true,
                                  w = w,
                                  fov = fov,
                                  perfect_fov = perfect_fov,
                                  n_cam = ncams,
                                  width = width,
                                  height = height,
                                  distribution = "uniform",
                                  cluster_radius = clus_rad,
                                  mean_cluster_size =  clus_size,
                                  binned = FALSE,
                                  breaks = NULL)

  #write_rds(res_unif_cl, "outputs/res_unif_cl_vlow.rds")
  summarise_density_sim(res_unif_cl)

  res_unif_lps <- run_density_lapse(n_rep = n_rep,
                                     D_true = D,
                                     sigma_true = sigma_true,
                                     w = w,
                                     fov = fov,
                                     perfect_fov = perfect_fov,
                                     n_cam = ncams,
                                     width = width,
                                     height = height,
                                     distribution = "uniform",
                                     cluster_radius = clus_rad,
                                     mean_cluster_size =  clus_size,
                                     binned = FALSE,
                                     breaks = NULL)

  #write_rds(res_unif_lps, "outputs/res_unif_lps_vlow.rds")

  print(bind_rows(
    summarise_density_sim(res_unif),
    summarise_density_sim(res_unif_cl),
    summarise_density_sim(res_unif_lps)
  ))

#  win.graph(10,10)
#  plot_density_sim(res_unif)

#  win.graph(10,10)
#  plot_density_sim(res_unif_cl)

  ##---------------------------------------
  ## Clustered density field
  ##---------------------------------------
  res_clus <- run_density_ctds(n_rep = n_rep,
                              D_true = D,
                              sigma_closest = sigma_closest,
                              sigma_true = sigma_true,
                              w = w,
                              fov = fov,
                              perfect_fov = perfect_fov,
                              n_cam = ncams,
                              width = width,
                              height = height,
                              distribution = "clustered",
                              cluster_radius = clus_rad,
                              mean_cluster_size =  clus_size,
                              binned = FALSE,
                              breaks = NULL)

  summarise_density_sim(res_clus)

  #write_rds(res_clus, "outputs/res_clus_ctds_vlow.rds")

  res_clus_cl <- run_density_closest(n_rep = n_rep,
                                     D_true = D,
                                     sigma_closest = sigma_closest,
                                     sigma_true = sigma_true,
                                     w = w,
                                     fov = fov,
                                     perfect_fov = perfect_fov,
                                     n_cam = ncams,
                                     width = width,
                                     height = height,
                                     distribution = "clustered",
                                     cluster_radius = clus_rad,
                                     mean_cluster_size =  clus_size,
                                     binned = FALSE,
                                     breaks = NULL)

  #write_rds(res_clus_cl, "outputs/res_clus_cl_vlow.rds")

  res_clus_lps <- run_density_lapse(n_rep = n_rep,
                                     D_true = D,
                                     sigma_true = sigma_true,
                                     w = w,
                                     fov = fov,
                                     perfect_fov = perfect_fov,
                                     n_cam = ncams,
                                     width = width,
                                     height = height,
                                     distribution = "clustered",
                                     cluster_radius = clus_rad,
                                     mean_cluster_size =  clus_size,
                                     binned = FALSE,
                                     breaks = NULL)

  #write_rds(res_clus_lps, "outputs/res_clus_lps_vlow.rds")

  print(bind_rows(
    summarise_density_sim(res_clus),
    summarise_density_sim(res_clus_cl),
    summarise_density_sim(res_clus_lps)
  ))



  win.graph(10,10)
  plot_density_sim(res_unif)

  win.graph(10,10)
  plot_density_sim(res_clus)


## Summarise results--------------------------------------------

models<- c("res_unif_ctds","res_unif_cl","res_unif_lps",
           "res_clus_ctds","res_clus_cl","res_clus_lps")
rdens<- c("vlow","low","med","high")
stub4<- ".rds"

fnames<- crossing(models, rdens, stub4)

fnames_list <- fnames |>
  select(-stub4) |>
  tidyr::unite(fname, everything(), sep = "_") |>
  mutate(fname = file.path("Outputs", paste0(fname, ".rds"))) |>
  pull(fname)

## read all result files and row-bind, keeping parameter columns from fnames
res_all <- fnames |>
  select(-stub4) |>
  mutate(fname = fnames_list) |>
  filter(file.exists(fname)) |>
  mutate(data = purrr::map(fname, function(f) {
    d <- readRDS(f)
    tr <- attr(d, "truth")
    d <- tibble::as_tibble(d)
    for (nm in names(tr)) d[[nm]] <- tr[[nm]]
    d
  })) |>
  tidyr::unnest(data) |>
    select(-fname, -rdens, -w, -fov, -sigma_true, -n_cam)

##---- Identify outliers-----

res_all |>
  mutate(rel_bias = (D - D_true) / D_true) |>
  filter(abs(rel_bias) > 1) |>
  count(model, distribution, D_true, name = "n_outliers")

## ---- Tabulate
res_summ<- res_all |>
      summarise(
        n        = n(),
        mean_n   = mean(mean_n),
        median_D   = median(D),
        rel_bias = median((D - first(D_true)) / first(D_true)),
        cv_enc   = median(cv_encounter),
        cv_det   = median(cv_detection),
        cv       = median(cv),
        coverage = mean(cover),
        .by = c(model, distribution, D_true)
      ) |>
      mutate(across(c(median_D, rel_bias, cv_det, cv_enc, cv, coverage), ~round(.x, 3))) |>
      mutate(mean_n = round(mean_n, 1)) |>
      arrange(model, distribution, D_true)

write_csv(filter(res_summ, distribution=="uniform"), "outputs/sim_summary_unif.csv")
write_csv(filter(res_summ, distribution=="clustered"), "outputs/sim_summary_clus.csv")

res_summ |> group_by(model) |> summarise(cv = mean(cv), cover = mean(coverage))

##---- Box plot of relative bias ----
win.graph(10,8)
res_all |> mutate(model = factor(model,
                        levels = c("CTDS", "Closest", "lapse"),
                        labels = c("CTDS", "Closest", "Time lapse")),
         distribution = factor(distribution,
                               levels = c("uniform", "clustered"),
                               labels = c("Uniform", "Clustered"))) |>
ggplot(aes(x = factor(D_true), y = (D - D_true) / D_true, fill = model)) +
  geom_boxplot(outlier.size = 0.5, outliers = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(-0.1,0.1), linetype = "dotted") +
  facet_grid(distribution ~ model) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = "True density (D)", y = expression(bold("Relative bias of " * hat(D)))) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size=14),
        axis.text = element_text(size=11),
        strip.text = element_text(face = "bold", size =12))

ggsave("outputs/rel_bias.png")

##---- Box plot of mean precision (CV) ----
win.graph(10,10)
res_all |>
  mutate(model = factor(model,
                        levels = c("CTDS", "Closest", "lapse"),
                        labels = c("CTDS", "Closest", "Time lapse")),
         distribution = factor(distribution,
                               levels = c("uniform", "clustered"),
                               labels = c("Uniform", "Clustered"))) |>
  ggplot(aes(x = factor(D_true), y = cv, fill = model)) +
  geom_boxplot(outlier.size = 0.5, outliers = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_grid(distribution ~ model) +
  coord_cartesian(ylim = c(0, 0.5)) +
  labs(x = "True density (D)", y = "Relative precision (CV)") +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(face = "bold", size=12),
        strip.text = element_text(face = "bold", size =12))


