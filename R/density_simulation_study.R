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
library(docstring)

source("R/CTDS_density_functions.R")

# For function help use docstring e.g.
docstring(run_density_ctds)
##---------------------------------------
## Uniform density field
##---------------------------------------
n_rep<- 10                  # increase for better inference
D<- 0.01                    # Density (m2)
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
clus_size<- 30              # average size of cluster sd. set to 30 for D = 0.005 and 0.01,
                            # or 60 for D=0.05 or 0.1

delta<- 2.5
bin_bp <- seq(0, w, by = delta)  # bin breakpoints for binned = TRUE

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

  print(bind_rows(
    summarise_density_sim(res_unif),
    summarise_density_sim(res_unif_cl),
    summarise_density_sim(res_unif_lps)
  ))

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

  print(bind_rows(
    summarise_density_sim(res_clus),
    summarise_density_sim(res_clus_cl),
    summarise_density_sim(res_clus_lps)
  ))


# Copy truth attributes (model, distribution, D_true) onto the results as columns
add_truth_cols <- function(res) {
  truth <- attr(res, "truth")
  res$model <- truth$model
  res$distribution <- truth$distribution
  res$D_true <- truth$D_true
  res
}

res_all<- bind_rows(lapply(list(res_unif,res_unif_cl,res_unif_lps,
                                res_clus,res_clus_cl,res_clus_lps),
                           add_truth_cols))

##---- Box plot of relative bias (minimum example) ----

win.graph(10,8)
res_all |> mutate(model = factor(model,
                        levels = c("CTDS", "Closest", "lapse")),
         distribution = factor(distribution,
                               levels = c("uniform", "clustered"),
                               labels = c("Uniform", "Clustered"))) |>
ggplot(aes(x = factor(D_true), y = (D - D_true) / D_true, fill = model)) +
  geom_boxplot(outlier.size = 0.5, outliers = FALSE,
               position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = c(-0.1,0.1), linetype = "dotted") +
  facet_grid(distribution ~ model) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(x = "True density (D)", y = expression(bold("Relative bias of " * hat(D))),
       fill = "Detection") +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.title = element_text(face = "bold", size=14),
        axis.text = element_text(size=11),
        strip.text = element_text(face = "bold", size =12))


