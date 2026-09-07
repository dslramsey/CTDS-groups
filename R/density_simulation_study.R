# Simulation study: bias and coverage of standard CTDS distance-sampling
# compared with closest individual detetcion, using multiple independent
# camera placements per replicate.
# Density estimated using conditional likelihood and an empirical
# encounter-rate variance estimator.
#
# For each replicate we simulate one animal field with known global density,
# then place n_cam cameras at random or systematic locations.

library(tidyverse)
library(patchwork)
library(Distance)


source("R/group_size_funcs.R")
source("R/CTDS_density_functions.R")


##---------------------------------------
## Uniform density field
##---------------------------------------
n_rep<- 500   # increase for better inference
D<- 0.1
ncams<- 200
fov<- 40
w<- 15
sigma_closest<- 4
sigma_true<- 12
perfect_fov<- FALSE
bearing <- 270
width<- 500
height<- 500
clus_rad <- 5
clus_size<- 30

delta<- 2.5
bin_bp <- seq(0, w, by = delta)  # bin breakpoints

  res_unif <- run_density_sim(n_rep = n_rep,
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

  win.graph(10,10)
  plot_density_sim(res_unif)

  win.graph(10,10)
  plot_density_sim(res_unif_cl)

  ##---------------------------------------
  ## Clustered density field
  ##---------------------------------------
  res_clus <- run_density_sim(n_rep = n_rep,
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

  win.graph(10,10)
  plot_density_sim(res_unif)

  win.graph(10,10)
  plot_density_sim(res_clus)

##-----------------------------------------------------------------
##  plot a single density realisation
##-----------------------------------------------------------------

  D<- 0.001
  ncams<- 200
  fov<- 40
  w<- 15
  bearing <- 270
  width<- 500
  height<- 500
  clus_rad <- 50
  clus_size<- 5

  animals<- generate_animals(width = width,
                             height = height,
                             density = D,
                             distribution = "uniform",
                             cluster_radius = 5,
                             mean_cluster_size = 20)

  cams<- generate_cam_locs(ncams, w, width, height, bearing, camera_layout="random")

  win.graph(10,10)
  plot_animals_sectors(animals, cams, radius = w, angle = fov,
                       width = width, height = height,
                       show_detected = TRUE)



