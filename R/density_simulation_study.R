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

source("R/group_size_funcs.R")
source("R/CTDS_density_functions.R")


##---------------------------------------
## Uniform density field
##---------------------------------------
n_rep<- 1000   # increase for better inference
D<- 0.01
ncams<- 200
fov<- 40
w<- 15
sigma_true<- 5
bearing <- 270
width<- 500
height<- 500
clus_rad <- 5
clus_size<- 10
cam_layout<- "grid"

delta<- 2.5
bin_bp <- seq(0, w, by = delta)  # bin breakpoints

  res_unif <- run_density_sim(n_rep = n_rep,
                              D_true = D,
                              sigma_true = sigma_true,
                              w = w,
                              fov = fov,
                              n_cam = ncams,
                              width = width,
                              height = height,
                              distribution = "uniform",
                              camera_layout =  cam_layout,
                              cluster_radius = clus_rad,
                              mean_cluster_size =  clus_size,
                              binned = FALSE,
                              breaks = NULL)

  res_unif_cl <- run_density_closest(n_rep = n_rep,
                                  D_true = D,
                                  sigma_true = sigma_true,
                                  w = w,
                                  fov = fov,
                                  n_cam = ncams,
                                  width = width,
                                  height = height,
                                  distribution = "uniform",
                                  camera_layout =  cam_layout,
                                  cluster_radius = clus_rad,
                                  mean_cluster_size =  clus_size,
                                  binned = FALSE,
                                  breaks = NULL)



  print(bind_rows(
    summarise_density_sim(res_unif),
    summarise_density_sim(res_unif_cl)
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
                              sigma_true = sigma_true,
                              w = w,
                              fov = fov,
                              n_cam = ncams,
                              width = width,
                              height = height,
                              distribution = "clustered",
                              camera_layout = cam_layout,
                              cluster_radius = clus_rad,
                              mean_cluster_size =  clus_size,
                              binned = FALSE,
                              breaks = NULL)

  res_clus_cl <- run_density_closest(n_rep = n_rep,
                              D_true = D,
                              sigma_true = sigma_true,
                              w = w,
                              fov = fov,
                              n_cam = ncams,
                              width = width,
                              height = height,
                              distribution = "clustered",
                              camera_layout = cam_layout,
                              cluster_radius = clus_rad,
                              mean_cluster_size =  clus_size,
                              binned = FALSE,
                              breaks = NULL)


  print(bind_rows(
    summarise_density_sim(res_clus),
    summarise_density_sim(res_clus_cl)
  ))

  win.graph(10,10)
  plot_density_sim(res_unif)

  win.graph(10,10)
  plot_density_sim(res_clus)

##-----------------------------------------------------------------
##  plot a single density realisation
##-----------------------------------------------------------------
  animals<- generate_animals(width = width,
                             height = height,
                             density = D,
                             distribution = "uniform",
                             cluster_radius = 5,
                             mean_cluster_size = 10)

  cams<- generate_cam_locs(ncams, w, width, height, bearing, cam_layout)

  win.graph(10,10)
  plot_animals_sectors(animals, cams, radius = w, angle = fov,
                       width = width, height = height,
                       show_detected = TRUE)


  dd<- seq(0,15,0.1)
  pp<- hn_func(dd, 13)
  win.graph(7,7)
  plot(dd,pp, type="l", ylim=c(0,1))

