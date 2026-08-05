# Simulation study: bias and coverage of the sector distance-sampling
# density estimator, using MULTIPLE independent camera placements per replicate
# and an EMPIRICAL encounter-rate variance estimator.
#
# For each replicate we simulate one animal field with known global density,
# then place n_cam cameras at random locations and bearings (each sector kept
# fully inside the region). Detections are pooled to fit one half-normal
# detection function; density is estimated from total detections over total
# effective area; and the variance is estimated empirically from the
# among-camera count variation (plus a detection-function delta-method term).
#
# Depends on simulate_animals(), sample_sector(), fit_detection_hn(),
# estimate_density_multi(), and hn_func() in R/group_size_funcs.R.

library(tidyverse)
library(patchwork)

source("R/group_size_funcs.R")
source("R/CTDS_density_functions.R")


##---------------------------------------
## no group effects
##---------------------------------------
D<- 0.1
ncams<- 50
fov<- 40
w<- 15
sigma_true<- 5
bearing <- 270
width<- 200
height<- 200
clus_rad <- 5
clus_size<- 100
cam_layout<- "grid"

  res_unif <- run_density_sim(n_rep = 1000,
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
                              mean_cluster_size =  clus_size)

  res_clus <- run_density_sim(n_rep = 1000,
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
                              mean_cluster_size =  clus_size)


  print(bind_rows(
    summarise_density_sim(res_unif),
    summarise_density_sim(res_clus)
  ))

  win.graph(10,10)
  plot_density_sim(res_unif)

  win.graph(10,10)
  plot_density_sim(res_clus)

  ##---------------------------------------
  ## group effects
  ##---------------------------------------
  D<- 0.1
  ncams<- 50
  fov<- 40
  w<- 15
  sigma_true<- 5
  bearing <- 270
  width<- 200
  height<- 200
  clus_rad <- 5
  clus_size<- 100
  cam_layout<- "grid"

  res_unif <- run_density_sim_groups(n_rep = 1000,
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
                              mean_cluster_size =  clus_size)

  summarise_density_sim(res_unif)

  res_clus <- run_density_sim_groups(n_rep = 10,
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
                              mean_cluster_size =  clus_size)


  print(bind_rows(
    summarise_density_sim(res_unif),
    summarise_density_sim(res_clus)
  ))

  win.graph(10,10)
  plot_density_sim(res_unif)

  win.graph(10,10)
  plot_density_sim(res_clus)

##-----------------------------------------------------------------

  animals<- generate_animals(width = width,
                             height = height,
                             density = D,
                             distribution = "clustered",
                             cluster_radius = 5,
                             mean_cluster_size = 100)

  cams<- generate_cam_locs(ncams, w, width, height, bearing, cam_layout)

  win.graph(10,10)
  plot_animals_sectors(animals, cams, radius = w, angle = fov,
                       width = width, height = height,
                       show_detected = TRUE)
