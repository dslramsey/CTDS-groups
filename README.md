A new method to account for the detection of multiple individuals in
camera traps for use in camera trap distance sampling (CTDS)
applications
================

## Overview

This repository contains code for peer review only:

- Camera trap distance sampling (CTDS) is a popular recent method used
  to estimate wildlife abundance from camera trap images that uses the
  distances of detected individuals from the camera and point distance
  sampling methods to estimate animal density. When multiple individuals
  are detected in the camera field of view at the same time, standard
  practice involves recording the distance to all observed individuals
  when estimating a distance sampling detection function.
- For camera traps that rely on heat-in-motion sensors to trigger the
  camera, the closest individual is most likely to trigger the camera
  sensor. This means that other individuals in the field of view could
  contribute observations that do not depend on the sensitivity of the
  camera sensor. This likely represents a source of detection
  heterogeneity.
- Theory suggests that this heterogeneity can be accounted for by
  fitting a single detection function to the distances of all observable
  individuals, and unbiased estimates of density should then be obtained
  based on the property of pooling robustness. This is the approach
  currently advocated for analysing CTDS data.
- An alternative approach in this situation is to explicitly model the
  detection of the closest individual, which should directly represent
  the detection probability of the camera sensor. To address this, we
  developed a new availability model based on the distance distribution
  of the nearest individual using order statistics; to better reflect
  how multiple individuals are detected by camera traps.
- This repository contains `R` code and functions to conduct the
  simulations documented in the manuscript comparing the performance of
  the standard CTDS approach with the new method based on recording only
  the distance of the closest individual from a detected group. We also
  compare a third CTDS approach where cameras are not triggered by a
  motion sensor but are programmed to record snapshot moments at set
  intervals (i.e., time lapse mode). All simulations were conducted in
  `R (v. 4.5.3)`

### File descriptions:

- `r/density_simulation_study.r` simulation code to generate random
  animal locations within a rectangular area that are then sampled with
  random camera locations (circle sector with given angular field of
  view). Includes options for both random uniform and clustered
  distributions of animals. Camera snapshot moments are generated
  assuming the closest individual triggers the camera. Density
  estimation compares the standard CTDS approach with our method based
  on recording only the distance to the closest individual in the group.
  We also compare results with an alternative CTDS method where cameras
  are not triggered by a sensor.
- `r/subpopulation_simulation.r` simulation code (largely contained in `sim_once()`) 
  to generate random animal location in two "subpopulations" of the same density (e.g. 0.01) 
  as per `r/density_simulation_study.r`, but with different levels of clustering, essentially,
  impacting only the average group size. However, instead of estimating abundance
  of each subpopulation with the use of a separate detection function, a single
  detection function is used for the following methods: CTDS, MCDS (subpopulation 
  covariate), 'closest' and time-lapse. Plots and tables are generated at the end of
  the script using the output of the replicate `sim_once()` saved as `outputs/res_subpopulation.rds`.
- `r/CTDS_density_functions.r` contains various functions required by
  the main script. Help for each function can viewed using the
  `docstring` package.

## Prerequisites

The script require packages `tidyverse (v. 2.0.0)`,
`Distance (v.2.0.1)`, `numDeriv (v. 2016.8-1.1)`,
`docstring (v. 1.0.0)`.
