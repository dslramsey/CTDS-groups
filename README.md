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
  when estimating a distance sampling detection function.\
- For camera traps that rely on heat-in-motion sensors to trigger the
  camera, the closest individual is most likely to trigger the camera
  sensor. This means that other individuals in the field of view could
  contribute observations that do not depend on the sensitivity of the
  camera sensor. This likely represents a source of detection
  heterogeneity.
- Theory suggests that this heterogeneity can be accounted for by
  fitting a single detection function to the distances of all observable
  individuals, and unbiased estimates of density should then be obtained
  based on the property of pooling robustness.
- An alternative approach in this situation is to explicitly model the
  detection of the closest individual, which should directly represent
  the detection probability of the camera sensor. To address this, we
  developed a new availability model based on the distance distribution
  of the nearest individual using order statistics; to better reflect
  how multiple individuals are detected by camera traps.
- Simulation results show that the standard CTDS approach has a small
  negative bias under this source of detection heterogeneity, whereas
  the adjusted availability gives approximately unbiased estimates with
  near optimal confidence interval coverage when only the closest
  detection in a group is recorded. Simulation of standard CTDS when
  cameras are not triggered by a sensor (i.e. time-lapse mode) were
  unbiased under all scenarios we considered.
- Our study demonstrates that the proposed method has some advantages
  over the standard CTDS analyses for camera traps triggered by a
  heat-in-motion sensor as only a single distance need be recorded in
  images containing multiple individuals.

### File descriptions:

- `r/density_simulation_stydy.r` simulation code to generate random
  animal locations within a rectangular area that are then sampled with
  random camera locations (circle sector with given angular field of
  view). Includes options for both random uniform and clustered
  distributions of animals. Camera snapshot moments are generated
  assuming the closest individual triggers the camera. Density
  estimation compares the standard CTDS approach with our method based
  on recording only the distance to the closest individual in the group.
  We also compare results with an alternative CTDS method where cameras
  are not triggered by a sensor.
- `r/CTDS_density_functions.r` contains various functions required by
  the main script.

## Prerequisites

The script require packages `tidyverse`, `Distance`, `numDeriv`.
