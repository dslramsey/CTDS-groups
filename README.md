A new method to account for the detection of multiple individuals in
camera traps for use in camera trap distance sampling (CTDS)
applications
================

## Overview

This repository contains code from:

Ramsey, D.S.L., and Cally, J.G. (2026). A new method to account for the
detection of multiple individuals in camera traps for use in camera trap
distance sampling (CTDS) applications *Methods in Ecology and Evolution*

- Camera trap distance sampling (CTDS) is popular recent method used to
  estimate wildlife abundance from camera trap images that uses the
  distances of detected individuals from the camera and point distance
  sampling methods to estimate animal density. When multiple individuals
  are detected in the camera field of view at the same time, standard
  practice suggests that distances to all individuals in the group be
  used when estimating a distance sampling detection function. However,
  this may be problematic in CTDS studies as individuals in a group are
  likely to occur at a range of distances from the camera.
- For camera traps that rely on heat-in-motion sensors to trigger the
  camera, the closest individual is most likely to trigger the camera
  sensor. This means that other group members in the image may not
  represent independent detections, causing bias in the estimated
  detection function.
- To address this, we developed a new availability model based on the
  distance distribution of the nearest individual in a group using order
  statistics; to better reflect how grouped animals are detected by
  camera traps.
- Simulation results show that the standard CTDS approach can produce
  substantial positive bias, especially as group size increases, whereas
  the adjusted availability gives approximately unbiased estimates when
  only the closest detection in a group is recorded.
- Our study demonstrates that the proposed method improves CTDS analyses
  for social or group-living species, while noting that its validity
  depends on assumptions about which animal triggers the sensor and how
  individuals are distributed within the camera field of view

### File descriptions:

- `r/CTDS_groups_sims.r` simulation code to generate random distances of
  individuals from camera traps and fit distance sampling detection
  functions using a new availability model for the closest individual in
  a detected group in CTDS snapshot moments.
- `r/group_size_functions.r` contains various functions required by the
  main script.

## Prerequisites

The script require packages `tidyverse`, `tidyr`.
