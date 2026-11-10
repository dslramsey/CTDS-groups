A new method to account for the detection of multiple individuals in
camera traps for use in camera trap distance sampling (CTDS)
applications
================

## Overview

This repository contains data and code from:

Ramsey, D.S.L., and Cally, J.G. (2027). A new method to account for the
detection of multiple individuals in camera traps for use in camera trap
distance sampling (CTDS) applications. *Methods in Ecology and
Evolution* (in press).

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
- `r/subpopulation_simulation.r` simulation code to generate random
  animal location in two “subpopulations” of the same density but with
  different levels of clustering, essentially, impacting only the
  average group size. A single detection function is then fitted to the
  pooled data using the standard CTDS approach and our closest
  individual approach. We also include a model that fits a covariate for
  the subpopulation identity using MCDS.
- `r/CTDS_density_functions.r` contains various functions required by
  the main scripts. Help for each function can viewed using the
  `docstring` package.

## Prerequisites

The script require packages `tidyverse (v. 2.0.0)`,
`Distance (v.2.0.1)`, `numDeriv (v. 2016.8-1.1)`,
`docstring (v. 1.0.0)`.
