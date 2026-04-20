# Complementarity-Based Friction Modeling for Tendon-Driven Continuum Robots

This repository contains MATLAB code for implementing a complementarity-based friction model for tendon-driven continuum robots (TDCR).

## File Locations

1. For the 2-disk TDCR example, see `TDCR_2disk\fig_vid_pull.m`. This example complements the experimental results shown in Figure 7 of our paper.

2. For general TDCRs, see the demo `TDCR_general\simu_multiTendon.m`.

## Setup

Before running the code, please add all folders in this repository to the MATLAB path.

```matlab
addpath(genpath(pwd));