# Complementarity-Based Friction Modeling for Tendon-Driven Continuum Robots

This repository provides MATLAB code for friction modeling of tendon-driven continuum robots (TDCR) using a complementarity-based formulation. The proposed method incorporates the Capstan friction model as a set of complementarity constraints, enabling the prediction of tendon–disk friction loss, transitions between sticking and sliding, and friction-induced hysteresis. The repository includes demo scripts for both a 2-disk TDCR example and more general TDCR simulations.

## File Locations

- `TDCR_2disk/fig_vid_pull.m`: 2-disk TDCR example corresponding to the experimental results shown in Fig. 7 of the paper.
- `TDCR_general/simu_multiTendon.m`: demo for general TDCR simulations.

## Setup

Before running the code, add all folders in this repository to the MATLAB path:

```matlab
addpath(genpath(pwd));
```

## Citation

If you use this code in your research, please cite:

J. Shen, B. Browne, J. Ha and Y. Chen, "Friction Modeling of Tendon-Driven Continuum Robots Through Linear Complementarity Problem," in *IEEE Robotics and Automation Letters*, vol. 11, no. 6, pp. 7468-7475, June 2026, doi: 10.1109/LRA.2026.3688094.

### BibTeX

```bibtex
@ARTICLE{11495397,
  author={Shen, Jia and Browne, Brendan and Ha, Junhyoung and Chen, Yue},
  journal={IEEE Robotics and Automation Letters},
  title={Friction Modeling of Tendon-Driven Continuum Robots Through Linear Complementarity Problem},
  year={2026},
  volume={11},
  number={6},
  pages={7468-7475},
  doi={10.1109/LRA.2026.3688094}
}
```
