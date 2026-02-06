# Satellite Optimal Guidance

Academic project implementing orbit design and optimal guidance strategies.

## Overview

This repository contains three simulation-based studies developed as part of a university satellite guidance project.

The work is organized in three parts:
1. Characterization and continuation of Earth-Moon L1 Halo orbits family
2. Optimization of impulsive guidance Aphophis asteroid impactor mission
3. Time-of-flight optimization of continuous guidance Earth-Venus transfer

The full methodology and results (including plots) are available in 'docs/report.pdf'.

## Repository structure

- 'src/' - MATLAB implementations of each study
- 'docs/' - Project report

## Reproducibility and external dependencies

This project relies on external SPICE kernels for ephemerides and constants, not included in this repository.
To run the simulations, the following NAIF kernels are required:
- 'naif0012.tls'
- 'pck00010.tpc'
- 'gm_de432.tpc'
- 'de432s.bsp'
- '20099942_Apophis.bsp'