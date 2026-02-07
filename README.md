# Satellite Optimal Guidance

Individual project implementing orbit design for Libration point missions and optimal guidance strategies for deep-space scenarios.

## Overview

This repository contains three simulation-based studies developed as part of a university satellite guidance project.

The work is organized in three parts:
1. Characterization and continuation of Earth-Moon L1 Halo orbits family
2. Optimization of impulsive guidance Aphophis asteroid impactor mission
3. Time-of-flight optimization of continuous guidance Earth-Venus transfer

## Results and Validation

Key results:
- Periodic Halo orbit computed through first-order STM iterative corrections. Orbit family identified through natural parameter continuation on initial conditions.
- Optimized impulsive impactor launch and deep-space spacecraft maneuvers, increasing Apophis close-approach distance by 2x factor.
- Computed time-optimal Earth-Moon transfer of ~ 141 days given 800mN thrust, with final heliocentric errors of ~ 0.04 km, ~ 0.01 km/s. Optimality affirmed through Hamiltonian testing. 500mN optimal transfer obtained through numerical continuation

Representative outputs:
- Periodic Halo orbits plots
- Impactor spacecraft optimal maneuver sequence
- Time-optimal Earth-Venus transfers

Representative figures are available in 'results/' (PNG format).
See 'results/results.txt' for figure-by-figure explanations.
The full methodology and results are documented in 'docs/report.pdf'.

## Repository structure

- 'src/' - MATLAB implementations of each study
- 'docs/' - Project report
- 'results/' - Key result figures (PNG)
- 'figures/' - Source figures (EPS)

## Development notes

This repository was uploaded after project completion.
Commit history does not reflect the original development timeline.

## Reproducibility and external dependencies

This project relies on external SPICE kernels for ephemerides and constants, not included in this repository.
To run the simulations, the following NAIF kernels are required:
- 'naif0012.tls'
- 'pck00010.tpc'
- 'gm_de432.tpc'
- 'de432s.bsp'
- '20099942_Apophis.bsp'