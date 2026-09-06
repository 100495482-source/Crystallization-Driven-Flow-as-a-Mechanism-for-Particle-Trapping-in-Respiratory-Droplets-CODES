# Crystallisation-Driven Flow as a Mechanism for Particle Trapping in Respiratory Droplets

MATLAB code accompanying the Bachelor's Thesis in Biomedical Engineering by
Águeda Vargas-Zúñiga Areilza, Universidad Carlos III de Madrid, 2026.

Supervisors: Javier Martínez Puig and Théophile Gaichies.
Coordinator: Francisco Javier Rodríguez-Rodríguez.

## Overview

This repository contains the full analysis pipeline used in the thesis. The
work investigates whether the growth of NaCl crystals during the evaporation
of a sessile saline droplet generates an internal flow capable of transporting
suspended tracer particles towards the crystals, with implications for the
survival of viruses in drying respiratory droplets.

The code covers:

- image preprocessing and background subtraction
- Eulerian velocimetry (PIV), from PIVlab exports
- Lagrangian particle tracking, from Fiji/TrackMate exports
- crystal segmentation, area and front velocity measurement
- particle accumulation and trapping analysis around each crystal
- a dimensional comparison of three candidate flow mechanisms: capillary
  replenishment, solutal Marangoni flow and natural convection

## Repository structure

```
runDropletPipeline.m        entry point: set the data root and the droplet, then run
preprocessing/              image preparation and background subtraction
analysis/                   the measurement stages called by the pipeline
analysis/accumulation/      particle trapping analysis, run in order accum1 to accum7
figures/                    scripts that produce the figures and tables of the thesis
utils/                      shared helper functions
```

Every script carries a header explaining what it does, what it reads and what
it writes.

## Requirements

- MATLAB R2025b
- Image Processing Toolbox

The statistical routines (Pearson, Spearman) are implemented locally, so the
Statistics Toolbox is not required.

Two external tools are used to produce inputs this code then reads:

- [PIVlab](https://pivlab.blogspot.com/) for the velocity fields
- [Fiji / TrackMate](https://imagej.net/plugins/trackmate/) for the particle
  trajectories

## Input data

Each droplet lives in its own folder and must contain:

- a folder of image frames
- a plain-text file named `especificaciones` with the acquisition metadata
  (frame rate, spatial calibration, temperature, humidity, the clock time of
  the first frame, and the nucleation frames)
- the PIVlab and TrackMate exports

The specifications file is the single authoritative source for the time axis:
nothing stored in a result file is allowed to override it.

Raw image sequences are not included in this repository due to their size.

## How to run

1. Set `DATA_ROOT` at the top of `runDropletPipeline.m` to the folder holding
   the acquisition sessions.
2. Choose the droplet and switch on the stages you need, then run the pipeline.
   Repeat for each droplet.
3. Run the scripts in `figures/`, setting `DATA_ROOT` in each one.

`runDropletPipeline.m` adds the subfolders to the MATLAB path automatically,
so no manual path setup is needed.

### Order of execution within `figures/`

Most scripts are independent, but four have dependencies:

- `crystalGrowthFigures` writes `summary_crystal_growth.csv`, which
  `interDropletFigures` reads. Run it first.
- `growthVsFlowFigures` writes `Table2_percrystal_growth_vs_flow.csv`, which
  `trappingVsGrowthFigures` reads. Run it first.

Within `analysis/accumulation/`, run the seven `accum` scripts in numerical
order; each consumes what the previous one wrote.

## License

MIT — see [LICENSE](LICENSE).
