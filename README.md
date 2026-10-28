# Low-cost NanoVNA-based stepped-frequency GPR for plastic landmine surrogates

MATLAB code and measured data for the paper:

> A. Alareqi, R. Szabolcsi, "Low-Cost NanoVNA-Based Stepped-Frequency GPR for Plastic Landmine Surrogates: Template-Based Detection in Air and Practical Limits in a Sandbox," submitted to *Radioengineering*, 2026.

## Hardware
- NanoVNA-H4 (firmware 1.2.44), SOLT + isolation calibration, 0.6–1.5 GHz, 201 points, IF bandwidth 100 Hz
- Two printed log-periodic dipole array (LPDA) antennas, 600–6000 MHz, bistatic, 24 cm apart
- Port 1 = transmit antenna, Port 2 = receive antenna, USB serial (COM port set in each script)

## Repository layout
| Folder | Content |
|---|---|
| `measure/` | Scripts that drive the NanoVNA over USB and save `radar_*.mat` / `bscan_*.mat` |
| `analysis/` | Offline analysis (no hardware needed): matched filter, surface-change removal, sparse deconvolution, B-scan processing |
| `data/in_air/` | In-air tests (tray, water, oil, empty bottle, empty runs, ranging) |
| `data/sandbox/` | Static sandbox session (surface targets, buried bottles, dig-and-refill controls) |
| `data/bscan/` | Sandbox B-scans and the kernel file `bscan_kernel.mat` |
| `figures/` | Figures of the paper (600 dpi) |

## Data format
Each `radar_<label>.mat` contains
- `f` – frequency vector (Hz), 1 × 201
- `Sref` – complex S21 of the empty-scene baseline
- `Stgt` – complex S21 of the scene with the target
- `info` (sandbox files) – label, note, time, baseline age

Each `bscan_<label>.mat` contains `f`, `S` (positions × frequencies, complex S21), `x_cm` (positions) and `info`.

## How to run the analysis
1. Open MATLAB (R2021a or newer) and add the repository to the path: `addpath(genpath(pwd))`
2. `cd` into the data folder you want to analyse.
3. Run, for example, `analyze_mf`, `analyze_sand3` or `bscan_analyze('water3', {'empty_a','empty_b'})`.
Settings (file names, windows) are at the top of each script.

## How to measure
Edit the COM port at the top of the script, then run e.g. `radar_measure` (label `baseline` first, then one label per target) or `sand_session`.

## License
Code: MIT License (see `LICENSE`). Data: CC BY 4.0.

## Citation
If you use this code or data, please cite the paper above (DOI to be added after publication).
