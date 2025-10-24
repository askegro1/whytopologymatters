# README — Figure 3 & Supplementary Figure 1

This repository contains the MATLAB script `Fig3_SuppFig1.m` and its helper subfunctions for generating Figure 3 and Supplementary Figure 1. 
The script reads a single precomputed dataset, produces several invisible plots, and then assembles two composite figures out of the individual plots. 
Final outputs are exported as vector PDFs (and `.fig` for reproducibility).

---

## 1) Overview

- Input (required): `results_EFC_2025_06.mat` containing a table variable named `resultsEFC_2025_06`.
- Outputs (auto-created in folder `PlotResults_Fig3_SuppFig1/`)  
  - `Figure_3_pdfformat.pdf` and `Figure_3_figformat.fig`  
  - `Supplementary_Figure_1_pdfformat.pdf` and `Supplementary_Figure_1_figformat.fig`
- Run command: open MATLAB in the folder and run:
  ```matlab
  Fig3_SuppFig1  ```

---

## 2) Software Environment

### Recommended MATLAB
- Tested on MATLAB R2023b.  

### Toolboxes
- Statistics and Machine Learning Toolbox (for `quantile`).

### OS
- Windows, macOS, or Linux. 

---

## 3) Input Data Contract

The MAT file `results_EFC_2025_06.mat` must define a variable named `resultsEFC_2025_06` that is a MATLAB table with (at least) the following columns:

| Column       | Type    | Units / Notes                                                                 |
|--------------|---------|--------------------------------------------------------------------------------|
| `Chemistry`  | string / char | Expected values include `CHEM_1`/`LFP` and `CHEM_2`/`NMC`. |
| `TC`         | double  | Manufacturing variability class. The script expects 1 or 2.            |
| `Ns`         | double  | Series cell count (e.g., 50, 100, 200).                                        |
| `Temp`       | double  | Mean temperature.                                                |
| `Tsig`       | double  | Temperature standard deviation.             |
| `Trest`      | double  | Rest time fraction.                |
| `meanEFC`    | double  | Output metric (mean), percent scale.                                          |
| `stdEFC`     | double  | Output metric (std), percent scale.                                           |

---

## 4) How to Run (Step‑by‑Step)

1. Launch **MATLAB** and set the working directory to the folder containing `Fig3_SuppFig1.m` and `results_EFC_2025_06.mat`.
2. Run the script (from the Command Window):
   ```matlab
   Fig3_SuppFig1
   ```
3. On success you should see the message:
   ```
   Figure 3 and Supplementary Figure 1 exported.
   ```
4. Open the outputs in `PlotResults_Fig3_SuppFig1/`:
   - `Figure_3_pdfformat.pdf`, `Figure_3_figformat.fig`
   - `Supplementary_Figure_1_pdfformat.pdf`, `Supplementary_Figure_1_figformat.fig`

> Component figures are created invisible and are not exported individually by default.


---

## 5) What the Script Does

1. Loads the table from `results_EFC_2025_06.mat` and validates required columns.
2. Splits data by Chemistry (LFP/NMC) and manufacturing variability class `TC` (1/2).
3. Generates invisible component figures:
   - Sensitivity to manufacturing variability (Lower vs Higher Manufacturing Tolerance)
   - Parameter sensitivities (Mean temperature, Temperature variation, Rest time, Number of cells) per chemistry
   - Chemistry × voltage group comparisons (low‑V / 400V / 800V)
   - Two‑chemistry comparison at nominal thermal conditions (LFP vs NMC)
4. Assembles composite figures.
5. Exports PDFs + MATLAB `.fig` for each composite figure.

---

## 6) Troubleshooting

- Figures appear off‑screen on multi‑monitor setups:  
  → Set `targetMonitor = 1` near the end of the script.


---
