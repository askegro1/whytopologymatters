# README — Figure 4 & Supplementary Table S2

This folder contains the MATLAB script `Fig4_SuppTableS2.m` and embedded helpers to:
1) Build an EV metadata file from a spreadsheet (`EV_Data.mat`), and
2) Load the results summary, prepare subsets, fit trends, overlay EVs, and export Figure 4 (PDF/FIG) plus Supplementary Table S2 (LaTeX).

---

## 1) Overview

Inputs (required in the working folder):
- `EVs_Voltages_Nominal.xlsx`  
- `results_EFC_2025_06.mat` (must contain a table variable named `resultsEFC_2025_06`)

Outputs (auto-created):
- `EV_Data.mat`
- `PlotResults_Fig4_SuppTable2/Figure4_pdfformat.pdf`
- `PlotResults_Fig4_SuppTable2/Figure4_figformat.fig`
- `PlotResults_Fig4_SuppTable2/SuppTable2_Data.tex`

Run (from MATLAB, with the working directory set to this folder):
```matlab
Fig4_SuppTableS2
```

---

## 2) Software Environment

- MATLAB: tested on R2023b.
- Toolboxes: None required.
- OS: Windows, macOS, or Linux.

---

## 3) What the Script Does

### Part 1 — Build `EV_Data.mat` from Excel
- Reads the first 6 columns of `EVs_Voltages_Nominal.xlsx` and renames them:
  `Brand, Model, Type, Year, VNom, ChemistryRaw`.
- Constructs `FullName = "<Brand> <Model>"`.
- Keeps only the columns needed downstream and saves:
  ```matlab
  EV_Data = table(FullName, Chemistry, VNom);
  save('EV_Data.mat','EV_Data');
  ```

### Part 2 — Load results, prepare data, plot, fit, overlay, export
- Loads `results_EFC_2025_06` from `results_EFC_2025_06.mat` and creates chemistry-specific subsets:
  - `preprocessChemistry` filters by chemistry key (`"CHEM_1"` → LFP, `"CHEM_2"` → NMC).
- Extracts TC = 2 slices with exact thermal settings:
  - LFP: `Temp=25`, `Tsig=0.42`, `Trest=0.95`
  - NMC: `Temp=25`, `Tsig=0.00`, `Trest=0.95`
- Computes NominalVoltage from chemistry and `Ns` via `convertVoltageVector`:
  - Rules:  
    - LFP, `Ns=4` → ~12 V (≈ 4 × 3.2 V)  
    - NMC, `Ns=4` → ~15 V (≈ 4 × 3.7 V)  
    - LFP, `Ns=16` → ~50 V (≈ 16 × 3.2 V)  
    - NMC, `Ns=14` → ~50 V (≈ 14 × 3.6 V)  
    - Fallback: `4 V × Ns` if no explicit rule matches.
- Sensitivity plot (`plotSensitivity`):
  - X = `NominalVoltage`; Y = `meanEFC` (std assumed in `stdEFC`).
  - Computes per-chemistry medians at each voltage grid point.
  - Selects the best trend among `{log, sqrt, power, poly2}` by lowest RSS using `fminsearch` and reports RMSE/MAE/R² in a summary table.
  - Draws the best-fit curve per chemistry and translucent median ± median(std) bands.
  - Overlays EV markers from `EV_Data.mat`: if data exists at the same X, uses the median; otherwise, evaluates the best-fit curve at that voltage.
- Exports:
  - `Figure4_pdfformat.pdf` and `Figure4_figformat.fig` to `PlotResults_Fig4_SuppTable2/`
  - `SuppTable2_Data.tex` summarizing model fits for both chemistries

---

## 4) Data Contract

### Input summary table (`results_EFC_2025_06`)
Must be a MATLAB table with (at least) the columns below:
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

### EV spreadsheet (`EVs_Voltages_Nominal.xlsx`)
- The script expects the first 6 columns to be `Brand, Model, Type, Year, VNom, ChemistryRaw` (in that order).
- `VNom` should be a numeric nominal pack voltage (in volts).

---

## 5) How to Run (Step-by-Step)

1. Place `Fig4_SuppTableS2.m`, `EVs_Voltages_Nominal.xlsx`, and `results_EFC_2025_06.mat` in the same folder.
2. Open MATLAB and set the working directory to that folder.
3. Run:
   ```matlab
   Fig4_SuppTableS2
   ```
4. Console should print progress, e.g.:
   ```
   Step 1/2: Building EV_Data.mat from "EVs_Voltages_Nominal.xlsx"...
     -> Saved EV_Data.mat (### rows)
   Step 2/2: Loading results and plotting sensitivity...
   All done. Results saved under "PlotResults_Fig4_SuppTable2/".
   ```
5. Open the outputs in `PlotResults_Fig4_SuppTable2/`:
   - `Figure4_pdfformat.pdf`, `Figure4_figformat.fig`
   - `SuppTable2_Data.tex`
   - (The helper `EV_Data.mat` is written in the working folder.)
