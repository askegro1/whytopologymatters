# README — Figure 5, Supplementary Figure 2 & Supplementary Figure 3

This folder contains the MATLAB script `Fig5_SuppFig2_SuppFig3.m`, which performs a global economic sensitivity analysis for reconfigurable vs conventional EV battery systems and generates:

- Figure 5 (multivariate sensitivity & cost breakdown),
- Supplementary Figure 2 (additional 1D sensitivity results),
- Supplementary Figure 3 (dominant parameter rectangles in 2D feature space).

Outputs are exported as publication-ready PDF (A4 size) and MATLAB .fig files.

---

## 1) Overview

1. Defines baseline economic parameters (`getParams`) and parameter ranges (`getParamRanges`).
2. Samples 100,000 configurations using Latin Hypercube Sampling:
   - Vehicle lifetime, driving distance, pack cost factor `ν`, ohmic losses `δ_loss`, discount rates, etc.
3. For every sampled case:
   - Computes Net Present Cost (NPC) for both CBP (conventional) and RBP (reconfigurable battery pack),
   - Calculates the NPC gain ΔNPC = NPC_CBP − NPC_RBP.
4. Identifies the best, baseline, and worst cases from ΔNPC.
5. Generates:
   - Figure 5: Cost breakdown (a), cumulative NPC curves (b–d), tornado plot (e), and 1D sensitivity panels (f–j).
   - Supp. Fig. 2: Sensitivity to remaining cost drivers (δ_loss, α difference, r, Y_EV, SOH ends).
   - Supp. Fig. 3: Regions in 2D parameter space where >99.7% of samples give ΔNPC > 0.

---

## 2) Required Inputs

This script is self-contained — no external `.mat` or `.xlsx` files are needed.

All parameters are generated internally using `getParams` and `getParamRanges`.

---

## 3) Software Environment

| Requirement | Details |
|-------------|---------|
| **MATLAB Version** | Tested on MATLAB R2023b.  |
| **Toolboxes Required** | Statistics & Machine Learning Toolbox (for `lhsdesign`) |
| **OS Compatibility** | Windows, macOS, Linux |

---

## 4) How to Run (Step-by-step)

1. Open MATLAB.
2. Set the working directory to the folder containing `Fig5_SuppFig2_SuppFig3.m`
3. Run the script:
   ```matlab
   Fig5_SuppFig2_SuppFig3
   ```

Note: With `n_samples = 100000`, runtime can be 5-19 minutes, depending on CPU.

---

## 5) Output Files

All results are saved inside:

```
PlotResults_Fig5_SuppFig2_SuppFig3/
```

| File | Description |
|------|-------------|
| `Figure5_pdfformat.pdf` | Main figure for paper (4×3 tiled layout) |
| `Figure5_figformat.fig` | MATLAB figure (editable, full resolution) |
| `SuppFig2_pdfformat.pdf` | Supplementary Figure 2: 1D sensitivity analyses |
| `SuppFig2_figformat.fig` | MATLAB figure for Supp. Fig. 2 |
| `SuppFig3_pdfformat.pdf` | Supplementary Figure 3: Dominant rectangles in parameter space |
| `SuppFig3_figformat.fig` | MATLAB figure for Supp. Fig. 3 |

---

## 6) Optional Configuration (inside the script)

You can modify:

| Variable | Purpose | Default |
|----------|---------|---------|
| `n_samples` | Number of LHS samples | `100000` |
| `targetMonitor` | Monitor index for figure placement | `2` |
| `paramRanges` | Parameter limits (via `getParamRanges`) | See script |
| `paperWidth`, `paperHeight` | A4 figure size for export | 21 × 29.7 cm |
| `successThreshold` | % of ΔNPC > 0 for dominant rectangle search | `0.997` (99.7%) |

---

## 7) Troubleshooting

| Issue | Solution |
|-------|----------|
| Figure windows open on wrong screen | Set `targetMonitor = 1` in script |
| Script runs slowly | Reduce `n_samples` to e.g. `10000` |

---

## 8) Reproducibility Notes

- Random sampling is reproducible (`rng(42)` fixed in code).
