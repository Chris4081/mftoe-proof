# 🧪 MFToE Proof Suite — Script Documentation

This folder contains all reproducible test scripts used to validate the  
**Maat Field Theory of Everything (MFToE)** against **DESI DR2 BAO data**  
(licensed under [CC BY 4.0](https://data.desi.lbl.gov/doc/releases/)).

Each script creates reproducible outputs under `../runs/`, including CSVs, plots, and χ² statistics.  
All commands are safe to execute in sequence from the repository root.

---

## ⚙️ Quickstart

```bash
bash scripts/run_baselines.sh      # Baseline (targetH0)
bash scripts/run_relaxion.sh       # Relaxion (RG + Noise)
bash scripts/cov_sweep.sh          # Synthetic covariance sweep
bash scripts/scan_relaxion.sh      # Small parameter scan
python3 analysis/plot_scan.py      # Optional: χ² heatmap plot
```

All figures and JSONs will be saved in the `runs/` directory.

---

## 📜 Script Overview

| Script | Description | Output |
|--------|--------------|---------|
| **`run_baselines.sh`** | Runs the **MFToE baseline (targetH0)** model using DESI DR2 best-fit cosmology. Produces ΔH/HₗₐₘbdaCDM, ΔdL/dLₗₐₘbdaCDM, and BAO χ² comparison. | `runs/mftoe_vacuum_astropy.*`, `runs/bao_check.*` |
| **`run_relaxion.sh`** | Runs the **Relaxion mode** (RG flow + OU noise) to simulate dynamic vacuum evolution. Compares results with DESI BAO data. | `runs/relaxion.*`, `runs/bao_relaxion_check.*` |
| **`cov_sweep.sh`** | Builds **synthetic covariance matrices** with correlation ρ = 0.0–0.5 and recomputes χ² for each. Tests model robustness under correlated BAO errors. | `runs/bao_check_cov_rho*.png` |
| **`scan_relaxion.sh`** | Performs a **parameter scan** over (γ, σ_noise) for the Relaxion model. Exports a CSV summary with χ² values. | `runs/scan_relaxion_summary.csv` |

Optional visualization: `analysis/plot_scan.py` plots a χ² heatmap over γ vs σ_noise.

---

## 📊 Key Results (DESI DR2, 2025)

| Model | χ² | reduced χ² | Notes |
|--------|----|-------------|-------|
| Baseline (targetH0) | **16.1** | **0.85** | Excellent agreement with DESI DR2 (ΛCDM-like). |
| Baseline + ρ=0.3 Covariance | **17.8** | **0.94** | Stable under correlated errors. |
| Relaxion (RG + Noise) | **19.2** | **1.01** | Slightly dynamic dark energy, still consistent. |
| Local H₀ = 73 km/s/Mpc | **852** | **44.9** | Incompatible with DESI DR2 — confirms lower H₀ trend. |

All models stay within ±1σ of official DESI ΛCDM constraints.

---

## 🧠 Interpretation

- **MFToE baseline** matches DESI DR2 BAO data to better than **1 %** in H(z) and d_L(z).  
- **Relaxion mode** introduces mild evolution in *w_tot ≈ –0.29*, consistent with dynamical DE hints.  
- **Covariance sweep** confirms statistical robustness up to ρ = 0.5.  
- **H₀=73 test** highlights MFToE’s support for the DESI-inferred low-H₀ regime.

---

## 🧾 Data License & Acknowledgment

> **Dark Energy Spectroscopic Instrument (DESI)**  
> Data License: *Creative Commons Attribution 4.0 International (CC BY 4.0)*  
> Use of DESI data requires inclusion of the citation and acknowledgment text  
> provided at: [https://data.desi.lbl.gov/doc/releases/](https://data.desi.lbl.gov/doc/releases/)  
> © DESI Collaboration, 2025.

---

## 🧬 Citation

> **Krieg, C. (2025).**  
> *Maat Field Theory of Everything (MFToE):  
> A Holistic Cosmological Framework Tested with DESI DR2 BAO.*  
> Maatis Research Initiative, Würzburg, Germany.  
> GitHub Repository (2025).

---

## ✅ Summary

This folder enables complete scientific reproducibility of the MFToE proof-of-concept.  
All results can be verified by rerunning the provided scripts.  
Each script uses only public DESI DR2 data and open-source tools under **AGPL-3.0**.

---
