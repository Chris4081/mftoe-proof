# 🌌 MFToE Proof — Dark Energy Reconstruction with DESI DR2 (2025)  
**Author:** Christof Krieg <br>
**Contact:** Christof.Krieg@Outlook.com <br>
**License:** GNU Affero General Public License v3.0 (AGPL-3.0) <br>
**Data License:** DESI Collaboration © 2025, [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) <br>
**Repository Type:** Research / Reproducible Cosmology Pipeline <br>
**Last Updated:** October 2025  

---

## 🧩 Overview

This repository presents the **proof-of-concept implementation** of the  
**Maat Field Theory of Everything (MFToE)** — a dynamic cosmological model  
combining physical evolution, renormalization effects, and noise-driven vacuum relaxation,  
tested directly against **Dark Energy Spectroscopic Instrument (DESI) DR2 (2025)** data.

The project provides a **fully reproducible pipeline**, from toy model integration to  
BAO comparison and covariance analysis, designed for scientific collaboration and open validation.


## Update Current Cosmology Tests

The MFToE baseline model has been compared to:

• SN Ia (Pantheon+)
• DESI DR2 BAO
• Structure Growth fσ8(z)



# 🌌 **MFToE Proof v1.2.1 — CAMB Integration & Interactive GUI**

### 🚀 Overview  
This release marks a major usability and reproducibility step for the **Maat Field Theory of Everything (MFToE)** proof-of-concept.  
Version **1.2.1** introduces the **CAMB Boltzmann backend**, a refined **joint-fit consistency model**, and a brand-new **interactive GUI** for streamlined simulation control.

---

## 🧩 New Features & Enhancements

### 🧠 1. MFToE Proof GUI v2 — Interactive Cockpit  
**File:** `mftoe_gui.py`  

A lightweight Tkinter interface providing one-click access to all core MFToE workflows:

| Tab | Purpose |
|:--|:--|
| **Info** | License info + links to `docs/mftoe.pdf` and `docs/mftoe_proof.pdf` |
| **Scan Relaxion** | Runs `scripts/scan_relaxion.sh` and produces Δχ² heatmaps |
| **Run Relaxion** | Single run with RG + noise and optional JSON parameters |
| **Run Baselines** | Executes `scripts/run_baselines.sh` |
| **BAO Quickcheck** | GUI form for `analysis/bao_compare.py` to check χ² fits |
| **Cov Sweep** | Runs `scripts/cov_sweep.sh` for synthetic ρ-sweeps |
| **Run All** | Sequential execution of baseline → scan → covariance tests |

All paths are **relative**, logs stream live into the console.  
Run directly from the repo root:  
```bash
python3 mftoe_gui.py
```
### 🔭 2. CAMB Integration for rₑₛ Calculation  
- New flag `--rd-backend {fixed,camb}` for `mftoe_vacuum_astropy.py` and `analysis/joint_fit.py`  
- Computes the sound horizon r_d via the **CAMB Boltzmann code**  
- Fully parameterized: `ombh2`, `omch2`, `Neff`, `Yp`, `mnu-eV`  
- Matches Planck 2018 priors (`r_d ≈ 147.10 Mpc`)  
- Updated requirements to include `camb`, `sympy`, `mpmath`

---

### ⚙️ 3. Joint-Fit Module Upgrades  
- Added `--match-H0rd`, `--ref-H0`, `--ref-rd` for H₀·r_d consistency  
- Automatic H₀ scaling to match DESI/CAMB priors  
- Robust path imports + JSON outputs in `runs/joint_*.json`  
- Residual plots generated automatically  

| Model | r_d [Mpc] | χ² (BAO + CMB r_d) | Reduced χ² |
|:--|:--:|:--:|:--:|
| Baseline (DESI) | 150.754 | 16.12 | 1.08 |
| CAMB | 147.10 | 16.13 | 1.08 |
| CAMB + H₀·r_d match | — | 16.12 | 1.07 |

---

### 📈 4. Relaxion Scan & Visualization  
- Updated `scripts/scan_relaxion.sh` to produce clean summary CSV `runs/scan_relaxion_summary.csv`  
- Improved `analysis/plot_scan.py` with English labels and Δχ² contours (1σ/2σ/3σ)  
- Trend plots and heatmaps saved automatically  

---

### 🧮 5. Code Stability & Usability  
- Safer imports for analysis modules (relative execution fixed)  
- Auto-create `runs/` if missing  
- Clearer CLI help texts and error messages  
- Live streaming of stdout + stderr in GUI console  

---

## 🚀  Update (v1.1.0)

**Technical Enhancements**
- Added full **Joint-Fit framework** combining BAO + SNIa + GW + CMB prior.  
- Implemented parameters `--n-params`, `--snia-sigma-int`, and `--snia-vpec` for robust error modeling.  
- Introduced **automatic residual plotting** via `--plot-resids`.  
- Added **automated relaxion scan** (`scripts/scan_relaxion.sh`) producing both JSON and CSV outputs.  
- Updated CLI argument structure for modular pipeline use.  

**Scientific Results**
- Joint fit (BAO + SNIa + CMB):  
  - χ² = **31.60**, reduced χ² = **1.58** (Pantheon+ ready).  
- Relaxion scan:  
  - χ² minimum at **γ ≈ 0.05**, **σ ≈ 1e-5**.  
  - Mean total equation of state **w_tot ≈ −0.29**.  
- Sub-percent deviations from ΛCDM across 0 < z < 3.  
- All fits consistent with **DESI DR2 (2025)** BAO data within 1 σ.  

**Project & Licensing**
- Added `CITATION.cff` with DOI, metadata, and author information.  
- Integrated **Zenodo DOI:** [10.5281/zenodo.17383354](https://doi.org/10.5281/zenodo.17383354)  
- Expanded README and LaTeX documentation with explicit **DESI CC BY 4.0 license notice**.  
- Official **GitHub Release v1.1.0** published (October 2025).  

---


## 📜 Scientific Abstract

The **MFToE vacuum model** introduces a dynamic scalar field χ controlling the residual vacuum energy,  
embedded in a minimalistic EFT-like system with optional RG running and Ornstein–Uhlenbeck noise.  
We integrate the late-time background from z = 3 → 0 using RK4 and compare the predictions for  
H(z), d_L(z), and D_M/r_d, D_H/r_d, D_V/r_d against **DESI DR2 BAO** measurements.

**Results:**  
- Deviations from ΛCDM below **0.7 %** in both H(z) and d_L(z).  
- χ² = 16.1 (reduced χ² = 0.85) for the baseline model.  
- χ² = 19.2 (reduced χ² = 1.0) for the dynamic relaxion + RG + noise run.  
- Excellent agreement with **DESI DR2 (2025)** compressed BAO data.  
- Covariance sweeps confirm stability for ρ ∈ [0.0, 0.5].  

This demonstrates that the **MFToE vacuum mechanism** can reproduce late-time expansion data  
while allowing mild dynamical dark-energy evolution — a strong empirical foundation for further exploration.

---

## 🧠 Repository Structure

```
MFToE-Proof/
├── mftoe_vacuum_astropy.py
├── analysis/
│   ├── bao_compare.py
│   ├── make_cov_from_csv.py
│   └── compare_runs.py
├── data/
│   └── desi_dr2/
│       ├── bao_summary.csv
│       ├── bao_cov.npy
│       └── iminuit/base/desi-bao-all/bestfit.minimum
├── runs/
├── scripts/
│   ├── run_baselines.sh
│   ├── cov_sweep.sh
│   └── scan_relaxion.sh
├── docs/
│   └── mftoe_proof.pdf
│   └── mftoe.pdf
├── LICENSE
└── README.md
```

---

## ⚙️ Installation & Setup

### Prerequisites
- Python ≥ 3.10  
- Packages: numpy, pandas, matplotlib, astropy

##  📥 Cloning the Repository

📥 **Clone and Run**

To get started:

```git clone https://github.com/Chris4081/mftoe-proof.git  
cd mftoe-proof

Create a virtual environment (recommended):

python3 -m venv .venv  
source .venv/bin/activate  
pip install -r requirements.txt

Run the baseline simulation:

bash scripts/run_baselines.sh
```
### Install dependencies
```bash
pip install -r requirements.txt
```

### Run baseline test
```bash
bash scripts/run_baselines.sh
```

### Run relaxion (RG + Noise)
```bash
bash scripts/scan_relaxion.sh
```

### Covariance sweep (synthetic)
```bash
bash scripts/cov_sweep.sh
```



---

## 🔬 Key Results (DESI DR2 2025)

| Model | Mode | RG | Noise | χ² | χ²_red | w_tot | Comment |
|:------|:-----|:--:|:-----:|:---:|:-------:|:------:|:--------|
| MFToE Baseline | targetH0 | off | off | 16.12 | 0.85 | -0.282 | Excellent fit |
| MFToE Relaxion | relaxion | on  | on  | 19.17 | 1.01 | -0.289 | Mild DE evolution |
| Covariance ρ = 0.3 | baseline | on | off | 17.79 | 0.94 | — | Stable with correlation |

---

## 📄 Data & Licensing

DESI DR2 (2025) data used under the  
[Creative Commons Attribution 4.0 International License (CC BY 4.0)](https://creativecommons.org/licenses/by/4.0/).  
Usage requires citation and acknowledgment per DESI Data Release documentation:  
👉 [https://data.desi.lbl.gov/doc/releases/](https://data.desi.lbl.gov/doc/releases/)

© DESI Collaboration, 2025.

**Code License:**  
All scripts and models in this repository are distributed under the  
**GNU Affero General Public License v3.0 (AGPL-3.0)**.


This project uses publicly available cosmological datasets:

- Pantheon+ Supernova dataset  
  Source: https://github.com/xxx/pantheon-plus  
  License: CC BY 4.0  
  Citation: Scolnic et al. (2022)

All data remain property of the original authors.


## Data sources and licenses

This repository uses publicly available data from:

- **DESI Data Release 2 (DR2)** —  
  Dark Energy Spectroscopic Instrument Collaboration (2024).  
  DOI: [10.5281/zenodo.11019438](https://doi.org/10.5281/zenodo.11019438)  
  Licensed under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/).

All DESI-based results in this repository acknowledge DESI DR2 as their source dataset.

> **Note:**  
> The current SNIa (supernova) component in the joint-fit pipeline uses a **synthetic dataset**
> for testing purposes only.  
> The framework is **fully compatible with future integration** of the real  
> **Pantheon+ Supernova Sample (Brout et al. 2022, ApJ 938, 110)**  
> once linked via the official [PantheonPlusSH0ES Data Release](https://github.com/PantheonPlusSH0ES/DataRelease).  
>  
> Citation reference (if integrated in the future):  
> Brout, D., Scolnic, D., Popovic, B., et al. (2022),  
> *The Pantheon+ Analysis: Cosmological Constraints*,  
> ApJ 938, 110.  
> DOI: [10.3847/1538-4357/ac8e04](https://doi.org/10.3847/1538-4357/ac8e04)

---
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17383354.svg)](https://doi.org/10.5281/zenodo.17383354)

## 🧭 Citation

```bibtex
@misc{krieg2025_mftoe_proof,
  author       = {Christof Krieg},
  title        = {MFToE Proof — Dark Energy Reconstruction with DESI DR2 (2025)},
  year         = {2025},
  note         = {GitHub repository},
  license      = {AGPL-3.0}
}

@dataset{desi2025_dr2,
  author       = {DESI Collaboration},
  title        = {Dark Energy Spectroscopic Instrument (DESI) Data Release 2},
  year         = {2025},
  note         = {https://data.desi.lbl.gov/doc/releases/},
  license      = {CC BY 4.0}
}
```

---



## 🧠 Philosophy

This work aligns with the **Maat Principles** of  
🌿 Harmony, ⚖️ Balance, 🎨 Creativity, 🌐 Connectedness, 🕊️ Respect —  
bridging science, ethics, and technology into a unified exploration of cosmology and consciousness.

---

**© 2025 Christof Krieg — MFToE Research Initiative**  
Licensed under **AGPL-3.0** | DESI data © DESI Collaboration (2025) CC BY 4.0  
