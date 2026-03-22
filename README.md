# 🌌 MFToE Proof — Dark Energy Reconstruction with DESI DR2 (2025)  
**Author:** Christof Krieg <br>
**Contact:** Christof.Krieg@Outlook.com <br>
**License:** GNU Affero General Public License v3.0 (AGPL-3.0) <br>
**Data License:** DESI Collaboration © 2025, [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) <br>
**Repository Type:** Research / Reproducible Cosmology Pipeline <br>
**DOI:** [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17383354.svg)](https://doi.org/10.5281/zenodo.17383354) <br>
**Last Updated:** March 2026

---

## 🧩 Overview

This repository presents the **proof-of-concept implementation** of the  
**Maat Field Theory of Everything (MFToE)** — a dynamic cosmological model  
combining physical evolution, renormalization effects, and noise-driven vacuum relaxation,  
tested directly against **Dark Energy Spectroscopic Instrument (DESI) DR2 (2025)** data.

The project provides a **fully reproducible pipeline**, from toy model integration to  
BAO comparison and covariance analysis, designed for scientific collaboration and open validation.

---

## 📊 Current Cosmology Tests

The MFToE baseline model has been compared to:

| Dataset | Status | Result |
|:--------|:------:|:-------|
| DESI DR2 BAO | ✅ active | χ² = 16.1, reduced χ² = 0.85 |
| SN Ia (Pantheon+) | ⚙️ synthetic | framework compatible |
| Structure Growth fσ₈(z) | ⚙️ planned | — |

---

## 🔬 New in v1.2.1 — CAMB Integration & Interactive GUI

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

```bash
python3 mftoe_gui.py
```

### 🔭 2. CAMB Integration for r_d Calculation

- New flag `--rd-backend {fixed,camb}` for `mftoe_vacuum_astropy.py` and `analysis/joint_fit.py`
- Computes the sound horizon r_d via the **CAMB Boltzmann code**
- Matches Planck 2018 priors (`r_d ≈ 147.10 Mpc`)

| Model | r_d [Mpc] | χ² (BAO + CMB r_d) | Reduced χ² |
|:--|:--:|:--:|:--:|
| Baseline (DESI) | 150.754 | 16.12 | 1.08 |
| CAMB | 147.10 | 16.13 | 1.08 |
| CAMB + H₀·r_d match | — | 16.12 | 1.07 |

### ⚙️ 3. Joint-Fit Module Upgrades

- Added `--match-H0rd`, `--ref-H0`, `--ref-rd` for H₀·r_d consistency
- Automatic H₀ scaling to match DESI/CAMB priors
- JSON outputs in `runs/joint_*.json` + automatic residual plots

---

## 🆕 New Tests (2026)

### MFToE Phase Analysis — Quartic Toy Model

A separate numerical analysis pipeline (`mftoe_phase_v2.py`) implements a
minimal phenomenological quartic free-energy model for structural selection:

```
F_red(λ, η; m, u, v) = −λ(m² + u²) + η·v² + m⁴ + u⁴ + v⁴
```

**Key results:**
- Exact analytic phase boundary: `η_c(λ) = 0` (independent of λ)
- Two-sector structure: aligned (`v ≈ 0` for η ≥ 0) and frustrated (`|v| > 0` for η < 0)
- Symmetry reduction v → |v| eliminates artificial branch doubling
- Local boundary scan: λ ∈ [0.5, 2.0], η ∈ [−0.5, 0.1] at fine resolution
- Polynomial fit of phase boundary: `η_c(λ) ≈ a₀ + a₁λ + a₂λ²`
- Numerical pipeline: multi-start root finding, branch tracking, phase diagram

**Companion paper:**  
*A Minimal State–Structure Model for Dynamical Structural Selection* (Krieg, 2026)

### Boundary Scan & Phase Diagram

`mftoe_phase_v2.py` implements a full 5-step pipeline:

1. **Branch Tracking** — multi-start root finding with symmetry reduction
2. **Transition Analysis** — switch-filter for real vs. numerical transitions
3. **Phase Diagram** — 4-class classification (aligned/frustrated × single/multi)
4. **Boundary Fit** — `η_c(λ)` extracted by linear interpolation + polynomial fit
5. **Overlay** — phase diagram with fitted boundary curve

All plots automatically exported to `mftoe_output/`.

---

## 📜 Scientific Abstract

The **MFToE vacuum model** introduces a dynamic scalar field χ controlling the residual vacuum energy,  
embedded in a minimalistic EFT-like system with optional RG running and Ornstein–Uhlenbeck noise.  
We integrate the late-time background from z = 3 → 0 using RK4 and compare the predictions for  
H(z), d_L(z), and D_M/r_d, D_H/r_d, D_V/r_d against **DESI DR2 BAO** measurements.

**Results:**  
- Deviations from ΛCDM below **0.7%** in both H(z) and d_L(z)
- χ² = 16.1 (reduced χ² = 0.85) for the baseline model
- χ² = 19.2 (reduced χ² = 1.0) for the dynamic relaxion + RG + noise run
- Excellent agreement with **DESI DR2 (2025)** compressed BAO data
- Covariance sweeps confirm stability for ρ ∈ [0.0, 0.5]

---

## 🔬 Key Results (DESI DR2 2025)

| Model | Mode | RG | Noise | χ² | χ²_red | w_tot | Comment |
|:------|:-----|:--:|:-----:|:---:|:-------:|:------:|:--------|
| MFToE Baseline | targetH0 | off | off | 16.12 | 0.85 | −0.282 | Excellent fit |
| MFToE Relaxion | relaxion | on  | on  | 19.17 | 1.01 | −0.289 | Mild DE evolution |
| Covariance ρ = 0.3 | baseline | on | off | 17.79 | 0.94 | — | Stable |
| CAMB + H₀·r_d | joint | off | off | 16.12 | 1.07 | — | CAMB backend |

---

## 🧠 Repository Structure

```
MFToE-Proof/
├── mftoe_vacuum_astropy.py       # Main cosmology simulation
├── mftoe_gui.py                  # Interactive GUI (v2)
├── mftoe_phase_v2.py             # Phase analysis pipeline (2026)
├── analysis/
│   ├── bao_compare.py
│   ├── joint_fit.py
│   ├── make_cov_from_csv.py
│   └── compare_runs.py
├── data/
│   └── desi_dr2/
│       ├── bao_summary.csv
│       ├── bao_cov.npy
│       └── iminuit/base/desi-bao-all/bestfit.minimum
├── mftoe_output/                 # Phase analysis plots (auto-generated)
├── runs/
├── scripts/
│   ├── run_baselines.sh
│   ├── cov_sweep.sh
│   └── scan_relaxion.sh
├── docs/
│   ├── mftoe_proof.pdf
│   └── mftoe.pdf
├── DATA_LICENSES.md              # Full data license documentation
├── LICENSE
└── README.md
```

---

## ⚙️ Installation & Setup

### Prerequisites
- Python ≥ 3.10
- Packages: `numpy`, `pandas`, `matplotlib`, `astropy`, `scipy`
- Optional: `camb`, `sympy`, `mpmath`

```bash
git clone https://github.com/Chris4081/mftoe-proof.git
cd mftoe-proof
python3 -m venv .venv
source .venv/bin/activate
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

### Run phase analysis (new)
```bash
python3 mftoe_phase_v2.py
```

---

## 📄 Data & Licensing

This project uses publicly available cosmological datasets.  
See [`DATA_LICENSES.md`](DATA_LICENSES.md) for full details.

### DESI Data Release 2 (DR2)
- **Source:** https://data.desi.lbl.gov/doc/releases/
- **DOI:** https://doi.org/10.5281/zenodo.11019438
- **License:** CC BY 4.0 — © DESI Collaboration

Usage requires citation of the DESI DR2 data release and the relevant cosmology papers (see below).

### Pantheon+ Supernova Dataset
- **Source:** https://github.com/PantheonPlusSH0ES/DataRelease
- **License:** Data release — citation required

> **Note:** Pantheon+ is currently not used in the numerical pipeline.  
> A synthetic dataset is used instead. Integration is planned.

### Code License
All code in this repository is licensed under **GNU Affero General Public License v3.0 (AGPL-3.0)**.

---

## 🧭 Citation

```bibtex
@misc{krieg2026_mftoe_proof,
  author  = {Christof Krieg},
  title   = {MFToE Proof — Dark Energy Reconstruction with DESI DR2},
  year    = {2026},
  doi     = {10.5281/zenodo.17383354},
  url     = {https://github.com/Chris4081/mftoe-proof},
  license = {AGPL-3.0}
}

@article{DESI2025BAO,
  author  = {{DESI Collaboration}},
  title   = {{DESI DR2 Results II: Measurements of Baryon Acoustic Oscillations
              and Cosmological Constraints}},
  journal = {Phys. Rev. D},
  volume  = {112},
  pages   = {083515},
  year    = {2025},
  doi     = {10.1103/tr6y-kpc6},
  eprint  = {2503.14738}
}

@article{DESI2025Lya,
  author  = {{DESI Collaboration}},
  title   = {{DESI DR2 Results I: Baryon Acoustic Oscillations
              from the Lyman Alpha Forest}},
  journal = {Phys. Rev. D},
  volume  = {112},
  pages   = {083514},
  year    = {2025},
  doi     = {10.1103/PhysRevD.112.083514},
  eprint  = {2503.14739}
}

@article{Brout2022,
  author  = {Brout, D. and others},
  title   = {{The Pantheon+ Analysis: Cosmological Constraints}},
  journal = {Astrophys. J.},
  volume  = {938},
  pages   = {110},
  year    = {2022},
  doi     = {10.3847/1538-4357/ac8e04},
  eprint  = {2202.04077}
}

@article{Scolnic2022,
  author  = {Scolnic, D. and others},
  title   = {{The Pantheon+ Analysis: The Full Data Set and Light-curve Release}},
  journal = {Astrophys. J.},
  volume  = {938},
  pages   = {113},
  year    = {2022},
  doi     = {10.3847/1538-4357/ac8b7a},
  eprint  = {2112.03863}
}
```

---

## 🧠 Philosophy

This work aligns with the **Maat Principles** of  
🌿 Harmony, ⚖️ Balance, 🎨 Creativity, 🌐 Connectedness, 🕊️ Respect —  
bridging science, ethics, and technology into a unified exploration of cosmology and consciousness.

---

**© 2026 Christof Krieg — MFToE Research Initiative**  
Licensed under **AGPL-3.0** | DESI data © DESI Collaboration (2025) CC BY 4.0
