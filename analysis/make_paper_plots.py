#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

C_KMS = 299_792.458  # km/s

# -----------------------------
# Config
# -----------------------------
MFTOE_RUN = "runs/mftoe_baseline.csv"
SN_DATA = "data/pantheon_plus.csv"
OUT_PREFIX = "runs/paper"

H0 = 70.0
Omega_m = 0.3

# -----------------------------
# Helpers
# -----------------------------
def cumtrapz(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    out = np.zeros_like(x)
    if len(x) > 1:
        dx = np.diff(x)
        mid = 0.5 * (y[1:] + y[:-1])
        out[1:] = np.cumsum(dx * mid)
    return out

def lcdm_curves(z, H0=70.0, Om0=0.3):
    cosmo = FlatLambdaCDM(H0=H0 * (u.km/u.s/u.Mpc), Om0=Om0)
    Hz = cosmo.H(z).value / H0
    dL_mpc = cosmo.luminosity_distance(z).value
    dL_dimless = dL_mpc * (H0 / C_KMS)
    return Hz, dL_dimless

def mu_from_dL_dimless(dL_dimless, H0=70.0):
    dL_mpc = np.asarray(dL_dimless) * (C_KMS / H0)
    return 5.0 * np.log10(np.clip(dL_mpc, 1e-12, None)) + 25.0

def load_snia_csv(path):
    df = pd.read_csv(path)
    for col in ["z", "mu", "mu_err"]:
        if col not in df.columns:
            raise ValueError(f"SN Ia CSV missing column '{col}'")
    return df.sort_values("z").reset_index(drop=True)

# -----------------------------
# Load MFToE
# -----------------------------
mf = pd.read_csv(MFTOE_RUN)
z = mf["z"].values
H_mf = mf["H_over_H0"].values
dL_mf = mf["dL"].values
w_phi = mf["w_phi"].values
w_tot = mf["w_tot"].values

# ΛCDM reference
H_lcdm, dL_lcdm = lcdm_curves(z, H0=H0, Om0=Omega_m)

# Residuals
dH_rel = (H_mf - H_lcdm) / H_lcdm
ddL_rel = (dL_mf - dL_lcdm) / np.clip(dL_lcdm, 1e-30, None)

# Distance modulus
mu_mf = mu_from_dL_dimless(dL_mf, H0=H0)
mu_lcdm = mu_from_dL_dimless(dL_lcdm, H0=H0)
dmu = mu_mf - mu_lcdm

# -----------------------------
# Plot 1: H(z) residuals
# -----------------------------
plt.figure(figsize=(8, 5))
plt.plot(z, dH_rel * 100.0, lw=2, label=r"$\Delta H/H_{\Lambda{\rm CDM}}$")
plt.axhline(0.0, color="black", ls="--", lw=1)
plt.xlabel("redshift z")
plt.ylabel("Residual [%]")
plt.title(r"MFToE vs. $\Lambda$CDM: Hubble-rate residuals")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT_PREFIX}_Hz_residuals.png", dpi=220)

# -----------------------------
# Plot 2: distance modulus residuals
# -----------------------------
plt.figure(figsize=(8, 5))
plt.plot(z, dmu, lw=2, label=r"$\Delta \mu(z)$")
plt.axhline(0.0, color="black", ls="--", lw=1)

# optional SN points
try:
    sn = load_snia_csv(SN_DATA)
    z_sn = sn["z"].values
    mu_sn = sn["mu"].values
    mu_err = sn["mu_err"].values

    # compare SN to ΛCDM baseline for visualization
    mu_lcdm_sn = np.interp(z_sn, z, mu_lcdm)
    sn_resid = mu_sn - mu_lcdm_sn

    plt.errorbar(
        z_sn, sn_resid, yerr=mu_err,
        fmt=".", alpha=0.15, markersize=3, lw=0.5,
        label="Pantheon+-style data"
    )
except Exception as e:
    print(f"[warn] Could not overlay SN points: {e}")

plt.xlabel("redshift z")
plt.ylabel(r"$\Delta \mu$ [mag]")
plt.title(r"MFToE vs. $\Lambda$CDM: distance-modulus residuals")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT_PREFIX}_mu_residuals.png", dpi=220)

# -----------------------------
# Plot 3: dark-energy dynamics proxy
# -----------------------------
plt.figure(figsize=(8, 5))
plt.plot(z, w_phi, lw=2, label=r"$w_\phi(z)$")
plt.plot(z, w_tot, lw=2, label=r"$w_{\rm tot}(z)$")
plt.axhline(-1.0, color="black", ls="--", lw=1, label=r"$w=-1$")
plt.xlabel("redshift z")
plt.ylabel("equation of state")
plt.title(r"MFToE background dynamics")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT_PREFIX}_wphi_wtot.png", dpi=220)

print("Saved:")
print(f"  {OUT_PREFIX}_Hz_residuals.png")
print(f"  {OUT_PREFIX}_mu_residuals.png")
print(f"  {OUT_PREFIX}_wphi_wtot.png")