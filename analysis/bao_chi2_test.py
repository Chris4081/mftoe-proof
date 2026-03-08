#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

C_KMS = 299_792.458

# -----------------------------
# Config
# -----------------------------
MFTOE_RUN = "runs/mftoe_baseline.csv"
BAO_SUMMARY = "data/desi_dr2/bao_summary.csv"
BAO_COV = "data/desi_dr2/bao_cov.npy"

H0 = 70.0
Om0 = 0.3
rd = 150.754

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



def interp_model(z_query, z_grid, y_grid):
    return np.interp(z_query, z_grid, y_grid)

# -----------------------------
# Load MFToE model
# -----------------------------
mf = pd.read_csv(MFTOE_RUN)
z = mf["z"].values
H_mf = mf["H_over_H0"].values

# model distances
DH_mf = C_KMS / (H_mf * H0)
DM_mf = cumtrapz(z, C_KMS / (H_mf * H0))

# isotropic BAO distance
DV_mf = ((z * DM_mf**2 * DH_mf)) ** (1.0 / 3.0)

DHrd_mf = DH_mf / rd
DMrd_mf = DM_mf / rd
DVrd_mf = DV_mf / rd

# -----------------------------
# Load LCDM reference
# -----------------------------
cosmo = FlatLambdaCDM(H0=H0 * u.km / u.s / u.Mpc, Om0=Om0)
H_l = cosmo.H(z).value / H0

DH_l = C_KMS / (H_l * H0)
DM_l = cumtrapz(z, C_KMS / (H_l * H0))

# isotropic BAO distance for LCDM
DV_l = ((z * DM_l**2 * DH_l)) ** (1.0 / 3.0)

DHrd_l = DH_l / rd
DMrd_l = DM_l / rd
DVrd_l = DV_l / rd

# -----------------------------
# Load DESI BAO summary
# -----------------------------
bao = pd.read_csv(BAO_SUMMARY)
bao.columns = bao.columns.str.strip()

bao["z"] = pd.to_numeric(bao["z"], errors="coerce")
bao["DM_over_rd"] = pd.to_numeric(bao["DM_over_rd"], errors="coerce")
bao["DM_err"] = pd.to_numeric(bao["DM_err"], errors="coerce")
bao["DH_over_rd"] = pd.to_numeric(bao["DH_over_rd"], errors="coerce")
bao["DH_err"] = pd.to_numeric(bao["DH_err"], errors="coerce")
bao["DV_over_rd"] = pd.to_numeric(bao["DV_over_rd"], errors="coerce")
bao["DV_err"] = pd.to_numeric(bao["DV_err"], errors="coerce")

# keep only rows with at least one BAO observable
bao = bao.dropna(subset=["DM_over_rd", "DH_over_rd", "DV_over_rd"], how="all").reset_index(drop=True)

# -----------------------------
# Build data vector in DESI order
# Convention here:
# [DM(z1), DH(z1), DM(z2), DH(z2), ...]
# but only include whichever are present
# -----------------------------
z_data = []
data_vec = []
labels = []

for _, row in bao.iterrows():
    zi = row["z"]

    if not np.isnan(row["DM_over_rd"]):
        z_data.append(zi)
        data_vec.append(row["DM_over_rd"])
        labels.append((zi, "DM"))

    if not np.isnan(row["DH_over_rd"]):
        z_data.append(zi)
        data_vec.append(row["DH_over_rd"])
        labels.append((zi, "DH"))

    if not np.isnan(row["DV_over_rd"]):
        z_data.append(zi)
        data_vec.append(row["DV_over_rd"])
        labels.append((zi, "DV"))

data_vec = np.array(data_vec, float)

# -----------------------------
# Build MFToE / LCDM prediction vectors
# -----------------------------
pred_mf = []
pred_l = []

for zi, kind in labels:
    if kind == "DM":
        pred_mf.append(interp_model(zi, z, DMrd_mf))
        pred_l.append(interp_model(zi, z, DMrd_l))
    elif kind == "DH":
        pred_mf.append(interp_model(zi, z, DHrd_mf))
        pred_l.append(interp_model(zi, z, DHrd_l))
    elif kind == "DV":
        pred_mf.append(interp_model(zi, z, DVrd_mf))
        pred_l.append(interp_model(zi, z, DVrd_l))

pred_mf = np.array(pred_mf, float)
pred_l = np.array(pred_l, float)

# -----------------------------
# Load covariance
# -----------------------------
cov = np.load(BAO_COV)
icov = np.linalg.inv(cov)

# -----------------------------
# Residuals
# -----------------------------
delta_mf = pred_mf - data_vec
delta_l = pred_l - data_vec

chi2_mf = float(delta_mf @ icov @ delta_mf)
chi2_l = float(delta_l @ icov @ delta_l)

# naive dof estimate
n_data = len(data_vec)

print("DESI DR2 BAO comparison")
print("-----------------------")
print("Number of BAO data points =", n_data)
print("MFToE chi2 =", chi2_mf)
print("LCDM  chi2 =", chi2_l)
print("Δchi2 =", chi2_mf - chi2_l)
print("len(data_vec) =", len(data_vec))
print("cov shape =", cov.shape)
print(bao[["tracer", "z"]])


print("\nData vector order:")
for i, (zi, kind) in enumerate(labels):
    print(f"{i:2d}: z={zi:.3f}  {kind}")