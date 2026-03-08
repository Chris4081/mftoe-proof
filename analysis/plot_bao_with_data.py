# analysis/plot_bao_with_data.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

C_KMS = 299_792.458

# ---------- config ----------
MFTOE_RUN = "runs/mftoe_baseline.csv"
OUT_DM = "runs/bao_DMrd_with_data.png"
OUT_DH = "runs/bao_DHrd_with_data.png"

H0 = 70.0
Om0 = 0.3
rd = 150.754

# ---------- helpers ----------
def cumtrapz(x, y):
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    out = np.zeros_like(x)
    if len(x) > 1:
        dx = np.diff(x)
        mid = 0.5 * (y[1:] + y[:-1])
        out[1:] = np.cumsum(dx * mid)
    return out

# ---------- load MFToE ----------
df = pd.read_csv(MFTOE_RUN)
z = df["z"].values
H_mf = df["H_over_H0"].values

# ---------- LCDM ----------
cosmo = FlatLambdaCDM(H0=H0 * u.km / u.s / u.Mpc, Om0=Om0)
H_l = cosmo.H(z).value / H0

DH_mf = C_KMS / (H_mf * H0)
DH_l = C_KMS / (H_l * H0)

DM_mf = cumtrapz(z, C_KMS / (H_mf * H0))
DM_l = cumtrapz(z, C_KMS / (H_l * H0))

DMrd_mf = DM_mf / rd
DMrd_l  = DM_l / rd
DHrd_mf = DH_mf / rd
DHrd_l  = DH_l / rd

# ---------- load real DESI DR2 BAO ----------
bao = pd.read_csv("data/desi_dr2/bao_summary.csv")
bao.columns = bao.columns.str.strip()

print("Columns in BAO file:", bao.columns)
print(bao[["tracer","z"]])
# numeric conversion
bao["z"] = pd.to_numeric(bao["z"], errors="coerce")
bao["DM_over_rd"] = pd.to_numeric(bao["DM_over_rd"], errors="coerce")
bao["DM_err"] = pd.to_numeric(bao["DM_err"], errors="coerce")
bao["DH_over_rd"] = pd.to_numeric(bao["DH_over_rd"], errors="coerce")
bao["DH_err"] = pd.to_numeric(bao["DH_err"], errors="coerce")



# only rows with actual measurements
bao_dm = bao.dropna(subset=["DM_over_rd", "DM_err"])
bao_dh = bao.dropna(subset=["DH_over_rd", "DH_err"])

# ---------- DM/rd ----------
plt.figure(figsize=(8, 5))
plt.plot(z, DMrd_mf, lw=2.5, label="MFToE baseline")
plt.plot(z, DMrd_l, "--", lw=2.0, label=r"$\Lambda$CDM")
plt.errorbar(
    bao_dm["z"],
    bao_dm["DM_over_rd"],
    yerr=bao_dm["DM_err"],
    fmt="o",
    capsize=3,
    alpha=0.9,
    label="DESI DR2 BAO"
)
plt.xlabel("redshift z")
plt.ylabel(r"$D_M/r_d$")
plt.title(r"Transverse BAO distance")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(OUT_DM, dpi=220)

# ---------- DH/rd ----------
plt.figure(figsize=(8, 5))
plt.plot(z, DHrd_mf, lw=2.5, label="MFToE baseline")
plt.plot(z, DHrd_l, "--", lw=2.0, label=r"$\Lambda$CDM")
plt.errorbar(
    bao_dh["z"],
    bao_dh["DH_over_rd"],
    yerr=bao_dh["DH_err"],
    fmt="o",
    capsize=3,
    alpha=0.9,
    label="DESI DR2 BAO"
)
plt.xlabel("redshift z")
plt.ylabel(r"$D_H/r_d$")
plt.title(r"Radial BAO distance")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(OUT_DH, dpi=220)

print(f"Saved: {OUT_DM}")
print(f"Saved: {OUT_DH}")