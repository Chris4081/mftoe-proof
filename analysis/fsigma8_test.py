#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp

# -----------------------------
# Config
# -----------------------------
MFTOE_RUN = "runs/mftoe_baseline.csv"
OUT_PREFIX = "runs/fsigma8"

Omega_m0 = 0.3
sigma8_0 = 0.8

# Optional MFToE effective gravity:
USE_MFTOE_MU = False
XI_A = 1.6e-3   # illustrative product from your first fit

# -----------------------------
# Load background
# -----------------------------
df = pd.read_csv(MFTOE_RUN)

z = df["z"].values
H = df["H_over_H0"].values

# sort in increasing scale factor / increasing N = ln a
a = 1.0 / (1.0 + z)
N = np.log(a)

order = np.argsort(N)
N = N[order]
a = a[order]
z = z[order]
H = H[order]

# Interpolators
H_of_N = interp1d(N, H, kind="cubic", fill_value="extrapolate")

# d ln H / dN numerically on the grid
lnH = np.log(H)
dlnH_dN_grid = np.gradient(lnH, N)
dlnH_of_N = interp1d(N, dlnH_dN_grid, kind="cubic", fill_value="extrapolate")

def Omega_m_of_N(Nval):
    aval = np.exp(Nval)
    Hval = H_of_N(Nval)
    return Omega_m0 * aval**(-3) / np.clip(Hval**2, 1e-30, None)

def mu_of_N(Nval):
    """
    Effective gravitational coupling mu(a)=G_eff/G.
    Conservative baseline: mu=1.
    Optional simple MFToE correction:
        mu = 1/(1 + xiA)
    """
    if not USE_MFTOE_MU:
        return 1.0
    return 1.0 / (1.0 + XI_A)

def growth_ode(Nval, y):
    """
    y[0] = D
    y[1] = dD/dN
    """
    D, Dp = y
    Om = Omega_m_of_N(Nval)
    dlnH = dlnH_of_N(Nval)
    mu = mu_of_N(Nval)

    Dpp = -(2.0 + dlnH) * Dp + 1.5 * Om * mu * D
    return [Dp, Dpp]

# -----------------------------
# Initial conditions
# -----------------------------
# Start deep enough in matter domination:
N_ini = N.min()
N_fin = N.max()

# In matter domination: D ~ a => D = exp(N), D' = D
D_ini = np.exp(N_ini)
Dp_ini = D_ini

sol = solve_ivp(
    growth_ode,
    (N_ini, N_fin),
    [D_ini, Dp_ini],
    t_eval=N,
    rtol=1e-8,
    atol=1e-10,
)

if not sol.success:
    raise RuntimeError("Growth integration failed.")

D = sol.y[0]
Dp = sol.y[1]

# Normalize so that D(z=0)=1
D0 = D[-1]
D = D / D0
Dp = Dp / D0

# Growth rate
f = Dp / np.clip(D, 1e-30, None)

# sigma8(z)
sigma8_z = sigma8_0 * D

# f sigma8
fs8 = f * sigma8_z

# -----------------------------
# Save results
# -----------------------------
out_df = pd.DataFrame({
    "z": z,
    "a": a,
    "N": N,
    "H_over_H0": H,
    "D": D,
    "f": f,
    "sigma8_z": sigma8_z,
    "fsigma8": fs8,
})
csv_path = f"{OUT_PREFIX}.csv"
out_df.to_csv(csv_path, index=False)

# -----------------------------
# Plots
# -----------------------------
plt.figure(figsize=(8,5))
plt.plot(z, D, lw=2, label="D(z)")
plt.xlabel("redshift z")
plt.ylabel("growth factor D")
plt.title("MFToE linear growth factor")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT_PREFIX}_D.png", dpi=220)

plt.figure(figsize=(8,5))
plt.plot(z, f, lw=2, label="f(z)")
plt.xlabel("redshift z")
plt.ylabel("growth rate f")
plt.title("MFToE linear growth rate")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT_PREFIX}_f.png", dpi=220)

plt.figure(figsize=(8,5))
plt.plot(z, fs8, lw=2, label=r"$f\sigma_8(z)$")
plt.xlabel("redshift z")
plt.ylabel(r"$f\sigma_8$")
plt.title(r"MFToE structure growth prediction")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(f"{OUT_PREFIX}_fs8.png", dpi=220)

print("Saved:")
print(f"  {csv_path}")
print(f"  {OUT_PREFIX}_D.png")
print(f"  {OUT_PREFIX}_f.png")
print(f"  {OUT_PREFIX}_fs8.png")

# Quick summary
for zq in [0.0, 0.5, 1.0, 2.0]:
    i = np.abs(z - zq).argmin()
    print(f"z={z[i]:.2f} | D={D[i]:.4f} | f={f[i]:.4f} | fσ8={fs8[i]:.4f}")