# analysis/plot_fs8_with_data.py
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ---------- files ----------
MFTOE_FS8 = "runs/fsigma8.csv"
OUT = "runs/fsigma8_with_data.png"

# ---------- load model ----------
df = pd.read_csv(MFTOE_FS8)
z = df["z"].values
fs8 = df["fsigma8"].values

# ---------- small reference RSD dataset ----------
# illustrative literature-style points for quick paper plotting
# columns: z, fs8, err
rsd = pd.DataFrame([
    (0.02, 0.428, 0.046),   # 6dF-like
    (0.15, 0.490, 0.145),   # SDSS-like
    (0.32, 0.427, 0.056),   # BOSS LOWZ-like
    (0.57, 0.426, 0.029),   # BOSS CMASS-like
    (0.77, 0.490, 0.180),   # VVDS/WiggleZ-like
    (1.05, 0.280, 0.080),   # VIPERS-like
    (1.40, 0.482, 0.116),   # FastSound-like
], columns=["z", "fs8", "err"])

# ---------- plot ----------
plt.figure(figsize=(8,5))
plt.plot(z, fs8, lw=2.5, label=r"MFToE baseline $f\sigma_8(z)$")
plt.errorbar(
    rsd["z"], rsd["fs8"], yerr=rsd["err"],
    fmt="o", capsize=3, alpha=0.9, label="RSD data (illustrative)"
)

plt.xlabel("redshift z")
plt.ylabel(r"$f\sigma_8$")
plt.title(r"MFToE structure growth vs. RSD measurements")
plt.grid(alpha=0.3)
plt.legend()
plt.tight_layout()
plt.savefig(OUT, dpi=220)
print(f"Saved: {OUT}")