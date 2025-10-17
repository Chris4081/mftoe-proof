#!/usr/bin/env python3
import sys, numpy as np, pandas as pd

"""
Build a synthetic BAO covariance from a CSV like bao_summary.csv.
- Diagonal: uses the quoted errors^2.
- Off-diagonal: adds a correlation rho within the same redshift bin
  between (DM, DH, DV) if those exist at that z.
"""

def build_cov(bao_csv, rho=0.3):
    df = pd.read_csv(bao_csv)
    rows = []
    # build a flat list of available observations with (z, name, sigma)
    for _, r in df.iterrows():
        z = float(r["z"])
        for name, ecol in [("DM_over_rd","DM_err"),
                           ("DH_over_rd","DH_err"),
                           ("DV_over_rd","DV_err")]:
            if name in r and ecol in r and pd.notna(r[name]) and pd.notna(r[ecol]):
                rows.append((z, name, float(r[ecol])))

    N = len(rows)
    C = np.zeros((N, N), float)

    # fill diagonal
    for i in range(N):
        C[i, i] = rows[i][2]**2

    # add intra-z correlations
    for i in range(N):
        zi, ni, si = rows[i]
        for j in range(i+1, N):
            zj, nj, sj = rows[j]
            if abs(zi - zj) < 1e-6:  # same effective redshift bin
                C[i, j] = C[j, i] = rho * si * sj

    return C, rows

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: make_cov_from_csv.py bao_summary.csv out.npy [rho]")
        sys.exit(1)
    bao_csv, out_npy = sys.argv[1], sys.argv[2]
    rho = float(sys.argv[3]) if len(sys.argv) > 3 else 0.3
    C, rows = build_cov(bao_csv, rho=rho)
    np.save(out_npy, C)
    print(f"→ Saved synthetic covariance to {out_npy}  (shape={C.shape}, rho={rho})")