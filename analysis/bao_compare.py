#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, json, os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# -----------------------------
# BAO model helpers
# -----------------------------
def model_bao_at_z(z, z_mod, H_over_H0, dL_dimless, H0phys, rd):
    """
    Inputs:
      - z:        scalar redshift where to evaluate
      - z_mod:    array of model redshifts (sorted)
      - H_over_H0: model H(z)/H0 (same length as z_mod)
      - dL_dimless: model luminosity distance in units of c/H0 (same length)
      - H0phys:   physical H0 in km/s/Mpc (for unit conversions)
      - rd:       sound horizon r_d in Mpc

    Returns:
      (DM/rd, DH/rd, DV/rd) at redshift z
    """
    # Interpolationen auf genaues z
    Hn = np.interp(z, z_mod, H_over_H0)   # H/H0 (dimensionless)
    dL = np.interp(z, z_mod, dL_dimless)  # dimensionless, in units of c/H0

    # D_M = d_L / (1+z), aber d_L ist in (c/H0)-Einheiten
    DM_dimless = dL / (1.0 + z)                      # still in c/H0
    DH_dimless = 1.0 / max(Hn, 1e-30)                # D_H = c/H = (c/H0)/(H/H0)

    # In Mpc: multiply by (c/H0phys), dann / r_d → /rd
    # Faktor (c/H0phys) kürzt sich weg wenn beide Seiten gleich behandelt werden,
    # aber für DV ist die Mischformel praktisch:
    DM_over_rd = DM_dimless * (299792.458 / H0phys) / rd
    DH_over_rd = DH_dimless * (299792.458 / H0phys) / rd

    # Isotrope Kombination:
    # D_V = [ z * D_M^2 * D_H ]^{1/3}
    DV_over_rd = ( z * (DM_over_rd**2) * DH_over_rd ) ** (1.0/3.0)

    return DM_over_rd, DH_over_rd, DV_over_rd


def build_diag_cov_from_df(df):
    """Diagonal-Kovarianz aus den Fehlern der gestackten Zeilen (DM/DH/DV)."""
    errs2 = (df["err"].values.astype(float) ** 2)
    return np.diag(errs2)


# -----------------------------
# Main
# -----------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-csv", required=True, help="MFToE output CSV (z,H_over_H0,dL,...)")
    ap.add_argument("--bao-csv",   required=True, help="DESI BAO summary CSV")
    ap.add_argument("--H0phys", type=float, required=True, help="Physical H0 [km/s/Mpc]")
    ap.add_argument("--rd",     type=float, required=True, help="Sound horizon r_d [Mpc]")
    ap.add_argument("--cov", default="", help="Optional .npy covariance matrix (stacked DM,DH,DV rows order)")
    ap.add_argument("--show-residuals", action="store_true")
    ap.add_argument("--export-json",    action="store_true")
    ap.add_argument("--out", default="runs/bao_check", help="Output prefix (PNG/JSON)")
    args = ap.parse_args()

    # --- CSVs laden ---
    mod = pd.read_csv(args.model_csv)  # erwartet Spalten: z,H_over_H0,dL,...
    bao = pd.read_csv(args.bao_csv)    # erwartet Spalten: tracer,z, DM_over_rd,DM_err, DH_over_rd,DH_err, DV_over_rd,DV_err

    # --- gestackte Vergleichstabelle bauen ---
    rows = []
    for _, r in bao.iterrows():
        z = float(r["z"])
        DMm, DHm, DVm = model_bao_at_z(
            z,
            mod["z"].values, mod["H_over_H0"].values, mod["dL"].values,
            args.H0phys, args.rd
        )
        for name, mval, dcol, ecol in [
            ("DM_over_rd", DMm, "DM_over_rd", "DM_err"),
            ("DH_over_rd", DHm, "DH_over_rd", "DH_err"),
            ("DV_over_rd", DVm, "DV_over_rd", "DV_err"),
        ]:
            dv = r.get(dcol)
            ev = r.get(ecol)
            if pd.notna(dv) and pd.notna(ev):
                rows.append(dict(
                    tracer=r.get("tracer",""), z=z, obs=name,
                    data=float(dv), err=float(ev), model=float(mval)
                ))

    df = pd.DataFrame(rows)
    if df.empty:
        raise RuntimeError("No comparable BAO rows found. Check your bao_summary.csv columns.")

    df["resid"] = df["data"] - df["model"]
    df["pull"]  = df["resid"] / df["err"]

    # --- Kovarianz laden oder Diagonal-Fallback ---
    if args.cov:
        try:
            cov = np.load(args.cov)
        except FileNotFoundError:
            print(f"⚠️  Cov file not found: {args.cov} — falling back to diagonal covariance from stacked rows.")
            cov = build_diag_cov_from_df(df)
    else:
        cov = build_diag_cov_from_df(df)

    n = len(df)
    if cov.shape != (n, n):
        print(f"⚠️  Cov shape {cov.shape} mismatches N={n} — using diagonal covariance.")
        cov = build_diag_cov_from_df(df)

    # --- χ² ---
    r = df["resid"].values.astype(float)
    try:
        icov = np.linalg.inv(cov)
        chi2 = float(r.T @ icov @ r)
    except np.linalg.LinAlgError:
        print("⚠️  Cov not invertible — using diagonal fallback.")
        cov = build_diag_cov_from_df(df)
        chi2 = float(np.sum((r / df["err"].values) ** 2))

    ndof = n  # ohne Parameterabzug (reiner Konsistenz-Check)
    red  = chi2 / max(ndof, 1)
    print(f"χ² = {chi2:.3f}  (reduced χ² = {red:.3f}, N={ndof})")

    # --- Plot ---
    fig, axs = plt.subplots(2 if args.show_residuals else 1, 3,
                            figsize=(15, 10 if args.show_residuals else 6))
    axs = np.atleast_2d(axs)

    for j, (name, title) in enumerate([
        ("DM_over_rd", "D_M / r_d"),
        ("DH_over_rd", "D_H / r_d"),
        ("DV_over_rd", "D_V / r_d"),
    ]):
        sub = df[df["obs"] == name]
        axs[0, j].set_title(title); axs[0, j].set_xlabel("z"); axs[0, j].grid(True, alpha=0.3)
        if len(sub) == 0:
            axs[0, j].text(0.5, 0.5, "(no data)", ha="center", va="center")
            continue
        zvals = sub["z"].values
        axs[0, j].errorbar(zvals, sub["data"].values, yerr=sub["err"].values,
                           fmt="o", capsize=3, label="DESI DR2")
        axs[0, j].plot(zvals, sub["model"].values, "-", label="MFToE")
        if j == 0:
            axs[0, j].legend()

        if args.show_residuals:
            axs[1, j].axhline(0, color="k", lw=1)
            axs[1, j].errorbar(zvals, sub["pull"].values, yerr=np.ones_like(zvals),
                               fmt="o", capsize=3)
            axs[1, j].set_title(f"{title} residuals (pull)")
            axs[1, j].set_xlabel("z"); axs[1, j].set_ylabel("(data - model)/σ")
            axs[1, j].grid(True, alpha=0.3)

    plt.tight_layout()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    out_png = f"{args.out}.png"
    plt.savefig(out_png, dpi=160)
    print(f"→ Plot saved: {out_png}")

    if args.export_json:
        out_json = f"{args.out}.json"
        with open(out_json, "w") as f:
            json.dump(dict(
                chi2=chi2, reduced_chi2=red, N=int(ndof),
                rows=df.to_dict(orient="records"),
            ), f, indent=2)
        print(f"→ JSON saved: {out_json}")


if __name__ == "__main__":
    main()