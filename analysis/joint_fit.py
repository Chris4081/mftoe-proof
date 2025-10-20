#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Joint chi^2: BAO (compressed) + optional SNIa + optional GW + optional CMB r_d prior
Uses a model CSV with columns: z, H_over_H0, dL (d_L in dimensionless c/H0 units).

Example:
  python3 analysis/joint_fit.py \
    --model-csv runs/mftoe_vacuum_astropy.csv \
    --bao-csv   data/desi_dr2/bao_summary.csv \
    --H0phys 67.36 --rd 150.754 \
    --cmb-rd 150.74 --cmb-rd-sigma 0.30 \
    --snia-csv  data/snia/pantheonplus_summary.csv \
    --gw-csv    data/gw/bright_sirens.csv \
    --out runs/joint_baseline --save-json --plot-resids
"""

import argparse, json, sys, pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# --- robust path bootstrap (run as script OR as module) ---
HERE = pathlib.Path(__file__).resolve().parent          # .../analysis
ROOT = HERE.parent                                      # .../ (repo root)
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

# Optional rd_engine import (only needed if rd-backend=camb)
try:
    from analysis.rd_engine import EarlyParams, rd_camb
except Exception:
    try:
        from rd_engine import EarlyParams, rd_camb
    except Exception:
        EarlyParams = None
        rd_camb = None

# ---------- Constants ----------
C_KMS = 299_792.458  # km/s

# ---------- Helpers ----------
def chi2_rd_prior(rd_used, rd_mean, rd_sigma):
    """Gaussian prior on r_d (e.g. from CMB)."""
    if rd_mean is None or rd_sigma is None:
        return 0.0
    if rd_sigma == 0.0:
        return 0.0 if rd_used == rd_mean else np.inf
    return ((rd_used - rd_mean) / rd_sigma) ** 2

def load_snia_csv(csv_path):
    """Load SNIa CSV: expects columns z, mu, mu_err (distance modulus)."""
    df = pd.read_csv(csv_path)
    need = {"z", "mu", "mu_err"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"SNIa CSV missing columns: {missing}")
    return df.sort_values("z")

def chi2_snia_profileM(z_snia, mu_snia, mu_err_snia, z_mod, dL_dimless, H0phys,
                       cov_path=None, sigma_int=0.0, vpec_kms=0.0):
    """
    χ² for SNe Ia with analytic profiling of nuisance 𝓜.
    μ_mod(z) = 5 log10(d_L/Mpc) + 25;  d_L(Mpc) = dL_dimless * (c/H0phys)
    """
    dL_interp = np.interp(z_snia, z_mod, dL_dimless)
    dL_Mpc = dL_interp * (C_KMS / H0phys)
    mu_mod = 5.0 * np.log10(np.clip(dL_Mpc, 1e-6, None)) + 25.0

    err = np.array(mu_err_snia, float)

    # Peculiar-velocity error (mag): σ_μ,pec = (5/ln10) * σ_v / (c z)
    if vpec_kms and np.any(z_snia > 0):
        sigma_mu_pec = (5.0 / np.log(10.0)) * (vpec_kms / (C_KMS * np.clip(z_snia, 1e-4, None)))
        err = np.sqrt(err**2 + sigma_mu_pec**2)

    # Intrinsic scatter (mag)
    if sigma_int and sigma_int > 0:
        err = np.sqrt(err**2 + sigma_int**2)

    w = 1.0 / np.clip(err, 1e-12, None)**2
    resid0 = mu_snia - mu_mod
    M_star = float(np.sum(w * resid0) / np.sum(w))
    resid = resid0 - M_star

    if cov_path:
        C = np.load(cov_path)
        invC = np.linalg.inv(C)
        chi2 = float(resid @ invC @ resid)
    else:
        chi2 = float(np.sum((resid / np.clip(err, 1e-12, None))**2))

    return chi2, M_star, err

def load_bright_sirens(csv_path):
    """Load GW CSV: expects columns z, dL_Mpc, dL_err_Mpc."""
    df = pd.read_csv(csv_path)
    need = {"z", "dL_Mpc", "dL_err_Mpc"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"GW CSV missing columns: {missing}")
    return df.sort_values("z")

def chi2_gw(z_gw, dL_gw_Mpc, dL_gw_err_Mpc, z_mod, dL_dimless, H0phys, cov_path=None):
    """χ² for GW bright sirens."""
    dL_interp = np.interp(z_gw, z_mod, dL_dimless)
    dL_mod_Mpc = dL_interp * (C_KMS / H0phys)
    resid = dL_gw_Mpc - dL_mod_Mpc

    if cov_path:
        C = np.load(cov_path)
        if C.shape != (len(resid), len(resid)):
            raise ValueError(f"GW covariance shape mismatch: {C.shape} vs {len(resid)}")
        invC = np.linalg.inv(C)
        chi2 = float(resid @ invC @ resid)
    else:
        chi2 = float(np.sum((resid / np.clip(dL_gw_err_Mpc, 1e-12, None))**2))
    return chi2

def load_model_csv(path):
    df = pd.read_csv(path)
    for col in ["z", "H_over_H0", "dL"]:
        if col not in df.columns:
            raise ValueError(f"Model CSV missing column '{col}'")
    return df.sort_values("z").reset_index(drop=True)

def model_bao_at_z(zq, z_mod, H_over_H0, dL_dimless, H0phys, rd):
    """
    Returns (DM/rd, DH/rd, DV/rd) at z=zq.
    dL_dimless is d_L in (c/H0) units; H0phys in km/s/Mpc; rd in Mpc.
    """
    if not np.isfinite(rd) or rd <= 0:
        raise ValueError(f"Invalid r_d in BAO projection: {rd}")
    Hn     = np.interp(zq, z_mod, H_over_H0)
    dL_dim = np.interp(zq, z_mod, dL_dimless)

    dL_Mpc = dL_dim * (C_KMS / H0phys)
    DM_Mpc = dL_Mpc / (1.0 + zq)
    H_km   = Hn * H0phys
    DH_Mpc = C_KMS / np.clip(H_km, 1e-30, None)
    DV_Mpc = (zq * DM_Mpc**2 * DH_Mpc) ** (1.0/3.0)

    return DM_Mpc/rd, DH_Mpc/rd, DV_Mpc/rd

def chi2_bao(bao_csv, model_df, H0phys, rd, cov_path=None, export_vector=None):
    bao = pd.read_csv(bao_csv)
    rows = []
    for _, r in bao.iterrows():
        z = float(r["z"])
        DMm, DHm, DVm = model_bao_at_z(
            z, model_df["z"].values, model_df["H_over_H0"].values,
            model_df["dL"].values, H0phys, rd
        )
        for name, mval, dcol, ecol in [
            ("DM_over_rd", DMm, "DM_over_rd", "DM_err"),
            ("DH_over_rd", DHm, "DH_over_rd", "DH_err"),
            ("DV_over_rd", DVm, "DV_over_rd", "DV_err"),
        ]:
            dval = r.get(dcol, np.nan)
            eval_ = r.get(ecol, np.nan)
            if pd.notna(dval) and pd.notna(eval_):
                rows.append(dict(z=z, obs=name, data=float(dval), err=float(eval_), model=float(mval)))
    df = pd.DataFrame(rows)
    if len(df) == 0:
        return 0.0, 0, df

    resid = df["data"].values - df["model"].values

    if cov_path:
        C = np.load(cov_path)
        if C.shape != (len(df), len(df)):
            raise ValueError(f"BAO covariance shape {C.shape} vs vector length {len(df)}")
        invC = np.linalg.inv(C)
        chi2 = float(resid @ invC @ resid)
    else:
        chi2 = float(np.sum((resid / np.clip(df["err"].values, 1e-12, None))**2))

    if export_vector is not None:
        export_vector["bao_vector"] = df.to_dict(orient="records")
    return chi2, len(df), df

def build_bao_dataframe_for_residuals(bao_csv_df, model_df, H0phys, rd):
    rows = []
    z_mod = model_df["z"].values
    Hn    = model_df["H_over_H0"].values
    dLdim = model_df["dL"].values

    for _, r in bao_csv_df.iterrows():
        z = float(r["z"])
        DMm, DHm, DVm = model_bao_at_z(z, z_mod, Hn, dLdim, H0phys, rd)
        for name, mval, dcol, ecol in [
            ("DM_over_rd", DMm, "DM_over_rd", "DM_err"),
            ("DH_over_rd", DHm, "DH_over_rd", "DH_err"),
            ("DV_over_rd", DVm, "DV_over_rd", "DV_err"),
        ]:
            dval = r.get(dcol, np.nan)
            eval_ = r.get(ecol, np.nan)
            if pd.notna(dval) and pd.notna(eval_):
                resid = float(dval) - float(mval)
                pull  = resid / float(eval_) if float(eval_) > 0 else np.nan
                rows.append(dict(
                    tracer=r.get("tracer",""),
                    z=z, obs=name, data=float(dval), err=float(eval_),
                    model=float(mval), resid=resid, pull=pull
                ))
    return pd.DataFrame(rows)

def plot_bao_with_residuals(df_long, out_base, title_suffix=""):
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    panels = [("DM_over_rd", r"$D_M/r_d$"),
              ("DH_over_rd", r"$D_H/r_d$"),
              ("DV_over_rd", r"$D_V/r_d$")]

    for j, (name, title) in enumerate(panels):
        sub = df_long[df_long["obs"] == name]
        ax_top = axs[0, j]; ax_bot = axs[1, j]

        if len(sub) == 0:
            ax_top.text(0.5, 0.5, "(no data)", ha="center", va="center")
            ax_top.set_title(title); ax_top.set_xlabel("z")
            ax_top.grid(True, alpha=0.3)
            ax_bot.axis("off")
            continue

        z = sub["z"].values
        ax_top.errorbar(z, sub["data"].values, yerr=sub["err"].values,
                        fmt="o", capsize=3, label="DESI DR2")
        ax_top.plot(z, sub["model"].values, "-", label="MFToE")
        ax_top.set_title(f"{title} {title_suffix}".strip())
        ax_top.set_xlabel("z"); ax_top.grid(True, alpha=0.3)
        if j == 0:
            ax_top.legend()

        ax_bot.axhline(0.0, lw=1)
        ax_bot.errorbar(z, sub["pull"].values, yerr=np.ones_like(z),
                        fmt="o", capsize=3)
        ax_bot.set_xlabel("z"); ax_bot.set_ylabel("pull")
        ax_bot.set_title(f"{title} residuals")
        ax_bot.grid(True, alpha=0.3)

    plt.tight_layout()
    out_png = f"{out_base}_bao_residuals.png"
    plt.savefig(out_png, dpi=160)
    plt.close(fig)
    print(f"→ Residual plot saved: {out_png}")

# ---------- CLI ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model-csv", required=True)
    ap.add_argument("--bao-csv", required=True)
    ap.add_argument("--cov", default="", help="optional BAO covariance .npy")
    ap.add_argument("--snia-csv", default="")
    ap.add_argument("--snia-cov", default="", help="optional SNIa covariance .npy")
    ap.add_argument("--gw-csv", default="")
    ap.add_argument("--gw-cov", default="", help="optional GW covariance .npy")

    # CMB r_d prior (optional)
    ap.add_argument("--cmb-rd", type=float, default=None, help="mean of r_d prior [Mpc]")
    ap.add_argument("--cmb-rd-sigma", type=float, default=None, help="sigma of r_d prior [Mpc]")
    ap.add_argument("--no-cmb-prior", action="store_true",
                    help="Disable CMB r_d Gaussian prior entirely")

    ap.add_argument("--H0phys", type=float, required=True, help="physical H0 [km/s/Mpc]")
    ap.add_argument("--rd", type=float, default=None,
                    help="sound horizon r_d [Mpc] for BAO compression (required if --rd-backend=fixed)")

    # r_d backend (fixed vs. CAMB)
    ap.add_argument("--rd-backend", choices=["fixed", "camb"], default="fixed",
                    help="Use fixed r_d (default) or compute from CAMB Boltzmann code")
    ap.add_argument("--ombh2", type=float, default=0.02237, help="Physical baryon density Ω_b h²")
    ap.add_argument("--omch2", type=float, default=0.1200, help="Physical cold dark matter density Ω_c h²")
    ap.add_argument("--Neff",  type=float, default=3.046,  help="Effective number of neutrino species")
    ap.add_argument("--Yp",    type=float, default=0.245,  help="Primordial helium fraction")
    ap.add_argument("--mnu-eV", type=float, default=0.06,  help="Total neutrino mass [eV]")

    # Option: keep H0*rd constant via rescaling of H0phys
    ap.add_argument("--match-H0rd", action="store_true",
                    help="Rescale H0phys to keep H0*rd constant w.r.t. a reference pair")
    ap.add_argument("--ref-H0", type=float, default=67.36,
                    help="Reference H0 [km/s/Mpc] for H0*rd matching (default: DESI 67.36)")
    ap.add_argument("--ref-rd", type=float, default=150.754,
                    help="Reference r_d [Mpc] for H0*rd matching (default: DESI 150.754)")

    ap.add_argument("--out", default="runs/joint_out")
    ap.add_argument("--save-json", action="store_true")
    ap.add_argument("--plot-resids", action="store_true", help="plot BAO residuals")
    ap.add_argument("--n-params", type=int, default=5,
                    help="Number of free model parameters for reduced chi^2")

    # SN robust errors
    ap.add_argument("--snia-sigma-int", type=float, default=0.0, help="intrinsic SN scatter in mag (added in quadrature)")
    ap.add_argument("--snia-vpec",      type=float, default=0.0, help="peculiar-velocity scatter in km/s (adds to mag error)")

    args = ap.parse_args()

    # --- determine r_d to use (fixed or CAMB) ---
    if args.rd_backend == "camb":
        if rd_camb is None or EarlyParams is None:
            raise RuntimeError("CAMB backend requested but rd_engine not available/importable.")
        ep = EarlyParams(
            H0=args.H0phys,
            ombh2=args.ombh2,
            omch2=args.omch2,
            Neff=args.Neff,
            Yp=args.Yp,
            mnu_eV=args.mnu_eV
        )
        rd_used = float(rd_camb(ep))
        print(f"→ Using r_d = {rd_used:.3f} Mpc  [camb]")
    else:
        if args.rd is None:
            raise ValueError("Provide --rd when using --rd-backend fixed.")
        rd_used = float(args.rd)
        print(f"→ Using r_d = {rd_used:.3f} Mpc  [fixed]")

    # --- optional: rescale H0 to keep H0*rd ≈ const w.r.t. (ref_H0, ref_rd) ---
    H0phys_eff = float(args.H0phys)
    if args.match_H0rd:
        if not np.isfinite(rd_used) or rd_used <= 0:
            raise ValueError("match-H0rd requested but r_d is not positive.")
        scale = (args.ref_H0 * args.ref_rd) / (args.H0phys * rd_used)
        H0phys_eff = args.H0phys * scale
        print(f"→ match-H0rd: using H0_eff = {H0phys_eff:.3f} km/s/Mpc "
              f"(scale={scale:.5f}) to keep H0*rd ≈ const vs. "
              f"(H0_ref={args.ref_H0}, rd_ref={args.ref_rd})")
    else:
        print(f"→ no match-H0rd: using H0phys = {H0phys_eff:.3f} km/s/Mpc")

    if not np.isfinite(rd_used) or rd_used <= 0:
        raise ValueError(f"Invalid r_d: {rd_used}. Provide --rd > 0 or use --rd-backend camb.")

    # Load model AFTER effective H0/rd are decided (for clarity)
    model = load_model_csv(args.model_csv)

    report = {}

    # BAO
    bao_export = {}
    chi2_b, n_b, bao_df = chi2_bao(
        args.bao_csv, model, H0phys_eff, rd_used,
        cov_path=(args.cov or None), export_vector=bao_export
    )
    report["chi2_bao"] = chi2_b
    report["N_bao"] = n_b

    # SN Ia
    chi2_sn, n_sn = 0.0, 0
    sn_Mstar = None
    if args.snia_csv:
        sn = load_snia_csv(args.snia_csv)
        result = chi2_snia_profileM(
            sn["z"].values, sn["mu"].values, sn["mu_err"].values,
            model["z"].values, model["dL"].values, H0phys_eff,
            cov_path=(args.snia_cov or None),
            sigma_int=args.snia_sigma_int,
            vpec_kms=args.snia_vpec
        )
        chi2_sn = result[0]
        if len(result) > 1:
            sn_Mstar = result[1]
        n_sn = len(sn)
    report["chi2_snia"] = chi2_sn
    report["N_snia"] = n_sn
    if sn_Mstar is not None:
        report["SN_Mstar_profiled_mag"] = sn_Mstar

    # GW bright sirens
    chi2_gw_val, n_gw = 0.0, 0
    if args.gw_csv:
        gw = load_bright_sirens(args.gw_csv)
        chi2_gw_val = chi2_gw(
            gw["z"].values, gw["dL_Mpc"].values, gw["dL_err_Mpc"].values,
            model["z"].values, model["dL"].values, H0phys_eff,
            cov_path=(args.gw_cov or None)
        )
        n_gw = len(gw)
    report["chi2_gw"] = chi2_gw_val
    report["N_gw"] = n_gw

    # CMB r_d prior
    chi2_cmb = 0.0
    use_cmb_prior = (not args.no_cmb_prior) and \
                    (args.cmb_rd is not None and args.cmb_rd_sigma is not None)
    if use_cmb_prior:
        chi2_cmb = chi2_rd_prior(rd_used, args.cmb_rd, args.cmb_rd_sigma)
    report["chi2_cmb_prior"] = chi2_cmb
    report["use_cmb_prior"] = bool(use_cmb_prior)

    # Totals & reduced chi^2
    chi2_tot = chi2_b + chi2_sn + chi2_gw_val + chi2_cmb
    N_tot = n_b + n_sn + n_gw + (1 if use_cmb_prior else 0)
    nuisance = 1 if args.snia_csv else 0  # profiled 𝓜
    dof = max(N_tot - args.n_params - nuisance, 1)
    red_chi2 = chi2_tot / dof

    report["chi2_total"] = chi2_tot
    report["N_total"] = N_tot
    report["reduced_chi2"] = red_chi2
    report["dof"] = dof
    report["H0phys_input"] = args.H0phys
    report["H0phys_eff"] = H0phys_eff
    report["rd_used"] = rd_used
    report["inputs"] = {
        "model_csv": args.model_csv,
        "bao_csv": args.bao_csv,
        "bao_cov": (args.cov or None),
        "snia_csv": (args.snia_csv or None),
        "snia_cov": (args.snia_cov or None),
        "gw_csv": (args.gw_csv or None),
        "gw_cov": (args.gw_cov or None),
        "cmb_rd": (args.cmb_rd if args.cmb_rd is not None else None),
        "cmb_rd_sigma": (args.cmb_rd_sigma if args.cmb_rd_sigma is not None else None),
        "rd_backend": args.rd_backend,
        "ombh2": args.ombh2,
        "omch2": args.omch2,
        "Neff": args.Neff,
        "Yp": args.Yp,
        "mnu_eV": args.mnu_eV,
        "match_H0rd": bool(args.match_H0rd),
        "ref_H0": args.ref_H0,
        "ref_rd": args.ref_rd,
        "n_params": args.n_params,
        "snia_sigma_int": args.snia_sigma_int,
        "snia_vpec": args.snia_vpec,
        "no_cmb_prior": bool(args.no_cmb_prior),
    }

    _parts = [
        f"BAO: {chi2_b:.3f} (N={n_b})",
        f"SN: {chi2_sn:.3f} (N={n_sn})",
        f"GW: {chi2_gw_val:.3f} (N={n_gw})",
    ]
    if use_cmb_prior:
        _parts.append(f"CMB r_d: {chi2_cmb:.3f}")
    print("χ² components → " + " | ".join(_parts))
    print(f"χ² total = {chi2_tot:.3f}   (N={N_tot}, dof={dof}, reduced χ² = {red_chi2:.3f})")

    # Save JSON
    if args.save_json:
        out_json = Path(f"{args.out}.json")
        payload = {"report": report, "bao_vector": (bao_export.get("bao_vector", []))}
        out_json.parent.mkdir(parents=True, exist_ok=True)
        with open(out_json, "w") as f:
            json.dump(payload, f, indent=2)
        print(f"→ JSON saved: {out_json}")

    # Optional BAO residual plots
    if args.plot_resids:
        bao_long = build_bao_dataframe_for_residuals(pd.read_csv(args.bao_csv), model, H0phys_eff, rd_used)
        plot_bao_with_residuals(bao_long, args.out, title_suffix="(Joint Fit)")

if __name__ == "__main__":
    main()
