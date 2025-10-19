#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Joint chi^2: BAO (compressed) + optional SNIa + optional GW + optional CMB r_d prior
Uses your existing model CSV (z, H_over_H0, dL).

Usage example:
  python3 analysis/joint_fit.py \
    --model-csv runs/mftoe_vacuum_astropy.csv \
    --bao-csv   data/desi_dr2/bao_summary.csv \
    --H0phys 67.36 --rd 150.754 \
    --cmb-rd 150.74 --cmb-rd-sigma 0.30 \
    --snia-csv  data/snia/pantheonplus_summary.csv \
    --gw-csv    data/gw/bright_sirens.csv \
    --out runs/joint_baseline
"""
import argparse, json
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# ---------- Inline Helpers (fehlende Modules ersetzt) ----------
C_KMS = 299_792.458  # km/s

def chi2_rd_prior(rd_used, rd_mean, rd_sigma):
    """Simple Gaussian prior on r_d from CMB (e.g., Planck)."""
    if rd_sigma == 0.0:
        return 0.0 if rd_used == rd_mean else np.inf
    return ((rd_used - rd_mean) / rd_sigma)**2

def load_snia_csv(csv_path):
    """Load SNIa CSV: expects columns z, mu, mu_err (distance modulus)."""
    df = pd.read_csv(csv_path)
    need = {"z", "mu", "mu_err"}
    if not need.issubset(df.columns):
        raise ValueError(f"SNIa CSV missing {need - set(df.columns)}")
    return df.sort_values("z")

def chi2_snia_profileM(z_snia, mu_snia, mu_err_snia, z_mod, dL_dimless, H0phys,
                       cov_path=None, sigma_int=0.0, vpec_kms=0.0):
    """
    χ² für SNe mit analytischer Profilierung von 𝓜.
    Optional: intrinsische Streuung sigma_int (mag) und Peculiar-Velocity-Fehler vpec_kms.

    μ_mod(z) = 5 log10(d_L/Mpc) + 25  (aus dL_dimless * c/H0phys)
    """
    C_KMS = 299_792.458

    # Modell-μ
    dL_interp = np.interp(z_snia, z_mod, dL_dimless)
    dL_Mpc = dL_interp * (C_KMS / H0phys)
    mu_mod = 5.0 * np.log10(np.clip(dL_Mpc, 1e-6, None)) + 25.0

    # Fehler auf μ: kombinieren
    err = np.array(mu_err_snia, float)

    # Peculiar-Velocity-Fehler (in mag): σ_μ,pec = (5/ln 10) * σ_v / (c z)
    if vpec_kms and np.any(z_snia > 0):
        sigma_mu_pec = (5.0/np.log(10.0)) * (vpec_kms / (C_KMS * np.clip(z_snia, 1e-4, None)))
        err = np.sqrt(err**2 + sigma_mu_pec**2)

    # Intrinsische Streuung σ_int (mag) additiv in Quadratur
    if sigma_int and sigma_int > 0:
        err = np.sqrt(err**2 + sigma_int**2)

    # 𝓜 analytisch profilieren:
    # resid_raw = mu_data - mu_mod - 𝓜;  M* = <mu_data - mu_mod>_w
    w = 1.0 / np.clip(err, 1e-12, None)**2
    resid0 = mu_snia - mu_mod
    M_star = np.sum(w * resid0) / np.sum(w)
    resid = resid0 - M_star

    if cov_path:
        # Falls du später die volle COV nutzt, musst du hier gemeinsam M profilieren (leicht aufwendiger).
        # Fürs Quick-Setup bleiben wir beim diagonalen Fall.
        C = np.load(cov_path)
        invC = np.linalg.inv(C)
        chi2 = float(resid @ invC @ resid)
    else:
        chi2 = float(np.sum((resid / np.clip(err, 1e-12, None))**2))

    # Zusätzlich den profilierten Offset zurückgeben (nützlich fürs Logging)
    return chi2, float(M_star), err

def load_bright_sirens(csv_path):
    """Load GW CSV: expects z, dL_Mpc, dL_err_Mpc."""
    df = pd.read_csv(csv_path)
    need = {"z", "dL_Mpc", "dL_err_Mpc"}
    if not need.issubset(df.columns):
        raise ValueError(f"GW CSV missing {need - set(df.columns)}")
    return df.sort_values("z")

def chi2_gw(z_gw, dL_gw_Mpc, dL_gw_err_Mpc, z_mod, dL_dimless, H0phys, cov_path=None):
    """χ² for GW bright sirens: direct dL comparison."""
    dL_interp = np.interp(z_gw, z_mod, dL_dimless)
    dL_mod_Mpc = dL_interp * (C_KMS / H0phys)
    resid = dL_gw_Mpc - dL_mod_Mpc

    if cov_path:
        C = np.load(cov_path)
        if C.shape != (len(resid), len(resid)):
            raise ValueError(f"GW cov shape mismatch")
        invC = np.linalg.inv(C)
        chi2 = float(resid @ invC @ resid)
    else:
        chi2 = float(np.sum((resid / np.clip(dL_gw_err_Mpc, 1e-12, None))**2))
    return chi2

# ---------- BAO helpers (consistent with your bao_compare.py logic) ----------
def load_model_csv(path):
    df = pd.read_csv(path)
    for col in ["z", "H_over_H0", "dL"]:
        if col not in df.columns:
            raise ValueError(f"Model CSV missing column '{col}'")
    return df.sort_values("z").reset_index(drop=True)

def model_bao_at_z(zq, z_mod, H_over_H0, dL_dimless, H0phys, rd):
    """
    Liefert (DM_over_rd, DH_over_rd, DV_over_rd) beim Ziel-z=zq.
    Erwartet:
      - z_mod: z-Achse deines Modell-CSV
      - H_over_H0: H/H0 (modell)
      - dL_dimless: d_L * H0 / c (dimensionlos!) aus deinem CSV
      - H0phys: km/s/Mpc
      - rd: Mpc
    """
    # Interpolation auf zq
    Hn      = np.interp(zq, z_mod, H_over_H0)
    dL_dim  = np.interp(zq, z_mod, dL_dimless)

    # Dimensionlos -> Mpc
    dL_Mpc  = dL_dim * (C_KMS / H0phys)         # d_L in Mpc
    DM_Mpc  = dL_Mpc / (1.0 + zq)               # D_M = dL/(1+z)
    H_km    = Hn * H0phys                       # H(z) in km/s/Mpc
    DH_Mpc  = C_KMS / H_km                      # D_H = c/H(z) in Mpc

    # Eisenstein-Hu D_V (hat Distanz-Dimension) -> dann / r_d
    DV_Mpc  = (zq * DM_Mpc**2 * DH_Mpc) ** (1.0/3.0)

    # Für DESI-Vergleich: alles / r_d
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
            if dcol in bao.columns and ecol in bao.columns and pd.notna(r.get(dcol)) and pd.notna(r.get(ecol)):
                rows.append(dict(z=z, obs=name, data=float(r[dcol]), err=float(r[ecol]), model=float(mval)))
    df = pd.DataFrame(rows)
    if len(df) == 0:
        return 0.0, 0, df

    data = df["data"].values
    model = df["model"].values
    resid = data - model

    if cov_path:
        C = np.load(cov_path)
        if C.shape != (len(df), len(df)):
            raise ValueError(f"BAO covariance shape {C.shape} does not match vector length {len(df)}")
        invC = np.linalg.inv(C)
        chi2 = float(resid @ invC @ resid)
    else:
        chi2 = float(np.sum((resid / np.clip(df["err"].values, 1e-12, None))**2))

    if export_vector is not None:
        export_vector["bao_vector"] = df.to_dict(orient="records")
    return chi2, len(df), df

def build_bao_dataframe_for_residuals(bao_csv_df, model_df, H0phys, rd):
    """
    Baut ein langes DataFrame mit Spalten:
      tracer, z, obs, data, err, model, resid, pull
    obs ∈ {"DM_over_rd","DH_over_rd","DV_over_rd"}
    """
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
                dval = float(dval); eval_ = float(eval_)
                resid = dval - float(mval)
                pull  = resid / eval_ if eval_ > 0 else np.nan
                rows.append(dict(
                    tracer=r.get("tracer",""),
                    z=z, obs=name,
                    data=dval, err=eval_,
                    model=float(mval), resid=resid, pull=pull
                ))
    return pd.DataFrame(rows)

def plot_bao_with_residuals(df_long, out_base, title_suffix=""):
    """
    Erzeugt 2×3-Figure:
      obere Reihe: Daten±σ vs. Modell
      untere Reihe: Pulls ( (data-model)/σ )
    Speichert nach f"{out_base}_bao_residuals.png"
    """
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
    ap.add_argument("--cmb-rd", type=float, default=None, help="mean of r_d prior [Mpc]")
    ap.add_argument("--cmb-rd-sigma", type=float, default=None, help="sigma of r_d prior [Mpc]")
    ap.add_argument("--H0phys", type=float, required=True, help="physical H0 [km/s/Mpc]")
    ap.add_argument("--rd", type=float, required=True, help="sound horizon r_d [Mpc] used in BAO compression")
    ap.add_argument("--out", default="runs/joint_out")
    ap.add_argument("--save-json", action="store_true")
    ap.add_argument("--plot-resids", action="store_true", help="plot BAO residuals")
    ap.add_argument("--n-params", type=int, default=5,
                help="Number of free model parameters for reduced chi^2")
    ap.add_argument("--snia-sigma-int", type=float, default=0.0, help="intrinsische SN-Streuung in mag (additiv)")
    ap.add_argument("--snia-vpec",      type=float, default=0.0, help="Peculiar-velocity Streuung in km/s (additiv)")
    args = ap.parse_args()

    model = load_model_csv(args.model_csv)

    report = {}

    # BAO
    bao_export = {}
    chi2_b, n_b, bao_df = chi2_bao(args.bao_csv, model, args.H0phys, args.rd, cov_path=(args.cov or None), export_vector=bao_export)
    report["chi2_bao"] = chi2_b
    report["N_bao"] = n_b

    # SN Ia
    chi2_sn, n_sn = 0.0, 0
    sn_Mstar = None
    if args.snia_csv:
        sn = load_snia_csv(args.snia_csv)
        chi2_sn, sn_Mstar, sn_err_eff = chi2_snia_profileM(
            sn["z"].values, sn["mu"].values, sn["mu_err"].values,
            model["z"].values, model["dL"].values, args.H0phys,
            cov_path=(args.snia_cov or None),
            sigma_int=args.snia_sigma_int,
            vpec_kms=args.snia_vpec
        )
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
            model["z"].values, model["dL"].values, args.H0phys,
            cov_path=(args.gw_cov or None)
        )
        n_gw = len(gw)
    report["chi2_gw"] = chi2_gw_val
    report["N_gw"] = n_gw

    # CMB r_d prior
    chi2_cmb = 0.0
    if args.cmb_rd is not None and args.cmb_rd_sigma is not None:
        chi2_cmb = chi2_rd_prior(args.rd, args.cmb_rd, args.cmb_rd_sigma)
    report["chi2_cmb_prior"] = chi2_cmb

    # Totals
    chi2_tot = chi2_b + chi2_sn + chi2_gw_val + chi2_cmb
    N_tot = n_b + n_sn + n_gw + (1 if (args.cmb_rd is not None and args.cmb_rd_sigma is not None) else 0)
    
    # zusätzliche genutzte Nuisance-Parameter
    nuisance = 1 if args.snia_csv else 0  # 𝓜 profiliert
    red_chi2 = chi2_tot / max(N_tot - args.n_params - nuisance, 1)
    
    report["chi2_total"] = chi2_tot
    report["N_total"] = N_tot
    report["reduced_chi2"] = red_chi2
    report["H0phys"] = args.H0phys
    report["rd_used"] = args.rd
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
    }

    print(f"χ² components → BAO: {chi2_b:.3f} (N={n_b}) | SN: {chi2_sn:.3f} (N={n_sn}) | GW: {chi2_gw_val:.3f} (N={n_gw}) | CMB r_d: {chi2_cmb:.3f}")
    print(f"χ² total = {chi2_tot:.3f}   (N={N_tot}, reduced χ² = {red_chi2:.3f})")

    # Save JSON summary
    if args.save_json:
        out_json = Path(f"{args.out}.json")
        payload = {"report": report, "bao_vector": bao_export.get("bao_vector", [])}
        out_json.parent.mkdir(parents=True, exist_ok=True)
        with open(out_json, "w") as f:
            json.dump(payload, f, indent=2)
        print(f"→ JSON saved: {out_json}")

    # Optional BAO residuals plot
    if args.plot_resids:
        bao_long = build_bao_dataframe_for_residuals(pd.read_csv(args.bao_csv), model, args.H0phys, args.rd)
        plot_bao_with_residuals(bao_long, args.out, title_suffix="(Joint Fit)")

if __name__ == "__main__":
    main()