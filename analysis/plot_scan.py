#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot scan results for the MFToE relaxion grid with Δχ² contours.

Reads a CSV produced by scripts/scan_relaxion.sh and generates:
  - Heatmap of χ² (or Δχ²) over (gamma, sigma)
  - Δχ² contour lines at 1σ/2σ/3σ for 2 parameters: 2.30/6.18/11.83
  - Best-fit marker
  - Optional trend plot vs. gamma (for each sigma)

Usage:
  python3 analysis/plot_scan.py \
      --csv runs/scan_relaxion_summary.csv \
      --out runs/scan_relaxion_heatmap \
      [--metric total|delta] [--no-trend]

Outputs:
  <out>.png                     # heatmap + contours
  <out>_trend.png               # trends (unless --no-trend)

CSV compatibility:
  - New format (recommended):
      tag,gamma,sigma,chi2_bao,chi2_sn,chi2_gw,chi2_cmb,chi2_total,reduced_chi2,N_total
  - Minimal format (fallback):
      metric,gamma,sigma,chi2
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

LEVELS_2PAR = [2.30, 6.18, 11.83]  # 1σ, 2σ, 3σ (for 2 params)

def load_scan(csv_path: str):
    df = pd.read_csv(csv_path)
    cols = set(df.columns.str.lower())
    # Normalize columns
    if "gamma" not in df.columns or "sigma" not in df.columns:
        raise ValueError("CSV needs at least columns: gamma, sigma")
    # Choose chi2 source
    if "chi2_total" in df.columns:
        df["chi2_used"] = df["chi2_total"].astype(float)
    elif "chi2" in df.columns:
        df["chi2_used"] = df["chi2"].astype(float)
    else:
        raise ValueError("CSV must contain chi2_total or chi2.")
    # Ensure numeric
    df["gamma"] = df["gamma"].astype(float)
    df["sigma"] = df["sigma"].astype(float)
    return df

def pivot_grid(df):
    gammas = np.sort(df["gamma"].unique())
    sigmas = np.sort(df["sigma"].unique())
    # Build Z matrix with NaN fill
    Z = np.full((len(sigmas), len(gammas)), np.nan)
    for i, s in enumerate(sigmas):
        for j, g in enumerate(gammas):
            sub = df[(df["gamma"]==g) & (df["sigma"]==s)]
            if len(sub) == 1:
                Z[i, j] = float(sub["chi2_used"].values[0])
    return gammas, sigmas, Z

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out", required=True, help="output base path (no extension)")
    ap.add_argument("--metric", choices=["total","delta"], default="total",
                    help="'total' plots χ², 'delta' plots Δχ² from best")
    ap.add_argument("--no-trend", action="store_true", help="skip trend panel export")
    args = ap.parse_args()

    df = load_scan(args.csv)
    gammas, sigmas, Z = pivot_grid(df)

    if np.all(np.isnan(Z)):
        raise RuntimeError("Grid is empty after pivot; check CSV content.")

    # Best fit
    min_idx = np.nanargmin(Z)
    i_best, j_best = np.unravel_index(min_idx, Z.shape)
    chi2_min = Z[i_best, j_best]
    g_best, s_best = gammas[j_best], sigmas[i_best]

    # Choose plotting matrix
    if args.metric == "delta":
        Z_plot = Z - chi2_min
        title_metric = r"$\Delta\chi^2$"
    else:
        Z_plot = Z.copy()
        title_metric = r"$\chi^2$"

    # --- Heatmap + contours ---
    fig, ax = plt.subplots(figsize=(8.5, 6.2))

    # pcolormesh expects bin edges; build simple edges from unique centers
    def centers_to_edges(x):
        x = np.asarray(x, float)
        if len(x) == 1:
            dx = 0.5*max(1e-6, abs(x[0])*0.1)
            return np.array([x[0]-dx, x[0]+dx])
        dx = np.diff(x)
        edges = np.concatenate(([x[0]-dx[0]/2], x[:-1]+dx/2, [x[-1]+dx[-1]/2]))
        return edges

    g_edges = centers_to_edges(gammas)
    s_edges = centers_to_edges(sigmas)

    # Heatmap
    pcm = ax.pcolormesh(g_edges, s_edges, Z_plot, shading="auto")
    cbar = plt.colorbar(pcm, ax=ax, label=title_metric)

    # Δχ² contours (only meaningful on Δχ²)
    if args.metric == "delta":
        CS = ax.contour(gammas, sigmas, Z_plot, levels=LEVELS_2PAR, linewidths=1.8)
        ax.clabel(CS, inline=True, fmt={lv: f"{lv:.2f}" for lv in LEVELS_2PAR}, fontsize=9)

    # Best-fit marker
    ax.plot(g_best, s_best, "o", ms=8, label=f"best: γ={g_best:g}, σ={s_best:g}, χ²={chi2_min:.3f}")

    ax.set_xlabel(r"$\gamma$")
    ax.set_ylabel(r"$\sigma_{\rm noise}$")
    ax.set_title(f"Relaxion scan heatmap ({title_metric}); best χ²={chi2_min:.3f}")
    ax.legend(loc="best", frameon=True)
    ax.grid(alpha=0.25)

    out_png = f"{args.out}.png"
    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(out_png, dpi=160)
    plt.close(fig)
    print(f"✅ Heatmap saved: {out_png}")

    # --- Trend plot (optional): χ² vs gamma for each sigma ---
    if not args.no_trend:
        fig2, ax2 = plt.subplots(figsize=(8.5, 5.6))
        for i, s in enumerate(sigmas):
            y = Z_plot[i, :]
            ax2.plot(gammas, y, marker="o", label=f"σ={s:g}")
        ax2.axhline(0 if args.metric=="delta" else chi2_min, lw=1, alpha=0.5)
        ax2.set_xlabel(r"$\gamma$")
        ax2.set_ylabel(title_metric)
        ax2.set_title(f"Trend vs. γ ({title_metric})")
        ax2.grid(alpha=0.3)
        ax2.legend(ncol=2)
        out_trend = f"{args.out}_trend.png"
        plt.tight_layout()
        plt.savefig(out_trend, dpi=160)
        plt.close(fig2)
        print(f"✅ Trend plot saved: {out_trend}")

    # Print a tiny summary to stdout
    print(f"Best-fit: gamma={g_best:g}, sigma={s_best:g}, chi2_min={chi2_min:.3f}")
    if args.metric != "delta":
        print("Tip: re-run with --metric delta to see Δχ² contours (1σ/2σ/3σ).")

if __name__ == "__main__":
    main()