#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plot scan results for the MFToE relaxion grid with Δχ² contours.

Generates:
  - Heatmap of χ² (or Δχ²) over (gamma, sigma)
  - Δχ² contours at 1σ/2σ/3σ (2 params): 2.30 / 6.18 / 11.83
  - Best-fit marker
  - Optional trend plot vs gamma per sigma

Usage:
  python3 analysis/plot_scan.py \
      --csv runs/scan_relaxion_summary.csv \
      --out runs/scan_relaxion_heatmap \
      [--metric total|delta] [--no-trend] \
      [--style publication|notebook] \
      [--formats png,pdf,svg]

Outputs:
  <out>.(png|pdf|svg)                  # heatmap + contours
  <out>_trend.(png|pdf|svg)            # trends (unless --no-trend)
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from pathlib import Path

# Δχ² thresholds for 2 free parameters (68.3%, 95.4%, 99.7%)
LEVELS_2PAR = [2.30, 6.18, 11.83]


def apply_style(style: str):
    """Set Matplotlib rcParams for publication or notebook styles."""
    if style == "publication":
        mpl.rcParams.update({
            "figure.dpi": 150,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 1.0,
            "axes.grid": True,
            "grid.alpha": 0.25,
            "grid.linestyle": "-",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "lines.linewidth": 1.8,
            "patch.linewidth": 0.6,
            "pdf.fonttype": 42,   # TrueType in PDF
            "ps.fonttype": 42,
            "text.usetex": False, # keep simple; avoids TeX dependency
        })
    else:
        # notebook defaults (a bit larger labels)
        mpl.rcParams.update({
            "figure.dpi": 120,
            "savefig.dpi": 160,
            "axes.grid": True,
            "grid.alpha": 0.3,
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 12,
            "legend.fontsize": 10,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "lines.linewidth": 1.6,
        })


def load_scan(csv_path: str):
    """Load scan CSV and extract relevant χ² metric."""
    df = pd.read_csv(csv_path)
    if "gamma" not in df.columns or "sigma" not in df.columns:
        raise ValueError("CSV must include columns: gamma, sigma")

    # Auto-select chi² source
    if "chi2_total" in df.columns:
        df["chi2_used"] = df["chi2_total"].astype(float)
    elif "chi2" in df.columns:
        df["chi2_used"] = df["chi2"].astype(float)
    else:
        raise ValueError("CSV must contain either 'chi2_total' or 'chi2'.")

    df["gamma"] = df["gamma"].astype(float)
    df["sigma"] = df["sigma"].astype(float)
    return df


def pivot_grid(df):
    """Build a 2D grid (sigma × gamma) from the scan table."""
    gammas = np.sort(df["gamma"].unique())
    sigmas = np.sort(df["sigma"].unique())
    Z = np.full((len(sigmas), len(gammas)), np.nan)

    for i, s in enumerate(sigmas):
        for j, g in enumerate(gammas):
            sub = df[(df["gamma"] == g) & (df["sigma"] == s)]
            if len(sub) == 1:
                Z[i, j] = float(sub["chi2_used"].values[0])
    return gammas, sigmas, Z


def centers_to_edges(x):
    """Convert value centers to bin edges for pcolormesh."""
    x = np.asarray(x, float)
    if len(x) == 1:
        dx = 0.5 * max(1e-6, abs(x[0]) * 0.1)
        return np.array([x[0] - dx, x[0] + dx])
    dx = np.diff(x)
    edges = np.concatenate(([x[0] - dx[0] / 2], x[:-1] + dx / 2, [x[-1] + dx[-1] / 2]))
    return edges


def save_fig(fig: mpl.figure.Figure, base: str, formats):
    Path(base).parent.mkdir(parents=True, exist_ok=True)
    for ext in formats:
        fig.savefig(f"{base}.{ext}", dpi=mpl.rcParams["savefig.dpi"], bbox_inches="tight")
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--out", required=True, help="Output base path (no extension)")
    ap.add_argument("--metric", choices=["total", "delta"], default="total",
                    help="'total' plots absolute χ²; 'delta' plots Δχ² from best-fit")
    ap.add_argument("--no-trend", action="store_true", help="Skip trend plot")
    ap.add_argument("--style", choices=["publication", "notebook"], default="notebook",
                    help="Figure styling preset")
    ap.add_argument("--formats", default="png",
                    help="Comma-separated list of formats to save: e.g. 'png,pdf,svg'")
    args = ap.parse_args()

    apply_style(args.style)
    formats = [f.strip().lower() for f in args.formats.split(",") if f.strip()]

    df = load_scan(args.csv)
    gammas, sigmas, Z = pivot_grid(df)
    if np.all(np.isnan(Z)):
        raise RuntimeError("Empty grid — check CSV content or parsing.")

    # --- Best fit ---
    min_idx = np.nanargmin(Z)
    i_best, j_best = np.unravel_index(min_idx, Z.shape)
    chi2_min = Z[i_best, j_best]
    g_best, s_best = gammas[j_best], sigmas[i_best]

    # --- Metric selection ---
    if args.metric == "delta":
        Z_plot = Z - chi2_min
        title_metric = r"$\Delta\chi^2$"
    else:
        Z_plot = Z.copy()
        title_metric = r"$\chi^2$"

    # --- Heatmap + Contours ---
    fig, ax = plt.subplots(figsize=(8.0, 5.8))
    g_edges = centers_to_edges(gammas)
    s_edges = centers_to_edges(sigmas)

    pcm = ax.pcolormesh(g_edges, s_edges, Z_plot, shading="auto", cmap="viridis")
    cbar = plt.colorbar(pcm, ax=ax)
    cbar.set_label(title_metric)

    if args.metric == "delta":
        CS = ax.contour(gammas, sigmas, Z_plot, levels=LEVELS_2PAR,
                        colors="white", linewidths=1.6)
        ax.clabel(CS, inline=True, fmt={lv: f"{lv:.2f}" for lv in LEVELS_2PAR},
                  fontsize=8, colors="white")

    ax.plot(g_best, s_best, "o", ms=7, color="red",
            label=f"Best: γ={g_best:g}, σ={s_best:g}, χ²={chi2_min:.3f}")

    ax.set_xlabel(r"$\gamma$ (feedback strength)")
    ax.set_ylabel(r"$\sigma_{\mathrm{noise}}$")
    ax.set_title(f"Relaxion parameter scan ({title_metric})")
    ax.legend(loc="best", frameon=True)
    ax.grid(alpha=0.25)

    save_fig(fig, args.out, formats)
    print(f"✅ Heatmap saved: {', '.join(f'{args.out}.{e}' for e in formats)}")

    # --- Trend plot ---
    if not args.no_trend:
        fig2, ax2 = plt.subplots(figsize=(8.0, 5.0))
        for i, s in enumerate(sigmas):
            y = Z_plot[i, :]
            ax2.plot(gammas, y, marker="o", label=f"σ={s:g}")
        ax2.axhline(0 if args.metric == "delta" else chi2_min, lw=1.0, alpha=0.5)
        ax2.set_xlabel(r"$\gamma$")
        ax2.set_ylabel(title_metric)
        ax2.set_title(f"Trends across γ ({title_metric})")
        ax2.grid(alpha=0.3)
        ax2.legend(ncol=2)
        save_fig(fig2, f"{args.out}_trend", formats)
        print(f"✅ Trend plot saved: {', '.join(f'{args.out}_trend.{e}' for e in formats)}")

    # --- Summary to stdout ---
    print("Best-fit parameters:")
    print(f"  gamma  = {g_best:g}")
    print(f"  sigma  = {s_best:g}")
    print(f"  chi²_min = {chi2_min:.3f}")
    if args.metric != "delta":
        print("Tip: re-run with --metric delta for 1σ/2σ/3σ contours.")


if __name__ == "__main__":
    main()