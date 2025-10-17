#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
cd "$repo_root"
mkdir -p runs

echo "→ Running MFToE (relaxion, RG+Noise)…"
python3 mftoe_vacuum_astropy.py \
  --mode relaxion --rg on --noise on \
  --astropy on \
  --desi-bestfit data/desi_dr2/iminuit/base/desi-bao-all/bestfit.minimum \
  --out runs/relaxion

echo "→ BAO check (no covariance)…"
python3 analysis/bao_compare.py \
  --model-csv runs/relaxion.csv \
  --bao-csv   data/desi_dr2/bao_summary.csv \
  --H0phys 67.36 --rd 150.754 \
  --out runs/bao_relaxion_check

echo "✅ Done. Outputs in runs/:"
ls -lh runs | sed 's/^/   /'
