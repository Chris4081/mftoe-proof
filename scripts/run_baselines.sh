#!/usr/bin/env bash
set -euo pipefail

# Go to repo root (one level above scripts/)
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
cd "$repo_root"

mkdir -p runs

echo "→ Running MFToE baseline (Astropy ΛCDM reference + DESI bestfit params)…"
python3 mftoe_vacuum_astropy.py \
  --astropy on \
  --desi-bestfit data/desi_dr2/iminuit/base/desi-bao-all/bestfit.minimum \
  --out runs/mftoe_vacuum_astropy

echo "→ BAO quick check (no covariance)…"
python3 analysis/bao_compare.py \
  --model-csv runs/mftoe_vacuum_astropy.csv \
  --bao-csv   data/desi_dr2/bao_summary.csv \
  --H0phys 67.36 --rd 150.754 \
  --out runs/bao_check

echo "✅ Done. Outputs in runs/:"
ls -lh runs | sed 's/^/   /'
