#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
cd "$repo_root"
mkdir -p runs

for rho in 0.0 0.2 0.3 0.4 0.5; do
  echo "→ Building synthetic covariance (rho=${rho})…"
  python3 analysis/make_cov_from_csv.py \
    data/desi_dr2/bao_summary.csv \
    data/desi_dr2/bao_cov.npy \
    "${rho}"

  echo "→ BAO compare with rho=${rho}…"
  python3 analysis/bao_compare.py \
    --model-csv runs/mftoe_vacuum_astropy.csv \
    --bao-csv   data/desi_dr2/bao_summary.csv \
    --H0phys 67.36 --rd 150.754 \
    --cov    data/desi_dr2/bao_cov.npy \
    --out    "runs/bao_check_cov_rho${rho}"
done

echo "✅ Sweep complete. Outputs in runs/:"
ls -lh runs | sed 's/^/   /'
