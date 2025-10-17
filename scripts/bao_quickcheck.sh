#!/usr/bin/env bash
set -euo pipefail

# Compare the baseline run to DESI BAO, no covariance (or add --cov path below)
python3 analysis/bao_compare.py \
  --model-csv runs/mftoe_vacuum_astropy.csv \
  --bao-csv   data/desi_dr2/bao_summary.csv \
  --H0phys 67.36 --rd 150.754 \
  --show-residuals --export-json \
  --out runs/bao_check