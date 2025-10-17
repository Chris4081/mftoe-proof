#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
cd "$repo_root"
mkdir -p runs

# Simple grid
gammas=(0.05 0.10 0.20)
sigmas=(1e-5 1e-4)

echo "metric,gamma,sigma,chi2" > runs/scan_relaxion_summary.csv

for g in "${gammas[@]}"; do
  for s in "${sigmas[@]}"; do
    tag="relaxion_g${g}_s${s}"
    echo "→ Running ${tag}…"
    # temporäres Param-JSON (nur Overrides)
    tmpjson="$(mktemp)"
    cat > "${tmpjson}" <<JSON
{"gamma": ${g}, "sigma_noise": ${s}}
JSON

    python3 mftoe_vacuum_astropy.py \
      --mode relaxion --rg on --noise on \
      --astropy on \
      --desi-bestfit data/desi_dr2/iminuit/base/desi-bao-all/bestfit.minimum \
      --params "${tmpjson}" \
      --out "runs/${tag}"

    # BAO compare (no cov for speed)
    out_json="runs/bao_${tag}.json"
    python3 analysis/bao_compare.py \
      --model-csv "runs/${tag}.csv" \
      --bao-csv   data/desi_dr2/bao_summary.csv \
      --H0phys 67.36 --rd 150.754 \
      --export-json \
      --out "runs/bao_${tag}" >/dev/null

    # Pull χ² aus JSON
    chi2=$(python3 - <<PY
import json,sys
d=json.load(open("${out_json}"))
print(d.get("chi2", "NA"))
PY
)
    echo "bao,${g},${s},${chi2}" >> runs/scan_relaxion_summary.csv
    rm -f "${tmpjson}"
  done
done

echo "✅ Scan finished. See runs/scan_relaxion_summary.csv"
