#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
cd "$repo_root"
mkdir -p runs

# Grid
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

    # 1) MFToE-Integration
    python3 mftoe_vacuum_astropy.py \
      --mode relaxion --rg on --noise on \
      --astropy on \
      --desi-bestfit data/desi_dr2/iminuit/base/desi-bao-all/bestfit.minimum \
      --params "${tmpjson}" \
      --out "runs/${tag}"

    rm -f "${tmpjson}"

    # 2) Joint-Fit (BAO + SNIa + CMB r_d prior); robustere SN-Fehler
    joint_json="runs/${tag}_joint.json"
    python3 analysis/joint_fit.py \
      --model-csv "runs/${tag}.csv" \
      --bao-csv   data/desi_dr2/bao_summary.csv \
      --snia-csv  data/snia/pantheonplus_summary.csv \
      --snia-sigma-int 0.10 --snia-vpec 250 \
      --H0phys 67.36 --rd 150.754 \
      --cmb-rd 150.74 --cmb-rd-sigma 0.30 \
      --n-params 5 \
      --out "runs/${tag}_joint" \
      --save-json

    # 3) χ² aus JOINT-JSON extrahieren (Total)
    if [[ -f "${joint_json}" ]]; then
      chi2=$(python3 - <<PY
import json
with open("${joint_json}") as f:
    d = json.load(f)
print(d["report"]["chi2_total"])
PY
)
      echo "joint,${g},${s},${chi2}" >> runs/scan_relaxion_summary.csv
    else
      echo "❌ Joint JSON fehlt: ${joint_json} — überspringe Eintrag"
      echo "joint,${g},${s},NA" >> runs/scan_relaxion_summary.csv
    fi

  done
done

echo "✅ Scan finished. See runs/scan_relaxion_summary.csv"