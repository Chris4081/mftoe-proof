#!/usr/bin/env bash
# scan_relaxion.sh — BAO-only Joint-Fit (CAMB r_d, H0·r_d matching) ohne CMB-/SN-/GW-Prior
set -euo pipefail

# --- repo root & env ---
repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
cd "$repo_root"
mkdir -p runs

# Python-Importpfad
export PYTHONPATH="${repo_root}:${PYTHONPATH:-}"

# Reproducibility (optional)
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export NUMEXPR_NUM_THREADS="${NUMEXPR_NUM_THREADS:-1}"
export PYTHONHASHSEED="${PYTHONHASHSEED:-0}"

# --- required inputs sanity check ---
need_files=(
  "data/desi_dr2/bao_summary.csv"
  "data/desi_dr2/iminuit/base/desi-bao-all/bestfit.minimum"
)
for f in "${need_files[@]}"; do
  [[ -f "$f" ]] || { echo "❌ Missing required file: $f" >&2; exit 1; }
done

# optional BAO-Cov (falls vorhanden, wird genutzt)
bao_csv="data/desi_dr2/bao_summary.csv"
bao_cov_opt=()
if [[ -f "data/desi_dr2/bao_cov.npy" ]]; then
  bao_cov_opt=(--cov "data/desi_dr2/bao_cov.npy")
  echo "→ Using BAO covariance: data/desi_dr2/bao_cov.npy"
fi

# --- grid ---
gammas=(0.05 0.10 0.20)
sigmas=(1e-5 1e-4)

summary="runs/scan_relaxion_summary.csv"
echo "metric,gamma,sigma,chi2" > "$summary"

# --- helper: safe mktemp cleanup ---
_tmpfiles=()
cleanup() {
  for t in "${_tmpfiles[@]:-}"; do
    [[ -n "${t:-}" && -f "$t" ]] && rm -f "$t" || true
  done
  echo "✅ Scan finished. See ${summary}"
}
trap cleanup EXIT

for g in "${gammas[@]}"; do
  for s in "${sigmas[@]}"; do
    tag="relaxion_g${g}_s${s}"
    echo "→ Running ${tag}…"

    # 1) MFToE-Integration (nur Overrides via temp JSON)
    tmpjson="$(mktemp -t mftoe_relaxion.XXXXXX.json)"
    _tmpfiles+=("$tmpjson")
    cat > "${tmpjson}" <<JSON
{"gamma": ${g}, "sigma_noise": ${s}}
JSON

    model_csv="runs/${tag}.csv"
    if [[ ! -f "${model_csv}" ]]; then
      python3 ./mftoe_vacuum_astropy.py \
        --mode relaxion --rg on --noise on \
        --astropy on \
        --desi-bestfit data/desi_dr2/iminuit/base/desi-bao-all/bestfit.minimum \
        --params "${tmpjson}" \
        --out "runs/${tag}"
    else
      echo "  ↪︎ model already exists: ${model_csv} (skip generate)"
    fi

    [[ -f "${model_csv}" ]] || { echo "❌ Missing model CSV: ${model_csv}" >&2; exit 1; }

    # 2) Joint-Fit (BAO only, CAMB r_d, H0·r_d matching; KEIN CMB-Prior, KEINE SN/GW)
    out_base="runs/${tag}_joint"
    joint_json="${out_base}.json"

    if [[ ! -f "${joint_json}" ]]; then
      python3 ./analysis/joint_fit.py \
        --model-csv "${model_csv}" \
        --bao-csv   "${bao_csv}" \
        "${bao_cov_opt[@]}" \
        --H0phys 67.36 \
        --rd-backend camb \
        --ombh2 0.02237 --omch2 0.1200 --Neff 3.046 --Yp 0.245 --mnu-eV 0.06 \
        --match-H0rd --ref-H0 67.36 --ref-rd 150.754 \
        --no-cmb-prior \
        --out "${out_base}" \
        --save-json
    else
      echo "  ↪︎ joint JSON already exists: ${joint_json} (skip fit)"
    fi

    # 3) χ² aus JOINT-JSON extrahieren (Total)
    if [[ -f "${joint_json}" ]]; then
      chi2="$(python3 - <<'PY'
import json, sys
p=sys.argv[1]
with open(p) as f:
    d=json.load(f)
v=d.get("report",{}).get("chi2_total", None)
if v is None:
    raise SystemExit("missing chi2_total")
print(v)
PY
"${joint_json}")"
      echo "joint,${g},${s},${chi2}" >> "$summary"
      echo "→ χ²_total(${tag}) = ${chi2}"
    else
      echo "❌ Joint JSON fehlt: ${joint_json} — trage NA ein"
      echo "joint,${g},${s},NA" >> "$summary"
    fi

  done
done
