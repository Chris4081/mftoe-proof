#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")"/.. && pwd)"
cd "$repo_root"
mkdir -p runs

# -------------------- Scan-Grid --------------------
gammas=(0.05 0.10 0.20)
sigmas=(1e-5 1e-4)

# -------------------- Toggles ----------------------
USE_CAMB_RD=true         # true → r_d aus CAMB + H0*r_d-Matching aktivieren
REF_H0=67.36             # Referenz H0 (DESI DR2 compressed)
REF_RD=150.754           # Referenz r_d (DESI DR2 compressed)

# CAMB Early-Params (Planck-ähnlich)
OMBH2=0.02237
OMCH2=0.1200
NEFF=3.046
YP=0.245
MNU=0.06

# -------------------- Data -------------------------
DESI_BESTFIT="data/desi_dr2/iminuit/base/desi-bao-all/bestfit.minimum"
BAO_CSV="data/desi_dr2/bao_summary.csv"
SN_CSV="data/snia/pantheonplus_summary.csv"   # "" falls nicht vorhanden
SN_SIGMA_INT=0.10
SN_VPEC=250

# -------------------- Summary-Header ---------------
summary="runs/scan_relaxion_summary.csv"
echo "tag,gamma,sigma,chi2_total,reduced_chi2,N_total,chi2_bao,chi2_sn,chi2_gw,chi2_cmb,H0_eff,rd_used" > "$summary"

for g in "${gammas[@]}"; do
  for s in "${sigmas[@]}"; do
    tag="relaxion_g${g}_s${s}"
    echo "→ Running ${tag}…"

    # 1) MFToE-Integration (Relaxion + RG + Noise)
    tmpjson="$(mktemp)"
    cat > "${tmpjson}" <<JSON
{"gamma": ${g}, "sigma_noise": ${s}}
JSON

    python3 mftoe_vacuum_astropy.py \
      --mode relaxion --rg on --noise on \
      --astropy on \
      --desi-bestfit "${DESI_BESTFIT}" \
      --params "${tmpjson}" \
      --out "runs/${tag}"

    rm -f "${tmpjson}"

    # 2) Joint-Fit (BAO + optional SNIa + CMB r_d prior)
    joint_json="runs/${tag}_joint.json"

    # Basis-Args
    joint_args=(
      --model-csv "runs/${tag}.csv"
      --bao-csv   "${BAO_CSV}"
      --H0phys "${REF_H0}"
      --cmb-rd 150.74
      --cmb-rd-sigma 0.30
      --n-params 5
      --out "runs/${tag}_joint"
      --save-json
    )

    # r_d-Backend + Matching
    if [[ "${USE_CAMB_RD}" == "true" ]]; then
      # rd wird intern aus CAMB geholt; Match H0*rd an (REF_H0, REF_RD)
      joint_args+=( --rd 0 --rd-backend camb --ombh2 "${OMBH2}" --omch2 "${OMCH2}" --Neff "${NEFF}" --Yp "${YP}" --mnu-eV "${MNU}" )
      joint_args+=( --match-H0rd --ref-H0 "${REF_H0}" --ref-rd "${REF_RD}" )
    else
      # fixed rd (DESI-kompatibel)
      joint_args+=( --rd "${REF_RD}" )
    fi

    # SNIa robust errors, falls Datei existiert
    if [[ -n "${SN_CSV}" && -f "${SN_CSV}" ]]; then
      joint_args+=( --snia-csv "${SN_CSV}" --snia-sigma-int "${SN_SIGMA_INT}" --snia-vpec "${SN_VPEC}" )
    fi

    python3 analysis/joint_fit.py "${joint_args[@]}"

    # 3) Werte aus Joint-JSON extrahieren
    if [[ -f "${joint_json}" ]]; then
      python3 - "$tag" "$summary" <<'PY'
import json, sys, pathlib
tag, summary = sys.argv[1], sys.argv[2]
j = pathlib.Path(f"runs/{tag}_joint.json")
rep = json.loads(j.read_text())["report"]
row = ",".join([
  tag,
  str(tag.split("_g")[-1].split("_s")[0]),   # gamma
  str(tag.split("_s")[-1]),                  # sigma
  f'{rep.get("chi2_total",0.0):.3f}',
  f'{rep.get("reduced_chi2",0.0):.3f}',
  str(rep.get("N_total",0)),
  f'{rep.get("chi2_bao",0.0):.3f}',
  f'{rep.get("chi2_snia",0.0):.3f}',
  f'{rep.get("chi2_gw",0.0):.3f}',
  f'{rep.get("chi2_cmb_prior",0.0):.3f}',
  f'{rep.get("H0phys_eff", rep.get("H0phys_input", 0.0)):.6f}',
  f'{rep.get("rd_used", 0.0):.6f}',
])
with open(summary, "a") as f:
    f.write(row + "\n")
PY
      echo "→ Saved ${joint_json} and updated summary."
    else
      echo "❌ Joint JSON fehlt: ${joint_json}"
    fi

  done
done

echo "✅ Scan finished. See ${summary}"