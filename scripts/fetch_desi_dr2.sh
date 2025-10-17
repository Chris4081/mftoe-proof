#!/usr/bin/env bash
set -euo pipefail

# Wurzel relativ zum Repo (script darf von überall aufgerufen werden)
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RAW_DIR="$REPO_ROOT/data/desi_dr2/raw"
PROC_DIR="$REPO_ROOT/data/desi_dr2/processed"
URLS_FILE="$REPO_ROOT/scripts/desi_dr2_urls.txt"

mkdir -p "$RAW_DIR" "$PROC_DIR"

echo "→ Downloading DESI DR2/Y3 BAO FITS files into: $RAW_DIR"

# 1) URLs aus externer Textdatei lesen (eine URL pro Zeile).
#    Vorteil: Du kannst leicht die *korrekten* DR2-Links eintragen, sobald du sie hast.
if [[ ! -f "$URLS_FILE" ]]; then
  cat > "$URLS_FILE" <<'EOF'
# Put official DESI DR2 BAO FITS URLs here, one per line, no quotes.
# Example placeholders (these 404 by design – replace them!):
https://data.desi.lbl.gov/public/BAO/DR2/BAO_all_catalogs_v1.fits
https://data.desi.lbl.gov/public/BAO/DR2/BAO_results_v1.fits
https://data.desi.lbl.gov/public/BAO/DR2/BAO_covariance_v1.fits
EOF
  echo "⚠️  Created placeholder URL list at $URLS_FILE"
  echo "   Please open it and paste the official DR2 links (one per line)."
fi

# 2) Dateien laden (fehlschlagende URLs werden nur gewarnt)
while IFS= read -r url; do
  # leere Zeilen & Kommentare überspringen
  [[ -z "${url// }" ]] && continue
  [[ "$url" =~ ^# ]] && continue

  fname="$(basename "$url")"
  echo "  - $fname"
  if curl -fL --retry 3 --retry-delay 2 -o "$RAW_DIR/$fname" "$url"; then
    :
  else
    echo "⚠️  Warning: could not fetch $fname"
    # falls leere/kaputte Datei entstanden ist, löschen
    [[ -f "$RAW_DIR/$fname" && ! -s "$RAW_DIR/$fname" ]] && rm -f "$RAW_DIR/$fname"
  fi
done < "$URLS_FILE"

echo "→ Converting FITS → CSV with astropy…"
python3 "$REPO_ROOT/scripts/fits2csv.py" \
  --in-dir "$RAW_DIR" \
  --out-dir "$PROC_DIR"

echo "✅ Done."
echo "Raw FITS:      $RAW_DIR"
echo "Processed CSV: $PROC_DIR"