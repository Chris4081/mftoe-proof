#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert FITS tables -> CSV (or ECSV if multi-dim columns are present).
- Auto-detects vector/array columns and falls back to ECSV.
- Optional: --flatten to expand vector columns into scalar columns col_0, col_1, ...
- Optional: --columns col1,col2,... to select a subset of columns.
Usage:
  python scripts/fits2csv.py <raw_fits_dir> <out_dir> [--flatten] [--columns z,D_H_over_r_d]
"""

import sys, os, argparse
from pathlib import Path
import numpy as np

from astropy.table import Table
from astropy.io import fits

def has_vector_cols(tab: Table) -> bool:
    for name in tab.colnames:
        shape = tab[name].shape
        # shape like (N,) is scalar per row; (N,k,...) means vector per row
        if len(shape) > 1:
            return True
    return False

def flatten_table(tab: Table) -> Table:
    """Expand vector columns into scalar columns col_0, col_1, ..."""
    out = Table()
    # copy scalar columns
    for name in tab.colnames:
        arr = tab[name]
        if len(arr.shape) == 1:
            out[name] = arr
    # expand vector columns
    for name in tab.colnames:
        arr = tab[name]
        if len(arr.shape) > 1:
            # per-row shape excluding the first axis is arr.shape[1:]
            # For a typical table, shape is (N, k) -> vector of length k
            # If higher dims, flatten them
            tail = int(np.prod(arr.shape[1:]))
            flat = arr.reshape(len(arr), tail)
            for j in range(tail):
                out[f"{name}_{j}"] = flat[:, j]
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw_dir", help="Directory with FITS files")
    ap.add_argument("out_dir", help="Directory for CSV/ECSV output")
    ap.add_argument("--flatten", action="store_true", help="Flatten vector columns into scalar columns")
    ap.add_argument("--columns", default="", help="Comma-separated subset of columns to keep (applied after flatten)")
    args = ap.parse_args()

    raw = Path(args.raw_dir)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)

    fits_files = sorted(list(raw.glob("*.fits")))
    if not fits_files:
        print(f"No FITS files found in {raw}")
        sys.exit(0)

    subset = [c.strip() for c in args.columns.split(",") if c.strip()]

    for fpath in fits_files:
        try:
            # Read first table HDU automatically; if specific HDU needed, pass hdu=1
            tab = Table.read(fpath, format="fits")
        except Exception as e:
            print(f"⚠️  Could not read {fpath.name}: {e}")
            continue

        # Optionally flatten
        if args.flatten and has_vector_cols(tab):
            tab = flatten_table(tab)

        # Optional column subset
        if subset:
            missing = [c for c in subset if c not in tab.colnames]
            if missing:
                print(f"⚠️  {fpath.name}: missing requested columns {missing} — keeping available ones.")
            keep = [c for c in subset if c in tab.colnames]
            if keep:
                tab = tab[keep]

        out_name_base = fpath.stem

        # Decide CSV vs ECSV
        vector = has_vector_cols(tab)
        if vector and not args.flatten:
            # Use ECSV to preserve structured/array columns
            out_path = out / f"{out_name_base}.ecsv"
            tab.write(out_path, format="ascii.ecsv", overwrite=True)
            print(f"→ Wrote ECSV: {out_path.name} (vector columns preserved)")
        else:
            # Safe to write CSV
            out_path = out / f"{out_name_base}.csv"
            try:
                tab.write(out_path, format="ascii.csv", overwrite=True, fast_writer=False)
                print(f"→ Wrote CSV:  {out_path.name}")
            except Exception as e:
                # Last resort: ECSV
                out_path = out / f"{out_name_base}.ecsv"
                tab.write(out_path, format="ascii.ecsv", overwrite=True)
                print(f"→ CSV failed, wrote ECSV instead: {out_path.name} ({e})")

if __name__ == "__main__":
    main()