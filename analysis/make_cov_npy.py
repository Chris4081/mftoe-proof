#!/usr/bin/env python3
import sys, numpy as np, pandas as pd

def load_matrix(path):
    # Versucht zuerst als blank-getrennte Zahlen, sonst CSV mit Header
    try:
        return np.loadtxt(path)
    except Exception:
        df = pd.read_csv(path)
        # Falls die Datei Spaltennamen hat, entferne nicht-numerische
        df = df.select_dtypes(include=["number"])
        return df.values

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: make_cov_npy.py <cov.txt|csv> <out.npy>")
        sys.exit(1)
    M = load_matrix(sys.argv[1])
    np.save(sys.argv[2], M)
    print(f"Saved {sys.argv[2]} with shape {M.shape}")