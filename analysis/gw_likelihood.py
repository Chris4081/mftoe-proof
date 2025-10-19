#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd

C_KMS = 299_792.458  # km/s

def load_bright_sirens(path):
    """
    CSV columns: name, z, dL_Mpc, dL_err_Mpc
    """
    df = pd.read_csv(path)
    for col in ["z", "dL_Mpc", "dL_err_Mpc"]:
        if col not in df.columns:
            raise ValueError(f"GW CSV missing column '{col}'")
    return df.sort_values("z").reset_index(drop=True)

def interp_dL_Mpc(zq, z_mod, dL_dimless_mod, H0phys):
    dL_dimless = np.interp(zq, z_mod, dL_dimless_mod)
    return dL_dimless * (C_KMS / H0phys)

def chi2_gw(z, dL_Mpc, dL_err, z_mod, dL_dimless_mod, H0phys):
    dL_model = interp_dL_Mpc(z, z_mod, dL_dimless_mod, H0phys)
    res = (dL_model - dL_Mpc) / np.clip(dL_err, 1e-12, None)
    return float(np.sum(res**2))