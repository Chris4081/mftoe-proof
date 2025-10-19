#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd

C_KMS = 299_792.458  # km/s

def load_snia_csv(path):
    """
    Expect CSV with columns: z, mu, mu_err
    Returns DataFrame sorted by z.
    """
    df = pd.read_csv(path)
    for col in ["z", "mu", "mu_err"]:
        if col not in df.columns:
            raise ValueError(f"SN Ia CSV missing column '{col}'")
    return df.sort_values("z").reset_index(drop=True)

def interp_mu_model(z_query, z_mod, dL_dimless, H0phys):
    """
    Convert model dimensionless dL (in H0/c units) to Mpc and to distance modulus.
    dL_Mpc = dL_dimless * c / H0
    mu = 5 log10(dL / 10 pc) = 5 log10(dL_Mpc) + 25
    """
    dL_q = np.interp(z_query, z_mod, dL_dimless)
    dL_mpc = dL_q * (C_KMS / H0phys)
    mu_mod = 5.0 * np.log10(np.clip(dL_mpc, 1e-12, None)) + 25.0
    return mu_mod

def chi2_snia_profileM(z, mu_data, mu_err, z_mod, dL_dimless_mod, H0phys, cov=None):
    """
    SN Ia chi2 with analytic profiling over absolute magnitude M (additive offset).
    If 'cov' is None, use diagonal from mu_err.
    """
    mu_mod = interp_mu_model(z, z_mod, dL_dimless_mod, H0phys)
    r0 = mu_data - mu_mod  # vector without M

    if cov is None:
        W = np.diag(1.0 / np.clip(mu_err, 1e-12, None)**2)
    else:
        W = np.linalg.inv(cov)

    ones = np.ones_like(r0)
    a = ones @ W @ ones
    b = ones @ W @ r0
    # profiled chi2 (see e.g. Conley+ 2011)
    chi2 = (r0 @ W @ r0) - (b * b) / max(a, 1e-30)
    return float(chi2)