import numpy as np
import pandas as pd

C_KMS = 299_792.458  # km/s

def load_snia_csv(path):
    df = pd.read_csv(path)
    for col in ["z", "mu", "mu_err"]:
        if col not in df.columns:
            raise ValueError(f"SN Ia CSV missing column '{col}'")
    return df.sort_values("z").reset_index(drop=True)

def interp_mu_model(z_query, z_mod, dL_dimless, H0phys):

    dL_q = np.interp(z_query, z_mod, dL_dimless)
    dL_mpc = dL_q * (C_KMS / H0phys)

    mu_mod = 5.0 * np.log10(np.clip(dL_mpc, 1e-12, None)) + 25.0
    return mu_mod

def chi2_snia_profileM(z, mu_data, mu_err, z_mod, dL_dimless_mod, H0phys, cov=None):

    mu_mod = interp_mu_model(z, z_mod, dL_dimless_mod, H0phys)
    r0 = mu_data - mu_mod

    if cov is None:
        W = np.diag(1.0 / np.clip(mu_err, 1e-12, None)**2)
    else:
        W = np.linalg.inv(cov)

    ones = np.ones_like(r0)
    a = ones @ W @ ones
    b = ones @ W @ r0

    chi2 = (r0 @ W @ r0) - (b*b) / max(a, 1e-30)

    return float(chi2)