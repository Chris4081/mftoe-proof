import pandas as pd
from snia_test import load_snia_csv, chi2_snia_profileM
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

def lcdm_dL_dimless(z, H0=70.0, Om=0.3):
    """
    Returns dimensionless luminosity distance dL in units of (c/H0)
    """
    cosmo = FlatLambdaCDM(H0=H0 * (u.km / u.s / u.Mpc), Om0=Om)
    dL_mpc = cosmo.luminosity_distance(z).value
    dL_dimless = dL_mpc * (H0 / 299792.458)
    return dL_dimless

# ---------- config ----------
MFTOE_RUN = "runs/mftoe_baseline.csv"
SN_DATA = "data/pantheon_plus.csv"
H0 = 70.0
Omega_m = 0.3

# ---------- load MFToE ----------
mf = pd.read_csv(MFTOE_RUN)
z_mod = mf["z"].values
dL_mod = mf["dL"].values

# ---------- load SN ----------
sn = load_snia_csv(SN_DATA)
z = sn["z"].values
mu = sn["mu"].values
mu_err = sn["mu_err"].values

# ---------- MFToE chi2 ----------
chi2_mftoe = chi2_snia_profileM(
    z,
    mu,
    mu_err,
    z_mod,
    dL_mod,
    H0
)

# ---------- ΛCDM reference ----------
dL_lcdm = lcdm_dL_dimless(z_mod, H0=H0, Om=Omega_m)

chi2_lcdm = chi2_snia_profileM(
    z,
    mu,
    mu_err,
    z_mod,
    dL_lcdm,
    H0
)

# ---------- output ----------
dof = len(z) - 1

print("SN Ia comparison")
print("MFToE chi2 =", chi2_mftoe)
print("LCDM  chi2 =", chi2_lcdm)
print("Δchi2 =", chi2_mftoe - chi2_lcdm)
print("dof =", dof)
print("MFToE chi2/dof =", chi2_mftoe / dof)
print("LCDM  chi2/dof =", chi2_lcdm / dof)