from dataclasses import dataclass
import camb

@dataclass
class EarlyParams:
    H0: float           # km/s/Mpc
    ombh2: float        # Ω_b h^2
    omch2: float        # Ω_c h^2
    Neff: float = 3.046
    Yp: float = 0.245
    mnu_eV: float = 0.06
    Tcmb_K: float = 2.7255

def rd_fixed(value_mpc: float) -> float:
    return float(value_mpc)

def rd_camb(p: EarlyParams) -> float:
    pars = camb.CAMBparams()
    pars.set_cosmology(
        H0=p.H0, ombh2=p.ombh2, omch2=p.omch2,
        mnu=p.mnu_eV, YHe=p.Yp, TCMB=p.Tcmb_K,
        num_massive_neutrinos=1
    )
    pars.InitPower.set_params()
    results = camb.get_results(pars)
    rd = results.get_derived_params()['rdrag']  # Mpc
    return float(rd)