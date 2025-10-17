#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
MFToE (Single-Field) + dynamic vacuum field χ
Extended with optional Astropy cosmology reference (ΛCDM baseline)

Modes: off | sequester | targetH0 | relaxion
Options: --rg on/off, --noise on/off, --astropy on/off/auto
Time variable: N = ln(a), integrated from N_start to N_end
"""

import argparse, json, numpy as np, matplotlib.pyplot as plt
from pathlib import Path
from analysis.desi_bestfit import load_desi_bestfit_minimum
# ---------- Astropy Cosmology (optional) ----------
USE_ASTROPY = True
try:
    import astropy.units as u
    from astropy.cosmology import FlatLambdaCDM
except Exception:
    USE_ASTROPY = False

# ---------- Helpers ----------
def cumtrapz(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float)
    out = np.zeros_like(x, float)
    if len(x) > 1:
        dx = np.diff(x); mid = 0.5*(y[1:] + y[:-1])
        out[1:] = np.cumsum(dx * mid)
    return out

def safe_pos(x, eps=1e-12):
    return np.clip(x, eps, None)

# ---------- ΛCDM via Astropy (or fallback) ----------
C_KMS = 299_792.458  # km/s

def lcdm_curves_with_astropy(z, H0_phys=70.0, Om0=0.3, Or0=0.0):
    """
    Returns (H_LCDM/H0_phys, dL_LCDM_dimless).
    dL_dimless = dL_Mpc * (H0_phys / c)  -> same units as MFToE (c/H0).
    """
    z = np.asarray(z, float)
    if USE_ASTROPY:
        cosmo = FlatLambdaCDM(H0=H0_phys * (u.km/u.s/u.Mpc), Om0=Om0, Tcmb0=2.725*u.K)
        Hz = cosmo.H(z).to_value(u.km/u.s/u.Mpc)          # [km/s/Mpc]
        Hnorm = Hz / H0_phys                              # H/H0
        dL_Mpc = cosmo.luminosity_distance(z).to_value(u.Mpc)
        dL_dimless = dL_Mpc * (H0_phys / C_KMS)           # dimensionless (c/H0 units)
        return Hnorm, dL_dimless
    else:
        # fallback: flat ΛCDM in dimensionless units
        Ol = 1.0 - Om0 - Or0
        Hnorm = np.sqrt(Or0*(1+z)**4 + Om0*(1+z)**3 + Ol)  # H/H0
        chi = cumtrapz(z, 1.0/safe_pos(Hnorm))
        dL_dimless = (1.0 + z) * chi
        return Hnorm, dL_dimless

# ---------- OU noise ----------
class OU:
    def __init__(self, tau=0.5, sigma=1e-4, x0=0.0):
        self.tau=tau; self.sigma=sigma; self.x=x0
    def step(self, h):
        drift = -self.x/self.tau
        diffusion = np.sqrt(max(2*self.sigma**2/self.tau, 0.0))
        self.x += drift*h + diffusion*np.sqrt(max(h,1e-12))*np.random.randn()
        return self.x

# ---------- Potential ----------
def Vphi(phi, m2, lam): return 0.5*m2*phi**2 + lam*phi**4
def dVphi_dphi(phi, m2, lam): return m2*phi + 4.0*lam*phi**3

# ---------- Simple RG flow ----------
def rg_flow(m2, lam, H, a_beta=0.0, b_beta=0.0, mu_ref=1.0):
    mu = max(H, 1e-6)
    t = np.log(mu/mu_ref)
    return m2 + a_beta*m2*t, lam + b_beta*lam**2*t

# ---------- RHS builder ----------
def make_rhs(P):
    Om, Or = P["Omega_m"], P["Omega_r"]
    V0, g = P["V0_bare"], P["g"]
    mode, gamma, Mpen, rho_target = P["mode"], P["gamma"], P["Mpen"], P["rho_target"]
    ou, rg_on = P.get("ou", None), P["rg_on"]
    a_beta, b_beta, mu_ref = P["a_beta"], P["b_beta"], P["mu_ref"]

    def rhs(N, y):
        phi, pi, chi = y
        a = np.exp(N)
        rho_m, rho_r = Om*a**-3, Or*a**-4

        # RG at previous H
        H_alt = P.get("H_last", 1.0)
        if rg_on:
            m2_eff, lam_eff = rg_flow(P["m2"], P["lam"], H_alt, a_beta, b_beta, mu_ref)
        else:
            m2_eff, lam_eff = P["m2"], P["lam"]

        V_field = Vphi(phi, m2_eff, lam_eff)
        V_noise = ou.step(P["h"]) if ou else 0.0
        V_eff = V0 + g*chi + V_field + V_noise

        denom = max(1.0 - 0.5*pi*pi, 1e-12)
        H2 = (rho_m + rho_r + V_eff) / denom
        H = np.sqrt(max(H2, 1e-30))
        P["H_last"] = H  # for next RG

        # epsilon for damping
        KE = 0.5 * H2 * pi*pi
        rho_phi, p_phi = KE + V_field + V0 + g*chi, KE - (V_field + V0 + g*chi)
        rho_tot, p_tot = rho_m + rho_r + rho_phi, (1/3)*rho_r + p_phi
        eps = 1.5*(1 + p_tot/max(rho_tot, 1e-30))

        # equations
        dphi = pi
        dpi  = -(3.0 - eps)*pi - dVphi_dphi(phi, m2_eff, lam_eff)/H2

        diff = V_eff - rho_target
        if mode == "off":
            dchi = 0.0
        elif mode in ("sequester","targetH0"):
            dchi = -gamma * diff * g / max(Mpen**4,1e-24) / max(H,1e-6)
        elif mode == "relaxion":
            dchi = -gamma * (diff*g/max(Mpen**4,1e-24) + 0.1*phi*dVphi_dphi(phi,m2_eff,lam_eff)) / max(H,1e-6)
        else:
            dchi = 0.0

        return np.array([dphi, dpi, dchi], float)
    return rhs

# ---------- RK4 ----------
def rk4_step(f, N, y, h):
    k1=f(N, y)
    k2=f(N+0.5*h, y+0.5*h*k1)
    k3=f(N+0.5*h, y+0.5*h*k2)
    k4=f(N+h, y+h*k3)
    return y + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4)

# ---------- Main ----------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", default="", help="optional JSON file")
    ap.add_argument("--mode", default="targetH0", choices=["off","sequester","targetH0","relaxion"])
    ap.add_argument("--rg", default="off", choices=["on","off"])
    ap.add_argument("--noise", default="off", choices=["on","off"])
    ap.add_argument("--astropy", default="auto", choices=["auto","on","off"])
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", default="mftoe_vacuum_astropy", help="output prefix (CSV & PNG)")
    ap.add_argument("--H0phys", type=float, default=70.0,
                    help="physical H0 [km/s/Mpc] for Astropy ΛCDM baseline")
    ap.add_argument(
        "--desi-bestfit",
        default="",
        help="Relative path to DESI iminuit bestfit.minimum(.txt), e.g. data/desi_dr2/iminuit/base/desi-bao-all/bestfit.minimum"
    )
    
    args = ap.parse_args()
    np.random.seed(args.seed)

    # Defaults
    P = dict(
        m2=1e-6, lam=1e-6, V0_bare=0.7, g=1.0, gamma=0.1, Mpen=1.0,
        rho_target=0.0, Omega_m=0.3, Omega_r=0.0, H0=1.0,
        phi0=3.0, pi0=1e-4, chi0=0.0,
        N_start=0.0, N_end=-np.log(4.0), N_steps=2000,
        a_beta=0.0, b_beta=0.0, mu_ref=1.0,
        tau_noise=0.5, sigma_noise=1e-4
    )

    if args.params:
        with open(args.params,'r') as f:
            P.update(json.load(f))

    # step size
    P["h"] = (P["N_end"] - P["N_start"]) / max(P["N_steps"] - 1, 1)

    # OU init (only RHS steps it; history only reads ou.x)
    ou = OU(P["tau_noise"], P["sigma_noise"], 0.0) if args.noise=="on" else None
    P["ou"] = ou

    # Mode & RG flags
    P["mode"] = args.mode
    P["rg_on"] = (args.rg == "on")

    # Auto rho_target for targetH0 (H(N=0) match)
    if P["mode"] == "targetH0":
        pi0 = P["pi0"]
        denom0 = max(1.0 - 0.5*pi0*pi0, 1e-12)
        a0 = 1.0
        rho_m0 = P["Omega_m"] * a0**-3
        rho_r0 = P["Omega_r"] * a0**-4
        V_field0 = Vphi(P["phi0"], P["m2"], P["lam"])
        target_sum = (P["H0"]*P["H0"]) * denom0
        Veff_target0 = max(target_sum - rho_m0 - rho_r0, 0.0)
        P["rho_target"] = Veff_target0

  
    # ---- DESI bestfit hook (new parser) ----
    if args.desi_bestfit:
        bf = load_desi_bestfit_minimum(args.desi_bestfit)
        P["Omega_m"] = bf["Omega_m"]
        args.H0phys  = bf["H0"]  # physikalisches H0 (km/s/Mpc)
        print(f"[DESI] Using Ωm={P['Omega_m']:.6f}, H0={args.H0phys:.3f} km/s/Mpc, r_d={bf['rdrag']:.3f} Mpc")

    # ---- DESI bestfit loader (inline fallback) ----
    def load_desi_bestfit(path):
        from pathlib import Path
        txt = Path(path).read_text().splitlines()
        header = values = None
        for line in txt:
            if line.strip().startswith("#"):
                header = line.lstrip("#").strip().split()
            elif line.strip():
                values = line.strip().split()
                break
        if not header or not values:
            raise ValueError(f"Could not parse DESI bestfit file: {path}")
        row = {name: float(val) for name, val in zip(header, values)}
        Om = row.get("omegam", row.get("omm"))
        H0r = row.get("H0rdrag"); rd = row.get("rdrag")
        if Om is None or H0r is None or rd is None:
            raise KeyError("Missing columns omegam/omm, H0rdrag, rdrag")
        return {"Omega_m": Om, "H0": H0r/rd, "rdrag": rd, "H0rdrag": H0r}

        bf = load_desi_bestfit(args.desi_bestfit)
        # set ΛCDM baseline params from DESI bestfit
        P["Omega_m"] = bf["Omega_m"]
        # This H0 feeds the astropy baseline (ΛCDM); MFToE still uses H0=1 internally.
        args.H0phys = bf["H0"]
        print(f"[DESI] Using Ωm={P['Omega_m']:.4f}, H0={args.H0phys:.2f} km/s/Mpc from {args.desi_bestfit}")

    # Astropy toggle
    global USE_ASTROPY
    if args.astropy == "on":
        USE_ASTROPY = True
    elif args.astropy == "off":
        USE_ASTROPY = False

    # RHS
    rhs = make_rhs(P)

    # Integration
    N = np.linspace(P["N_start"], P["N_end"], P["N_steps"])
    y = np.array([P["phi0"], P["pi0"], P["chi0"]], float)
    h = P["h"]

    phi_hist = np.zeros(len(N))
    pi_hist  = np.zeros(len(N))
    chi_hist = np.zeros(len(N))
    H_hist   = np.zeros(len(N))
    wphi_hist= np.zeros(len(N))
    wtot_hist= np.zeros(len(N))

    for i in range(len(N)):
        phi_hist[i], pi_hist[i], chi_hist[i] = y

        # background
        a = np.exp(N[i])
        rho_m = P["Omega_m"] * a**-3
        rho_r = P["Omega_r"] * a**-4

        # same effective params as RHS (no second OU step!)
        if P["rg_on"]:
            Hold = H_hist[i-1] if i > 0 else 1.0
            m2_eff, lam_eff = rg_flow(P["m2"], P["lam"], Hold, a_beta=P["a_beta"], b_beta=P["b_beta"], mu_ref=P["mu_ref"])
        else:
            m2_eff, lam_eff = P["m2"], P["lam"]

        V_field = Vphi(y[0], m2_eff, lam_eff)
        V_noise = (P["ou"].x if P["ou"] is not None else 0.0)  # read-only
        V_eff   = P["V0_bare"] + P["g"]*y[2] + V_field + V_noise

        denom = max(1.0 - 0.5*y[1]**2, 1e-12)
        H2 = (rho_m + rho_r + V_eff) / denom
        H  = np.sqrt(max(H2, 1e-30))
        H_hist[i] = H
        P["H_last"] = H  # keep RG in sync

        KE = 0.5 * H2 * y[1]**2
        # w_phi: for the scalar field potential only
        wphi_hist[i] = (KE - V_field) / (KE + V_field) if (KE + V_field) > 1e-30 else -1.0
        # w_tot: includes vacuum contributions
        rho_phi_tot  = KE + (V_field + P["V0_bare"] + P["g"]*y[2])
        p_phi_tot    = KE - (V_field + P["V0_bare"] + P["g"]*y[2])
        rho_tot = rho_m + rho_r + rho_phi_tot
        p_tot  = (1/3.0)*rho_r + p_phi_tot
        wtot_hist[i] = p_tot / max(rho_tot, 1e-30)   # ← fixed

        if i < len(N)-1:
            y = rk4_step(rhs, N[i], y, h)

    # Observables
    z = np.exp(-N) - 1.0
    order = np.argsort(z)
    z = z[order]
    Hn = H_hist[order] / max(P["H0"], 1e-30)
    wphi = wphi_hist[order]
    wtot = wtot_hist[order]

    chi_int = cumtrapz(z, 1.0/safe_pos(Hn))
    dL  = (1.0 + z) * chi_int

    H_l, dL_l = lcdm_curves_with_astropy(z, H0_phys=args.H0phys, Om0=P["Omega_m"], Or0=P["Omega_r"])

    def nearest(x, y, xq): return y[np.abs(x - xq).argmin()]

    print("z | ΔH/H_LCDM | ΔdL/dL_LCDM")
    for zq in [0.2, 0.5, 1.0, 2.0, 3.0]:
        dH = (nearest(z,Hn,zq) - nearest(z,H_l,zq)) / nearest(z,H_l,zq)
        dd = (nearest(z,dL,zq) - nearest(z,dL_l,zq)) / max(nearest(z,dL_l,zq), 1e-30)
        print(f"{zq:.2f} | {dH*100:+.3f}% | {dd*100:+.3f}%")

    print(f"\nMean w_phi={np.mean(wphi):+.3f} | Mean w_tot={np.mean(wtot):+.3f}")

    # Save CSV (compatible with analysis/compare_runs.py)
    csv_path = f"{args.out}.csv"
    with open(csv_path, "w") as f:
        f.write("z,H_over_H0,dL,w_phi,w_tot\n")
        for zi, Hi, dLi, wpi, wti in zip(z, Hn, dL, wphi, wtot):
            f.write(f"{zi:.9e},{Hi:.9e},{dLi:.9e},{wpi:.9e},{wti:.9e}\n")

    # Plots
    fig, axs = plt.subplots(2, 3, figsize=(15, 10))
    axs[0,0].plot(z, wphi); axs[0,0].set_title("w_φ(z)")
    axs[0,1].plot(z, Hn, label="MFToE"); axs[0,1].plot(z, H_l, "--", label="ΛCDM"); axs[0,1].set_title("H/H0"); axs[0,1].legend()
    axs[0,2].plot(z, dL, label="MFToE"); axs[0,2].plot(z, dL_l, "--", label="ΛCDM"); axs[0,2].set_title("d_L (c/H0 units)"); axs[0,2].legend()
    axs[1,0].plot(z, (Hn - H_l)/safe_pos(H_l) * 100); axs[1,0].set_title("ΔH/H_ΛCDM [%]"); axs[1,0].set_xlabel("z")
    axs[1,1].plot(z, (dL - dL_l)/safe_pos(dL_l) * 100); axs[1,1].set_title("Δd_L/d_L_ΛCDM [%]"); axs[1,1].set_xlabel("z")
    kev = (1.0 + wphi) / np.clip(1.0 - wphi, 1e-30, None)
    axs[1,2].plot(z, np.log10(np.clip(kev, 1e-30, None))); axs[1,2].set_title("log10(KE/V)"); axs[1,2].set_xlabel("z")
    for ax in axs.ravel(): ax.grid(True, alpha=0.3)
    plt.tight_layout()
    png_path = f"{args.out}.png"
    plt.savefig(png_path, dpi=160)
    plt.show()

    print(f"→ Saved: {csv_path}, {png_path}")

# ---------- DESI bestfit parser (iminuit bestfit.minimum) ----------
def load_desi_bestfit_minimum(path):
    """
    Parse DESI iminuit 'bestfit.minimum' text file and return:
      {
        "Omega_m": float,
        "rdrag":   float,   # [Mpc]
        "H0":      float,   # [km/s/Mpc]
        "H0rdrag": float    # [km/s]
      }
    """
    import re
    from pathlib import Path

    # Regex: optional leading index, then float, then token name
    num = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?'
    line_re = re.compile(r'^\s*(?:\d+)?\s*(' + num + r')\s+([A-Za-z0-9_][A-Za-z0-9_\-]*)\b')

    vals = {}
    for raw in Path(path).read_text().splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        m = line_re.match(s)
        if not m:
            continue
        value = float(m.group(1))
        key   = m.group(2).lower()   # normalize
        # keep last occurrence (files often list some params twice)
        vals[key] = value

    # Map aliases & derive missing pieces
    Om     = vals.get("omegam", vals.get("omm"))
    rdrag  = vals.get("rdrag")
    H0     = vals.get("h0")
    H0r    = vals.get("h0rdrag")

    if H0 is None and (H0r is not None) and (rdrag is not None) and rdrag != 0.0:
        H0 = H0r / rdrag

    missing = [k for k,v in {"omegam/omm":Om, "rdrag":rdrag, "H0/H0rdrag":(H0 or H0r)}.items() if v is None]
    if missing:
        have = ", ".join(sorted(vals.keys()))
        raise ValueError(f"Missing keys {missing} in DESI bestfit file. Parsed keys: {have}")

    if H0r is None:
        H0r = H0 * rdrag

    return {"Omega_m": float(Om), "rdrag": float(rdrag), "H0": float(H0), "H0rdrag": float(H0r)}

if __name__ == "__main__":
    main()