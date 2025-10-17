# ---------- DESI bestfit parser ----------
def load_desi_bestfit_minimum(path):
    """
    Parse DESI iminuit bestfit.minimum text file.

    Returns dict with:
      Omega_m  : float
      rdrag    : float  [Mpc]
      H0       : float  [km/s/Mpc]
      H0rdrag  : float  [km/s]
    """
    import re
    from pathlib import Path

    txt = Path(path).read_text().splitlines()
    vals = {}
    line_re = re.compile(r"^\s*\d+\s+([+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)\s+([A-Za-z0-9_]+)")

    for line in txt:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        m = line_re.match(line)
        if not m:
            continue
        value_str, name = m.groups()
        try:
            value = float(value_str)
        except ValueError:
            continue
        vals[name.lower()] = value

    Om = vals.get("omegam", vals.get("omm"))
    if Om is None:
        raise KeyError("Could not find 'omegam' or 'omm' in DESI bestfit file.")
    if "rdrag" not in vals:
        raise KeyError("Could not find 'rdrag' in DESI bestfit file.")
    rdrag = vals["rdrag"]

    if "h0" in vals:
        H0 = vals["h0"]
    elif "h0rdrag" in vals:
        H0 = vals["h0rdrag"] / rdrag
    else:
        raise KeyError("Could not find 'H0' nor 'H0rdrag' in DESI bestfit file.")

    H0r = vals.get("h0rdrag", H0 * rdrag)
    return {"Omega_m": Om, "rdrag": rdrag, "H0": H0, "H0rdrag": H0r}