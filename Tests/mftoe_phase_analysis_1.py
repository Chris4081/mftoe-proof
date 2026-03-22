"""
MFToE / SSP Phasenanalyse — vollständiges Script
Fixes: all_results Initialisierung, automatische Plot-Speicherung
Christof Krieg, 2026
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# ============================================================
# Ausgabe-Ordner
# ============================================================
SAVE_DIR = "mftoe_output"
os.makedirs(SAVE_DIR, exist_ok=True)


# ============================================================
# Branch-Tracking (MUSS vor Analyse definiert sein)
# ============================================================
def compute_F_reduced(lam, eta, sol):
    """
    Reduzierte freie Energie F_red(lambda, eta, m, u, v).
    ANPASSEN an dein konkretes Modell.
    Standardansatz: F = -lam*(m^2 + u^2) + eta*v^2 + m^4 + u^4 + v^4
    """
    m, u, v = sol
    return -lam * (m**2 + u**2) + eta * v**2 + m**4 + u**4 + v**4


def find_saddle_solutions(lam, eta, n_init=80, tol=1e-9):
    """
    Findet stationäre Punkte via Gradient-Descent-Sampling.
    Sucht Nullstellen von dF/d(m,u,v).
    """
    from scipy.optimize import fsolve

    def grad_F(sol):
        m, u, v = sol
        dFdm = -2 * lam * m + 4 * m**3
        dFdu = -2 * lam * u + 4 * u**3
        dFdv = 2 * eta * v + 4 * v**3
        return [dFdm, dFdu, dFdv]

    solutions = []
    rng = np.random.default_rng(42)

    for _ in range(n_init):
        x0 = rng.uniform(-2, 2, 3)
        try:
            sol = fsolve(grad_F, x0, full_output=True, xtol=tol)
            x, info, ier, msg = sol
            if ier == 1:
                residual = np.max(np.abs(grad_F(x)))
                if residual < 1e-7:
                    # Deduplizieren
                    is_new = True
                    for prev in solutions:
                        if np.max(np.abs(x - prev)) < 1e-5:
                            is_new = False
                            break
                    if is_new:
                        solutions.append(x)
        except Exception:
            pass

    return solutions


def track_branches_for_lambda(lam, eta_vals, n_init=80):
    """
    Verfolgt Lösungszweige entlang eta für gegebenes lambda.
    Gibt Liste von Branch-Dicts zurück:
      { "eta": [...], "sol": [...], "F": [...] }
    """
    # Einfache Implementierung: alle Lösungen pro eta separat finden
    # (für Branch-Kontinuität müsste man Continuation nutzen)

    # Sammle alle Lösungen pro eta
    all_eta = []
    all_sol = []
    all_F = []

    for eta in eta_vals:
        sols = find_saddle_solutions(lam, eta, n_init=n_init)
        for s in sols:
            F = compute_F_reduced(lam, eta, s)
            all_eta.append(eta)
            all_sol.append(s)
            all_F.append(F)

    if not all_eta:
        return []

    # Einfaches Branch-Clustering nach Lösung-Nähe
    branches = []
    used = [False] * len(all_eta)

    for i in range(len(all_eta)):
        if used[i]:
            continue
        branch = {
            "eta": [all_eta[i]],
            "sol": [all_sol[i]],
            "F": [all_F[i]]
        }
        used[i] = True
        # Verbinde mit nahen Punkten bei benachbarten eta
        for j in range(i + 1, len(all_eta)):
            if used[j]:
                continue
            if abs(all_eta[j] - all_eta[i]) < 2 * (eta_vals[1] - eta_vals[0]) + 1e-10:
                if np.max(np.abs(all_sol[j] - all_sol[i])) < 0.15:
                    branch["eta"].append(all_eta[j])
                    branch["sol"].append(all_sol[j])
                    branch["F"].append(all_F[j])
                    used[j] = True
        if len(branch["eta"]) >= 2:
            # Sortiere nach eta
            order = np.argsort(branch["eta"])
            branch["eta"] = [branch["eta"][k] for k in order]
            branch["sol"] = [branch["sol"][k] for k in order]
            branch["F"] = [branch["F"][k] for k in order]
            branches.append(branch)

    return branches


# ============================================================
# Hilfsfunktionen
# ============================================================
def build_eta_candidate_map(branches, eta_tol=1e-10):
    eta_map = {}
    for bid, br in enumerate(branches):
        for eta, sol, Fv in zip(br["eta"], br["sol"], br["F"]):
            eta_key = float(np.round(eta, 10))
            eta_map.setdefault(eta_key, []).append((bid, np.array(sol), float(Fv)))
    return eta_map


def global_minimum_curve(branches):
    eta_map = build_eta_candidate_map(branches)
    eta_sorted = sorted(eta_map.keys())

    out = {"eta": [], "branch_id": [], "m": [], "u": [], "v": [], "F": [], "n_candidates": []}

    for eta in eta_sorted:
        candidates = eta_map[eta]
        best = min(candidates, key=lambda x: x[2])
        bid, sol, Fv = best
        out["eta"].append(eta)
        out["branch_id"].append(bid)
        out["m"].append(sol[0])
        out["u"].append(sol[1])
        out["v"].append(sol[2])
        out["F"].append(Fv)
        out["n_candidates"].append(len(candidates))

    for k in out:
        out[k] = np.array(out[k])

    return out


def find_branch_switch_points(global_curve):
    eta = global_curve["eta"]
    bid = global_curve["branch_id"]
    switch_points = []
    for i in range(1, len(eta)):
        if bid[i] != bid[i - 1]:
            switch_points.append({
                "eta_left": eta[i - 1],
                "eta_right": eta[i],
                "branch_left": int(bid[i - 1]),
                "branch_right": int(bid[i]),
            })
    return switch_points


def interpolate_branch_F_at_eta(branch, eta_query):
    et = np.array(branch["eta"], dtype=float)
    Fv = np.array(branch["F"], dtype=float)
    sols = np.array(branch["sol"], dtype=float)

    if eta_query < et.min() or eta_query > et.max():
        return None

    return {
        "F": np.interp(eta_query, et, Fv),
        "m": np.interp(eta_query, et, sols[:, 0]),
        "u": np.interp(eta_query, et, sols[:, 1]),
        "v": np.interp(eta_query, et, sols[:, 2]),
    }


def find_pairwise_crossings(branches, eta_grid_dense=None, F_tol=5e-3):
    if eta_grid_dense is None:
        all_eta = []
        for br in branches:
            all_eta.extend(br["eta"])
        eta_min = float(np.min(all_eta))
        eta_max = float(np.max(all_eta))
        eta_grid_dense = np.linspace(eta_min, eta_max, 600)

    crossings = []

    for i in range(len(branches)):
        for j in range(i + 1, len(branches)):
            br_i = branches[i]
            br_j = branches[j]

            et_i = np.array(br_i["eta"], dtype=float)
            et_j = np.array(br_j["eta"], dtype=float)

            overlap_min = max(et_i.min(), et_j.min())
            overlap_max = min(et_i.max(), et_j.max())

            if overlap_max <= overlap_min:
                continue

            mask = (eta_grid_dense >= overlap_min) & (eta_grid_dense <= overlap_max)
            etas = eta_grid_dense[mask]
            if len(etas) < 3:
                continue

            diff = []
            valid_etas = []
            for e in etas:
                vi = interpolate_branch_F_at_eta(br_i, e)
                vj = interpolate_branch_F_at_eta(br_j, e)
                if vi is not None and vj is not None:
                    valid_etas.append(e)
                    diff.append(vi["F"] - vj["F"])

            if len(valid_etas) < 3:
                continue

            valid_etas = np.array(valid_etas)
            diff = np.array(diff)

            for k in range(1, len(valid_etas)):
                d1 = diff[k - 1]
                d2 = diff[k]

                if abs(d1) < F_tol:
                    crossings.append({
                        "eta": float(valid_etas[k - 1]),
                        "branch_i": i, "branch_j": j, "type": "near_zero"
                    })

                if d1 != 0 and d1 * d2 < 0:
                    eta1, eta2 = valid_etas[k - 1], valid_etas[k]
                    eta_cross = eta1 - d1 * (eta2 - eta1) / (d2 - d1)
                    crossings.append({
                        "eta": float(eta_cross),
                        "branch_i": i, "branch_j": j, "type": "sign_change"
                    })

    dedup = []
    for c in crossings:
        if not any(
            abs(c["eta"] - d["eta"]) < 1e-2 and
            {c["branch_i"], c["branch_j"]} == {d["branch_i"], d["branch_j"]}
            for d in dedup
        ):
            dedup.append(c)

    return dedup


def classify_point(candidates, eps_v=1e-3, F_degen=1e-3):
    Fs = np.array([c[2] for c in candidates], dtype=float)
    Fmin = Fs.min()
    near_min = [c for c in candidates if abs(c[2] - Fmin) < F_degen]
    any_frustrated = any(abs(sol[2]) > eps_v for _, sol, _ in near_min)

    if len(near_min) == 1:
        return 1 if any_frustrated else 0
    else:
        return 3 if any_frustrated else 2


def build_phase_map(all_results, lambda_vals_sorted=None, eps_v=1e-3, F_degen=1e-3):
    if lambda_vals_sorted is None:
        lambda_vals_sorted = sorted(all_results.keys())

    eta_set = set()
    for lam in lambda_vals_sorted:
        branches = all_results[lam]
        eta_map = build_eta_candidate_map(branches)
        eta_set.update(eta_map.keys())

    eta_vals_sorted = np.array(sorted(eta_set), dtype=float)
    phase_map = np.full((len(eta_vals_sorted), len(lambda_vals_sorted)), np.nan)

    for j, lam in enumerate(lambda_vals_sorted):
        branches = all_results[lam]
        eta_map = build_eta_candidate_map(branches)

        for i, eta in enumerate(eta_vals_sorted):
            eta_key = float(np.round(eta, 10))
            if eta_key not in eta_map:
                continue
            candidates = eta_map[eta_key]
            phase_map[i, j] = classify_point(candidates, eps_v=eps_v, F_degen=F_degen)

    return eta_vals_sorted, np.array(lambda_vals_sorted, dtype=float), phase_map


def extract_v_boundary(all_results, lambda_vals_sorted=None, eps_v=1e-3):
    if lambda_vals_sorted is None:
        lambda_vals_sorted = sorted(all_results.keys())

    boundary = []
    for lam in lambda_vals_sorted:
        gc = global_minimum_curve(all_results[lam])
        eta = gc["eta"]
        vabs = np.abs(gc["v"])
        is_fr = vabs > eps_v
        idx = np.where(is_fr[1:] != is_fr[:-1])[0]
        if len(idx) == 0:
            continue
        k = idx[0]
        eta_b = 0.5 * (eta[k] + eta[k + 1])
        boundary.append((lam, eta_b))

    return np.array(boundary, dtype=float) if boundary else np.empty((0, 2))


def extract_switch_boundary(all_results, lambda_vals_sorted=None):
    if lambda_vals_sorted is None:
        lambda_vals_sorted = sorted(all_results.keys())

    pts = []
    for lam in lambda_vals_sorted:
        gc = global_minimum_curve(all_results[lam])
        sw = find_branch_switch_points(gc)
        for s in sw:
            eta_b = 0.5 * (s["eta_left"] + s["eta_right"])
            pts.append((lam, eta_b))

    return np.array(pts, dtype=float) if pts else np.empty((0, 2))


# ============================================================
# Plot-Funktionen (mit automatischer Speicherung)
# ============================================================
def plot_transition_summary_for_lambda(lam, branches):
    gc = global_minimum_curve(branches)
    switches = find_branch_switch_points(gc)
    pair_cross = find_pairwise_crossings(branches)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    ax_m, ax_u, ax_v, ax_F = axes.flat

    ax_m.plot(gc["eta"], gc["m"], "k-", lw=2)
    ax_u.plot(gc["eta"], gc["u"], "k-", lw=2)
    ax_v.plot(gc["eta"], gc["v"], "k-", lw=2)
    ax_F.plot(gc["eta"], gc["F"], "k-", lw=2)

    for s in switches:
        eta_s = 0.5 * (s["eta_left"] + s["eta_right"])
        for ax in axes.flat:
            ax.axvline(eta_s, color="red", ls="--", alpha=0.8)

    for c in pair_cross:
        ax_F.axvline(c["eta"], color="blue", ls=":", alpha=0.7)

    ax_m.set_title(rf"Global minimum: $m(\eta)$, $\lambda={lam:.2f}$")
    ax_u.set_title(rf"Global minimum: $u(\eta)$, $\lambda={lam:.2f}$")
    ax_v.set_title(rf"Global minimum: $v(\eta)$, $\lambda={lam:.2f}$")
    ax_F.set_title(rf"Global minimum: $F_\mathrm{{red}}(\eta)$, $\lambda={lam:.2f}$")

    for ax, ylabel in zip(axes.flat, ["$m$", "$u$", "$v$", "$F_\\mathrm{red}$"]):
        ax.set_ylabel(ylabel)
        ax.grid(alpha=0.25)

    axes.flat[2].set_xlabel(r"$\eta$")
    axes.flat[3].set_xlabel(r"$\eta$")

    plt.tight_layout()
    fname = os.path.join(SAVE_DIR, f"transitions_lambda_{lam:.3f}.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Gespeichert: {fname}")


def plot_global_minimum_for_lambda(lam, branches):
    gc = global_minimum_curve(branches)

    fig, axes = plt.subplots(1, 3, figsize=(12, 4))

    axes[0].plot(gc["eta"], gc["m"], "b-", lw=2, label="$m$")
    axes[0].plot(gc["eta"], gc["u"], "r--", lw=2, label="$u$")
    axes[0].legend(); axes[0].set_xlabel(r"$\eta$"); axes[0].grid(alpha=0.25)
    axes[0].set_title(rf"Order params, $\lambda={lam:.2f}$")

    axes[1].plot(gc["eta"], gc["v"], "g-", lw=2)
    axes[1].set_xlabel(r"$\eta$"); axes[1].set_ylabel("$v$"); axes[1].grid(alpha=0.25)
    axes[1].set_title("v (frustration)")

    axes[2].plot(gc["eta"], gc["F"], "k-", lw=2)
    axes[2].set_xlabel(r"$\eta$"); axes[2].set_ylabel("$F_\\mathrm{red}$"); axes[2].grid(alpha=0.25)
    axes[2].set_title("Free energy")

    plt.tight_layout()
    fname = os.path.join(SAVE_DIR, f"global_lambda_{lam:.3f}.png")
    plt.savefig(fname, dpi=150)
    plt.close()
    print(f"  Gespeichert: {fname}")


def plot_clean_phase_diagram(all_results, lambda_vals_sorted=None, eps_v=1e-3, F_degen=1e-3):
    eta_vals, lambda_vals, phase_map = build_phase_map(
        all_results, lambda_vals_sorted=lambda_vals_sorted,
        eps_v=eps_v, F_degen=F_degen)

    switch_pts = extract_switch_boundary(all_results, lambda_vals_sorted=lambda_vals)
    v_boundary = extract_v_boundary(all_results, lambda_vals_sorted=lambda_vals, eps_v=eps_v)

    fig, ax = plt.subplots(figsize=(9, 7))

    cmap = plt.cm.get_cmap("tab10", 4)
    im = ax.imshow(
        phase_map, origin="lower", aspect="auto",
        extent=[lambda_vals.min(), lambda_vals.max(), eta_vals.min(), eta_vals.max()],
        cmap=cmap, vmin=0, vmax=3)

    if len(v_boundary) > 0:
        ax.plot(v_boundary[:, 0], v_boundary[:, 1], color="white", lw=2.5,
                label=r"$|v|=0$ boundary")

    if len(switch_pts) > 0:
        ax.scatter(switch_pts[:, 0], switch_pts[:, 1], s=18, c="black",
                   marker="o", label="global branch switches")

    ax.set_xlabel(r"$\lambda$")
    ax.set_ylabel(r"$\eta$")
    ax.set_title("Phase Diagram — MFToE/SSP")

    cbar = plt.colorbar(im, ax=ax, ticks=[0.5, 1.5, 2.5, 3.5])
    cbar.ax.set_yticklabels(["aligned single", "frustrated single",
                              "aligned multi", "frustrated multi"])
    ax.legend(loc="best")
    plt.tight_layout()

    fname = os.path.join(SAVE_DIR, "phase_diagram.png")
    plt.savefig(fname, dpi=200)
    plt.close()
    print(f"  Gespeichert: {fname}")


# ============================================================
# Hauptaufruf
# ============================================================
if __name__ == "__main__":

    print("=== MFToE Phasenanalyse ===\n")

    # ---- Parameter (hier anpassen) ----
    lambda_cuts = [0.0, 0.5, 1.0, 1.5, 2.0]
    eta_vals    = np.linspace(-2.0, 2.0, 60)
    # -----------------------------------

    # 1) Branch Tracking
    # all_results initialisieren — das war der Fehler vorher!
    all_results = {}

    print("--- Schritt 1: Branch Tracking ---")
    for lam in lambda_cuts:
        print(f"\n  λ = {lam:.2f} ...")
        branches = track_branches_for_lambda(lam, eta_vals)
        all_results[lam] = branches
        print(f"  Gefundene Zweige: {len(branches)}")

        plot_global_minimum_for_lambda(lam, branches)

    # 2) Übergangsanalyse pro lambda-Schnitt
    print("\n--- Schritt 2: Übergangsanalyse ---")
    for lam in sorted(all_results.keys()):
        print(f"\n  λ = {lam:.2f}")
        branches = all_results[lam]
        gc = global_minimum_curve(branches)
        switches = find_branch_switch_points(gc)
        crossings = find_pairwise_crossings(branches)
        print(f"  Globale Switches:    {len(switches)}")
        print(f"  Paarweise Kreuzungen: {len(crossings)}")
        plot_transition_summary_for_lambda(lam, branches)

    # 3) Phasendiagramm
    print("\n--- Schritt 3: Phasendiagramm ---")
    plot_clean_phase_diagram(all_results)

    print(f"\n✓ Fertig. Alle Plots in: {SAVE_DIR}/")
    print("\nDateien:")
    for f in sorted(os.listdir(SAVE_DIR)):
        print(f"  {SAVE_DIR}/{f}")
