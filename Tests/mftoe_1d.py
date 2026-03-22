import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root


# ============================================================
# Modellparameter
# ============================================================
kB = 1.0
T_eff = 0.8
omega = 1.0
mu_A = 1.2
g_A = 1.0

# 1D-Schnitte bei festen lambda-Werten
lambda_cuts = [0.0, 0.5, 1.0, 1.5]
eta_vals = np.linspace(0.0, 2.0, 121)

ROOT_TOL = 1e-9
MAXFEV = 6000

# Toleranzen fürs Branch-Matching
DEDUP_ATOL = 1e-5
MATCH_DIST = 0.18


# ============================================================
# Hilfsfunktionen
# ============================================================
def entropy_binary(m: float, kB_val: float = 1.0) -> float:
    eps = 1e-14
    m_clamped = np.clip(m, 0.0, 1.0 - eps)
    p_plus = 0.5 * (1.0 + m_clamped)
    p_minus = 0.5 * (1.0 - m_clamped)
    return -kB_val * (p_plus * np.log(p_plus) + p_minus * np.log(p_minus))


def stationarity_equations(x: np.ndarray, lam: float, eta: float) -> np.ndarray:
    """
    Reduziertes stationäres System:
    x = [m, u, v]
    """
    m, u, v = x
    eps = 1e-12
    m_reg = np.clip(m, eps, 1.0 - eps)

    eq1 = (
        -0.5 * omega
        + 0.5 * kB * T_eff * np.log((1.0 + m_reg) / (1.0 - m_reg))
        + 4.0 * lam * (v ** 2) * m_reg
        - eta * u
    )

    eq2 = -eta * m - mu_A * u + g_A * (u ** 2 + v ** 2) * u
    eq3 = v * (4.0 * lam * (m ** 2) - mu_A + g_A * (u ** 2 + v ** 2))

    return np.array([eq1, eq2, eq3], dtype=float)


def free_energy_reduced(m: float, u: float, v: float, lam: float, eta: float) -> float:
    if not (0.0 <= m < 1.0):
        return np.inf

    return (
        -0.5 * omega * m
        - T_eff * entropy_binary(m, kB)
        + 2.0 * lam * (v ** 2) * (m ** 2)
        - eta * u * m
        - 0.5 * mu_A * (u ** 2 + v ** 2)
        + 0.25 * g_A * (u ** 2 + v ** 2) ** 2
    )


def valid_solution(sol: np.ndarray, lam: float, eta: float, tol: float = 1e-7) -> bool:
    m, u, v = sol
    if not np.isfinite(sol).all():
        return False
    if m < -1e-8 or m >= 1.0:
        return False
    residual = np.linalg.norm(stationarity_equations(sol, lam, eta), ord=2)
    return residual < tol


def deduplicate_solutions(solutions: list[np.ndarray], atol: float = DEDUP_ATOL) -> list[np.ndarray]:
    unique = []
    for s in solutions:
        if not any(np.allclose(s, t, atol=atol, rtol=0.0) for t in unique):
            unique.append(s)
    return unique


def initial_guesses(lam: float, eta: float) -> list[np.ndarray]:
    """
    Multi-Start-Seeds.
    """
    m0 = np.tanh(omega / (2.0 * kB * T_eff))
    seeds = [
        np.array([m0, 0.0, 0.0]),
        np.array([0.15, 0.0, 0.0]),
        np.array([0.30, 0.2, 0.0]),
        np.array([0.30, -0.2, 0.0]),
        np.array([0.30, 0.2, 0.4]),
        np.array([0.30, 0.2, -0.4]),
        np.array([0.55, 0.4, 0.6]),
        np.array([0.55, -0.4, 0.6]),
        np.array([0.55, 0.4, -0.6]),
        np.array([0.75, 0.0, 0.5]),
    ]

    if mu_A > 1e-12:
        u_est = eta * m0 / mu_A
        seeds.append(np.array([m0, u_est, 0.0]))

    return seeds


def solve_point_all(lam: float, eta: float, extra_guesses: list[np.ndarray] | None = None) -> list[np.ndarray]:
    """
    Finde alle (numerisch auffindbaren) stationären Lösungen an einem Punkt.
    """
    guesses = initial_guesses(lam, eta)
    if extra_guesses:
        guesses = list(extra_guesses) + guesses

    found = []
    for guess in guesses:
        result = root(
            stationarity_equations,
            guess,
            args=(lam, eta),
            method="hybr",
            options={"maxfev": MAXFEV}
        )
        if result.success and valid_solution(result.x, lam, eta):
            found.append(result.x)

    return deduplicate_solutions(found)


def solution_distance(a: np.ndarray, b: np.ndarray) -> float:
    """
    Gewichtete Distanz zum Verknüpfen von Zweigen.
    """
    dm = 1.2 * abs(a[0] - b[0])
    du = 1.0 * abs(a[1] - b[1])
    dv = 1.0 * abs(a[2] - b[2])
    return np.sqrt(dm * dm + du * du + dv * dv)


# ============================================================
# Branch Tracking
# ============================================================
def track_branches_for_lambda(lam: float, eta_grid: np.ndarray):
    """
    Verfolge alle gefundenen Zweige entlang eines festen lambda-Schnitts.
    """
    branches = []  # Liste von Dicts: {"eta": [...], "sol": [...], "F": [...]}

    prev_solutions = []

    for i, eta in enumerate(eta_grid):
        extra_guesses = []

        # Fortsetzung: nutze vorige Lösungen als Seeds
        for s in prev_solutions:
            extra_guesses.append(s.copy())

        # zusätzlich Spiegelung v -> -v helfen
        for s in prev_solutions:
            mirrored = s.copy()
            mirrored[2] *= -1.0
            extra_guesses.append(mirrored)

        current_solutions = solve_point_all(lam, eta, extra_guesses=extra_guesses)

        # Falls erster Punkt: neue Zweige starten
        if i == 0:
            for sol in current_solutions:
                branches.append({
                    "eta": [eta],
                    "sol": [sol],
                    "F": [free_energy_reduced(sol[0], sol[1], sol[2], lam, eta)]
                })
            prev_solutions = current_solutions
            continue

        # Tracken: jede neue Lösung bestmöglich an existierende offene Zweige anhängen
        assigned_branch_ids = set()
        unmatched_solutions = []

        # Kandidaten: nur Zweige, die am vorigen eta-Wert endeten
        open_branch_ids = []
        for bid, br in enumerate(branches):
            if len(br["eta"]) > 0 and np.isclose(br["eta"][-1], eta_grid[i - 1]):
                open_branch_ids.append(bid)

        for sol in current_solutions:
            best_bid = None
            best_dist = np.inf

            for bid in open_branch_ids:
                if bid in assigned_branch_ids:
                    continue
                last_sol = branches[bid]["sol"][-1]
                d = solution_distance(sol, last_sol)
                if d < best_dist:
                    best_dist = d
                    best_bid = bid

            if best_bid is not None and best_dist < MATCH_DIST:
                branches[best_bid]["eta"].append(eta)
                branches[best_bid]["sol"].append(sol)
                branches[best_bid]["F"].append(
                    free_energy_reduced(sol[0], sol[1], sol[2], lam, eta)
                )
                assigned_branch_ids.add(best_bid)
            else:
                unmatched_solutions.append(sol)

        # Nicht gematchte Lösungen starten neue Zweige
        for sol in unmatched_solutions:
            branches.append({
                "eta": [eta],
                "sol": [sol],
                "F": [free_energy_reduced(sol[0], sol[1], sol[2], lam, eta)]
            })

        prev_solutions = current_solutions

    # Kurze Zweige wegfiltern
    filtered = []
    for br in branches:
        if len(br["eta"]) >= 4:
            filtered.append(br)

    return filtered


# ============================================================
# Plot-Funktion
# ============================================================
def plot_branches_for_lambda(lam: float, eta_grid: np.ndarray, branches: list[dict]):
    """
    Erzeuge 4er-Panel: m(eta), u(eta), v(eta), F(eta)
    """
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    ax_m, ax_u, ax_v, ax_F = axes.flat

    cmap = plt.cm.tab10
    for i, br in enumerate(branches):
        color = cmap(i % 10)
        et = np.array(br["eta"])
        sols = np.array(br["sol"])
        Fvals = np.array(br["F"])

        mvals = sols[:, 0]
        uvals = sols[:, 1]
        vvals = sols[:, 2]

        label = f"Zweig {i+1}"

        ax_m.plot(et, mvals, color=color, lw=2, label=label)
        ax_u.plot(et, uvals, color=color, lw=2)
        ax_v.plot(et, vvals, color=color, lw=2)
        ax_F.plot(et, Fvals, color=color, lw=2)

    ax_m.set_title(rf"$m(\eta)$ bei festem $\lambda={lam}$")
    ax_u.set_title(rf"$u(\eta)$ bei festem $\lambda={lam}$")
    ax_v.set_title(rf"$v(\eta)$ bei festem $\lambda={lam}$")
    ax_F.set_title(rf"$F_{{\mathrm{{red}}}}(\eta)$ bei festem $\lambda={lam}$")

    ax_m.set_ylabel(r"$m$")
    ax_u.set_ylabel(r"$u$")
    ax_v.set_ylabel(r"$v$")
    ax_F.set_ylabel(r"$F_{\mathrm{red}}$")

    ax_v.set_xlabel(r"$\eta$")
    ax_F.set_xlabel(r"$\eta$")

    ax_m.grid(alpha=0.25)
    ax_u.grid(alpha=0.25)
    ax_v.grid(alpha=0.25)
    ax_F.grid(alpha=0.25)

    ax_m.legend(fontsize=9, loc="best")
    plt.tight_layout()
    plt.show()


# ============================================================
# Optional: Plot des globalen Minimums entlang des Schnitts
# ============================================================
def plot_global_minimum_for_lambda(lam: float, eta_grid: np.ndarray, branches: list[dict]):
    """
    Zeigt den global energetisch bevorzugten Zweig entlang eta.
    """
    eta_to_candidates = {}

    for bid, br in enumerate(branches):
        for eta, sol, Fv in zip(br["eta"], br["sol"], br["F"]):
            eta_to_candidates.setdefault(float(eta), []).append((bid, sol, Fv))

    eta_sorted = sorted(eta_to_candidates.keys())
    m_best, u_best, v_best, F_best = [], [], [], []

    for eta in eta_sorted:
        candidates = eta_to_candidates[eta]
        best = min(candidates, key=lambda x: x[2])
        _, sol, Fv = best
        m_best.append(sol[0])
        u_best.append(sol[1])
        v_best.append(sol[2])
        F_best.append(Fv)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
    ax_m, ax_u, ax_v, ax_F = axes.flat

    ax_m.plot(eta_sorted, m_best, "k-", lw=2.5)
    ax_u.plot(eta_sorted, u_best, "k-", lw=2.5)
    ax_v.plot(eta_sorted, v_best, "k-", lw=2.5)
    ax_F.plot(eta_sorted, F_best, "k-", lw=2.5)

    ax_m.set_title(rf"Globales Minimum: $m(\eta)$ bei $\lambda={lam}$")
    ax_u.set_title(rf"Globales Minimum: $u(\eta)$ bei $\lambda={lam}$")
    ax_v.set_title(rf"Globales Minimum: $v(\eta)$ bei $\lambda={lam}$")
    ax_F.set_title(rf"Globales Minimum: $F_{{\mathrm{{red}}}}(\eta)$ bei $\lambda={lam}$")

    ax_m.set_ylabel(r"$m$")
    ax_u.set_ylabel(r"$u$")
    ax_v.set_ylabel(r"$v$")
    ax_F.set_ylabel(r"$F_{\mathrm{red}}$")

    ax_v.set_xlabel(r"$\eta$")
    ax_F.set_xlabel(r"$\eta$")

    for ax in axes.flat:
        ax.grid(alpha=0.25)

    plt.tight_layout()
    plt.show()


# ============================================================
# Hauptlauf
# ============================================================
if __name__ == "__main__":
    all_results = {}

    for lam in lambda_cuts:
        print(f"\n=== Tracking branches for lambda = {lam:.2f} ===")
        branches = track_branches_for_lambda(lam, eta_vals)
        all_results[lam] = branches
        print(f"Gefundene Zweige: {len(branches)}")

        plot_branches_for_lambda(lam, eta_vals, branches)
        plot_global_minimum_for_lambda(lam, eta_vals, branches)