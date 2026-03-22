import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root


# -----------------------------
# Model parameters
# -----------------------------
kB = 1.0
T_eff = 0.8
omega = 1.0
mu_A = 1.2
g_A = 1.0

# Parameter grid
lambda_vals = np.linspace(0.0, 2.0, 81)
eta_vals = np.linspace(0.0, 2.0, 81)

# Numerical tolerances
ROOT_TOL = 1e-10
MAXFEV = 5000


# -----------------------------
# Helper functions
# -----------------------------
def entropy_binary(m: float, kB_val: float = 1.0) -> float:
    """
    Binary von Neumann entropy for rho = 1/2 (I - m sigma_z), 0 <= m < 1.
    """
    eps = 1e-14
    m_clamped = np.clip(m, 0.0, 1.0 - eps)
    p_plus = 0.5 * (1.0 + m_clamped)
    p_minus = 0.5 * (1.0 - m_clamped)
    return -kB_val * (p_plus * np.log(p_plus) + p_minus * np.log(p_minus))


def stationarity_equations(x: np.ndarray, lam: float, eta: float) -> np.ndarray:
    """
    Reduced stationary equations for variables x = [m, u, v].
    """
    m, u, v = x

    # Mild regularization to avoid log singularities during iteration
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
    """
    Reduced free energy F_red(m,u,v).
    Returns +inf for invalid states.
    """
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
    """
    Check if a candidate solution is physically and numerically acceptable.
    """
    m, u, v = sol
    if not np.isfinite(sol).all():
        return False
    if m < -1e-8 or m >= 1.0:
        return False
    residual = np.linalg.norm(stationarity_equations(sol, lam, eta), ord=2)
    return residual < tol


def deduplicate_solutions(solutions: list[np.ndarray], atol: float = 1e-5) -> list[np.ndarray]:
    """
    Remove near-duplicate solutions.
    """
    unique = []
    for s in solutions:
        if not any(np.allclose(s, t, atol=atol, rtol=0.0) for t in unique):
            unique.append(s)
    return unique


def initial_guesses(lam: float, eta: float) -> list[np.ndarray]:
    """
    Generate a small bank of initial guesses.
    Includes Gibbs start and several symmetry-related seeds.
    """
    m0 = np.tanh(omega / (2.0 * kB * T_eff))

    guesses = [
        np.array([m0, 0.0, 0.0]),
        np.array([max(0.05, 0.8 * m0), 0.2, 0.0]),
        np.array([max(0.05, 0.8 * m0), -0.2, 0.0]),
        np.array([max(0.05, 0.8 * m0), 0.2, 0.3]),
        np.array([max(0.05, 0.8 * m0), 0.2, -0.3]),
        np.array([0.2, 0.0, 0.5]),
        np.array([0.5, 0.5, 0.5]),
        np.array([0.5, -0.5, 0.5]),
    ]

    # A rough aligned estimate for weak v:
    if mu_A > 1e-12:
        u_est = eta * m0 / mu_A
        guesses.append(np.array([m0, u_est, 0.0]))

    return guesses


def solve_point(lam: float, eta: float, prev_sol: np.ndarray | None = None) -> tuple[np.ndarray | None, list[np.ndarray]]:
    """
    Solve the stationary system at one grid point.
    Returns the minimum-free-energy solution and the list of valid solutions found.
    """
    guesses = initial_guesses(lam, eta)
    if prev_sol is not None:
        guesses = [prev_sol.copy()] + guesses

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

    found = deduplicate_solutions(found)

    if not found:
        return None, []

    energies = [free_energy_reduced(sol[0], sol[1], sol[2], lam, eta) for sol in found]
    best_idx = int(np.argmin(energies))
    return found[best_idx], found


# -----------------------------
# Continuation over parameter grid
# -----------------------------
M = np.full((len(eta_vals), len(lambda_vals)), np.nan)
U = np.full_like(M, np.nan)
V = np.full_like(M, np.nan)
F = np.full_like(M, np.nan)
BRANCH_COUNT = np.zeros_like(M)

# Scan eta row by row, lambda column by column
for i, eta in enumerate(eta_vals):
    prev = None
    for j, lam in enumerate(lambda_vals):
        best_sol, all_solutions = solve_point(lam, eta, prev_sol=prev)

        if best_sol is not None:
            m, u, v = best_sol
            M[i, j] = m
            U[i, j] = u
            V[i, j] = v
            F[i, j] = free_energy_reduced(m, u, v, lam, eta)
            BRANCH_COUNT[i, j] = len(all_solutions)
            prev = best_sol
        else:
            prev = None


# -----------------------------
# Plot 1: heatmap m(lambda, eta)
# -----------------------------
fig, ax = plt.subplots(figsize=(8, 6))
im = ax.imshow(
    M,
    origin="lower",
    aspect="auto",
    extent=[lambda_vals.min(), lambda_vals.max(), eta_vals.min(), eta_vals.max()],
    cmap="viridis"
)
cbar = plt.colorbar(im, ax=ax)
cbar.set_label(r"$m(\lambda,\eta)$")

ax.set_xlabel(r"$\lambda$")
ax.set_ylabel(r"$\eta$")
ax.set_title(r"Longitudinal order parameter $m(\lambda,\eta)$")
plt.tight_layout()
plt.show()


# -----------------------------
# Plot 2: transverse branch amplitude |v|
# -----------------------------
fig, ax = plt.subplots(figsize=(8, 6))
im = ax.imshow(
    np.abs(V),
    origin="lower",
    aspect="auto",
    extent=[lambda_vals.min(), lambda_vals.max(), eta_vals.min(), eta_vals.max()],
    cmap="magma"
)
cbar = plt.colorbar(im, ax=ax)
cbar.set_label(r"$|v(\lambda,\eta)|$")

ax.set_xlabel(r"$\lambda$")
ax.set_ylabel(r"$\eta$")
ax.set_title(r"Transverse structural amplitude $|v|$")
plt.tight_layout()
plt.show()


# -----------------------------
# Plot 3: number of solutions found
# -----------------------------
fig, ax = plt.subplots(figsize=(8, 6))
im = ax.imshow(
    BRANCH_COUNT,
    origin="lower",
    aspect="auto",
    extent=[lambda_vals.min(), lambda_vals.max(), eta_vals.min(), eta_vals.max()],
    cmap="cividis"
)
cbar = plt.colorbar(im, ax=ax)
cbar.set_label("number of stationary branches found")

ax.set_xlabel(r"$\lambda$")
ax.set_ylabel(r"$\eta$")
ax.set_title("Multiplicity of stationary solutions")
plt.tight_layout()
plt.show()