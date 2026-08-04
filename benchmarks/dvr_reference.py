"""Independent reference spectra for the adversarial cases.

Every reference here is computed by a method that shares no machinery with the
pipelines under test, so agreement is evidence rather than a tautology:

* **sinc-DVR** (Colbert & Miller, J. Chem. Phys. 96:1982, 1992) --- a spectral
  discrete variable representation, diagonalised densely.  Exponential
  convergence in the number of points; nothing in common with second-order
  finite differences plus Richardson.
* **Airy zeros** --- exact closed form for V(x) = |x|.
* the **double-well calibration** sweep, used to pick a barrier height whose
  tunnelling doublet splitting falls below the solver's bisection half-width.

Run::

    python3 dvr_reference.py
"""

from __future__ import annotations

import numpy as np
from scipy.special import ai_zeros


def sinc_dvr(V, half_width: float, n: int) -> np.ndarray:
    """Eigenvalues of -1/2 d^2/dx^2 + V on [-L, L] by sinc-DVR.

    Uniform grid, Colbert--Miller kinetic matrix, dense symmetric
    diagonalisation.  Returns all eigenvalues in ascending order.
    """
    x = np.linspace(-half_width, half_width, n + 2)[1:-1]
    h = x[1] - x[0]
    i = np.arange(n)
    d = i[:, None] - i[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        T = np.where(d == 0, 0.0, 2.0 * (-1.0) ** d / np.maximum(d ** 2, 1)) / (2 * h ** 2)
    np.fill_diagonal(T, np.pi ** 2 / (6 * h ** 2))
    return np.linalg.eigvalsh(T + np.diag(V(x)))


def airy_reference(k: int = 6) -> np.ndarray:
    """Exact spectrum of -1/2 psi'' + |x| psi = E psi.

    For x > 0 the substitution t = 2^(1/3)(x - E) gives psi = Ai(t), so the odd
    states (psi(0) = 0) sit at the Airy zeros a_n and the even states
    (psi'(0) = 0) at the zeros a'_n of Ai'.
    """
    a, ap, _, _ = ai_zeros(k)
    return np.sort(np.concatenate([-a / 2 ** (1 / 3), -ap / 2 ** (1 / 3)]))[:k]


def double_well_calibration():
    """Find a barrier for which the doublet splitting is below the bisection width.

    The engine bisects to max(0.01*tol, 2e-12) = 2e-10 at tol = 2e-8.  A
    splitting well below that is the interesting stress case.
    """
    rows = []
    for lam, a2 in ((1.0, 4.0), (1.0, 9.0), (2.0, 9.0)):
        E = sinc_dvr(lambda x, lam=lam, a2=a2: lam * (x * x - a2) ** 2, 8.0, 4000)[:6]
        rows.append((lam, a2, E, E[1] - E[0], E[3] - E[2]))
    return rows


if __name__ == "__main__":
    print("V(x) = |x|  -- exact Airy reference")
    print("  ", np.round(airy_reference(6), 10))
    print()
    print("V(x) = lam (x^2 - a2)^2  -- sinc-DVR, N = 4000, L = 8")
    for lam, a2, E, s01, s23 in double_well_calibration():
        print(f"  lam={lam} a2={a2}: E={np.round(E[:4], 10)}")
        print(f"      doublet splittings: {s01:.3e}, {s23:.3e}")
    print()
    print("high-accuracy reference used in the paper (N = 6000):")
    print("  ", np.round(sinc_dvr(lambda x: (x * x - 9.0) ** 2, 8.0, 6000)[:4], 10))
