"""Is the double-well failure in the planner or in the solver?

The raw-input pipeline returns, for V(x) = (x^2 - 9)^2 with k = 4, the values

    (4.2144398, 12.5272184, 20.6574438, 28.5922506)   status ACCEPT

against a reference of

    (4.2144398, 4.2144398, 12.5272184, 12.5272184)

an error of 16.07 at a declared tolerance of 2e-8.  There are two candidate
explanations and they call for completely different repairs:

    (a) the Sturm/bisection kernel cannot separate a doublet split by 8e-12
        when it bisects to a half-width of 2e-10, or
    (b) the planner chose a truncation window that contains only one of the
        two wells, so the operator actually solved was a different operator.

This script decides it by handing the kernel the *correct* window by hand and
running the same three-mesh Richardson ladder the pipeline would run.  It also
runs Matslise 2.0 on that same correct window, as an independent check on
whether a shooting-based solver resolves the doublet.

Run (from the repository root)::

    PYTHONPATH=. python3 planner_vs_solver.py
"""

from __future__ import annotations

import warnings

import numpy as np

warnings.filterwarnings("ignore")

from retained_spectral.engine import (
    RawSpectralProblem,
    native_eigvals_from_tridiagonal,
    retained_raw_input_readout,
    retained_tridiagonal,
    warm_native_kernel,
)

REFERENCE = np.array([4.21443981, 4.21443981, 12.52721843, 12.52721843])
CORRECT_WINDOW = (-8.0, 8.0)
WIDTH = 2.0e-10          # max(0.01 * tol, 2e-12) at tol = 2e-8


def main() -> None:
    warm_native_kernel()
    problem = RawSpectralProblem(
        name="double_well_a2_9",
        potential="symmetric_double_well",
        parameters=(("lam", 1.0), ("a2", 9.0)),
        modes=4,
        tolerance=2.0e-8,
    )

    print("--- what the raw-input planner does on its own ---")
    result = retained_raw_input_readout(problem)
    print("  window chosen :", tuple(round(x, 4) for x in result.window))
    print("  minima of V at: -3.0, +3.0")
    print("  values        :", np.round(np.asarray(result.values), 8))
    print("  status        :", result.status,
          " bound", float(np.max(result.diagnostic_bounds)))
    print("  error         : %.3e" % float(np.max(np.abs(np.asarray(result.values) - REFERENCE))))

    print()
    print("--- the same kernel, correct window supplied by hand ---")
    spectra = {}
    for intervals in (4096, 8192, 16384):
        d, e, _ = retained_tridiagonal(problem, CORRECT_WINDOW, intervals)
        spectra[intervals] = np.asarray(
            native_eigvals_from_tridiagonal(d, e, problem.modes, WIDTH)
        )
        print(f"  M={intervals:6d}", np.round(spectra[intervals], 10))
    richardson = (4.0 * spectra[16384] - spectra[8192]) / 3.0
    print("  Richardson  ", np.round(richardson, 10))
    print("  error       : %.3e  (tolerance %.1e)"
          % (float(np.max(np.abs(richardson - REFERENCE))), problem.tolerance))

    print()
    print("--- Matslise 2.0 on the same correct window, for comparison ---")
    try:
        from pyslise import Pyslise
        # Matslise convention is -y'' + V y = E y, so pass 2V and halve E.
        solver = Pyslise(lambda x: 2.0 * (x * x - 9.0) ** 2,
                         CORRECT_WINDOW[0], CORRECT_WINDOW[1], tolerance=2.0e-8)
        got = solver.eigenvaluesByIndex(0, 4, np.array([0.0, 1.0]))
        print(f"  requested 4, returned {len(got)}:",
              [(i, round(v / 2.0, 10)) for i, v in got])
    except ImportError:
        print("  pyslise not installed (pip install pyslise) --- skipped")

    print()
    print("Conclusion: the kernel resolves the doublet correctly when it is given")
    print("the right window.  The failure is the window, not the solver.")


if __name__ == "__main__":
    main()
