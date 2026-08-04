"""SciPy pipeline WITH retention — the missing arm of the end-to-end comparison.

The shipped comparator (`scipy_pipeline.scipy_raw_input_readout`) rebuilds the
whole Richardson ladder after the truncation window is widened.  RMS does not:
it keeps the mesh correction already purchased on the narrow window and
transports it across a single same-spacing solve on the wide window.

If the comparator is not given the same option, an end-to-end timing ratio
conflates two different things:

  * how fast the native Sturm kernel is, and
  * how much the *idea* of retention saves,

and retention is not native-specific — it is ten lines in any pipeline.  This
module ports retention into the SciPy pipeline and changes NOTHING else: same
window search, same decay gate, same mesh ladder, same DSTEBZ calls at the same
matched tolerance.  The difference between this arm and the shipped arm is the
value of the idea; the difference between this arm and RMS is the value of the
kernel.
"""

from __future__ import annotations

import math
import time

import numpy as np
from scipy.linalg import eigh_tridiagonal

from retained_spectral.engine import RawSpectralProblem, RawSpectralResult
from retained_spectral.competition.scipy_pipeline import (
    _expanded_window,
    _initial_window,
    _scipy_tridiagonal,
)


def _matched_tol(problem: RawSpectralProblem) -> float:
    return max(problem.tolerance * 0.01, 2.0e-12)


def _solve(problem: RawSpectralProblem, window, intervals) -> np.ndarray:
    diagonal, off_diagonal = _scipy_tridiagonal(problem, window, intervals)
    return np.asarray(
        eigh_tridiagonal(
            diagonal,
            off_diagonal,
            select="i",
            select_range=(0, problem.modes - 1),
            eigvals_only=True,
            check_finite=False,
            tol=_matched_tol(problem),
            lapack_driver="stebz",
        )
    )


def _richardson_mesh_raw(problem, window, initial_intervals, max_intervals):
    """As scipy_pipeline._richardson_mesh, but also returns the RAW finest spectrum.

    The raw finest value is what retention needs: the window shift is a
    difference of raw solves at identical spacing, so the extrapolated value
    cannot be used for it.
    """
    spectra: list[np.ndarray] = []
    intervals = max(64, initial_intervals)
    solves = 0
    while intervals <= max_intervals:
        spectra.append(_solve(problem, window, intervals))
        solves += 1
        if len(spectra) >= 3:
            coarse, fine, finest = spectra[-3:]
            r_coarse = (4.0 * fine - coarse) / 3.0
            r_fine = (4.0 * finest - fine) / 3.0
            shift = np.abs(r_fine - r_coarse)
            if float(np.max(shift)) <= 0.45 * problem.tolerance:
                return r_fine, spectra[-1], shift, intervals, solves
        intervals *= 2
    if len(spectra) < 3:
        raise RuntimeError("SciPy mesh cap did not admit three grids")
    coarse, fine, finest = spectra[-3:]
    r_coarse = (4.0 * fine - coarse) / 3.0
    r_fine = (4.0 * finest - fine) / 3.0
    return r_fine, spectra[-1], np.abs(r_fine - r_coarse), intervals // 2, solves


def scipy_retained_readout(
    problem: RawSpectralProblem,
    *,
    max_intervals: int = 1_048_576,
    max_window_rounds: int = 10,
) -> RawSpectralResult:
    started = time.perf_counter()
    scale, window = _initial_window(problem)
    initial_intervals = max(
        64,
        2 ** math.ceil(
            math.log2(max((window[1] - window[0]) / scale * 8.0, 64.0))
        ),
    )
    total_solves = 0
    values = None
    mesh_shift = np.full(problem.modes, math.inf)
    window_shift = np.full(problem.modes, math.inf)
    finest = initial_intervals
    reason = "window rounds exhausted"

    for _round in range(max_window_rounds):
        r_fine, raw_finest, mesh_shift, finest, solves = _richardson_mesh_raw(
            problem, window, initial_intervals, max_intervals
        )
        total_solves += solves

        expanded, tail_pass = _expanded_window(
            problem, window, float(r_fine[-1]), scale
        )

        # --- RETENTION: one same-spacing solve on the wide window, and the
        #     mesh correction already paid for is transported, not recomputed.
        ratio = (expanded[1] - expanded[0]) / (window[1] - window[0])
        expanded_intervals = max(64, int(round(finest * ratio)))
        raw_expanded = _solve(problem, expanded, expanded_intervals)
        total_solves += 1

        window_shift = np.abs(raw_expanded - raw_finest)
        values = r_fine + (raw_expanded - raw_finest)
        diagnostic = mesh_shift + window_shift
        window = expanded
        finest = expanded_intervals

        if tail_pass and float(np.max(diagnostic)) <= problem.tolerance:
            reason = "retained: mesh, decay-boundary, and window gates passed"
            break
        initial_intervals = max(64, expanded_intervals // 8)

    diagnostic = mesh_shift + window_shift
    return RawSpectralResult(
        values=np.asarray(values),
        window=window,
        finest_intervals=finest,
        solve_count=total_solves,
        mesh_shift=np.asarray(mesh_shift),
        window_shift=np.asarray(window_shift),
        diagnostic_bounds=np.asarray(diagnostic),
        status="ACCEPT" if float(np.max(diagnostic)) <= problem.tolerance else "HOLD",
        reason=reason,
        elapsed_seconds=time.perf_counter() - started,
        method="SciPy Richardson-stability pipeline WITH retention",
        recurrence_updates=0,
        peak_working_bytes=0,
        reference_used_for_planning=False,
        per_case_schedule=False,
        tier="finite_diagnostic",
    )
