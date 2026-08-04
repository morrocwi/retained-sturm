"""Adversarial Field-A suite: the matrices DSTEBZ was hardened for.

The native kernel clamps a vanishing pivot at 1e-300 and does nothing else.
LAPACK's DSTEBZ carries scaling and block splitting developed precisely for
clustered spectra, badly scaled matrices, and negligible off-diagonals.  None of
the seven declared physics instances exercises any of that: they are all smooth
single- or double-well potentials whose discretisations have a constant
off-diagonal -1/(2h^2) and a modest condition number.

This suite goes straight at the gap, using the standard hard cases from the
symmetric-tridiagonal literature.  Correctness is judged against dense LAPACK
(a different algorithm: implicit QR on the full matrix), timing second --- a
speed win on a matrix whose eigenvalues are wrong is worth nothing.
"""

from __future__ import annotations

import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import json
import random
import statistics
import sys
import time

import numpy as np
import scipy.linalg as sla

from pathlib import Path

_DEFAULT_IDM_REPO = Path(__file__).resolve().parent.parent.parent / "information-discrete-math"
sys.path.insert(0, os.environ.get("IDM_REPO_PATH", str(_DEFAULT_IDM_REPO)))
from retained_spectral.engine import native_eigvals_from_tridiagonal, warm_native_kernel


def wilkinson(n_half: int = 10):
    """W+_(2m+1): the classic near-degenerate-pair matrix."""
    n = 2 * n_half + 1
    d = np.abs(np.arange(n, dtype=float) - n_half)
    e = np.ones(n - 1)
    return d, e


def glued_wilkinson(copies: int = 8, glue: float = 1.0e-14):
    """Dhillon's glued Wilkinson: blocks joined by a negligible off-diagonal.

    This is the canonical block-splitting test.  DSTEBZ detects the tiny
    coupling, splits into independent blocks, and solves each; a solver that
    does not split sees eigenvalues of multiplicity `copies`.
    """
    d0, e0 = wilkinson(10)
    d = np.tile(d0, copies)
    e = np.empty(len(d) - 1)
    m = len(d0)
    for c in range(copies):
        lo = c * m
        e[lo:lo + m - 1] = e0
        if c < copies - 1:
            e[lo + m - 1] = glue
    return d, e


def toeplitz_121(n: int = 4000):
    """d=2, e=-1: exact spectrum 2 - 2cos(k*pi/(n+1)), crowded at both ends."""
    return np.full(n, 2.0), np.full(n - 1, -1.0)


def toeplitz_exact(n: int, k: int):
    idx = np.arange(1, k + 1)
    return 2.0 - 2.0 * np.cos(idx * np.pi / (n + 1))


def scaled_extremes(n: int = 2000, span: float = 1.0e10):
    """Diagonal entries spanning twenty orders of magnitude."""
    rng = np.random.default_rng(7)
    expo = np.linspace(-np.log10(span), np.log10(span), n)
    d = 10.0 ** expo * rng.choice([-1.0, 1.0], size=n)
    e = 10.0 ** (0.5 * (expo[:-1] + expo[1:])) * 0.5
    return d, e


def denormal_split(n: int = 1000):
    """Two blocks joined by an off-diagonal at the pivot floor itself."""
    rng = np.random.default_rng(11)
    d = rng.uniform(-1.0, 1.0, size=n)
    e = rng.uniform(0.5, 1.5, size=n - 1)
    e[n // 2] = 1.0e-300
    return d, e


def zero_offdiag_cluster(n: int = 600):
    """Exactly repeated diagonal blocks with zero coupling: multiplicity n/3."""
    block = np.array([1.0, 2.0, 3.0])
    d = np.tile(block, n // 3)
    e = np.zeros(n - 1)
    e[np.arange(n - 1) % 3 != 2] = 1.0e-8
    return d, e


CASES = {
    "wilkinson_w21": wilkinson(10),
    "glued_wilkinson_8x21": glued_wilkinson(8),
    "toeplitz_121_n4000": toeplitz_121(4000),
    "scaled_extremes_1e10": scaled_extremes(2000),
    "denormal_split_1e-300": denormal_split(1000),
    "zero_offdiag_cluster": zero_offdiag_cluster(600),
}

WIDTH = 2.0e-10          # the same declared bisection half-width as Run D
K = 4


def hot(fn, reps, rng, other=None):
    fn()
    if other is not None:
        other()
    a, b = [], []
    for _ in range(reps):
        order = [0, 1] if rng.random() < 0.5 else [1, 0]
        for which in order:
            f = fn if which == 0 else other
            t = time.perf_counter()
            f()
            (a if which == 0 else b).append(time.perf_counter() - t)
    return a, b


def main():
    warm_native_kernel()
    rng = random.Random(20260804)
    out = {}
    print(f"{'case':24s} {'n':>6} {'native err':>12} {'dstebz err':>12} "
          f"{'nat ms':>8} {'ste ms':>8} {'ratio':>6}")
    for name, (d, e) in CASES.items():
        n = len(d)
        dense = np.diag(d) + np.diag(e, 1) + np.diag(e, -1)
        ref = np.sort(np.linalg.eigvalsh(dense))[:K]

        try:
            nat = np.asarray(native_eigvals_from_tridiagonal(d, e, K, WIDTH))
            nat_err = float(np.max(np.abs(np.sort(nat) - ref)))
        except Exception as exc:                       # noqa: BLE001
            nat, nat_err = None, float("inf")
            print(f"  native raised: {type(exc).__name__}: {exc}")

        try:
            ste = np.asarray(sla.eigh_tridiagonal(
                d, e, select="i", select_range=(0, K - 1), eigvals_only=True,
                check_finite=False, tol=WIDTH, lapack_driver="stebz"))
            ste_err = float(np.max(np.abs(np.sort(ste) - ref)))
        except Exception as exc:                       # noqa: BLE001
            ste, ste_err = None, float("inf")
            print(f"  dstebz raised: {type(exc).__name__}: {exc}")

        reps = 9 if n < 3000 else 5
        ta, tb = hot(lambda: native_eigvals_from_tridiagonal(d, e, K, WIDTH), reps, rng,
                     other=lambda: sla.eigh_tridiagonal(
                         d, e, select="i", select_range=(0, K - 1), eigvals_only=True,
                         check_finite=False, tol=WIDTH, lapack_driver="stebz"))
        mn, ms = statistics.median(ta) * 1e3, statistics.median(tb) * 1e3
        print(f"{name:24s} {n:>6d} {nat_err:>12.3e} {ste_err:>12.3e} "
              f"{mn:>8.3f} {ms:>8.3f} {ms/mn:>6.3f}")
        out[name] = {"n": n, "native_err": nat_err, "dstebz_err": ste_err,
                     "native_ms": mn, "dstebz_ms": ms, "ratio": ms / mn,
                     "reference_first_k": ref.tolist(),
                     "native_values": None if nat is None else np.sort(nat).tolist(),
                     "dstebz_values": None if ste is None else np.sort(ste).tolist()}
        sys.stdout.flush()
    _out_path = Path(__file__).resolve().parent / "results" / "adversarial_kernel.json"
    _out_path.parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(_out_path, "w"), indent=2)


if __name__ == "__main__":
    main()
