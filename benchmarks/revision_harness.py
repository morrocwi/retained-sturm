"""Revision harness: paired timing, three end-to-end arms, one process launch.

Differences from `competition.run`, all of them responses to referee points:

* **Paired timing in Field A.**  The shipped kernel audit times the native
  kernel in one loop and each competitor in a later loop, so any thermal or
  frequency drift between the loops is charged systematically to one arm.  Here
  the arms are interleaved within every repeat under a recorded seed, and the
  bootstrap resamples *repeat indices*, not each arm independently, so the
  interval respects the pairing.

* **More samples.**  50 kernel repeats and 15 end-to-end repeats per launch,
  against 11 and 21 before.

* **Multiple process launches.**  This script does one launch and writes one
  JSON; run it several times and aggregate.  Median-of-samples inside a single
  process cannot see between-process variance in allocator state and code
  placement.

* **Three end-to-end arms.**  RMS, the shipped SciPy pipeline, and the SciPy
  pipeline with retention ported in.  Without the third arm an end-to-end ratio
  cannot separate "our kernel is faster" from "our comparator was not given the
  idea under test".

Usage::

    PYTHONPATH=. python3 -m retained_spectral.competition.revision_harness \
        --launch 1 --json /tmp/launch1.json
"""

from __future__ import annotations

import os

for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
             "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import argparse
import json
import platform
import random
import statistics
import sys
import time
from pathlib import Path

import numpy as np
import scipy.linalg as sla

from retained_spectral.engine import (
    KERNEL_FIELD,
    NATIVE_KERNEL_COMPILED,
    native_eigvals_from_tridiagonal,
    raw_benchmark_targets,
    retained_raw_input_readout,
    warm_native_kernel,
)
from retained_spectral.competition.executor_audit import retained_tridiagonal
from retained_spectral.competition.scipy_pipeline import scipy_raw_input_readout
from retained_spectral.competition.scipy_retained import scipy_retained_readout

BOOTSTRAP_RESAMPLES = 4000


def matched_width(problem) -> float:
    return max(problem.tolerance * 0.01, 2.0e-12)


def paired_ratio_ci(a, b, *, seed: int, resamples: int = BOOTSTRAP_RESAMPLES):
    """95% CI for median(b)/median(a) resampling PAIRED repeat indices."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    point = float(np.median(b) / np.median(a))
    n = len(a)
    if n < 2:
        return {"ratio_median": point, "ci95_low": point, "ci95_high": point,
                "verdict": "insufficient-samples", "n_pairs": n}
    rng = np.random.default_rng(seed)
    boot = np.empty(resamples)
    for i in range(resamples):
        idx = rng.integers(0, n, size=n)
        boot[i] = np.median(b[idx]) / np.median(a[idx])
    lo, hi = float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))
    verdict = "native_faster" if lo > 1.0 else ("competitor_faster" if hi < 1.0 else "tie")
    return {"ratio_median": point, "ci95_low": lo, "ci95_high": hi,
            "verdict": verdict, "n_pairs": n, "resamples": resamples}


def interleaved(calls: dict, repeats: int, rng: random.Random) -> dict:
    """Time every arm once per repeat, in a shuffled order each repeat."""
    for fn in calls.values():          # warm / compile, untimed
        fn()
    samples = {name: [] for name in calls}
    names = list(calls)
    for _ in range(repeats):
        rng.shuffle(names)
        for name in names:
            started = time.perf_counter()
            calls[name]()
            samples[name].append(time.perf_counter() - started)
    return samples


def kernel_field(target, *, intervals: int, repeats: int, seed: int) -> dict:
    problem = target.problem
    window = retained_raw_input_readout(problem).window
    diagonal, off_diagonal, _ = retained_tridiagonal(problem, window, intervals)
    k, w = problem.modes, matched_width(problem)
    rng = random.Random(f"{seed}:kernel:{problem.name}")
    samples = interleaved(
        {
            "native": lambda: native_eigvals_from_tridiagonal(diagonal, off_diagonal, k, w),
            "dstebz": lambda: sla.eigh_tridiagonal(
                diagonal, off_diagonal, select="i", select_range=(0, k - 1),
                eigvals_only=True, check_finite=False, tol=w, lapack_driver="stebz"),
        },
        repeats, rng,
    )
    out = {"intervals": intervals, "matched_width": w,
           "median_ms": {n: statistics.median(s) * 1e3 for n, s in samples.items()},
           "samples_ms": {n: [x * 1e3 for x in s] for n, s in samples.items()}}
    out["dstebz_over_native"] = paired_ratio_ci(samples["native"], samples["dstebz"], seed=seed)
    return out


def end_to_end_field(target, *, repeats: int, seed: int) -> dict:
    problem = target.problem
    reference = np.asarray(target.reference, dtype=float)
    rng = random.Random(f"{seed}:e2e:{problem.name}")
    arms = {
        "rms": lambda: retained_raw_input_readout(problem),
        "scipy_shipped": lambda: scipy_raw_input_readout(problem),
        "scipy_retained": lambda: scipy_retained_readout(problem),
    }
    samples = interleaved(arms, repeats, rng)
    results = {name: fn() for name, fn in arms.items()}
    out = {
        "median_ms": {n: statistics.median(s) * 1e3 for n, s in samples.items()},
        "samples_ms": {n: [x * 1e3 for x in s] for n, s in samples.items()},
        "solve_count": {n: int(r.solve_count) for n, r in results.items()},
        "max_reference_abs_error": {
            n: float(np.max(np.abs(np.asarray(r.values) - reference)))
            for n, r in results.items()},
        "status": {n: r.status for n, r in results.items()},
        "tolerance": problem.tolerance,
    }
    out["shipped_over_rms"] = paired_ratio_ci(samples["rms"], samples["scipy_shipped"], seed=seed)
    out["retained_over_rms"] = paired_ratio_ci(samples["rms"], samples["scipy_retained"], seed=seed)
    out["shipped_over_retained"] = paired_ratio_ci(
        samples["scipy_retained"], samples["scipy_shipped"], seed=seed)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--launch", type=int, required=True)
    ap.add_argument("--kernel-repeats", type=int, default=50)
    ap.add_argument("--e2e-repeats", type=int, default=15)
    ap.add_argument("--intervals", type=int, default=768)
    ap.add_argument("--seed", type=int, default=20260804)
    ap.add_argument("--json", type=Path, required=True)
    args = ap.parse_args()

    warm_native_kernel()
    seed = args.seed + args.launch
    record = {
        "schema": "idm.retained-spectral-revision.v1",
        "launch": args.launch,
        "seed": seed,
        "design": "arms interleaved within each repeat; bootstrap resamples paired repeat indices",
        "environment": {
            "python": sys.version.split()[0],
            "platform": platform.platform(),
            "numpy": np.__version__,
            "native_kernel_compiled": bool(NATIVE_KERNEL_COMPILED),
            "kernel_field": KERNEL_FIELD,
        },
        "kernel": {},
        "end_to_end": {},
    }
    for target in raw_benchmark_targets():
        name = target.problem.name
        record["kernel"][name] = kernel_field(
            target, intervals=args.intervals, repeats=args.kernel_repeats, seed=seed)
        record["end_to_end"][name] = end_to_end_field(
            target, repeats=args.e2e_repeats, seed=seed)
        print(f"launch {args.launch}  {name} done", flush=True)

    args.json.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("wrote", args.json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
