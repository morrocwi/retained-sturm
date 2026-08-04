"""Aggregate the revision-harness launches into the paper's tables.

`revision_harness.py` writes one JSON per process launch.  This pools them and
produces:

* Field A (kernel-only): the paired ratio DSTEBZ/native with a 95% bootstrap
  interval over pooled repeat indices, plus the per-launch point estimates so a
  reader can see whether a verdict is stable across launches;
* Field B (end-to-end): three arms --- RMS, the shipped SciPy pipeline, and the
  SciPy pipeline with retention ported in --- reported as shipped/RMS,
  retained/RMS, and shipped/retained.  The third column isolates the value of
  the retention idea from the value of the kernel;
* the diagnostic-bound soundness table: released bound against audited error.

The bootstrap resamples *paired* indices, because the harness interleaves the
arms within every repeat.  Resampling the arms independently, as the shipped
statistics module does, throws that pairing away.

Run::

    python3 aggregate_launches.py launch1.json launch2.json launch3.json
    python3 aggregate_launches.py --results competition_results.json   # bound table
"""

from __future__ import annotations

import argparse
import json
import math

import numpy as np

RESAMPLES = 4000
SEED = 20260804


def paired_ci(a, b, seed: int = SEED, resamples: int = RESAMPLES):
    a, b = np.asarray(a, float), np.asarray(b, float)
    n = len(a)
    point = float(np.median(b) / np.median(a))
    rng = np.random.default_rng(seed)
    boot = np.empty(resamples)
    for i in range(resamples):
        idx = rng.integers(0, n, size=n)
        boot[i] = np.median(b[idx]) / np.median(a[idx])
    lo, hi = float(np.percentile(boot, 2.5)), float(np.percentile(boot, 97.5))
    verdict = "win" if lo > 1.0 else ("LOSS" if hi < 1.0 else "TIE")
    return point, lo, hi, verdict


def geomean(xs):
    return math.exp(sum(math.log(x) for x in xs) / len(xs))


def pooled(launches, field, case, arm):
    return [x for r in launches for x in r[field][case]["samples_ms"][arm]]


def kernel_table(launches):
    cases = list(launches[0]["kernel"])
    print("=== Field A (kernel-only), paired, "
          f"{len(launches)} launches x {len(launches[0]['kernel'][cases[0]]['samples_ms']['native'])} "
          "repeats ===")
    ratios = []
    for c in cases:
        pt, lo, hi, v = paired_ci(pooled(launches, "kernel", c, "native"),
                                  pooled(launches, "kernel", c, "dstebz"))
        per = [r["kernel"][c]["dstebz_over_native"]["ratio_median"] for r in launches]
        ratios.append(pt)
        print(f"{c:32s} {pt:6.3f} [{lo:5.3f},{hi:5.3f}] {v:4s} "
              f"per-launch {['%.3f' % x for x in per]}")
    print(f"geometric mean: {geomean(ratios):.4f}")
    losses = sum(1 for c in cases
                 if paired_ci(pooled(launches, "kernel", c, "native"),
                              pooled(launches, "kernel", c, "dstebz"))[3] == "LOSS")
    ties = sum(1 for c in cases
               if paired_ci(pooled(launches, "kernel", c, "native"),
                            pooled(launches, "kernel", c, "dstebz"))[3] == "TIE")
    verdict = "HOLD" if losses else ("TIE" if ties else "ACCEPT")
    print(f"declared speed verdict: {verdict} "
          f"({len(cases) - losses - ties} wins, {ties} ties, {losses} losses)")


def end_to_end_table(launches):
    cases = list(launches[0]["end_to_end"])
    print()
    print("=== Field B (end-to-end), three arms ===")
    sh, rt, sr = [], [], []
    for c in cases:
        rms = pooled(launches, "end_to_end", c, "rms")
        shipped = pooled(launches, "end_to_end", c, "scipy_shipped")
        retained = pooled(launches, "end_to_end", c, "scipy_retained")
        p1, *_ = paired_ci(rms, shipped)
        p2, l2, h2, v2 = paired_ci(rms, retained)
        p3, *_ = paired_ci(retained, shipped)
        sh.append(p1); rt.append(p2); sr.append(p3)
        counts = launches[0]["end_to_end"][c]["solve_count"]
        print(f"{c:32s} solves {counts['rms']:3d}/{counts['scipy_shipped']:3d}/"
              f"{counts['scipy_retained']:3d} | shipped/RMS {p1:5.3f} | "
              f"retained/RMS {p2:5.3f} [{l2:5.3f},{h2:5.3f}] {v2:4s} | "
              f"shipped/retained {p3:5.3f}")
    print(f"geomean shipped/RMS      : {geomean(sh):.4f}   <- overstated headline")
    print(f"geomean retained/RMS     : {geomean(rt):.4f}   <- honest headline")
    print(f"geomean shipped/retained : {geomean(sr):.4f}   <- value of retention alone")


def bound_table(results_path):
    d = json.load(open(results_path))
    print()
    print("=== diagnostic bound vs audited error ===")
    violations = []
    for name, case in sorted(d["end_to_end"]["cases"].items()):
        native = case["native"]
        bound = max(native["diagnostic_bounds"])
        err = native["max_reference_abs_error"]
        print(f"{name:32s} B={bound:10.3e}  err={err:10.3e}  ratio={bound / err:6.2f}")
        if bound < err:
            violations.append(name)
    print("bound below the audited error on:", violations or "no case")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("launches", nargs="*", help="launch*.json from revision_harness.py")
    ap.add_argument("--results", help="competition_results.json for the bound table")
    args = ap.parse_args()
    if args.launches:
        launches = [json.load(open(p)) for p in args.launches]
        kernel_table(launches)
        end_to_end_table(launches)
    if args.results:
        bound_table(args.results)


if __name__ == "__main__":
    main()
