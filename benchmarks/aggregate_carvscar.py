"""Aggregate the RMS vs Matslise 2.0 launches.

`carvscar.py` writes one JSON per launch with interleaved paired samples.  This
pools them and reports the ratio RMS/Matslise, so a value above 1 means RMS is
slower.

Two asymmetries must be read alongside the numbers, and the script prints them:

* **Matslise is not paying for domain discovery.**  It takes (V, min, max, tol);
  the window handed to it here is the one RMS derived for itself.  RMS's timing
  therefore includes well search, gate evaluation, and the window witness, and
  Matslise's does not.  The comparison answers "how long does the whole task
  take with each tool", not "which solver core is faster".
* **Rows where the two arms saw different domains are invalid for timing** and
  are flagged rather than silently averaged.

Run::

    python3 aggregate_carvscar.py carvscar1.json carvscar2.json carvscar3.json
"""

from __future__ import annotations

import json
import sys

import numpy as np

RESAMPLES = 3000


def paired_ci(a, b, seed=1, resamples=RESAMPLES):
    a, b = np.asarray(a, float), np.asarray(b, float)
    n = len(a)
    rng = np.random.default_rng(seed)
    boot = np.empty(resamples)
    for i in range(resamples):
        idx = rng.integers(0, n, size=n)
        boot[i] = np.median(a[idx]) / np.median(b[idx])
    return (float(np.median(a) / np.median(b)),
            float(np.percentile(boot, 2.5)),
            float(np.percentile(boot, 97.5)))


def main(paths):
    launches = [json.load(open(p)) for p in paths]
    print(f"{'case':30s} {'RMS ms':>8} {'MAT ms':>8} {'RMS/MAT':>8} "
          f"{'CI':>16} {'RMS err':>10} {'MAT err':>10}")
    for name in launches[0]:
        first = launches[0][name]
        if "matslise_error" in first:
            print(f"{name[:30]:30s}   matslise did not run: {first['matslise_error'][:44]}")
            continue
        a = np.array([x for r in launches for x in r[name]["rms_samples"]])
        b = np.array([x for r in launches for x in r[name]["mat_samples"]])
        pt, lo, hi = paired_ci(a, b)
        flag = ""
        if first.get("hand_window"):
            flag = "   INVALID for timing: arms saw different domains"
        print(f"{name[:30]:30s} {np.median(a):8.2f} {np.median(b):8.2f} {pt:8.2f} "
              f"[{lo:6.2f},{hi:6.2f}] {first['rms_err']:10.2e} "
              f"{first['mat_err']:10.2e}{flag}")
    print()
    print("Ratio > 1 means RMS is slower.  Matslise is given the window; RMS derives it.")


if __name__ == "__main__":
    main(sys.argv[1:] or ["carvscar1.json", "carvscar2.json", "carvscar3.json"])
