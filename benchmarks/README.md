# RMS benchmark scripts

## Dependency: this repo does not vendor the engine

These benchmark scripts import `retained_spectral`, a package that lives in the
sibling repository **`information-discrete-math`**, not in this repo. Nothing
under `retained_spectral/` is copied here — the scripts resolve it at import
time.

`information-discrete-math` ships a `pyproject.toml` (no `setup.py`), so it
installs as an editable package with `pip`.

### Setup

Clone the sibling repo next to this one (the default path resolution below
assumes exactly this layout: `information-discrete-math/` and
`retained-sturm/` as siblings under the same parent directory):

```bash
cd ~/your-workspace   # or wherever this repo's parent lives
git clone <information-discrete-math-repo-url> information-discrete-math
```

Install it (editable, with the benchmark extras — `numpy`, `scipy`, `numba`,
`matplotlib`, `jax`):

```bash
cd information-discrete-math
pip install -e ".[spectral-bench]"
```

If you clone `information-discrete-math` somewhere else (not as a sibling
directory), point the scripts at it explicitly:

```bash
export IDM_REPO_PATH=/path/to/information-discrete-math
```

`adv_kernel.py` and `adv_kernel2.py` resolve the engine path as
`os.environ.get("IDM_REPO_PATH", <sibling-default-via-pathlib>)` — the env var
wins if set; otherwise they default to
`Path(__file__).resolve().parent.parent.parent / "information-discrete-math"`,
i.e. this repo's own parent directory, sibling-layout assumed.

```bash
export PYTHONPATH=.
```

Requirements beyond the repository's own: `numba` (the compiled kernel —
without it the timings are meaningless), and `pyslise` for the Matslise
comparison (`pip install pyslise`).

---

## Science status (read before citing any number below)

- **Overall benchmark verdict is HOLD**, not ACCEPT — an ARPACK comparator
  failed to converge on the Morse case. Nothing below overrides that.
- **Speed vs Matslise is MIXED**, not a win: 1.33–13x slower on smooth
  potentials (harmonic family, `|x|`, Pöschl–Teller), 8x faster only on
  `pure_quartic`. Do not state "faster" without that qualification.
- **Known live bug**: the domain-truncation planner can return ACCEPT while
  silently missing a second well on a double-well potential
  (`V=(x^2-9)^2, k=4`), error 16.07 against tolerance 2e-8.
  `planner_vs_solver.py` isolates this as a planner/window-selection bug, not
  a solver/pivot bug — the Sturm kernel itself is correct given the right
  window (see that script's section below).
- **Headline claim is auditable completeness**, not speed: Sturm/inertia
  counting structurally cannot silently skip a root the way shooting-method
  solvers (Matslise/SLEDGE/SLEIGN2) can. This connects to an already
  machine-checked (`[Th_coqc]`) theorem elsewhere in
  `information-discrete-math` — `ResolvedInertia` / `IDM_ResolvedCount.v`, a
  forced three-value certain+/certain-/unresolved-⊥ readout. That theorem is
  cited here, not reproved here.
- **Self-correction history is a feature, not something to hide** — see "Note
  on an earlier error" below.
- Tier tags used throughout: `[Th_coqc]` machine-checked elsewhere,
  `[finite_diagnostic]` measured on a finite test set, `[Open]` known
  unresolved question.

---

Two of the scripts (`scipy_retained.py`, `revision_harness.py`) are modules
that belong inside `retained_spectral/competition/` (in the
`information-discrete-math` checkout); the rest are standalone and run from
this `benchmarks/` directory.

---

## Modules to install into the package

| file | destination | what it is |
|---|---|---|
| `scipy_retained.py` | `retained_spectral/competition/` | The SciPy pipeline **with retention ported in** — the third end-to-end arm. Same window search, same decay gate, same ladder, same DSTEBZ calls at the same matched tolerance; the only change is that the widened window costs one same-spacing solve plus transport of the correction instead of a rebuilt ladder. Without this arm an end-to-end ratio cannot separate "our kernel is fast" from "our comparator was denied the idea under test". |
| `revision_harness.py` | `retained_spectral/competition/` | Paired, interleaved, multi-launch timing for both fields, plus the three-arm end-to-end comparison. Replaces the shipped design, which timed each arm in a separate loop and bootstrapped them independently. |

```bash
python3 -m retained_spectral.competition.revision_harness --launch 1 --json launch1.json
python3 -m retained_spectral.competition.revision_harness --launch 2 --json launch2.json
python3 -m retained_spectral.competition.revision_harness --launch 3 --json launch3.json
```

---

## Standalone scripts

### `aggregate_launches.py`
Pools the launches and produces the paper's three main tables. The bootstrap
resamples **paired** repeat indices, which the shipped statistics module does
not.

```bash
python3 aggregate_launches.py launch1.json launch2.json launch3.json \
        --results retained_spectral/results/competition_results.json
```

Produces `[finite_diagnostic]`: kernel geometric mean 1.867 and the speed
verdict (5 wins, 1 tie, 1 **loss** on `factorized_sextic_ground`); end-to-end
geomeans 3.279 shipped / **2.411 retained** / 1.360 retention-alone; and the
bound-versus-error table (bound exceeds the audited error on all seven, ratios
7.7–21.4).

### `kernel_tolerance_sweep.py`
The tolerance-matching diagnosis. Times the native kernel at the instance
tolerance and at the engine's true working half-width, against DSTEBZ at its
default and at both matched values, on all seven instances.

Produces `[finite_diagnostic]`: kernel geomean 2.560 as shipped → 1.908
matched. Establishes that the comparator was being charged for work the
native side never did, and that the correct matched quantity is
`max(0.01*tol, 2e-12)`, not `tol`.

### `adv_scaling.py`
Field A on `harmonic_low4` at matched width from n = 767 to n = 786,431.

Produces `[finite_diagnostic]` the flat band 1.85–2.56 that refutes the "it is
only wrapper overhead" objection.

### `adv_strawman.py`
Measures, solve by solve, how much of the shipped comparator's time is the
Richardson ladder it rebuilds after widening the window.

Produces `[finite_diagnostic]`: 26.9% of comparator solve time, and the
projected end-to-end ratio 4.565 → 3.382 on `harmonic_low4`. This was the
estimate later confirmed by the direct third-arm measurement (1.360 geomean,
i.e. 26.5%).

### `adv_kernel.py`, `adv_kernel2.py`
The adversarial tridiagonal suite: Wilkinson W⁺₂₁ and W⁺₄₀₁, glued Wilkinson
(eight blocks joined by 1e-14 — the canonical block-splitting test), Toeplitz
[-1,2,-1], a matrix spanning 10^±10, a denormal split at the pivot floor
itself, and a zero-coupling cluster. Correctness is judged against dense
`eigvalsh`, a different algorithm.

`adv_kernel.py` runs k = 4; `adv_kernel2.py` pushes k up to all eigenvalues.

Output: `results/adversarial_kernel.json` and `results/adversarial_kernel2.json`
(written relative to the script's own directory, not the working directory).

Produces `[finite_diagnostic]`: accuracy is never the problem, including on
the cases that target the 1e-300 pivot floor. Speed is: on glued Wilkinson
the ratio falls to **0.478 at k = 64 and 0.326 at k = 168** — the native
kernel is three times slower, because DSTEBZ splits the matrix into eight
blocks and it cannot.

### `dvr_reference.py`
Independent reference spectra: sinc-DVR (Colbert–Miller), exact Airy zeros
for V = |x|, and the double-well calibration sweep used to pick a barrier
whose doublet splitting (7.9e-11) falls below the bisection half-width.

### `planner_vs_solver.py`
The decisive isolation test for the double-well failure. Hands the kernel the
correct window by hand and runs the same Richardson ladder.

Produces `[finite_diagnostic]`: with the correct window the kernel returns
the degenerate pairs correctly at error 2.05e-09 < tolerance. **The failure
is the window, not the solver.** Also runs Matslise on the same window, which
returns only 2 of the 4 requested eigenvalues — a shooting method can step
over a doublet; a Sturm count cannot.

### `carvscar.py`, `aggregate_carvscar.py`
RMS against Matslise 2.0 on the same domain and tolerance, paired and
interleaved, three launches.

```bash
for L in 1 2 3; do python3 carvscar.py $L; done
python3 aggregate_carvscar.py carvscar1.json carvscar2.json carvscar3.json
```

Produces `[finite_diagnostic]`: RMS is 1.33–1.39x slower on the harmonic
family, 7.6x slower on `|x|`, 13x slower on Pöschl–Teller, and **8x faster**
on `pure_quartic`. Two cases did not run under this author's use of the
pyslise API and are reported as such rather than scored as failures.

**Read this comparison with its asymmetry in view:** Matslise takes
`(V, min, max, tolerance)` and is handed the window RMS derived, so RMS's
timing includes well search, gate evaluation and the window witness while
Matslise's does not. The `double_well_a2_9` row is flagged invalid for timing
because the two arms saw different domains.

---

## Note on an earlier error

A first pass at the Matslise comparison used three unpaired samples from a
single launch and reported RMS as "2–10x slower" across the board. That is
the design defect criticised elsewhere in this work, committed in the
measurement that criticised it. The paired result above supersedes it: the
picture is mixed, not uniformly unfavourable. This correction stands as
evidence the tier discipline works, not as something to hide — it is
recorded here rather than quietly overwritten.
