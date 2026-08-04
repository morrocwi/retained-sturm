# Reproducibility

## Science status (read before running or citing anything)

- **Overall benchmark verdict is HOLD**, not ACCEPT — an ARPACK comparator
  (`eigsh`) fails to converge on `morse_lambda5_all_bound` in Runs A–C
  (fairness gate); Run D's fairness gate would pass with the ARPACK-as-primary
  repair, but the speed gate then fails on its own merits because
  `factorized_sextic_ground` is a measured **loss**. See `docs/paper-map.md`
  row for `sec:fairness`.
- **Speed vs. MATSLISE is MIXED**, not a win: 1.33–13x slower on smooth
  potentials (harmonic family, `\|x\|`, Pöschl–Teller), 8x faster only on
  `pure_quartic`. Never state "faster" unqualified.
- **Headline claim is auditable completeness**, not speed: Sturm/inertia
  counting structurally cannot silently skip a root the way shooting-method
  solvers (Matslise/SLEDGE/SLEIGN2) can. This connects to an already
  machine-checked (`[Th_coqc]`) theorem elsewhere in
  `information-discrete-math` — `ResolvedInertia` / `IDM_ResolvedCount.v`, a
  forced three-value certain+/certain-/unresolved-⊥ readout — cited here, not
  reproved here.
- **Known live bug**: the domain-truncation planner can return ACCEPT while
  silently missing a second well on the double-well potential
  `V=(x^2-9)^2, k=4`, error 16.07 against tolerance 2e-8.
  `benchmarks/planner_vs_solver.py` isolates this as a planner/window-selection
  bug, not a solver/pivot bug — the Sturm kernel is correct given the right
  window. `benchmarks/weyl_gate.py` is a reference/diagnostic implementation of
  the necessary-condition repair the paper describes but did not implement; it
  is `[Open]` — not wired into the engine's planner, and its own docstring
  states a residual gap (a component narrower than the scan grid can still be
  missed).
- **Self-correction is on the record, not hidden**: an early Matslise
  comparison used 3 unpaired samples from a single launch and reported
  "2–10x slower" across the board; the paired/interleaved result (see
  `carvscar.py`, `aggregate_carvscar.py`) superseded it with the mixed picture
  above. Tier discipline caught its own error; that is treated as evidence the
  discipline works, not as something to remove from the record.
- Legitimacy here comes from the reproducibility map below and the disclosed
  tier/patch/correction history, not from any outside sign-off.
- Tier tags: `[Th_coqc]` machine-checked elsewhere, `[finite_diagnostic]`
  measured on a finite test set, `[Open]` known unresolved question.

## Engine dependency

This repo does not vendor the solver. All benchmark scripts import
`retained_spectral`, which lives in the sibling repository
**`information-discrete-math`**.

Pinned commit for the numbers in `paper/v1/rms_harmonic.tex` and the results
committed under `benchmarks/results/`:

```
ed51c63be84c915ebdfc643a195322b0bf5feed5
```

**Updated for `paper/v2`'s `sec:v2repair` additions**: `information-discrete-math`'s
`main` branch has since advanced to

```
a05841dd68153e73cd76ea6b1f1d98c54efac6a3
```

(merge commit for PR #109, the ARPACK shift-invert fix — includes `72d5eec`
and `acf4b51`). The `k>1` scope correction and the `k1_discrete_readout`
instruction-count instrument are on branch `fix/scope-speed-gate-to-retained-k-gt-1`,
tip commit

```
3823245720c6f2e31b2760c649b7f80b543768c1
```

**not yet merged to `main`** at the time of writing.

(obtained by running `git -C information-discrete-math rev-parse HEAD` in
this workspace at doc-writing time — 2026-08-04. Re-verify this value
yourself before relying on it; it is not guaranteed to stay HEAD as that
sibling repo continues to change.)

To reproduce against the same commit:

```bash
cd ~/your-workspace   # or wherever this repo's parent lives
git clone <information-discrete-math-repo-url> information-discrete-math
git -C information-discrete-math checkout ed51c63be84c915ebdfc643a195322b0bf5feed5
```

`information-discrete-math` ships a `pyproject.toml` (no `setup.py`), so it
installs as an editable package:

```bash
cd information-discrete-math
pip install -e ".[spectral-bench]"
```

The `spectral-bench` extra (per `pyproject.toml`) pins:

```
numpy>=1.24
scipy>=1.10
numba>=0.58
matplotlib>=3.6
jax>=0.4
```

`numba` is used inside the engine itself (`retained_spectral/engine.py`,
`retained_spectral/retained_mode.py`, and the `competition/` chart/run
modules) — none of the scripts in *this* repo's `benchmarks/` directory
import `numba` directly, but the timings are meaningless without it since the
native kernel it JIT-compiles is what is being timed.

### Path resolution — two conventions coexist, check per script

The benchmark scripts in this repo do **not** all resolve the engine path the
same way; check which convention a given script uses before assuming an env
var takes effect:

- `adv_kernel.py`, `adv_kernel2.py` read **`IDM_REPO_PATH`**, defaulting to
  `Path(__file__).resolve().parent.parent.parent / "information-discrete-math"`
  — i.e. this repo's own parent directory, sibling-layout assumed.
- `kernel_tolerance_sweep.py`, `carvscar.py`, `adv_scaling.py`,
  `adv_strawman.py` read **`IDM_REPO`** (a different variable name), defaulting
  to the current working directory (`os.getcwd()`) if unset.
- `planner_vs_solver.py`, `dvr_reference.py`, `aggregate_carvscar.py`,
  `aggregate_launches.py` do not set `sys.path` themselves — they rely on
  `retained_spectral` already being importable (installed editable, or
  `PYTHONPATH` already set).
- `revision_harness.py`, `scipy_retained.py` are meant to be installed
  **into** the `information-discrete-math` checkout (see "Setup" below), so
  they import `retained_spectral` as a sibling module, not via `sys.path`
  insertion.

If you run scripts from both groups in the same session, set **both**
`IDM_REPO_PATH` and `IDM_REPO` to the same path, plus `PYTHONPATH`, to avoid a
script silently falling back to `os.getcwd()`:

```bash
export IDM_REPO_PATH=~/your-workspace/information-discrete-math
export IDM_REPO=~/your-workspace/information-discrete-math
export PYTHONPATH=.
```

### Additional per-script requirements

- `numba` — required for the compiled kernel; without it, all timing numbers
  in `docs/paper-map.md` are meaningless (correctness is unaffected).
- `pyslise` — required only by `benchmarks/carvscar.py` and
  `benchmarks/planner_vs_solver.py` (the Matslise comparison arm). Install
  with `pip install pyslise`. Not required for any other script.

## Setup: installing the two modules that belong inside the engine checkout

Two scripts are not standalone — they are modules meant to live inside
`retained_spectral/competition/` in the `information-discrete-math` checkout:

| file (this repo) | destination (`information-discrete-math` checkout) |
|---|---|
| `benchmarks/scipy_retained.py` | `retained_spectral/competition/scipy_retained.py` |
| `benchmarks/revision_harness.py` | `retained_spectral/competition/revision_harness.py` |

Copy (or symlink) both into place before running Run C/D reproduction:

```bash
cp benchmarks/scipy_retained.py \
   ~/your-workspace/information-discrete-math/retained_spectral/competition/
cp benchmarks/revision_harness.py \
   ~/your-workspace/information-discrete-math/retained_spectral/competition/
```

## Applying the Run C patch

Run C is Run B plus the single change in `patches/revision_major.patch`
(matched-tolerance fix, see `patches/README.md` for the full rationale).
Apply it inside the `information-discrete-math` checkout at the pinned
commit:

```bash
cd ~/your-workspace/information-discrete-math
git apply ~/your-workspace/retained-sturm/patches/revision_major.patch
```

## Run instructions

All commands below assume the engine is installed per "Setup" above and (for
`revision_harness.py`) the patch is applied.

### Run D — paired/interleaved timing (Field A + Field B, three end-to-end arms)

```bash
export PYTHONPATH=.
python3 -m retained_spectral.competition.revision_harness --launch 1 --json launch1.json
python3 -m retained_spectral.competition.revision_harness --launch 2 --json launch2.json
python3 -m retained_spectral.competition.revision_harness --launch 3 --json launch3.json
```

### Pool the launches into the paper's three main tables

```bash
cd retained-sturm/benchmarks
python3 aggregate_launches.py launch1.json launch2.json launch3.json \
        --results retained_spectral/results/competition_results.json
```

Path to `--results` should point at the pinned engine checkout's committed
results file, e.g.
`~/your-workspace/information-discrete-math/retained_spectral/results/competition_results.json`.

### Tolerance-matching diagnosis

```bash
export IDM_REPO=~/your-workspace/information-discrete-math
python3 benchmarks/kernel_tolerance_sweep.py
```

### Scaling (kernel advantage vs. problem size)

```bash
export IDM_REPO=~/your-workspace/information-discrete-math
python3 benchmarks/adv_scaling.py
```

### Comparator ladder-cost attribution

```bash
export IDM_REPO=~/your-workspace/information-discrete-math
python3 benchmarks/adv_strawman.py
```

### Adversarial tridiagonal suite

```bash
export IDM_REPO_PATH=~/your-workspace/information-discrete-math
export PYTHONPATH=.
python3 benchmarks/adv_kernel.py     # k = 4; writes benchmarks/results/adversarial_kernel.json
python3 benchmarks/adv_kernel2.py    # k up to all eigenvalues; writes benchmarks/results/adversarial_kernel2.json
```

### Independent double-well / non-smooth reference spectra

```bash
export PYTHONPATH=.:~/your-workspace/information-discrete-math
python3 benchmarks/dvr_reference.py
```

### Double-well planner-vs-solver isolation test

```bash
export PYTHONPATH=.:~/your-workspace/information-discrete-math
pip install pyslise   # required for the Matslise comparison in this script
python3 benchmarks/planner_vs_solver.py
```

### Matslise paired comparison (not currently cited in `rms_harmonic.tex` v1 — see `docs/paper-map.md`)

```bash
export IDM_REPO=~/your-workspace/information-discrete-math
pip install pyslise
for L in 1 2 3; do python3 benchmarks/carvscar.py $L; done
python3 benchmarks/aggregate_carvscar.py carvscar1.json carvscar2.json carvscar3.json
```

### `weyl_gate.py` — reference diagnostic for the unimplemented repair

```bash
python3 benchmarks/weyl_gate.py
```

Standalone; no dependency on `retained_spectral`. `[Open]` — see
`docs/paper-map.md` and the science-status section above: this is a
diagnostic-grade global-coverage scan, not a certified fix wired into the
engine's planner.

## What is and is not committed in this repo

- `benchmarks/results/adversarial_kernel.json`,
  `benchmarks/results/adversarial_kernel2.json` — committed here.
- `competition_results.json` (Run A) — lives in the `information-discrete-math`
  checkout at `retained_spectral/results/competition_results.json`, not
  duplicated in this repo.
- `launch1.json` / `launch2.json` / `launch3.json` (Runs D) and
  `carvscar1.json` / `carvscar2.json` / `carvscar3.json` — not committed in
  either repo; regenerate them with the commands above.
