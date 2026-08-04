# Retained Multilevel Sturm (RMS)

A raw-input-planning, diagnostic-acceptance-gated architecture for the lowest
$k$ eigenvalues of a one-dimensional Schrödinger operator — and a benchmark
against itself that was twice corrected against its own first pass.

**Author:** Yaoharee Lahtee, Open Civil Science Initiative, Bangkok, Thailand
(ORCID: 0009-0005-3861-0626)

---

## Abstract

> We describe and evaluate *Retained Multilevel Sturm* (RMS), a computational
> architecture for the lowest $k$ eigenvalues of a one-dimensional
> Schrödinger operator that accepts only raw problem data — the potential
> callable, its parameters, the number of levels requested, and a
> tolerance — and derives the truncation window, the mesh ladder, and the
> eigenvalue brackets internally. The architecture composes five classical
> ingredients (central finite differences, Sturm sequences on a symmetric
> tridiagonal matrix, bisection, Richardson extrapolation, and a WKB-style
> boundary diagnostic) under two design rules that are not standard practice:
> search brackets and mesh corrections computed at one level are *retained
> and revalidated* at the next rather than recomputed, and mesh convergence
> is separated from domain-truncation stability into two independent
> acceptance gates whose sum forms the reported diagnostic bound.
> On seven declared spectra — harmonic, displaced, squeezed, quartic,
> factorised sextic, Pöschl–Teller and Morse — all released values fall
> within their declared tolerances against analytic, factorisation-derived,
> or published references consulted only after both pipelines returned; on
> the harmonic instance the four lowest levels are reproduced to
> $1.465\times10^{-10}$ at a tolerance of $2\times10^{-8}$. Retention reduces
> the eigensolve count from $12$ to $7$ on the four-level cases.
> Our main results concern the benchmark rather than the method. Three
> defects are diagnosed and repaired, and each repair costs the method
> something. (i) The shipped comparator was left at its LAPACK default
> tolerance while the native kernel bisected to a much looser declared width,
> so the two sides were not doing the same work; we also show that the
> obvious repair — matching the comparator to the instance tolerance
> $\varepsilon$ — makes the pipeline's own mesh gate at $0.45\varepsilon$
> unsatisfiable and the refinement loop non-terminating, and we identify the
> correct matched quantity. (ii) The kernel field timed each arm in a
> separate loop and bootstrapped them independently; we interleave the arms
> within every repeat, resample paired indices, and take $150$ pairs across
> three process launches. Under this design one instance becomes a measured
> *loss* (ratio $0.986$, interval upper bound $0.999$) and a second a tie, so
> the project's own declared speed verdict falls from ACCEPT to HOLD; the
> kernel geometric mean falls from $2.560$ to $1.867$. (iii) The end-to-end
> comparator did not implement the retention that is the architecture's own
> contribution, so the comparison conflated kernel speed with the value of an
> idea any pipeline could adopt. We port retention into the comparator and
> add it as a third arm: the honest end-to-end geometric mean is $2.411$, not
> $3.279$, with $26.5\%$ of the original margin attributable to the missing
> idea rather than to the method. RMS still wins all seven end-to-end cases
> by paired interval against this fairer opponent. Two further attacks fail:
> the kernel advantage does not decay with problem size (ratio $1.85$–$2.56$
> from $n = 767$ to $n = 786{,}431$), and the reported diagnostic bound
> exceeds the audited error on all seven instances by factors of $7.7$ to
> $21.4$.
> Finally we evaluate outside the declared suite, on the matrices and
> potentials the design is least equipped for, and report two losses. On a
> glued Wilkinson matrix — eight blocks joined by an off-diagonal of
> $10^{-14}$ — accuracy is unaffected, but the native kernel is $2.1\times$
> slower than DSTEBZ at $k = 64$ and $3.1\times$ slower with all eigenvalues
> requested, because DSTEBZ splits the matrix into blocks and the native
> kernel cannot. And on a symmetric double well with a barrier above the
> requested levels, all three pipelines return the spectrum of a single well
> with status ACCEPT and a passing diagnostic bound, an error of $16.07$
> against a tolerance of $2\times10^{-8}$: the truncation gate tests whether
> the boundary lies in a classically forbidden region, and a barrier between
> two wells passes that test exactly as a decaying tail does. We give the
> repair and restrict the planner's honest scope to single-well potentials
> with monotone tails.

**Editorial note on the abstract's kernel geomean (not an edit to the
abstract itself — quoted verbatim above):** "kernel geometric mean falls
from $2.560$ to $1.867$" is a **two-step** chain, not attributable to fix
(ii) alone: fix (i) (tolerance-matching) accounts for $2.560\to1.908$, and
fix (ii) (paired/interleaved timing) accounts for the remaining
$1.908\to1.867$ — both visible in the paper's own table (`sec:matched`,
`benchmarks/kernel_tolerance_sweep.py`, and `sec:bench`/Run C,
`benchmarks/revision_harness.py`; see `docs/paper-map.md` rows for
`sec:matched` and `sec:bench`) `[finite_diagnostic]`. Flagged by this
repo's own independent adversarial review — noted here rather than
silently reworded into the abstract.

---

## STATUS (read this before the numbers below)

**Overall benchmark verdict is HOLD, not ACCEPT.** The original cause — the
SciPy/ARPACK comparator (`eigsh`) failing to converge on
`morse_lambda5_all_bound`, which failed the fairness gate outright in Runs
A–C — has been diagnosed and fixed: `information-discrete-math` PR #109
(branch `fix/arpack-morse-shift-invert`, commit `72d5eec`, **open, not yet
merged**) switches the ARPACK comparator to shift-invert mode with a
Gershgorin-derived `sigma` computed from the operator alone (never from
native's own eigenvalues, so the comparator stays independent). Verified:
all seven declared cases now converge and cross-check, zero regression on
the six that already worked `[finite_diagnostic]`.

That closes the fairness gate deterministically — it does **not** close the
overall HOLD. The remaining and now sole cause is `factorized_sextic_ground`:
its speed vs. `SciPy eigh_tridiagonal` sits at genuine measurement parity.
Across 5 independent local runs (`audit_repeats=5`, the audit's default),
overall verdict was 3/5 ACCEPT, 2/5 HOLD; raising to `audit_repeats=30`
improves this to 4/5 ACCEPT but does not eliminate the flip
`[finite_diagnostic]`. CPU instruction-count measurement (`perf stat`,
repeated trials, median reported) shows why: native uses ~1.9–2.3x *more*
instructions per call than SciPy on this exact case, yet is marginally
faster in wall-clock — the edge is a microarchitectural
execution-efficiency effect, not fewer operations, which is exactly why it
sits at the noise floor. See `docs/paper-map.md`'s closing section for the
fuller reading (**corrected once** — the mechanism at play is narrower than
first stated): this architecture retains information two distinct ways —
mesh-level bracket carry-over across refinement levels (runs at any `k`,
including `k=1`) and cross-mode batching, evaluating multiple requested
eigenvalue indices together (only this needs `k>1` to have a second index to
batch with). `factorized_sextic_ground` and `pure_quartic_ground` both
request `k=1`, so cross-mode batching has nothing to act on there — and the
mesh-level retention that does still run does not visibly pay for itself
(native uses more instructions, not fewer). Near-parity here is consistent
with that reading, not an anomaly. No result in this repo overrides the
HOLD, and no sentence here claims the `k=1` cases as evidence for
cross-mode batching, in either direction.

**Known live bug — domain-truncation gate can silently miss a well.** On the
symmetric double well $V=(x^2-9)^2$, $k=4$, the planner locates one well,
builds a truncation window that excludes the other, and all three pipelines
return status ACCEPT with a passing diagnostic bound ($1.58\times10^{-9}$)
against an actual error of **16.07**. `benchmarks/planner_vs_solver.py`
isolated this as a planner/window-selection bug, not a solver/pivot bug — the
Sturm kernel is correct once handed the right window. The paper states a
repair (require the classically-forbidden region to extend to infinity, not
just pass a local boundary test) but did not implement it.
`benchmarks/weyl_gate.py` implements that repair as a **reference/diagnostic
scan only** — it is not wired into the `retained_spectral` engine's planner,
and its own docstring discloses a residual gap: a classically-allowed
component narrower than the scan's grid spacing can still be missed. The
production fix inside the engine itself remains **[Open]**.

---

## Tier legend

| Tag | Meaning |
|---|---|
| `[Th_coqc]` | Machine-checked theorem, proved elsewhere in `information-discrete-math` — cited here, not reproved here. |
| `[finite_diagnostic]` | Measured on a finite, disclosed test set — a readout, not a universal guarantee. |
| `[Open]` | Known unresolved question or an identified fix not yet implemented. |
| `HOLD` | The benchmark's own declared gate verdict — a fairness or speed gate failed; distinct from any external judgement. |

Every numeric or status claim below carries one of these tags inline.

---

## What this repo claims

The headline claim is **auditable completeness**, not speed. Sturm/inertia
counting on a symmetric tridiagonal matrix is a global sign-count over the
whole domain; it structurally cannot silently skip an eigenvalue the way a
shooting-method solver (Matslise, SLEDGE, SLEIGN2) can step over a close
doublet. This connects to — and cites, without re-deriving — an already
machine-checked `[Th_coqc]` result elsewhere in `information-discrete-math`:
`ResolvedInertia` / `IDM_ResolvedCount.v`, a forced three-value
certain+/certain-/unresolved-⊥ readout. `benchmarks/planner_vs_solver.py`
demonstrates the practical side of this directly `[finite_diagnostic]`: fed
the correct window by hand on the double-well case, Matslise returns only 2
of the 4 requested eigenvalues, while the Sturm kernel returns all 4 correctly
at error $2.05\times10^{-9}$.

On the seven declared spectra, retention reduces the eigensolve count from 12
to 7 on the four-level cases `[finite_diagnostic]`, and the reported
diagnostic bound exceeds the audited error on all seven instances by factors
of 7.7–21.4 `[finite_diagnostic]`.

**The `12 → 7` cross-mode-batching claim is scoped to `k>1`** (corrected once
— see `docs/paper-map.md`'s closing section for the full correction record).
This architecture retains information two distinct ways: mesh-level bracket
carry-over across refinement levels (runs at any `k`, including `k=1`), and
cross-mode batching, evaluating multiple requested eigenvalue indices
together (the `12 → 7` reduction above) — which needs `k>1` to have a second
index to batch with. At `k=1` (`factorized_sextic_ground`,
`pure_quartic_ground`) cross-mode batching has nothing to act on, and
instruction-count measurement shows the mesh-level retention that *does*
still run there does not visibly pay for itself (native uses *more* CPU
instructions than SciPy on both `k=1` cases). Neither `k=1` case is cited
anywhere in this repo as evidence for cross-mode batching, in either
direction (see `docs/paper-map.md`'s closing section — untiered, an
interpretive reading, not a `[Th_coqc]` proof).

## What this repo does not claim

**Speed vs. Matslise is mixed, not a win.** `[finite_diagnostic]`: RMS is
1.33–1.39x slower on the harmonic family, 7.6x slower on $|x|$, and 13x
slower on Pöschl–Teller; RMS is 8x faster only on `pure_quartic`. No sentence
in this repo states "faster" without that qualification. Note also the
asymmetry in that comparison: Matslise is handed the window RMS itself
derived, so RMS's own timing includes well search, gate evaluation, and the
window witness that Matslise's does not.

Against the (fairer, retention-ported) SciPy end-to-end comparator, the
honest end-to-end geometric mean is 2.411x, not the original 3.279x
`[finite_diagnostic]`, with one instance (`factorized_sextic_ground`) a
measured loss (ratio 0.986) and a second a tie — which is exactly why the
project's own speed verdict is HOLD, not ACCEPT. (This 0.986 is the paper's
own **end-to-end** Run D number, a different measurement field from the
**kernel-only** executor-audit CI discussed in the STATUS section above and
in `docs/paper-map.md` — both point at `factorized_sextic_ground`, but they
are not the same measurement; do not conflate them.) On the adversarial
tridiagonal suite, the native kernel is measurably slower on glued Wilkinson
matrices (ratio 0.478 at $k=64$, 0.326 at $k=168$) because DSTEBZ can
block-split and the native kernel cannot `[finite_diagnostic]`.

## Correction history

An early pass at the Matslise comparison used three unpaired samples from a
single launch and reported RMS as "2–10x slower" across the board — the
exact design defect (unpaired, non-interleaved timing) that this same work
criticises in the original SciPy comparator, committed in the measurement
that criticised it. A later paired/interleaved comparison
(`benchmarks/carvscar.py`, `benchmarks/aggregate_carvscar.py`) superseded it
with the mixed picture stated above: 1.33–13x slower on smooth potentials, 8x
faster on `pure_quartic`. This correction is recorded here, not overwritten
or hidden, and is treated as evidence that the tier discipline works — it
caught its own earlier error the same way it caught the comparator's tolerance
mismatch and the missing-retention arm.

---

## Repo layout

- [`paper/latest/`](paper/latest/) — symlink to the current paper revision.
  Currently `v2` (`rms_harmonic.tex` only — not yet compiled to PDF on a
  machine with the full LaTeX toolchain; see `docs/paper-map.md`). `v2` adds
  a `sec:v2repair` subsection reporting the ARPACK fix and the corrected
  `k>1` speed-claim scope on top of `v1`'s unchanged text; `v1`
  (`rms_harmonic.tex` + `rms_harmonic.pdf`) remains available at
  [`paper/v1/`](paper/v1/) unchanged.
- [`benchmarks/`](benchmarks/) — all benchmark scripts, their own `README.md` with per-script science status, and committed adversarial-suite result JSON under `benchmarks/results/`.
- [`patches/`](patches/) — the single Run C fix (`revision_major.patch`) applied to the sibling engine checkout, with its own `README.md` rationale.
- [`docs/`](docs/) — [`paper-map.md`](docs/paper-map.md) (every paper claim traced to its source script/tier) and [`reproducibility.md`](docs/reproducibility.md) (full setup/run instructions, pinned engine commit).

## How to reproduce

See [`docs/reproducibility.md`](docs/reproducibility.md) for the full
setup — the sibling `information-discrete-math` engine dependency, the pinned
commit the paper's numbers were measured against, per-script path-resolution
conventions, and exact run commands for every result cited above.

## License

This repo splits license by content type: code (`benchmarks/`, other
scripts) is MIT; the paper and prose documentation (`paper/`, and `docs/`/
`environment/` non-script files) is CC-BY-4.0. See [`LICENSE`](LICENSE) for
the exact split and full text.

The sibling dependency `information-discrete-math` (which provides
`retained_spectral`, the engine these benchmarks import — see
`docs/reproducibility.md`) is separately licensed **MIT**, same author. No
license conflict: both code trees are MIT. This repo does not vendor or
redistribute that code, only depends on it at a pinned commit.

## Citation

See [`CITATION.cff`](CITATION.cff).
