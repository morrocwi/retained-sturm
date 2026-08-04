# `revision_major.patch` — what it actually fixed, and why it stays here

## What was broken

The kernel-level speed comparator (`executor_audit.py`) and the end-to-end
pipeline comparator (`scipy_pipeline.py`) were each handing the two solvers
under comparison a *different amount of numerical work*, then reporting the
resulting timing ratio as if it were a fair "same problem, same accuracy
target" comparison:

- The native Sturm/bisection kernel was called with `problem.tolerance` as
  its bisection half-width.
- The SciPy/LAPACK comparator (`sla.eigh_tridiagonal(..., select="i", ...)`)
  was called with `tol` left unset, which lets DSTEBZ fall back to its own
  machine-precision-based default half-width — a different, and generally
  much tighter, target than `problem.tolerance`.

Two solvers converging to two different half-widths are not doing the same
amount of work, so any speed ratio measured between them conflates "which
kernel is faster" with "which kernel was asked to try less hard." That is a
tolerance/bisection-half-width mismatch bug in the benchmark harness itself,
not in either solver.

## What the patch does (read from the diff)

1. Introduces one shared function, `_matched_width(problem)` /
   `matched_width(problem)`, defined identically in both
   `executor_audit.py` and the new `revision_harness.py`:
   `max(problem.tolerance * 0.01, 2e-12)`. This is the same half-width value
   the end-to-end engine already uses internally — the patch does not invent
   a new number, it makes the *comparator* honor the number the engine
   already commits to.
2. Passes that matched width to **both** sides of every timed comparison:
   the native kernel call (already was) and now also the
   `sla.eigh_tridiagonal(..., tol=_matched_width(problem),
   lapack_driver="stebz", ...)` call in `executor_audit.py` and in
   `scipy_pipeline.py`'s `_richardson_mesh`. Pinning `lapack_driver="stebz"`
   also removes driver-choice as a second silent variable.
3. Adds `revision_harness.py`, a new paired-timing driver that goes further
   than the tolerance fix alone: it interleaves the timed arms within every
   repeat (instead of timing one solver in a whole loop, then the other, so
   thermal/frequency drift cannot land unevenly on one arm), bootstraps
   confidence intervals over *paired* repeat indices, and — the part
   relevant to the "faster/slower" claim in the paper — adds a **third**
   end-to-end arm: `scipy_retained.py`, which ports the same window-retention
   idea RMS uses into the SciPy Richardson pipeline, changing nothing else
   about it. Without that third arm, a two-arm ratio cannot distinguish "the
   native kernel is faster" from "the comparator was never given the idea
   under test," which is the distinction the paper's speed claim depends on.
4. Adds two new potential cases to `engine.py` (`abs_linear`, a non-smooth
   `|x|` kink where the O(h^2) Richardson assumption is not justified, and
   `symmetric_double_well`, a clustered-spectrum case with exponentially
   small tunnelling splittings) — harder cases for the corrected harness to
   be re-run against, not softer ones.

## Why this stays here as a first-class artifact instead of being squashed

This repo's tier discipline (`[Th_coqc]` / `[finite_diagnostic]` / `[Open]`,
HOLD-not-ACCEPT where warranted) only means anything if the history of a
caught measurement bug is visible, not folded invisibly into a "clean"
benchmark script and forgotten. Keeping the raw diff:

- lets anyone re-derive, from the diff itself, exactly what changed between
  the pre-fix and post-fix timing numbers, instead of trusting a prose
  summary;
- is the same self-correction pattern already documented for the "2-10x
  slower vs Matslise" claim (an early unpaired-sample error later corrected
  by a fair paired/interleaved comparison, see `carvscar.py`,
  `aggregate_carvscar.py`, and the "Note on an earlier error" in the
  `rms_inspect14` README) — evidence the tier discipline catches and
  corrects its own errors, which is the actual source of this work's
  legitimacy, not any outside sign-off;
- keeps the mixed, qualified speed verdict honest: even after this fix,
  the native kernel is 1.33-13x *slower* than MATSLISE on smooth potentials
  (harmonic family, `|x|`, Poschl-Teller) and only 8x faster on
  `pure_quartic` — a tolerance-matching bug fix is not the same claim as
  "our kernel is fast," and the diff is the record that separates the two.
