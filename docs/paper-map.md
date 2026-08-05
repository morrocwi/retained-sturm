# Paper map — `rms_harmonic.tex` → source scripts and data

Readout, not a certification: this table cross-references what the paper
(`paper/v1/rms_harmonic.tex`, unchanged; `paper/v2/rms_harmonic.tex` adds
`sec:v2repair` on top of the same v1 text — see `paper/latest`, currently
`v2`) claims to the actual script(s)/file(s) that produced each number, so a
claim can be re-derived rather than taken on faith. "Paper section" gives
the `\label{}` where one exists. Tier tags follow `benchmarks/README.md`:
`[Th_coqc]` machine-checked elsewhere, `[finite_diagnostic]` measured on a
finite test set, `[Open]` known unresolved question. HOLD/ACCEPT is the
benchmark's own gate verdict, not an external judgement.

**Overall status, updated:** v1's declared verdict was **HOLD**. Both
underlying causes have since been addressed in `information-discrete-math`:
1. The ARPACK/Morse non-convergence (`sec:fairness`) — shift-invert with a
   Gershgorin-derived `sigma` (computed from the operator alone, never from
   native's own eigenvalues — stays non-circular), verified on all seven
   cases, **merged to `main`** (PR #109, commits `72d5eec`/`acf4b51`).
2. The remaining `factorized_sextic_ground` speed parity — traced to `k=1`
   (see `sec:v2repair` in `paper/v2` and the closing section of this file);
   the strict speed gate is now scoped to `k>1` cases, with `k=1` reported
   separately via a CPU instruction-count instrument. PR #110 (commits
   `43c0ac4`/`3823245`) passed two rounds of independent review, all CI
   checks green, and **is now merged to `main`** (fast-forward, no conflicts).

Under the corrected `k>1` scope, the overall verdict was **ACCEPT in five
independent runs out of five** in local testing (previously 3/5 ACCEPT, 2/5
HOLD under the original all-seven-case scope) — a narrower claim than v1's
original all-seven-case ACCEPT, not the same claim reasserted. `k=1` results
are not overturned; they are excluded from a claim they cannot support.

| Paper section | Claim/Table | Source script(s) | Data file | Tier |
|---|---|---|---|---|
| Abstract / `sec:intro` Contributions | Headline claim: auditable completeness (Sturm/inertia counting cannot silently skip a root) | not reproved in this repo — cited from `information-discrete-math`: `ResolvedInertia` / `IDM_ResolvedCount.v` | — | `[Th_coqc]` (cited, not reproduced here) |
| `sec:reference` | Exact harmonic spectrum $E_n = n+\tfrac12$ used only as post-hoc audit, quarantined from the planner | analytic derivation in the text itself; no script | — | analytic (exact, not measured) |
| `sec:well`, `sec:gate`, `sec:discretisation`, `sec:sturm`, `sec:brackets`, `sec:mesh`, `sec:retain`, `sec:gates` | Method description: raw-input window/mesh/bracket planning, retained batched bisection, two-gate acceptance | engine implementation, not this repo — `information-discrete-math/retained_spectral/engine.py` (see `docs/reproducibility.md` for the pinned commit) | — | method definition, not a measurement |
| `sec:results`, Table `tab:results` | Released `harmonic_low4` eigenvalues, post-hoc audit against $E_n=n+\tfrac12$, $1.465\times10^{-10}$ at $\eps=2\times10^{-8}$ | Run A: recorded committed results (`aggregate_launches.py --results …` reads this file, does not regenerate Run A) | `information-discrete-math/retained_spectral/results/competition_results.json` | `[finite_diagnostic]` (recorded run) |
| `sec:boundsound`, Table `tab:bound` | Reported diagnostic bound exceeds audited error on all seven instances, factors 7.7–21.4 | `benchmarks/aggregate_launches.py` | `competition_results.json` (+ launch JSONs it pools) | `[finite_diagnostic]` |
| `sec:suite`, Table `tab:suite` | The seven declared spectra; solve counts, RMS vs. independent | engine `raw_benchmark_targets` (in `information-discrete-math`); pooled by `benchmarks/aggregate_launches.py` | `competition_results.json` | `[finite_diagnostic]` |
| `sec:work` | Retention reduces 12 solves to 7 on the four-level cases | same source as `tab:suite`; described, not a separate script | `competition_results.json` | `[finite_diagnostic]` |
| `sec:bench` "Three runs" / "Two fields" | Run design (A recorded, B/C/D executed for this paper) | `patches/revision_major.patch` (Run C's single change); `benchmarks/revision_harness.py` (Run D paired/interleaved harness) | — | design description |
| `sec:matched` | Same-work tolerance-mismatch diagnosis; matched quantity `max(0.01*tol, 2e-12)`; kernel geomean $2.560\to1.908$ | `benchmarks/kernel_tolerance_sweep.py`; the fix itself is `patches/revision_major.patch` | script output (not committed as JSON in this repo) | `[finite_diagnostic]` |
| Run C, Table `tab:kernel` | Field A kernel-only ratios on the replication host, matched width | `benchmarks/revision_harness.py` (Field A path) | `launch{1,2,3}.json` (produced by running the harness; not pre-committed in this repo) | `[finite_diagnostic]` |
| `sec:scaling`, Table `tab:scaling` | Kernel advantage flat 1.85–2.56 from $n=767$ to $n=786{,}431$, refutes wrapper-overhead objection | `benchmarks/adv_scaling.py` | script output | `[finite_diagnostic]` |
| `sec:paired` Run D, `sec:thirdarm`, Table `tab:e2e` | Paired/interleaved timing, three end-to-end arms (shipped / retained / retention-alone); one measured **loss** on `factorized_sextic_ground` (ratio 0.986); end-to-end geomean $3.279\to2.411$; retention alone accounts for 26.5% of the original margin | `benchmarks/revision_harness.py` (Field B, three arms; the retained arm is `benchmarks/scipy_retained.py`, which must be installed into `retained_spectral/competition/`); pooled by `benchmarks/aggregate_launches.py` | `launch{1,2,3}.json` → `--results competition_results.json` | `[finite_diagnostic]` |
| `sec:thirdarm` (supporting estimate) | Solve-by-solve estimate: 26.9% of comparator time is the rebuilt Richardson ladder; projected ratio $4.565\to3.382$ on `harmonic_low4` (later confirmed by the direct Run D measurement, 26.5%) | `benchmarks/adv_strawman.py` | script output | `[finite_diagnostic]` |
| `sec:fairness` | HOLD verdicts traced to cause: Runs A–C fail on the fairness gate (ARPACK `ArpackNoConvergence` on `morse_lambda5_all_bound`); Run D fairness would pass but speed then fails on `factorized_sextic_ground` | `benchmarks/revision_harness.py` / `aggregate_launches.py` (verdict fields in pooled output) | `competition_results.json` | `[finite_diagnostic]` |
| `sec:fairness` (update, **now merged**) | ARPACK non-convergence root-caused (plain `which="SA"`, no shift-invert, fails on a wide window / 5 requested modes) and fixed via shift-invert with a Gershgorin-derived `sigma`; all seven cases now converge and cross-check, zero regression on the six that already worked. Does NOT by itself flip the all-seven-case verdict — under that original scope `factorized_sextic_ground`'s speed vs. `SciPy eigh_tridiagonal` is at genuine measurement parity (3/5 ACCEPT, 2/5 HOLD across 5 independent runs at `audit_repeats=5`; raising to `audit_repeats=30` improves to 4/5 ACCEPT but does not eliminate the flip); PR #110 (row below) resolves this by rescoping the speed gate to `k>1`, giving 5/5 ACCEPT | `information-discrete-math` PR #109, branch `fix/arpack-morse-shift-invert`, commit `72d5eec`, merge commit `a05841d` — **merged to `main`**, edits `retained_spectral/competition/executor_audit.py` in the sibling repo, not this one | manual verification only (not yet a committed JSON in either repo) | `[finite_diagnostic]` |
| — (diagnosis of the sextic parity, IDM-grounded) | Root cause is *not* an algorithmic deficiency: CPU instruction-count measurement (`perf stat -e instructions:u`, two independent slope estimates) shows native uses ~1.86–2.01x MORE instructions per call than `SciPy eigh_tridiagonal` on this exact case, yet is slightly faster in wall-clock (~5–10%) — the wall-clock edge is a microarchitectural execution-efficiency effect (higher instructions-per-cycle), not fewer operations, which is why it sits at the noise floor and flips under repeated measurement | ad hoc instrumentation (`/tmp/instr_native.py` / `/tmp/instr_scipy.py` + `perf stat`, ephemeral, not committed to either repo) | perf counter output, recorded in conversation only | `[finite_diagnostic]` |
| `sec:advmatrix`, Table `tab:adversarial` | Adversarial tridiagonal suite (Wilkinson $W^+_{21}$, $W^+_{401}$, glued Wilkinson, Toeplitz, extreme-scale, denormal pivot floor, zero-coupling cluster) vs. dense `eigvalsh`; glued Wilkinson ratio falls to 0.478 at $k=64$, 0.326 at $k=168$ (native kernel ~3x slower — DSTEBZ block-splits, native kernel cannot) | `benchmarks/adv_kernel.py` (k=4), `benchmarks/adv_kernel2.py` (k up to all eigenvalues) | `benchmarks/results/adversarial_kernel.json`, `benchmarks/results/adversarial_kernel2.json` | `[finite_diagnostic]` |
| `sec:advnonsmooth` | Non-smooth $V=\lvert x\rvert$, Airy-zero reference; correct but costly (22 solves vs 7, 21.2ms vs 4.4ms) | reference spectrum: `benchmarks/dvr_reference.py` (Airy zeros); case added to engine by `patches/revision_major.patch` (`abs_linear`); run via `revision_harness.py` | `competition_results.json` / launch JSONs | `[finite_diagnostic]` |
| `sec:advfailure` | **Silent failure**: symmetric double well $V=(x^2-9)^2$, $k=4$ — planner window excludes second well; all three pipelines return status ACCEPT, bound $1.58\times10^{-9}\le\eps$, actual error **16.07** | independent reference: `benchmarks/dvr_reference.py` (sinc-DVR, Colbert–Miller, $N=6000$, double-well calibration sweep); isolation test: `benchmarks/planner_vs_solver.py` (hand-fed correct window, kernel returns error $2.05\times10^{-9}$ — confirms the bug is the window, not the solver; Matslise on the same window returns only 2 of 4 requested eigenvalues) | script output; also referenced in `information-discrete-math`'s `results/competition_results.json` for the ACCEPT-with-error-16.07 instance | `[finite_diagnostic]` |
| `sec:advfailure` (repair, **not implemented** in the engine) | Paper states the repair "is not difficult to state, though we have not implemented it": require the classically-forbidden region to extend to infinity (global scan / WKB action past $D_{\min}$), not just a local boundary test | `benchmarks/weyl_gate.py` — a **standalone reference/diagnostic** in *this* repo, written after the paper, implementing exactly the necessary-condition global scan the paper describes as unimplemented. It is NOT wired into `retained_spectral`'s planner and is NOT a certified Dirichlet-bracketing proof (its own docstring states the residual gap: a component of the classically-allowed region narrower than the scan grid spacing can still be missed). | `benchmarks/weyl_gate.py` (no persisted results file; run directly) | `[Open]` — production fix not yet in the engine; reference diagnostic only in this repo |
| `sec:novelty`, `sec:threats` | Scope of novelty / threats to validity discussion | prose only, cross-references the rows above | — | discussion |
| Reproducibility statement / Data & code availability | Run A read from committed results file; Runs B–D executed on the host of `tab:host`; single-line Run C change in the appendix | `patches/revision_major.patch`, `patches/README.md` | — | provenance statement |
| Appendix "Pipeline in pseudocode", "Recorded configuration", "Replication host" | Declared constants (`tab:config`), host environment for Run A (`tab:hostA`) and Runs B/C (`tab:host`) | none (recorded metadata, not a script) | — | recorded configuration |

## Material referenced in `benchmarks/README.md` but not (yet) cited in `rms_harmonic.tex` v1

These scripts exist in this repo and are documented in `benchmarks/README.md`,
but the paper text itself (checked by grep for "Matslise", "carvscar",
"pyslise") does not cite their numbers in v1 — they are additional benchmark
material, not paper claims:

| Script(s) | What it produces | Tier |
|---|---|---|
| `benchmarks/carvscar.py`, `benchmarks/aggregate_carvscar.py` | RMS vs. Matslise 2.0, paired/interleaved: 1.33–1.39x slower (harmonic family), 7.6x slower (`\|x\|`), 13x slower (Pöschl–Teller), **8x faster** (`pure_quartic`) — MIXED, not a win | `[finite_diagnostic]` |
| (superseded) earlier Matslise pass | "2–10x slower" from 3 unpaired samples, single launch — corrected by the paired result above; kept as a documented self-correction, not hidden | `[finite_diagnostic]` (superseded) |

If a future paper revision cites these numbers, add a row above pointing at
the same scripts.

## `k=1` and cross-mode batching (interpretive, not machine-checked — corrected once)

`factorized_sextic_ground` and `pure_quartic_ground` both request `k=1`
(a single eigenvalue). This section is an *interpretive argument*, tagged
distinctly from the empirical rows above — it is not a measurement and not a
`[Th_coqc]` proof, so it carries neither tier; treat it as a proposed reading,
open to revision.

**Correction record, stated first (tier-honesty, not hidden):** an earlier
version of this section claimed "the retention operator degenerates to the
identity on a singleton" — i.e., that *all* retention has nothing to act on
at `k=1`. Independent review read `information-discrete-math/retained_spectral/engine.py`
directly and found that is not accurate: `_native_mesh_readout` /
`_validated_brackets` carry a bracket found at one mesh-refinement level
forward as a hint to the next *regardless of how many modes are requested*
— that retention runs identically at `k=1`. What genuinely has nothing to
act on at `k=1` is a **different, narrower mechanism**: batching multiple
*requested eigenvalue indices* together in one traversal, so index 1's
bracket informs index 2's (`sec:retain`; the paper's own `12 → 7`
solve-count reduction is stated only for the four-level case, `sec:work` —
that reduction is this cross-mode batching, not mesh-level retention). The
corrected claim below reflects this; the over-broad original claim is kept
here, struck through in spirit, so the correction itself is traceable.

`information-discrete-math`'s own definition of retained information (δ_R,
`textbook/INFORMATION_DISCRETE_MATHEMATICS.md`) requires at least two related
readouts across which a distinction is retained and reused. Cross-mode
batching fits that definition and needs `k>1` to have a second member to
retain across; mesh-level bracket retention is a *different* application of
the same δ_R idea (across refinement levels, not across requested indices)
and applies at any `k`, including `k=1`.

Read this way, the near-parity measured on `factorized_sextic_ground` (and
the fact that `pure_quartic_ground`'s 8x margin, mixed as it is against
Matslise, is never attributed anywhere in this repo to cross-mode batching —
see "What this repo claims" in `README.md`) is not simply "nothing happened."
It is a more specific and more concerning finding: cross-mode batching has
nothing to act on at `k=1` as its own definition predicts, **and** the
mesh-level retention that *does* still run at `k=1` does not visibly pay for
itself there — instruction-count measurement shows native uses *more* CPU
instructions than SciPy on both `k=1` cases, not fewer, so the small
wall-clock edge that remains is a microarchitectural efficiency effect, not
an algorithmic one. Whether mesh-level retention's benefit itself scales
with the number of requested modes (plausible, structurally similar to
cross-mode batching's own scope limit) or is unrelated overhead is an open
question this repo does not resolve.

The corresponding claim-scope statement, honest under this reading: **the
cross-mode-batching claim (the `12 → 7` solve-count reduction) is scoped to
`k>1`**; at `k=1`, any speed result — win, loss, or tie — should be read as
a kernel-implementation/microarchitecture data point, not as evidence for or
against cross-mode batching, and says nothing about mesh-level retention's
value at `k>1` either way. Must not be captioned as if it were.
