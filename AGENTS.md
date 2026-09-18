# AGENTS.md - retained-sturm

## What this repository is

Retained Multilevel Sturm (RMS) is a raw-input-planning, diagnostic-acceptance-gated architecture for
the lowest k eigenvalues of a one-dimensional Schroedinger operator, and a benchmark against itself
that was twice corrected against its own first pass. It is a single-author showcase repository
accompanying the RMS paper and its benchmark/reproducibility artifacts. The headline claim is
auditable completeness, not speed: speed vs. Matslise is mixed, not a win. The README records a known
live bug (the domain-truncation gate can silently miss a well; the production fix remains `[Open]`),
and the paper restricts the planner's honest scope to single-well potentials with monotone tails.

## Read first

1. `README.md` - read its STATUS section before the numbers; then "What this repo claims",
   "What this repo does not claim" and "Correction history".
2. `docs/paper-map.md` - every paper claim traced to its source script and tier.
3. `docs/reproducibility.md` - full setup and run instructions, pinned engine commit.
4. `CITATION.cff` - citation metadata.

## Rules

- Tier legend per `README.md`: `[Th_coqc]` (cited here, not reproved here), `[finite_diagnostic]`
  (a readout on a finite, disclosed test set, not a universal guarantee), `[Open]`, and `HOLD`
  (the benchmark's own declared gate verdict). Every numeric or status claim carries a tag inline.
- No sentence in this repo states "faster" without the qualification given in `README.md`.
- Corrections are recorded, not overwritten or hidden (`README.md`, "Correction history").
- Do not conflate the end-to-end Run D field with the kernel-only executor-audit verdict.
- Contributions: open an issue first to discuss a fix or improvement (`CONTRIBUTING.md`).
- Licence split per `LICENSE`: code (`benchmarks/`, other scripts) is MIT; the paper and prose
  documentation (`paper/`, and `docs/` / `environment/` non-script files) is CC-BY-4.0.

## Programme map

This repository is one node of the Human-AI Readout Programme. Which repository answers which kind of
question, what to read first and which gate applies is kept in one place, the routing hub:
<https://github.com/morrocwi/main.hub> (start at its `AGENTS.md`, then `ROUTES.md`).
The hub holds pointers and pinned links only. It is a readout of one moment: when the hub and this
repository disagree, this repository wins.
