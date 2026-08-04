"""weyl_gate.py -- a Dirichlet-monotonicity / component-coverage gate for the
domain-truncation acceptance decision, plus a semiclassical (Weyl) count kept
strictly as a logged, non-deciding diagnostic.

WHY THIS FILE EXISTS
---------------------
``planner_vs_solver.py`` in this same directory isolated a real, live bug in
the ``retained_spectral`` raw-input planner: on the symmetric double well
V(x) = lambda * (x^2 - a^2)^2 with lambda = 1, a^2 = 9, k = 4 requested modes,
the planner locates the well near x = -3 and returns the truncation window
(-6.43, +0.43) [source: rms_harmonic.tex, Sec. "A silent failure: the barrier
that looks like a tail"].  That window contains ONE of the two wells and
excludes the other entirely.  The existing decay gate and window-witness gate
both test a SUFFICIENT, LOCAL condition -- "is the boundary point deep in a
classically forbidden region?" -- and a central barrier between two wells
answers that question exactly as a decaying tail does.  Neither gate ever asks
the NECESSARY, GLOBAL question: does the classically allowed region
{x : V(x) < mu} have any component that pokes outside the chosen domain at
all?  The result: status ACCEPT, reported diagnostic bound 1.58e-9 <= 2e-8,
actual error 16.07 -- nine orders of magnitude above the declared tolerance.

WHAT THIS FILE IS, AND IS NOT
------------------------------
This is a REFERENCE / DIAGNOSTIC implementation of the necessary-condition
check only: it scans a finite, declared grid for connected components of
{x : V(x) < mu} and checks whether every component found lies inside the
declared domain (a, b).  That is useful and, as far as it goes, correct: it
would have caught the double-well failure above (see main() below).

It is NOT a certified Dirichlet-bracketing proof.  A genuine certificate would
require an actual Dirichlet eigenvalue solve on the restricted domain (a, b)
to prove monotonicity of the truncated spectrum with respect to further
domain widening -- that capability does not exist in retained_spectral today
and building it is a separate, not-yet-done engineering task in the engine
itself, out of scope for this standalone script.  What is implemented here is
the cheaper, weaker, but honest thing: "did we find a component of the
classically allowed region outside the declared domain, on a grid of this
declared resolution, under this declared (attempted) Lipschitz bound".  That
is a REAL improvement over the status quo (which asks nothing global at all)
but it is a diagnostic-grade, not proof-grade, certificate. Tier: this file's
own scan is [finite_diagnostic] research code, not [Th_coqc]. The theorem
that motivates why "sound Sturm/inertia counting cannot silently skip a root"
is a live goal for shooting-method solvers is a DIFFERENT, ALREADY
machine-checked result living elsewhere in information-discrete-math
(ResolvedInertia / IDM_ResolvedCount.v -- a forced three-value
certain+/certain-/unresolved-bottom readout). This file does not reprove that
theorem; it cites it as motivation for why global coverage, not just local
boundary decay, is the right thing to certify.

RESIDUAL GAP, STATED PLAINLY
------------------------------
The grid scan below has a real, undischarged gap: a component of
{x : V(x) < mu} that is narrower than the declared grid spacing `delta`, and
that falls entirely between two adjacent grid samples, can be missed. The
"Lipschitz-ish bound" L this file estimates is a finite-difference empirical
estimate of max|V'| over the scan range, NOT a certified analytic bound -- for
smooth, simple potentials it is a reasonable diagnostic, but it is not a
proof that no such narrow component exists between samples. Any (R, delta, L)
triple returned by ACCEPT should be read as "no component was found missing
under this declared, finite resolution", not as "no component can exist".

TIER SUMMARY
------------
- PRIMARY certificate (this file's core contribution): dirichlet_coverage_
  certificate() -- [finite_diagnostic], decides ACCEPT/HOLD, always reports
  the (R, delta, L) it was accepted or held under.
- SECONDARY diagnostic: weyl_semiclassical_count() -- [finite_diagnostic],
  logged only, NEVER overrides or confirms the primary verdict.
- Cited-not-reproved: [Th_coqc] ResolvedInertia / IDM_ResolvedCount.v in
  information-discrete-math (sibling repo), motivating why a structural
  completeness argument -- not a local boundary heuristic -- is the right
  shape of certificate to aim for.

Run (from the repository root)::

    PYTHONPATH=. python3 benchmarks/weyl_gate.py
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum
from typing import Callable, Optional, Sequence

import numpy as np

try:
    from scipy import integrate as _scipy_integrate
except ImportError:  # pragma: no cover - exercised on dependency-light CI
    _scipy_integrate = None


# ---------------------------------------------------------------------------
# PRIMARY certificate: Dirichlet-monotonicity / component-outside-domain gate
# ---------------------------------------------------------------------------


class GateVerdict(Enum):
    """The only two outcomes the primary certificate may return.

    No bare bool anywhere in this module represents an acceptance decision --
    every decision point returns one of these two tagged values, carrying the
    evidence it was decided under.
    """

    ACCEPT = "ACCEPT"
    HOLD = "HOLD"


@dataclass(frozen=True)
class Component:
    """One connected component of {x : V(x) < mu} found on the scan grid."""

    left: float
    right: float

    @property
    def width(self) -> float:
        return self.right - self.left

    def outside_amount(self, domain: tuple[float, float]) -> float:
        """How far this component pokes outside `domain`, 0.0 if fully inside."""

        a, b = domain
        return max(0.0, a - self.left, self.right - b)


@dataclass(frozen=True)
class DirichletCoverageCertificate:
    """The typed verdict of the primary gate, with its full evidence trail.

    `verdict` is the only field a caller should branch on. Every other field
    exists so the verdict can be audited independently of this script -- the
    declared (R, delta, L) triple this ACCEPT or HOLD was issued under.
    """

    verdict: GateVerdict
    domain: tuple[float, float]
    mu: float
    scan_range: tuple[float, float]       # the actual (a - R, b + R) scanned
    R: float                              # declared scan half-width margin
    delta: float                          # declared grid spacing
    lipschitz_bound: Optional[float]      # empirical max|V'| estimate, or None
    components: tuple[Component, ...]     # every component found on the scan
    offending: tuple[Component, ...]      # components that poke outside domain
    reason: str

    def report(self) -> str:
        lines = [
            f"verdict       : {self.verdict.value}",
            f"domain (a,b)  : {self.domain}",
            f"mu            : {self.mu:.10g}",
            f"scan range    : {self.scan_range}  (R = {self.R:.6g})",
            f"delta         : {self.delta:.6g}",
            f"lipschitz L   : {self.lipschitz_bound!r} "
            f"(empirical estimate, NOT a certified bound)",
            f"components    : {[(round(c.left, 6), round(c.right, 6)) for c in self.components]}",
        ]
        if self.offending:
            lines.append(
                "offending     : "
                + ", ".join(
                    f"[{c.left:.6g},{c.right:.6g}] "
                    f"({c.outside_amount(self.domain):.6g} outside domain)"
                    for c in self.offending
                )
            )
        lines.append(f"reason        : {self.reason}")
        return "\n".join(lines)


def find_sublevel_components(
    V: Callable[[np.ndarray], np.ndarray],
    scan_range: tuple[float, float],
    mu: float,
    delta: float,
) -> tuple[Component, ...]:
    """Grid-scan {x in scan_range : V(x) < mu} for connected components.

    This is the finite, declared-resolution readout the whole certificate
    rests on: `delta` is the actual spacing used, stated explicitly rather
    than left implicit. A component that starts and ends strictly between two
    adjacent grid points can be missed -- see the module docstring's
    "residual gap" section.
    """

    a, b = scan_range
    if not b > a:
        raise ValueError("scan_range must have positive length")
    n_points = max(int(math.ceil((b - a) / delta)) + 1, 2)
    xs = np.linspace(a, b, n_points)
    vs = np.asarray(V(xs), dtype=float)
    mask = vs < mu

    components: list[Component] = []
    start: Optional[float] = None
    for i, inside in enumerate(mask):
        if inside and start is None:
            start = float(xs[i])
        if not inside and start is not None:
            components.append(Component(start, float(xs[i - 1])))
            start = None
    if start is not None:
        components.append(Component(start, float(xs[-1])))
    return tuple(components)


def _empirical_lipschitz_bound(
    V: Callable[[np.ndarray], np.ndarray],
    scan_range: tuple[float, float],
    *,
    probe_points: int = 4001,
) -> float:
    """Finite-difference estimate of max|V'| over scan_range.

    Explicitly NOT a certified analytic bound -- see module docstring. This is
    a generic, callable-only estimate (no symbolic differentiation of V is
    assumed) with a small safety factor, offered only as a diagnostic value
    to log alongside the certificate, never as a proof input.
    """

    a, b = scan_range
    xs = np.linspace(a, b, probe_points)
    vs = np.asarray(V(xs), dtype=float)
    h = xs[1] - xs[0]
    if h <= 0:
        return float("nan")
    slopes = np.abs(np.diff(vs) / h)
    if slopes.size == 0:
        return 0.0
    return float(1.25 * np.max(slopes))  # 25% safety margin, empirical only


def dirichlet_coverage_certificate(
    V: Callable[[np.ndarray], np.ndarray],
    domain: tuple[float, float],
    mu: float,
    *,
    R: Optional[float] = None,
    delta: Optional[float] = None,
    estimate_lipschitz: bool = True,
) -> DirichletCoverageCertificate:
    """The primary ACCEPT/HOLD certificate.

    Scans {x : V(x) < mu} over (a - R, b + R) at grid spacing `delta` and
    checks whether every component found lies inside the declared domain
    (a, b). Returns HOLD, naming the offending component and how far outside
    it extends, the instant one is found; ACCEPT only if none is found on the
    declared scan -- carrying the (R, delta, L) triple it was accepted under.

    Defaults for R and delta are declared explicitly, not silently picked:
      R     defaults to max(4.0, 0.75 * (b - a)) -- a scan margin comparable
            to the domain width itself, wide enough in practice to reach a
            second well sitting just beyond a mis-placed boundary (as in the
            double-well case this file was written to catch), but still an
            arbitrary finite choice, not a proof it reaches every well.
      delta defaults to (b - a) / 4096 -- four times the mesh density used in
            the reference harmonic case in rms_harmonic.tex's own working
            window, chosen for readable component boundaries in the demo,
            not derived from any error-bound argument.
    """

    a, b = domain
    if not b > a:
        raise ValueError("domain must have positive length")
    width = b - a
    if R is None:
        R = max(4.0, 0.75 * width)
    if delta is None:
        delta = width / 4096.0

    scan_range = (a - R, b + R)
    components = find_sublevel_components(V, scan_range, mu, delta)
    offending = tuple(c for c in components if c.outside_amount(domain) > 0.0)

    lipschitz_bound = (
        _empirical_lipschitz_bound(V, scan_range) if estimate_lipschitz else None
    )

    if offending:
        detail = "; ".join(
            f"component [{c.left:.6g}, {c.right:.6g}] extends "
            f"{c.outside_amount(domain):.6g} outside domain {domain}"
            for c in offending
        )
        reason = (
            f"HOLD: {len(offending)} of {len(components)} component(s) of "
            f"{{x : V(x) < {mu:.6g}}} found on the scan range {scan_range} "
            f"(delta={delta:.6g}) extend outside the declared domain "
            f"{domain}. {detail}."
        )
        verdict = GateVerdict.HOLD
    else:
        reason = (
            f"ACCEPT: all {len(components)} component(s) of "
            f"{{x : V(x) < {mu:.6g}}} found on the scan range {scan_range} "
            f"(R={R:.6g}, delta={delta:.6g}) fall inside the declared domain "
            f"{domain}. This is a necessary-condition check on a finite grid, "
            f"not a certified Dirichlet-bracketing proof -- see module "
            f"docstring 'residual gap'."
        )
        verdict = GateVerdict.ACCEPT

    return DirichletCoverageCertificate(
        verdict=verdict,
        domain=domain,
        mu=mu,
        scan_range=scan_range,
        R=R,
        delta=delta,
        lipschitz_bound=lipschitz_bound,
        components=components,
        offending=offending,
        reason=reason,
    )


# ---------------------------------------------------------------------------
# SECONDARY diagnostic: semiclassical (Weyl) phase-space count -- LOGGED ONLY
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class WeylDiagnostic:
    """[finite_diagnostic] semiclassical cross-check, not a certificate.

    This value is computed over the SAME wide scan range as the primary
    certificate (never the operator's own, narrower window) and is logged
    next to the primary verdict for a human to read. It must never be used to
    override or confirm dirichlet_coverage_certificate()'s ACCEPT/HOLD.
    """

    n_weyl: float
    scan_range: tuple[float, float]
    mu: float
    n_reference: Optional[float] = None
    difference: Optional[float] = None
    note: str = "[finite_diagnostic] semiclassical cross-check, not a certificate."


def weyl_semiclassical_count(
    V: Callable[[np.ndarray], np.ndarray],
    scan_range: tuple[float, float],
    mu: float,
) -> float:
    """N_Weyl(mu) = (1/pi) * integral_scan_range sqrt(2 * max(mu - V(x), 0)) dx.

    Computed by numerical quadrature (scipy.integrate.quad if available,
    else a plain composite Simpson's rule) over `scan_range` -- the same wide
    range the primary certificate scans, not the operator's narrower window.
    """

    a, b = scan_range

    def integrand(x: float) -> float:
        return math.sqrt(2.0 * max(mu - float(V(np.array([x]))[0]), 0.0))

    if _scipy_integrate is not None:
        value, _ = _scipy_integrate.quad(integrand, a, b, limit=400)
    else:  # pragma: no cover - exercised only without scipy
        n = 20001
        xs = np.linspace(a, b, n)
        ys = np.array([integrand(x) for x in xs])
        h = (b - a) / (n - 1)
        # composite Simpson's rule, n-1 even required
        value = (h / 3.0) * (
            ys[0] + ys[-1] + 4.0 * ys[1:-1:2].sum() + 2.0 * ys[2:-1:2].sum()
        )
    return value / math.pi


def run_weyl_diagnostic(
    V: Callable[[np.ndarray], np.ndarray],
    scan_range: tuple[float, float],
    mu: float,
    *,
    n_reference: Optional[float] = None,
) -> WeylDiagnostic:
    n_weyl = weyl_semiclassical_count(V, scan_range, mu)
    diff = None if n_reference is None else n_weyl - n_reference
    return WeylDiagnostic(
        n_weyl=n_weyl,
        scan_range=scan_range,
        mu=mu,
        n_reference=n_reference,
        difference=diff,
    )


# ---------------------------------------------------------------------------
# Potentials -- mirrors the formulas in retained_spectral.engine.potential_values
# (sibling repo, information-discrete-math/retained_spectral/engine.py).
# Reference only: re-typed here, not imported, so this repo has no runtime
# dependency on the engine repo. If the two ever diverge, engine.py is the
# one that matters for retained_spectral's actual behaviour; these are only
# for demonstrating this gate on the same declared test potentials.
# ---------------------------------------------------------------------------


def harmonic(x: np.ndarray, *, omega: float = 1.0, center: float = 0.0) -> np.ndarray:
    return 0.5 * omega**2 * (x - center) ** 2


def pure_quartic(x: np.ndarray, *, coupling: float = 0.25) -> np.ndarray:
    return coupling * x**4


def symmetric_double_well(x: np.ndarray, *, lam: float = 1.0, a2: float = 9.0) -> np.ndarray:
    """V(x) = lambda * (x^2 - a^2)^2 -- the case that broke the planner.

    Source: rms_harmonic.tex, Sec. "A silent failure: the barrier that looks
    like a tail"; lam=1, a2=9 gives a doublet splitting of 8.3e-12, below the
    bisection half-width, and a barrier of V(0) = lam*a2^2 = 81.
    """

    return lam * (x**2 - a2) ** 2


# ---------------------------------------------------------------------------
# Demonstration cases
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class DemoCase:
    name: str
    V: Callable[[np.ndarray], np.ndarray]
    domain: tuple[float, float]
    mu: float
    expected: GateVerdict
    note: str
    R: Optional[float] = None
    delta: Optional[float] = None
    n_reference: Optional[float] = None


def _double_well_bad_window_case() -> DemoCase:
    # V = (x^2 - 9)^2, k = 4 requested modes.
    # planner-chosen window (-6.43, +0.43), reference spectrum
    # (4.2144398, 4.2144398, 12.5272184, 12.5272184); the planner's own
    # 4-mode readout on this bad window returns (4.2144398, 12.5272184,
    # 20.6574438, 28.5922506) at status ACCEPT -- error 16.07 against
    # tolerance 2e-8.  Numbers: rms_harmonic.tex Sec. "A silent failure:
    # the barrier that looks like a tail"; also referenced in
    # planner_vs_solver.py (REFERENCE array, CORRECT_WINDOW = (-8, 8)).
    return DemoCase(
        name="double_well_a2_9 (KNOWN BAD narrow window)",
        V=lambda x: symmetric_double_well(x, lam=1.0, a2=9.0),
        domain=(-6.43, 0.43),
        mu=28.5922506,  # the planner's own (wrong) 4th-mode estimate, used
                         # as the target level per this gate's design: mu is
                         # the highest REQUESTED level, exactly as the
                         # existing decay gate uses E_{k-1}.
        expected=GateVerdict.HOLD,
        note=(
            "The known live bug this file targets: the second well "
            "(minimum near x=+3) is entirely outside this domain."
        ),
        R=6.0,  # wide enough to reach the second well's component fully;
                 # see module docstring default-R note for why this is a
                 # declared, not derived, choice.
        delta=0.01,
    )


def _harmonic_working_window_case() -> DemoCase:
    # V = 0.5 * x^2, k = 4 requested modes, mu = E_3 = 3.5 (n=3, n+0.5).
    # Working (verification) window [-10, 10] and mu = E_3 are exactly the
    # numbers derived in rms_harmonic.tex Sec. "Well localisation and the
    # initial window" / "The decay gate and the boundary witness"
    # (eq:window1: [a_1, b_1] = [-10, 10]; E_3 ~ 3.5).
    return DemoCase(
        name="harmonic_low4 (working window from rms_harmonic.tex)",
        V=lambda x: harmonic(x, omega=1.0, center=0.0),
        domain=(-10.0, 10.0),
        mu=3.5,
        expected=GateVerdict.ACCEPT,
        note="Working case from the paper's own worked example; single well.",
        n_reference=4,
    )


def _quartic_working_window_case() -> DemoCase:
    # V = 0.25 * x^4, k = 1 requested mode, mu = reference ground state
    # 0.420804974478 (raw_benchmark_targets() in retained_spectral/engine.py,
    # "published numerical comparator"). The exact planner-selected window
    # for this case could not be independently re-derived from the source
    # material available to this script (it depends on the internal
    # curvature-based well-localisation rule, which this file does not
    # reproduce -- see module docstring on not copying engine internals).
    # The domain below is a manually chosen, generously wide window, picked
    # only to demonstrate ACCEPT on a single-well case with a non-quadratic
    # potential; it is NOT claimed to be the planner's actual window.
    return DemoCase(
        name="pure_quartic_ground (manually chosen wide window, NOT the "
        "planner's actual window -- see comment in source)",
        V=lambda x: pure_quartic(x, coupling=0.25),
        domain=(-6.0, 6.0),
        mu=0.420804974478,
        expected=GateVerdict.ACCEPT,
        note=(
            "Single well; domain picked generously wide by hand, not taken "
            "from a planner run (not independently reproducible from the "
            "source material this script had access to)."
        ),
        n_reference=1,
    )


DEMO_CASES: tuple[DemoCase, ...] = (
    _double_well_bad_window_case(),
    _harmonic_working_window_case(),
    _quartic_working_window_case(),
)


def main() -> int:
    all_ok = True
    for case in DEMO_CASES:
        print("=" * 78)
        print(case.name)
        print("-" * 78)
        print(f"note: {case.note}")

        cert = dirichlet_coverage_certificate(
            case.V, case.domain, case.mu, R=case.R, delta=case.delta
        )
        print()
        print("[PRIMARY certificate -- decides ACCEPT/HOLD]")
        print(cert.report())

        diag = run_weyl_diagnostic(
            case.V, cert.scan_range, case.mu, n_reference=case.n_reference
        )
        print()
        print("[SECONDARY diagnostic -- logged only, does not decide anything]")
        print(f"{diag.note}")
        print(f"N_Weyl(mu) = {diag.n_weyl:.6f}  over scan range {diag.scan_range}")
        if diag.n_reference is not None:
            print(
                f"N_reference (requested modes) = {diag.n_reference}  "
                f"difference = {diag.difference:.6f}"
            )

        status = "OK" if cert.verdict is case.expected else "MISMATCH"
        if status == "MISMATCH":
            all_ok = False
        print()
        print(
            f"result: verdict={cert.verdict.value} expected={case.expected.value} "
            f"[{status}]"
        )
        print()

    print("=" * 78)
    print("ALL DEMO CASES MATCHED EXPECTED VERDICT" if all_ok else "SOME CASES MISMATCHED -- SEE ABOVE")
    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
