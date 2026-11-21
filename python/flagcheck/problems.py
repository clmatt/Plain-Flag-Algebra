"""Standard problems as data, so every backend (Python reference, legacy C++, new C++)
can run exactly the same specification.

Every problem proves `objective . s >= d`. To bound max G, use objective = -G; the
proven upper bound is then -d.
"""

from dataclasses import dataclass, field
from fractions import Fraction
from typing import List, Optional, Tuple

from .algebra import Combo, combo
from .graphs import Flag, from_edges
from .ramsey import (
    GREEN,
    NUM_COLORS,
    degree_polynomial,
    green_clique_constraint,
    monochromatic_clique,
    ramsey_zeros,
)


@dataclass
class Problem:
    name: str
    n: int
    colors: int
    zeros: Tuple[Flag, ...]
    objective: Combo
    eqs: List[Tuple[Combo, Fraction]] = field(default_factory=list)
    ineqs: List[Tuple[Combo, Fraction]] = field(default_factory=list)
    include_empty_type: bool = False
    normalize: bool = False          # use for feasibility (infeasibility-proof) problems
    # Localising PSD blocks (m, g): [[ g v v^T ]] >= 0. Only valid when g >= 0 pointwise.
    localizing: List[Tuple[int, object]] = field(default_factory=list)
    # Facts from outside flag algebra that the constraints rely on (printed by verify).
    assumptions: List[str] = field(default_factory=list)
    # N > 0: eliminate the clone variables of an N-vertex blow-up (flagcheck.blowup); the
    # SDP is then over the green-free graphs only.
    blowup: int = 0
    blowup_levels: bool = False      # with blowup: the sparse "levels" form (flagcheck.blowup)


K3 = from_edges(3, [(0, 1, 1), (0, 2, 1), (1, 2, 1)])
E3 = from_edges(3, [])
EDGE = from_edges(2, [(0, 1, 1)])
P3 = from_edges(3, [(0, 1, 1), (1, 2, 1)])
P4 = from_edges(4, [(0, 1, 1), (1, 2, 1), (2, 3, 1)])
C5 = from_edges(5, [(0, 1, 1), (1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 0, 1)])


def goodman(n: int = 3) -> Problem:
    """min K3 + E3 = 1/4 (Goodman)."""
    return Problem("goodman", n, 2, (), combo((K3, 1), (E3, 1)))


def mantel(n: int = 3) -> Problem:
    """Triangle-free => edge density <= 1/2. Proves -edge >= d, so the bound is -d."""
    return Problem("mantel", n, 2, (K3,), combo((EDGE, -1)))


def inducibility(G: Flag, n: int, name: str = "") -> Problem:
    """max p(G); the bound is -d."""
    return Problem(name or f"ind{G.k}", n, 2, (), combo((G, -1)))


def ramsey_objective(k: int, l: int, n: int) -> Problem:
    """Lidicky-Pfender / Appendix A: min p(green K2) over blow-ups of (k,l)-graphs.
    A (k,l)-graph on N vertices has a balanced blow-up with green density 1/N, so
    d > 1/N proves R(k,l) <= N, i.e. R(k,l) <= ceil(1/d) when 1/d is not an integer."""
    return Problem(f"R{k}{l}-obj", n, NUM_COLORS, tuple(ramsey_zeros(k, l)),
                   combo((monochromatic_clique(2, GREEN), 1)))


def ramsey_feasibility(k: int, l: int, N: int, n: int, green_powers: Optional[int] = None,
                       degree_window: Optional[Tuple[int, int]] = None) -> Problem:
    """Is there a (k,l)-graph on exactly N vertices? Constraints: p(green K_j) = 1/N^(j-1)
    for j = 2..green_powers (default n), optionally with every red degree in
    [lo, hi] (as [[prod_{d=lo..hi} (R - d/N)]] = 0, when it fits in n vertices).
    A certificate with d > 0 for objective 0 proves R(k,l) <= N."""
    zeros = tuple(ramsey_zeros(k, l))
    top = n if green_powers is None else green_powers
    eqs = [green_clique_constraint(j, N) for j in range(2, top + 1)]
    if degree_window is not None:
        lo, hi = degree_window
        if hi - lo + 2 > n:
            raise ValueError("degree polynomial needs more vertices than n")
        eqs.append((degree_polynomial(range(lo, hi + 1), N, zeros), Fraction(0)))
    return Problem(f"R{k}{l}-N{N}", n, NUM_COLORS, zeros, {}, eqs=eqs, normalize=True)
