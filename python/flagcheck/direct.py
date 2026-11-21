"""Evaluate flag-algebra expressions in one finite graph G, straight from the definitions.

This is a third, independent implementation used by `verify` to spot-check constraint
coefficients. It never forms products or averages symbolically. For an expression e with
at most |G| vertices, the coefficient of G in lift(e, |G|) is

  * p(H, G) for a graph H: the probability that a uniform |H|-subset of V(G) induces H;
  * for [[R]], R a rooted expression on s roots: the mean over all injective
    theta: [s] -> V(G) of R(G, theta);
  * for a product F_1 ... F_t of theta-rooted flags: the probability that uniform
    disjoint subsets X_1, ..., X_t of V(G) minus theta, of sizes |F_i| - s, make every
    (G[theta + X_i], theta) isomorphic to F_i (the chain rule plus Thm 3).

Sums and scalings are linear; a product of sums is expanded.
"""

from fractions import Fraction
from itertools import combinations, permutations
from math import comb, perm
from typing import Dict, List, Tuple

from .graphs import Flag, canonical, restrict
from .spec import Add, Avg, Lit, Mul, Scale

Monomial = Tuple[Fraction, Tuple[Flag, ...]]


def expand(e) -> List[Monomial]:
    """A rooted expression as a list of (coefficient, product of flags)."""
    if isinstance(e, dict):
        return [(Fraction(c), (F,)) for F, c in e.items()]
    if isinstance(e, Lit):
        return [(Fraction(1), (e.flag,))]
    if isinstance(e, Add):
        return [m for a in e.args for m in expand(a)]
    if isinstance(e, Scale):
        return [(Fraction(e.q) * c, fl) for c, fl in expand(e.arg)]
    if isinstance(e, Mul):
        out = [(Fraction(1), ())]
        for a in e.args:
            out = [(c1 * c2, f1 + f2) for c1, f1 in out for c2, f2 in expand(a)]
        return out
    if isinstance(e, Avg):
        raise ValueError("averaging inside a rooted expression is not supported")
    raise TypeError(f"not an expression: {e!r}")


class GraphEvaluator:
    def __init__(self, G: Flag):
        if G.s != 0:
            raise ValueError("G must be unlabelled")
        self.G = G
        self.n = G.k
        self._canon: Dict[tuple, Flag] = {}

    def _induced(self, theta: tuple, X: tuple) -> Flag:
        key = (theta, X)
        F = self._canon.get(key)
        if F is None:
            F = canonical(restrict(self.G, theta + X, len(theta)))
            self._canon[key] = F
        return F

    def _monomial(self, theta: tuple, flags: Tuple[Flag, ...]) -> Fraction:
        s = len(theta)
        targets = [canonical(F) for F in flags]
        sizes = [F.k - s for F in flags]
        if any(F.s != s for F in flags):
            raise ValueError("flags in a product must have the same number of roots")
        rest = tuple(v for v in range(self.n) if v not in theta)
        if sum(sizes) > len(rest):
            raise ValueError("expression does not fit in G")
        total = 1
        left = len(rest)
        for a in sizes:
            total *= comb(left, a)
            left -= a

        def count(i: int, avail: tuple) -> int:
            if i == len(flags):
                return 1
            hits = 0
            for X in combinations(avail, sizes[i]):
                if self._induced(theta, X) == targets[i]:
                    remaining = tuple(v for v in avail if v not in X)
                    hits += count(i + 1, remaining)
            return hits

        return Fraction(count(0, rest), total)

    def rooted(self, e, theta: tuple) -> Fraction:
        return sum((c * self._monomial(theta, fl) for c, fl in expand(e)), Fraction(0))

    def value(self, e) -> Fraction:
        """Coefficient of G in lift(e, |G|) for an unlabelled expression e."""
        if isinstance(e, Add):
            return sum((self.value(a) for a in e.args), Fraction(0))
        if isinstance(e, Scale):
            return Fraction(e.q) * self.value(e.arg)
        if isinstance(e, Avg):
            monos = expand(e.arg)
            s = {F.s for _, fl in monos for F in fl}
            if len(s) != 1:
                raise ValueError("mixed numbers of roots under [[.]]")
            s = s.pop()
            total = Fraction(0)
            for theta in permutations(range(self.n), s):
                total += sum((c * self._monomial(theta, fl) for c, fl in monos), Fraction(0))
            return total / perm(self.n, s)
        # unlabelled flags, combinations or products of them: zero roots
        return self.rooted(e, ())
