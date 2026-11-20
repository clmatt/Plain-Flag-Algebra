"""Exact limit densities of step graphons (weighted blow-ups), in rational arithmetic.

A step graphon has parts 0..p-1 with rational weights summing to 1, and a symmetric
colour matrix C (p x p) that includes the diagonal: C[a][a] is the colour between two
distinct points of the same part. Sampling k points i.i.d. from the weights and colouring
each pair with C gives the induced density exactly in the limit, with no o(1) terms.
This is the ground truth the flag-algebra identities must reproduce exactly.

Ramsey blow-ups: parts are the vertices of a finite graph, C[a][a] = 0 (green/clone),
C[a][b] = 1 (red edge) or 2 (blue non-edge).
"""

from fractions import Fraction
from functools import lru_cache
from itertools import product
from typing import Dict, Sequence

from .graphs import Flag, canonical, from_color_function


class StepGraphon:
    def __init__(self, weights: Sequence, colors: Sequence[Sequence[int]]):
        self.w = tuple(Fraction(x) for x in weights)
        if sum(self.w) != 1 or any(x < 0 for x in self.w):
            raise ValueError("weights must be non-negative and sum to 1")
        p = len(self.w)
        self.C = tuple(tuple(row) for row in colors)
        if len(self.C) != p or any(len(r) != p for r in self.C):
            raise ValueError("colour matrix has the wrong shape")
        for a in range(p):
            for b in range(p):
                if self.C[a][b] != self.C[b][a]:
                    raise ValueError("colour matrix must be symmetric")
        self.p = p
        self._graph_cache: Dict[int, dict] = {}
        self._flag_cache: dict = {}

    @classmethod
    def blowup(cls, graph: Flag, weights=None, clone_color: int = 0) -> "StepGraphon":
        """Blow up a finite coloured graph: each vertex becomes a part whose internal colour
        is `clone_color`. With weights=None the parts are balanced."""
        k = graph.k
        if weights is None:
            weights = [Fraction(1, k)] * k
        C = [[clone_color if a == b else graph.color(a, b) for b in range(k)] for a in range(k)]
        return cls(weights, C)

    def _sample_flag(self, parts: tuple, s: int) -> Flag:
        return from_color_function(len(parts), lambda i, j: self.C[parts[i]][parts[j]], s=s)

    def graph_distribution(self, k: int) -> dict:
        """{canonical graph on k vertices: limit density}."""
        if k not in self._graph_cache:
            dist: dict = {}
            for parts in product(range(self.p), repeat=k):
                weight = Fraction(1)
                for a in parts:
                    weight *= self.w[a]
                if weight == 0:
                    continue
                G = canonical(self._sample_flag(parts, 0))
                dist[G] = dist.get(G, Fraction(0)) + weight
            self._graph_cache[k] = dist
        return self._graph_cache[k]

    def density(self, G: Flag) -> Fraction:
        return self.graph_distribution(G.k).get(canonical(Flag(G.k, 0, G.cols)), Fraction(0))

    def flag_distribution(self, root_parts: tuple, m: int) -> dict:
        """Roots placed in the given parts (one root per entry). Returns
        {canonical flag on m vertices: limit density} for flags rooted at those points."""
        key = (root_parts, m)
        if key not in self._flag_cache:
            s = len(root_parts)
            dist: dict = {}
            for extra in product(range(self.p), repeat=m - s):
                weight = Fraction(1)
                for a in extra:
                    weight *= self.w[a]
                if weight == 0:
                    continue
                F = canonical(self._sample_flag(tuple(root_parts) + extra, s))
                dist[F] = dist.get(F, Fraction(0)) + weight
            self._flag_cache[key] = dist
        return self._flag_cache[key]

    def flag_density(self, F: Flag, root_parts: tuple) -> Fraction:
        return self.flag_distribution(tuple(root_parts), F.k).get(canonical(F), Fraction(0))

    def evaluate(self, combo: dict, root_parts: tuple = ()) -> Fraction:
        """Value of a linear combination of flags (all with s = len(root_parts))."""
        total = Fraction(0)
        for F, v in combo.items():
            if F.s != len(root_parts):
                raise ValueError("root count mismatch")
            total += v * self.flag_density(F, root_parts)
        return total

    def root_placements(self, s: int):
        """Yield (root_parts, probability) for s i.i.d. roots."""
        for parts in product(range(self.p), repeat=s):
            weight = Fraction(1)
            for a in parts:
                weight *= self.w[a]
            if weight:
                yield parts, weight
