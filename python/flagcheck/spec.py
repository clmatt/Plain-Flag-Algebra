"""Flag-algebra expressions and the spec format read by the C++ `flagsdp` tool.

An expression is a tree of Lit (flag literal), Add, Scale, Mul and Avg nodes. The same
tree can be evaluated by the Python reference (`evaluate`) or serialised for C++
(`sexpr`), so constraint construction can be cross-checked between implementations.
A plain Combo (dict) is accepted wherever an expression is expected.
"""

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable, Tuple, Union

from .algebra import Combo, add, average, combo, multiply, scale
from .graphs import Flag


@dataclass(frozen=True)
class Lit:
    flag: Flag


@dataclass(frozen=True)
class Add:
    args: Tuple


@dataclass(frozen=True)
class Scale:
    q: Fraction
    arg: object


@dataclass(frozen=True)
class Mul:
    args: Tuple


@dataclass(frozen=True)
class Avg:
    arg: object


Expr = Union[Lit, Add, Scale, Mul, Avg, dict]


def _flag_sexpr(F: Flag) -> str:
    return "(flag " + " ".join(str(x) for x in (F.k, F.s) + tuple(F.cols)) + ")"


def _rat(q) -> str:
    q = Fraction(q)
    return f"{q.numerator}/{q.denominator}" if q.denominator != 1 else str(q.numerator)


def sexpr(e: Expr) -> str:
    if isinstance(e, dict):
        if not e:
            raise ValueError("cannot serialise an empty combination")
        return "(add " + " ".join(f"(scale {_rat(v)} {_flag_sexpr(F)})" for F, v in e.items()) + ")"
    if isinstance(e, Lit):
        return _flag_sexpr(e.flag)
    if isinstance(e, Add):
        return "(add " + " ".join(sexpr(a) for a in e.args) + ")"
    if isinstance(e, Scale):
        return f"(scale {_rat(e.q)} {sexpr(e.arg)})"
    if isinstance(e, Mul):
        return "(mul " + " ".join(sexpr(a) for a in e.args) + ")"
    if isinstance(e, Avg):
        return f"(avg {sexpr(e.arg)})"
    raise TypeError(f"not an expression: {e!r}")


def evaluate(e: Expr, colors: int, zeros: Iterable = ()) -> Combo:
    """Evaluate with the Python reference (products keep only admissible flags)."""
    zeros = tuple(zeros)
    if isinstance(e, dict):
        return dict(e)
    if isinstance(e, Lit):
        return combo((e.flag, 1))
    if isinstance(e, Add):
        out: Combo = {}
        for a in e.args:
            out = add(out, evaluate(a, colors, zeros))
        return out
    if isinstance(e, Scale):
        return scale(evaluate(e.arg, colors, zeros), e.q)
    if isinstance(e, Mul):
        out = evaluate(e.args[0], colors, zeros)
        for a in e.args[1:]:
            out = multiply(out, evaluate(a, colors, zeros), colors, zeros)
        return out
    if isinstance(e, Avg):
        return average(evaluate(e.arg, colors, zeros))
    raise TypeError(f"not an expression: {e!r}")


def write_spec(problem, path: str) -> None:
    lines = [f"colors {problem.colors}", f"n {problem.n}",
             f"empty_type {1 if problem.include_empty_type else 0}"]
    for Z in problem.zeros:
        lines.append(f"zero {_flag_sexpr(Flag(Z.k, 0, Z.cols))}")
    if problem.objective:
        lines.append(f"objective {sexpr(problem.objective)}")
    for e, rhs in problem.eqs:
        lines.append(f"eq {_rat(rhs)} {sexpr(e)}")
    for e, rhs in problem.ineqs:
        lines.append(f"le {_rat(rhs)} {sexpr(e)}")
    for m, g in getattr(problem, "localizing", ()):
        lines.append(f"psd {m} {sexpr(g)}")
    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")
