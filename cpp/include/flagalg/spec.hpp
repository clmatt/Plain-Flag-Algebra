// Problem specifications (written by python/flagcheck/spec.py).
//
// Line-based; '#' starts a comment. Directives:
//   colors <c>
//   n <n>
//   empty_type <0|1>
//   zero <expr>                    (a single flag literal; roots ignored)
//   objective <expr>               (prove objective >= d)
//   eq <rational> <expr>           (expr = rational)
//   le <rational> <expr>           (expr <= rational)
// Expressions (S-expressions):
//   (flag k s c0 c1 ...)           flag literal, colours in lexicographic pair order
//   (add e1 e2 ...)  (scale q e)  (mul e1 e2 ...)  (avg e)
// Objective and constraints must evaluate to graphs (s = 0); they are lifted to n vertices.
#pragma once

#include <gmpxx.h>

#include <memory>
#include <string>
#include <vector>

#include "flagalg/theory.hpp"

namespace flagalg {

struct Expr {
    std::string op;                          // "flag", "add", "scale", "mul", "avg"
    Flag flag;                               // for "flag"
    mpq_class scalar;                        // for "scale"
    std::vector<std::shared_ptr<Expr>> args;
};

struct Constraint {
    enum Kind { Eq, Le } kind;
    mpq_class rhs;
    std::shared_ptr<Expr> expr;
};

struct Spec {
    int colors = 2;
    int n = 0;
    bool emptyType = false;
    std::vector<Flag> zeros;
    std::shared_ptr<Expr> objective;  // may be null (objective 0)
    std::vector<Constraint> constraints;
};

Spec parseSpec(const std::string& text);
mpq_class parseRational(const std::string& s);
Combo evaluate(const Theory& T, const Expr& e);

}  // namespace flagalg
