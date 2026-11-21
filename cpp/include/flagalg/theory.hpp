// A "theory": the number of edge colours plus a list of forbidden induced subgraphs
// ("zeros"). Generates admissible graphs and flags up to isomorphism, and implements the
// flag algebra operations (chain rule, product, averaging) with exact rationals.
#pragma once

#include <gmpxx.h>

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "flagalg/flag.hpp"

namespace flagalg {

using Combo = std::map<Flag, mpq_class>;  // canonical flags -> coefficient

class Theory {
public:
    Theory(int numColors, std::vector<Flag> zeros);

    int numColors() const { return numColors_; }
    const std::vector<Flag>& zeros() const { return zeros_; }

    Flag canon(const Flag& f) const { return canonical(f, numColors_); }

    // Does g contain a forbidden configuration? Only subsets containing the last vertex
    // are checked, which decides admissibility when g minus its last vertex is admissible.
    bool admissibleWithLastVertex(const Flag& g) const;
    bool admissible(const Flag& g) const;

    // All admissible graphs on n vertices (canonical, sorted). Cached.
    const std::vector<Flag>& graphs(int n) const;
    // All admissible sigma-flags on m vertices (sigma: labelled type, s = 0). Cached.
    const std::vector<Flag>& flags(const Flag& sigma, int m) const;

    // p(F, H): F and H flags of the same labelled type.
    mpq_class density(const Flag& F, const Flag& H) const;
    // c(F1, F2, H) of the product formula.
    mpq_class productCoefficient(const Flag& F1, const Flag& F2, const Flag& H) const;
    // q_sigma(F) of the averaging formula.
    mpq_class averageCoefficient(const Flag& F) const;

    Combo lift(const Combo& c, int m) const;
    Combo multiply(const Combo& a, const Combo& b) const;
    Combo average(const Combo& c) const;

private:
    int numColors_;
    std::vector<Flag> zeros_;
    std::map<int, std::unordered_set<Flag, FlagHash>> zerosBySize_;
    mutable std::map<int, std::vector<Flag>> graphCache_;
    mutable std::map<std::pair<Flag, int>, std::vector<Flag>> flagCache_;
    mutable std::map<std::pair<Flag, int>, Combo> liftCache_;  // single flag -> lifted combination

    // For unrooted lifts: (k, m) -> canonical k-vertex graph F -> [(index into graphs(m),
    // number of k-subsets of that graph inducing F)]. Built once from the m-vertex side.
    using SubIndex = std::unordered_map<Flag, std::vector<std::pair<std::uint32_t, std::uint32_t>>, FlagHash>;
    const SubIndex& subgraphIndex(int k, int m) const;
    mutable std::map<std::pair<int, int>, SubIndex> subIndexCache_;
};

void addTerm(Combo& c, const Flag& f, const mpq_class& v);

}  // namespace flagalg
