// Plain flag algebra SDP data (thesis Algorithm 2) with exact integer counts.
//
// For a graph G on n vertices and a labelled type sigma on s vertices (n - s even,
// r = (n - s)/2), Y_{G,sigma}[i][j] = count / denominator, where count is the number
// of pairs (theta, X) with theta: [s] -> V(G) injective, G[theta] exactly sigma, X an
// r-subset of the remaining vertices, (G[theta + X], theta) ~ flags[i] and
// (G[theta + rest], theta) ~ flags[j]; and denominator = n!/(n-s)! * C(n-s, r).
// Only i <= j is stored; Y is symmetric.
#pragma once

#include <cstdint>
#include <ostream>
#include <unordered_map>
#include <vector>

#include "flagalg/theory.hpp"

namespace flagalg {

struct Block {
    Flag sigma;                 // labelled type (s = 0)
    int m = 0;                  // flag size
    std::vector<Flag> flags;    // canonical sigma-flags on m vertices
    std::unordered_map<Flag, int, FlagHash> index;
    std::uint64_t denominator = 0;
};

struct YEntry {
    std::uint32_t block, i, j;
    std::uint64_t count;
};

// One representative labelled type per isomorphism class, for every s < n with
// s = n (mod 2); s = 0 only if includeEmptyType.
std::vector<Block> makeBlocks(const Theory& T, int n, bool includeEmptyType);

// Y entries for one graph (sorted by block, i, j).
std::vector<YEntry> computeY(const Theory& T, const std::vector<Block>& blocks, const Flag& G);

}  // namespace flagalg
