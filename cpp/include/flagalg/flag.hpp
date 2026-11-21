// Complete edge-coloured graphs with labelled roots ("flags").
//
// A Flag has k vertices 0..k-1, the first s of which are the roots (the labelled type),
// and a colour for every pair i < j, stored in lexicographic pair order. Plain graphs
// have s = 0. Two flags are isomorphic iff some permutation fixing every root maps one
// onto the other; `canonical` returns a representative that is equal for isomorphic
// flags and different otherwise.
#pragma once

#include <array>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace flagalg {

constexpr int kMaxVertices = 10;
constexpr int kMaxPairs = kMaxVertices * (kMaxVertices - 1) / 2;

// Index of pair (i, j), i != j, in lexicographic order on k vertices.
int pairIndex(int k, int i, int j);
inline int numPairs(int k) { return k * (k - 1) / 2; }

struct Flag {
    std::uint8_t k = 0;
    std::uint8_t s = 0;
    std::array<std::uint8_t, kMaxPairs> cols{};  // unused tail entries are always 0

    int color(int i, int j) const { return cols[pairIndex(k, i, j)]; }
    void setColor(int i, int j, int c) { cols[pairIndex(k, i, j)] = static_cast<std::uint8_t>(c); }

    bool operator==(const Flag& o) const { return k == o.k && s == o.s && cols == o.cols; }
    bool operator!=(const Flag& o) const { return !(*this == o); }
    bool operator<(const Flag& o) const;

    // The labelled type: the subgraph induced on the roots (as a plain graph, s = 0).
    Flag type() const;
    // The underlying graph (roots forgotten, not canonicalised).
    Flag forgetRoots() const { Flag f = *this; f.s = 0; return f; }

    std::string str() const;  // "k s c0 c1 ..." (the exchange format)
};

struct FlagHash {
    std::size_t operator()(const Flag& f) const;
};

Flag makeFlag(int k, int s, const std::vector<int>& cols);
// Every pair gets `defaultColor`, except the listed (i, j, colour) triples.
Flag fromEdges(int k, const std::vector<std::array<int, 3>>& edges, int defaultColor = 0, int s = 0);

// The flag induced on `vertices` (in that order); the first `s` of them become the roots.
Flag restrict(const Flag& f, const int* vertices, int count, int s);
inline Flag restrict(const Flag& f, const std::vector<int>& v, int s) {
    return restrict(f, v.data(), static_cast<int>(v.size()), s);
}

// Add a vertex k with colour colors[i] on the pair (i, k).
Flag extend(const Flag& f, const std::uint8_t* colors);

// Canonical form. Uses brute force over permutations of the non-roots when there are
// at most kBruteForceNonRoots of them, and nauty otherwise. The choice depends only on
// (k, s), so representatives are consistent.
constexpr int kBruteForceNonRoots = 4;
Flag canonical(const Flag& f, int numColors);
Flag canonicalBrute(const Flag& f);
Flag canonicalNauty(const Flag& f, int numColors, double* groupSize = nullptr);
// Number of automorphisms (permutations fixing every root and every colour); nauty.
std::uint64_t automorphismCount(const Flag& f, int numColors);
std::uint64_t automorphismCountBrute(const Flag& f);   // for tests

}  // namespace flagalg
