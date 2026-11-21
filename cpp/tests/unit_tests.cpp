// C++ unit tests for the core. The heavier cross-checks against the independent Python
// reference live in python/tests/test_cpp_backend.py.
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>
#include <string>
#include <tuple>

#include "flagalg/sdp.hpp"
#include "flagalg/theory.hpp"
#include "flagalg/util.hpp"

using namespace flagalg;

static int failures = 0;
#define CHECK(cond)                                                                   \
    do {                                                                              \
        if (!(cond)) {                                                                \
            std::cerr << __FILE__ << ":" << __LINE__ << ": CHECK failed: " #cond "\n"; \
            ++failures;                                                               \
        }                                                                             \
    } while (0)

static Flag randomFlag(std::mt19937& rng, int k, int s, int c) {
    std::vector<int> cols(numPairs(k));
    for (int& x : cols) x = static_cast<int>(rng() % c);
    return makeFlag(k, s, cols);
}

// Brute force and nauty canonical forms must induce the same isomorphism classes.
static void testCanonicalAgreement() {
    std::mt19937 rng(12345);
    for (int trial = 0; trial < 3000; ++trial) {
        int k = 2 + static_cast<int>(rng() % 6);  // 2..7
        int s = static_cast<int>(rng() % (k + 1));
        int c = 2 + static_cast<int>(rng() % 3);
        Flag f = randomFlag(rng, k, s, c);
        // random relabelling of the non-roots
        std::vector<int> perm(k);
        for (int i = 0; i < k; ++i) perm[i] = i;
        std::shuffle(perm.begin() + s, perm.end(), rng);
        Flag g = restrict(f, perm, s);
        CHECK(canonicalBrute(f) == canonicalBrute(g));
        CHECK(canonicalNauty(f, c) == canonicalNauty(g, c));
        // nauty and brute force classes agree on a second random flag
        Flag h = randomFlag(rng, k, s, c);
        CHECK((canonicalBrute(f) == canonicalBrute(h)) == (canonicalNauty(f, c) == canonicalNauty(h, c)));
        // roots are never moved
        CHECK(canonicalNauty(f, c).type() == f.type());
    }
    // a non-root permutation that swaps a root is NOT an isomorphism of flags
    Flag a = fromEdges(3, {{0, 1, 1}}, 0, 1);  // root 0 adjacent to vertex 1
    Flag b = fromEdges(3, {{1, 2, 1}}, 0, 1);  // root 0 isolated
    CHECK(canonical(a, 2) != canonical(b, 2));
    CHECK(canonicalNauty(a, 2) != canonicalNauty(b, 2));
}

static void testCounts() {
    Theory two(2, {});
    const int graphs[] = {1, 1, 2, 4, 11, 34, 156, 1044};
    for (int n = 0; n <= 7; ++n) CHECK(static_cast<int>(two.graphs(n).size()) == graphs[n]);
    Flag K3 = fromEdges(3, {}, 1);
    Theory triFree(2, {K3});
    const int tf[] = {1, 1, 2, 3, 7, 14, 38, 107};
    for (int n = 0; n <= 7; ++n) CHECK(static_cast<int>(triFree.graphs(n).size()) == tf[n]);
    Theory three(3, {});
    CHECK(three.graphs(4).size() == 66);
}

static void testNextSubsetAndBinom() {
    std::vector<int> sub = {0, 1};
    int count = 1;
    while (nextSubset(sub, 5)) ++count;
    CHECK(count == 10);
    std::vector<int> empty;
    CHECK(!nextSubset(empty, 3));
    CHECK(binom(8, 4) == 70);
    CHECK(fallingFactorial(8, 3) == 336);
}

static void testYRowSums() {
    // For every graph and block, sum_{i,j} Y[i][j] (full matrix) = P[G[theta] = sigma].
    Theory T(3, {});
    const int n = 5;
    auto blocks = makeBlocks(T, n, false);
    for (const Flag& G : T.graphs(n)) {
        auto Y = computeY(T, blocks, G);
        std::vector<std::uint64_t> total(blocks.size(), 0);
        for (const auto& e : Y) total[e.block] += (e.i == e.j ? 1 : 2) * e.count;
        for (std::size_t b = 0; b < blocks.size(); ++b) {
            // count ordered theta with G[theta] == sigma, times C(n-s, r)
            const int s = blocks[b].sigma.k;
            std::uint64_t matches = 0;
            std::vector<int> set(s);
            for (int i = 0; i < s; ++i) set[i] = i;
            do {
                std::vector<int> th = set;
                do {
                    if (restrict(G, th, 0) == blocks[b].sigma) ++matches;
                } while (std::next_permutation(th.begin(), th.end()));
            } while (nextSubset(set, n));
            CHECK(total[b] == matches * binom(n - s, (n - s) / 2));
        }
    }
}

// The indexed unrooted lift must equal the chain rule computed from density() directly,
// including for inadmissible F (empty lift) and for F = H (k = m).
static void testIndexedLift() {
    std::mt19937 rng(7);
    Flag redTriangle = makeFlag(3, 0, {1, 1, 1});
    for (auto [c, zeros, mmax] : {std::make_tuple(2, std::vector<Flag>{}, 6),
                                  std::make_tuple(3, std::vector<Flag>{redTriangle}, 5)}) {
        Theory T(c, zeros), ref(c, zeros);
        for (int m = 1; m <= mmax; ++m)
            for (int k = 1; k <= m; ++k) {
                std::vector<Flag> Fs = T.graphs(k);
                for (int r = 0; r < 5; ++r) Fs.push_back(T.canon(randomFlag(rng, k, 0, c)));
                for (const Flag& F : Fs) {
                    Combo fast = T.lift(Combo{{F, 1}}, m), slow;
                    for (const Flag& H : ref.graphs(m)) {
                        mpq_class d = ref.density(F, H);
                        if (d != 0) slow.emplace(H, d);
                    }
                    CHECK(fast == slow);
                }
            }
    }
}

// nauty's group size must equal the brute-force automorphism count, rooted or not, and
// for flags with many symmetries (monochromatic, two colour classes) as well as random.
static void testAutomorphismCount() {
    std::mt19937 rng(11);
    for (int k = 1; k <= 7; ++k)
        for (int s = 0; s <= std::min(k, 2); ++s)
            for (int c : {2, 3})
                for (int t = 0; t < 60; ++t) {
                    Flag f = randomFlag(rng, k, s, t < 10 ? 1 : (t < 20 ? 2 : c));
                    CHECK(automorphismCount(f, c) == automorphismCountBrute(f));
                }
    CHECK(automorphismCount(makeFlag(5, 0, std::vector<int>(10, 1)), 3) == 120);
}

int main() {
    testAutomorphismCount();
    testNextSubsetAndBinom();
    testCanonicalAgreement();
    testCounts();
    testYRowSums();
    testIndexedLift();
    if (failures) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "all unit tests passed\n";
    return 0;
}
