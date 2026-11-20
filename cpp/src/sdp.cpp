#include "flagalg/sdp.hpp"

#include <algorithm>
#include <map>
#include <numeric>
#include <stdexcept>

#include "flagalg/util.hpp"

namespace flagalg {

std::vector<Block> makeBlocks(const Theory& T, int n, bool includeEmptyType) {
    std::vector<Block> blocks;
    for (int s = n % 2; s < n; s += 2) {
        if (s == 0 && !includeEmptyType) continue;
        const int m = (n + s) / 2;
        for (const Flag& sigma : T.graphs(s)) {
            Block b;
            b.sigma = sigma;
            b.m = m;
            b.flags = T.flags(sigma, m);
            if (b.flags.empty()) continue;
            for (int i = 0; i < static_cast<int>(b.flags.size()); ++i) b.index.emplace(b.flags[i], i);
            b.denominator = fallingFactorial(n, s) * binom(n - s, (n - s) / 2);
            blocks.push_back(std::move(b));
        }
    }
    return blocks;
}

std::vector<YEntry> computeY(const Theory& T, const std::vector<Block>& blocks, const Flag& G) {
    const int n = G.k;
    // labelled type (exact colours) -> block, per type size
    std::map<int, std::unordered_map<Flag, int, FlagHash>> typeIndex;
    for (int b = 0; b < static_cast<int>(blocks.size()); ++b) typeIndex[blocks[b].sigma.k].emplace(blocks[b].sigma, b);

    std::map<std::tuple<std::uint32_t, std::uint32_t, std::uint32_t>, std::uint64_t> acc;
    for (const auto& [s, types] : typeIndex) {
        const int r = (n - s) / 2;
        std::vector<int> set(s);
        std::iota(set.begin(), set.end(), 0);
        std::array<int, kMaxVertices> va{}, vb{};
        do {  // root sets
            std::vector<int> theta = set;
            do {  // orderings of the roots
                Flag t = restrict(G, theta, 0);
                auto it = types.find(t);
                if (it == types.end()) continue;
                const Block& B = blocks[it->second];
                std::vector<int> rest;
                std::vector<bool> isRoot(n, false);
                for (int x : theta) isRoot[x] = true;
                for (int x = 0; x < n; ++x)
                    if (!isRoot[x]) rest.push_back(x);
                for (int i = 0; i < s; ++i) va[i] = vb[i] = theta[i];
                std::vector<int> X(r);
                std::iota(X.begin(), X.end(), 0);
                do {
                    std::vector<bool> inX(rest.size(), false);
                    for (int i = 0; i < r; ++i) {
                        inX[X[i]] = true;
                        va[s + i] = rest[X[i]];
                    }
                    int q = s;
                    for (std::size_t i = 0; i < rest.size(); ++i)
                        if (!inX[i]) vb[q++] = rest[i];
                    const Flag fa = T.canon(restrict(G, va.data(), s + r, s));
                    const Flag fb = T.canon(restrict(G, vb.data(), s + r, s));
                    auto ia = B.index.find(fa);
                    auto ib = B.index.find(fb);
                    if (ia == B.index.end() || ib == B.index.end())
                        throw std::logic_error("computeY: subflag missing from block (inadmissible graph?)");
                    // Count the ordered pair (i, j) only when i <= j: Y[i][j] for i <= j
                    // is then exactly the number of ordered occurrences.
                    if (ia->second <= ib->second) ++acc[{static_cast<std::uint32_t>(it->second),
                                                         static_cast<std::uint32_t>(ia->second),
                                                         static_cast<std::uint32_t>(ib->second)}];
                } while (nextSubset(X, static_cast<int>(rest.size())));
            } while (std::next_permutation(theta.begin(), theta.end()));
        } while (nextSubset(set, n));
    }
    std::vector<YEntry> out;
    out.reserve(acc.size());
    for (const auto& [key, count] : acc) out.push_back({std::get<0>(key), std::get<1>(key), std::get<2>(key), count});
    return out;
}

}  // namespace flagalg
