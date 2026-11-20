#include "flagalg/flag.hpp"

#include <algorithm>
#include <numeric>
#include <sstream>
#include <stdexcept>

extern "C" {
#include <nauty.h>
}

namespace flagalg {

namespace {

struct PairTable {
    // idx[k][i][j]
    std::array<std::array<std::array<std::int8_t, kMaxVertices>, kMaxVertices>, kMaxVertices + 1> idx{};
    PairTable() {
        for (int k = 0; k <= kMaxVertices; ++k) {
            int p = 0;
            for (int i = 0; i < k; ++i)
                for (int j = i + 1; j < k; ++j) {
                    idx[k][i][j] = static_cast<std::int8_t>(p);
                    idx[k][j][i] = static_cast<std::int8_t>(p);
                    ++p;
                }
        }
    }
};

const PairTable& pairTable() {
    static const PairTable t;
    return t;
}

}  // namespace

int pairIndex(int k, int i, int j) { return pairTable().idx[k][i][j]; }

bool Flag::operator<(const Flag& o) const {
    if (k != o.k) return k < o.k;
    if (s != o.s) return s < o.s;
    return cols < o.cols;
}

Flag Flag::type() const {
    std::vector<int> roots(s);
    std::iota(roots.begin(), roots.end(), 0);
    return restrict(*this, roots, 0);
}

std::string Flag::str() const {
    std::ostringstream os;
    os << int(k) << ' ' << int(s);
    for (int p = 0; p < numPairs(k); ++p) os << ' ' << int(cols[p]);
    return os.str();
}

std::size_t FlagHash::operator()(const Flag& f) const {
    std::uint64_t h = 1469598103934665603ULL ^ (std::uint64_t(f.k) << 8) ^ f.s;
    const int np = numPairs(f.k);
    for (int p = 0; p < np; ++p) {
        h ^= f.cols[p];
        h *= 1099511628211ULL;
    }
    return static_cast<std::size_t>(h);
}

Flag makeFlag(int k, int s, const std::vector<int>& cols) {
    if (k < 0 || k > kMaxVertices || s < 0 || s > k) throw std::invalid_argument("makeFlag: bad k/s");
    if (static_cast<int>(cols.size()) != numPairs(k)) throw std::invalid_argument("makeFlag: wrong number of colours");
    Flag f;
    f.k = static_cast<std::uint8_t>(k);
    f.s = static_cast<std::uint8_t>(s);
    for (int p = 0; p < numPairs(k); ++p) {
        if (cols[p] < 0 || cols[p] > 255) throw std::invalid_argument("makeFlag: bad colour");
        f.cols[p] = static_cast<std::uint8_t>(cols[p]);
    }
    return f;
}

Flag fromEdges(int k, const std::vector<std::array<int, 3>>& edges, int defaultColor, int s) {
    std::vector<int> cols(numPairs(k), defaultColor);
    for (const auto& e : edges) {
        if (e[0] == e[1] || e[0] < 0 || e[1] < 0 || e[0] >= k || e[1] >= k)
            throw std::invalid_argument("fromEdges: bad edge");
        cols[pairIndex(k, e[0], e[1])] = e[2];
    }
    return makeFlag(k, s, cols);
}

Flag restrict(const Flag& f, const int* v, int count, int s) {
    Flag g;
    g.k = static_cast<std::uint8_t>(count);
    g.s = static_cast<std::uint8_t>(s);
    int p = 0;
    for (int a = 0; a < count; ++a)
        for (int b = a + 1; b < count; ++b) g.cols[p++] = f.cols[pairIndex(f.k, v[a], v[b])];
    return g;
}

Flag extend(const Flag& f, const std::uint8_t* colors) {
    Flag g;
    const int k = f.k;
    if (k + 1 > kMaxVertices) throw std::length_error("extend: more than kMaxVertices vertices");
    g.k = static_cast<std::uint8_t>(k + 1);
    g.s = f.s;
    for (int i = 0; i < k; ++i)
        for (int j = i + 1; j < k; ++j) g.setColor(i, j, f.color(i, j));
    for (int i = 0; i < k; ++i) g.setColor(i, k, colors[i]);
    return g;
}

Flag canonicalBrute(const Flag& f) {
    const int k = f.k, s = f.s;
    std::array<int, kMaxVertices> perm{};
    std::iota(perm.begin(), perm.begin() + k, 0);
    Flag best = f;
    bool first = true;
    do {
        Flag g = restrict(f, perm.data(), k, s);
        if (first || g.cols < best.cols) {
            best = g;
            first = false;
        }
    } while (std::next_permutation(perm.begin() + s, perm.begin() + k));
    return best;
}

Flag canonicalNauty(const Flag& f, int numColors) {
    const int k = f.k, s = f.s;
    if (k <= 1) return f;
    int layers = 1;
    while ((1 << layers) < numColors) ++layers;
    const int nv = k * layers;
    const int m = SETWORDSNEEDED(nv);
    std::vector<setword> g(static_cast<std::size_t>(m) * nv, 0), cg(static_cast<std::size_t>(m) * nv, 0);
    std::vector<int> lab(nv), ptn(nv), orbits(nv);
    auto vid = [k](int v, int layer) { return layer * k + v; };
    for (int L = 0; L < layers; ++L) {
        for (int i = 0; i < k; ++i)
            for (int j = i + 1; j < k; ++j)
                if ((f.color(i, j) >> L) & 1) {
                    ADDONEEDGE(g.data(), vid(i, L), vid(j, L), m);
                }
        if (L + 1 < layers)
            for (int v = 0; v < k; ++v) ADDONEEDGE(g.data(), vid(v, L), vid(v, L + 1), m);
    }
    // Cells, layer by layer: each root alone (in order), then all non-roots.
    int pos = 0;
    for (int L = 0; L < layers; ++L) {
        for (int r = 0; r < s; ++r) {
            lab[pos] = vid(r, L);
            ptn[pos] = 0;
            ++pos;
        }
        for (int v = s; v < k; ++v) {
            lab[pos] = vid(v, L);
            ptn[pos] = 1;
            ++pos;
        }
        if (k > s) ptn[pos - 1] = 0;
    }
    DEFAULTOPTIONS_GRAPH(options);
    options.getcanon = TRUE;
    options.defaultptn = FALSE;
    statsblk stats;
    densenauty(g.data(), lab.data(), ptn.data(), orbits.data(), &options, &stats, m, nv, cg.data());
    // Positions 0..k-1 hold the layer-0 vertices (roots first, in order).
    std::array<int, kMaxVertices> order{};
    for (int p = 0; p < k; ++p) {
        if (lab[p] >= k) throw std::logic_error("canonicalNauty: layer-0 cell order violated");
        order[p] = lab[p];
    }
    for (int r = 0; r < s; ++r)
        if (order[r] != r) throw std::logic_error("canonicalNauty: roots moved");
    return restrict(f, order.data(), k, s);
}

Flag canonical(const Flag& f, int numColors) {
    if (f.k - f.s <= kBruteForceNonRoots) return canonicalBrute(f);
    return canonicalNauty(f, numColors);
}

}  // namespace flagalg
