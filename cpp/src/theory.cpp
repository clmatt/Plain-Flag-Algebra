#include "flagalg/theory.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>

#include "flagalg/util.hpp"

namespace flagalg {

void addTerm(Combo& c, const Flag& f, const mpq_class& v) {
    if (v == 0) return;
    auto it = c.find(f);
    if (it == c.end()) {
        c.emplace(f, v);
        return;
    }
    it->second += v;
    if (it->second == 0) c.erase(it);
}

Theory::Theory(int numColors, std::vector<Flag> zeros) : numColors_(numColors), zeros_(std::move(zeros)) {
    if (numColors < 2 || numColors > 16) throw std::invalid_argument("Theory: numColors must be in [2, 16]");
    for (const Flag& z : zeros_) {
        for (int p = 0; p < numPairs(z.k); ++p)
            if (z.cols[p] >= numColors) throw std::invalid_argument("Theory: zero uses a colour out of range");
        zerosBySize_[z.k].insert(canon(z.forgetRoots()));
    }
}

bool Theory::admissibleWithLastVertex(const Flag& g) const {
    const int last = g.k - 1;
    std::array<int, kMaxVertices> v{};
    for (const auto& [size, targets] : zerosBySize_) {
        if (size > g.k) continue;
        if (size == 0) return false;
        // all (size-1)-subsets of {0..last-1}, plus `last`
        std::vector<int> sub(size - 1);
        std::iota(sub.begin(), sub.end(), 0);
        do {
            for (int i = 0; i < size - 1; ++i) v[i] = sub[i];
            v[size - 1] = last;
            if (targets.count(canon(restrict(g, v.data(), size, 0)))) return false;
        } while (nextSubset(sub, last));
    }
    return true;
}

bool Theory::admissible(const Flag& g) const {
    for (const auto& [size, targets] : zerosBySize_) {
        if (size > g.k) continue;
        std::vector<int> sub(size);
        std::iota(sub.begin(), sub.end(), 0);
        do {
            if (targets.count(canon(restrict(g, sub, 0)))) return false;
        } while (nextSubset(sub, g.k));
    }
    return true;
}

namespace {

// Extend every parent by one vertex in all colourings; keep admissible ones, deduplicated.
std::vector<Flag> extendAll(const Theory& T, const std::vector<Flag>& parents) {
    const int c = T.numColors();
    std::vector<std::unordered_set<Flag, FlagHash>> local;
#pragma omp parallel
    {
#pragma omp single
        local.resize(omp_threads());
        auto& mine = local[omp_thread()];
#pragma omp for schedule(dynamic, 16)
        for (std::size_t pi = 0; pi < parents.size(); ++pi) {
            const Flag& P = parents[pi];
            const int k = P.k;
            std::array<std::uint8_t, kMaxVertices> col{};
            // iterate over all colour vectors of length k (base c counter)
            while (true) {
                Flag g = extend(P, col.data());
                if (T.admissibleWithLastVertex(g)) mine.insert(T.canon(g));
                int i = 0;
                while (i < k && ++col[i] == c) col[i++] = 0;
                if (i == k) break;
            }
        }
    }
    std::unordered_set<Flag, FlagHash> all;
    for (auto& s : local) all.insert(s.begin(), s.end());
    std::vector<Flag> out(all.begin(), all.end());
    std::sort(out.begin(), out.end());
    return out;
}

}  // namespace

const std::vector<Flag>& Theory::graphs(int n) const {
    if (n < 0 || n > kMaxVertices) throw std::length_error("graphs: n out of range");
    auto it = graphCache_.find(n);
    if (it != graphCache_.end()) return it->second;
    std::vector<Flag> out;
    if (n == 0) {
        out.push_back(Flag{});
    } else {
        out = extendAll(*this, graphs(n - 1));
    }
    return graphCache_[n] = std::move(out);
}

const std::vector<Flag>& Theory::flags(const Flag& sigmaIn, int m) const {
    Flag sigma = sigmaIn.forgetRoots();
    if (m < sigma.k) throw std::invalid_argument("flags: m smaller than the type");
    auto key = std::make_pair(sigma, m);
    auto it = flagCache_.find(key);
    if (it != flagCache_.end()) return it->second;
    std::vector<Flag> out;
    if (m == sigma.k) {
        Flag id = sigma;
        id.s = sigma.k;
        if (admissible(sigma)) out.push_back(id);
    } else {
        out = extendAll(*this, flags(sigma, m - 1));
    }
    return flagCache_[key] = std::move(out);
}

mpq_class Theory::density(const Flag& F, const Flag& H) const {
    if (F.s != H.s || F.type() != H.type()) throw std::invalid_argument("density: different labelled types");
    if (F.k > H.k) throw std::invalid_argument("density: F larger than H");
    const int s = F.s, r = F.k - s;
    const Flag target = canon(F);
    std::array<int, kMaxVertices> v{};
    for (int i = 0; i < s; ++i) v[i] = i;
    std::vector<int> X(r);
    std::iota(X.begin(), X.end(), 0);
    long hits = 0, total = 0;
    do {
        for (int i = 0; i < r; ++i) v[s + i] = s + X[i];
        ++total;
        if (canon(restrict(H, v.data(), F.k, s)) == target) ++hits;
    } while (nextSubset(X, H.k - s));
    mpq_class q(hits, total);
    q.canonicalize();
    return q;
}

mpq_class Theory::productCoefficient(const Flag& F1, const Flag& F2, const Flag& H) const {
    const int s = F1.s, a = F1.k - s, b = F2.k - s;
    if (H.k != s + a + b) throw std::invalid_argument("productCoefficient: wrong size");
    const Flag t1 = canon(F1), t2 = canon(F2);
    std::array<int, kMaxVertices> v1{}, v2{};
    for (int i = 0; i < s; ++i) v1[i] = v2[i] = i;
    std::vector<int> X(a);
    std::iota(X.begin(), X.end(), 0);
    long hits = 0, total = 0;
    do {
        std::vector<bool> inX(a + b, false);
        for (int i = 0; i < a; ++i) {
            inX[X[i]] = true;
            v1[s + i] = s + X[i];
        }
        int q = 0;
        for (int i = 0; i < a + b; ++i)
            if (!inX[i]) v2[s + q++] = s + i;
        ++total;
        if (canon(restrict(H, v1.data(), s + a, s)) == t1 && canon(restrict(H, v2.data(), s + b, s)) == t2) ++hits;
    } while (nextSubset(X, a + b));
    mpq_class q(hits, total);
    q.canonicalize();
    return q;
}

mpq_class Theory::averageCoefficient(const Flag& F) const {
    const int k = F.k, s = F.s;
    const Flag target = canon(F);
    std::vector<int> theta(s);
    std::iota(theta.begin(), theta.end(), 0);
    long hits = 0, total = 0;
    std::array<int, kMaxVertices> v{};
    do {  // subsets
        std::vector<int> perm = theta;
        do {  // orderings
            std::vector<bool> used(k, false);
            for (int i = 0; i < s; ++i) {
                v[i] = perm[i];
                used[perm[i]] = true;
            }
            int q = s;
            for (int x = 0; x < k; ++x)
                if (!used[x]) v[q++] = x;
            ++total;
            if (canon(restrict(F, v.data(), k, s)) == target) ++hits;
        } while (std::next_permutation(perm.begin(), perm.end()));
    } while (nextSubset(theta, k));
    mpq_class q(hits, total);
    q.canonicalize();
    return q;
}

const Theory::SubIndex& Theory::subgraphIndex(int k, int m) const {
    auto key = std::make_pair(k, m);
    auto it = subIndexCache_.find(key);
    if (it != subIndexCache_.end()) return it->second;
    const auto& Hs = graphs(m);
    // per graph H: the canonical k-vertex induced subgraphs with multiplicities
    std::vector<std::vector<std::pair<Flag, std::uint32_t>>> perH(Hs.size());
#pragma omp parallel for schedule(dynamic, 64)
    for (std::size_t h = 0; h < Hs.size(); ++h) {
        std::unordered_map<Flag, std::uint32_t, FlagHash> counts;
        std::vector<int> sub(k);
        std::iota(sub.begin(), sub.end(), 0);
        do {
            ++counts[canon(restrict(Hs[h], sub, 0))];
        } while (nextSubset(sub, m));
        perH[h].assign(counts.begin(), counts.end());
    }
    SubIndex index;
    for (std::size_t h = 0; h < Hs.size(); ++h)  // in graph order, so each list is sorted by h
        for (const auto& [F, n] : perH[h]) index[F].emplace_back(static_cast<std::uint32_t>(h), n);
    return subIndexCache_[key] = std::move(index);
}

Combo Theory::lift(const Combo& c, int m) const {
    Combo out;
    for (const auto& [F, coeff] : c) {
        if (F.k > m) throw std::invalid_argument("lift: cannot lift to fewer vertices");
        auto key = std::make_pair(F, m);
        auto it = liftCache_.find(key);
        if (it == liftCache_.end()) {
            Combo single;
            if (F.s == 0 && F.k > 0) {
                // p(F, H) = #{k-subsets of H inducing F} / C(m, k), read off the index
                long subsets = 1;
                for (int i = 0; i < F.k; ++i) subsets = subsets * (m - i) / (i + 1);
                const auto& index = subgraphIndex(F.k, m);
                auto f = index.find(canon(F));
                if (f != index.end()) {
                    const auto& Hs = graphs(m);
                    for (const auto& [h, n] : f->second) {
                        mpq_class d(static_cast<long>(n), subsets);
                        d.canonicalize();
                        single.emplace(Hs[h], d);
                    }
                }
            } else {
                for (const Flag& H : flags(F.type(), m)) {
                    mpq_class d = density(F, H);
                    if (d != 0) single.emplace(H, d);
                }
            }
            it = liftCache_.emplace(key, std::move(single)).first;
        }
        for (const auto& [H, d] : it->second) addTerm(out, H, coeff * d);
    }
    return out;
}

Combo Theory::multiply(const Combo& a, const Combo& b) const {
    Combo out;
    for (const auto& [F1, v1] : a)
        for (const auto& [F2, v2] : b) {
            if (F1.s != F2.s || F1.type() != F2.type())
                throw std::invalid_argument("multiply: flags of different labelled types");
            const int m = F1.k + F2.k - F1.s;
            for (const Flag& H : flags(F1.type(), m)) {
                mpq_class c = productCoefficient(F1, F2, H);
                if (c != 0) addTerm(out, H, v1 * v2 * c);
            }
        }
    return out;
}

Combo Theory::average(const Combo& c) const {
    Combo out;
    for (const auto& [F, v] : c) {
        Flag G = canon(F.forgetRoots());
        addTerm(out, G, v * averageCoefficient(F));
    }
    return out;
}

}  // namespace flagalg
