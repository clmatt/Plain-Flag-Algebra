// flagsdp: build the exact data of a plain flag algebra SDP from a problem spec.
//
// Usage: flagsdp <spec-file> <output-file> [--threads N] [--sdpa <file> [--normalize]]
//   --sdpa writes the CSDP input directly (same layout as python/flagcheck/sdp.py write_sdpa),
//   so large problems never need to be loaded into Python.
//
// Output (text, read by python/flagcheck/cppbackend.py):
//   FLAGSDP 1
//   n <n>
//   colors <c>
//   graphs <G>                 followed by G lines "k s c0 c1 ..."
//   blocks <B>                 then per block:
//     block <b> <m> <num flags> <denominator> <plain|localizing> <sigma: k s c...>
//     <num flags> lines "k s c..."
//   objective <nnz>            followed by "<graph> <rational>" lines
//   constraint <eq|le> <rhs> <nnz>   followed by "<graph> <rational>" lines (repeated)
//   Y <num entries>            followed by "<graph> <block> <i> <j> <value>" lines (i <= j);
//                              entry = value / denominator (integer counts for plain blocks,
//                              exact rationals with denominator 1 for localizing blocks)
//   end
#include <chrono>
#include <cstdio>
#include <fstream>
#include <map>
#include <tuple>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "flagalg/sdp.hpp"
#include "flagalg/spec.hpp"
#include "flagalg/util.hpp"

using namespace flagalg;

namespace {

double seconds(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

std::vector<std::pair<int, mpq_class>> toVector(const Theory& T, const Combo& c, int n,
                                                const std::unordered_map<Flag, int, FlagHash>& gi,
                                                const char* what) {
    for (const auto& [f, v] : c)
        if (f.s != 0) throw std::invalid_argument(std::string(what) + " has rooted flags; average it first");
    Combo lifted = T.lift(c, n);
    std::vector<std::pair<int, mpq_class>> out;
    for (const auto& [f, v] : lifted) {
        auto it = gi.find(f);
        if (it == gi.end()) throw std::logic_error("lifted graph not in the graph list");
        out.emplace_back(it->second, v);
    }
    std::sort(out.begin(), out.end(), [](auto& a, auto& b) { return a.first < b.first; });
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    try {
        if (argc < 3) {
            std::cerr << "usage: flagsdp <spec> <out> [--threads N]\n";
            return 2;
        }
        std::string sdpaPath;
        bool normalize = false;
        for (int a = 3; a < argc; ++a) {
            const std::string arg = argv[a];
            if (arg == "--threads" && a + 1 < argc) {
#ifdef _OPENMP
                omp_set_num_threads(std::stoi(argv[a + 1]));
#endif
                ++a;
            } else if (arg == "--sdpa" && a + 1 < argc) {
                sdpaPath = argv[++a];
            } else if (arg == "--normalize") {
                normalize = true;
            } else {
                throw std::invalid_argument("unknown option " + arg);
            }
        }
        auto t0 = std::chrono::steady_clock::now();
        std::ifstream in(argv[1]);
        if (!in) throw std::runtime_error("cannot read spec");
        std::stringstream buf;
        buf << in.rdbuf();
        Spec spec = parseSpec(buf.str());
        Theory T(spec.colors, spec.zeros);
        const int n = spec.n;

        const auto& graphs = T.graphs(n);
        std::unordered_map<Flag, int, FlagHash> gi;
        for (int i = 0; i < static_cast<int>(graphs.size()); ++i) gi.emplace(graphs[i], i);
        std::cerr << "graphs on " << n << " vertices: " << graphs.size() << " (" << seconds(t0) << "s)\n";

        // Clone elimination (python/flagcheck/blowup.py): y_G = sum_Q M[G][Q] z_Q over the
        // green-free graphs Q, with D M integral for D = N^n n!:
        //   D M[G][Q] = n!/aut(G) (N)_r aut(Q_r) (n-r)! c(Q_r, Q),
        // c(Q_r, Q) = number of r-subsets of Q inducing the quotient Q_r of G.
        // Levels mode keeps Y sparse: y_G = lambda_G w(Q_r) with w(Q) the injective density of
        // the green-free Q on |Q| <= n vertices, so each G feeds one variable. The variable
        // vector is v = w / n (each level sums to 1, so v is a probability vector), with
        //   D n lambda_G = n * n!/aut(G) (N)_r aut(Q_r) n!/r!,
        // and the levels are tied by w(Q') = sum_Q d(Q', Q) w(Q), |Q| = |Q'| + 1 (written as
        // equality constraints after the spec's own).
        const bool reduced = spec.blowup > 0;
        const bool levels = reduced && spec.blowupLevels;
        std::vector<std::vector<std::pair<int, mpq_class>>> levelEqs;
        std::vector<Flag> reducedGraphs;
        std::vector<int> toReduced(graphs.size(), -1);
        std::vector<std::vector<std::pair<int, mpz_class>>> rowsM(reduced ? graphs.size() : 0), colsM;
        mpz_class D = 1;
        if (reduced) {
            const int N = spec.blowup;
            if (spec.colors != 3) throw std::invalid_argument("blowup needs 3 colours (0 = clone)");
            // levels: rows 0..g_n-1 are the green-free n-vertex graphs (graph order), then n-1, ..., 1
            std::vector<std::unordered_map<Flag, int, FlagHash>> levelIndex(n + 1);
            for (std::size_t g = 0; g < graphs.size(); ++g) {
                bool free = true;
                for (int i = 0; i < n && free; ++i)
                    for (int j = i + 1; j < n; ++j)
                        if (graphs[g].color(i, j) == 0) {
                            free = false;
                            break;
                        }
                if (free) {
                    toReduced[g] = static_cast<int>(reducedGraphs.size());
                    if (levels) levelIndex[n].emplace(graphs[g], toReduced[g]);
                    reducedGraphs.push_back(graphs[g]);
                }
            }
            for (int r = n - 1; levels && r >= 1; --r)
                for (const Flag& Q : T.graphs(r)) {
                    bool free = true;
                    for (int i = 0; i < r && free; ++i)
                        for (int j = i + 1; j < r; ++j)
                            if (Q.color(i, j) == 0) {
                                free = false;
                                break;
                            }
                    if (!free) continue;
                    levelIndex[r].emplace(Q, static_cast<int>(reducedGraphs.size()));
                    reducedGraphs.push_back(Q);
                }
            mpz_class nFact = 1;
            for (int i = 2; i <= n; ++i) nFact *= i;
            mpz_pow_ui(D.get_mpz_t(), mpz_class(N).get_mpz_t(), n);
            D *= nFact;
            for (std::size_t g = 0; g < graphs.size(); ++g) {
                const Flag& G = graphs[g];
                std::vector<int> classOf(n, -1), reps;
                for (int v = 0; v < n; ++v) {
                    if (classOf[v] >= 0) continue;
                    classOf[v] = static_cast<int>(reps.size());
                    for (int u = v + 1; u < n; ++u)
                        if (G.color(u, v) == 0) classOf[u] = classOf[v];
                    reps.push_back(v);
                }
                for (int u = 0; u < n; ++u)       // blow-up structure: green classes, uniform colours
                    for (int w = u + 1; w < n; ++w) {
                        const bool same = classOf[u] == classOf[w];
                        if (same != (G.color(u, w) == 0) ||
                            (!same && G.color(u, w) != G.color(reps[classOf[u]], reps[classOf[w]])))
                            throw std::logic_error("blowup: graph is not a blow-up");
                    }
                const int r = static_cast<int>(reps.size());
                const Flag Qr = T.canon(restrict(G, reps, 0));
                mpz_class A = nFact;
                const unsigned long autG = static_cast<unsigned long>(automorphismCount(G, 3));
                if (A % autG != 0) throw std::logic_error("blowup: aut(G) does not divide n factorial");
                A /= autG;
                for (int i = 0; i < r; ++i) A *= (N - i);
                A *= static_cast<unsigned long>(automorphismCount(Qr, 3));
                if (levels) {      // A = n!/aut(G) (N)_r aut(Q_r); times n * n!/r!
                    for (int i = r + 1; i <= n; ++i) A *= i;
                    A *= n;
                    auto it = levelIndex[r].find(Qr);
                    if (it == levelIndex[r].end()) throw std::logic_error("blowup: quotient not in its level");
                    rowsM[g].emplace_back(it->second, A);
                    continue;
                }
                for (int i = 2; i <= n - r; ++i) A *= i;
                const long subsets = static_cast<long>(binom(n, r));
                for (const auto& [H, d] : T.lift(Combo{{Qr, 1}}, n)) {
                    auto it = gi.find(H);
                    if (it == gi.end()) throw std::logic_error("blowup: lifted graph not in the list");
                    if (toReduced[it->second] < 0) continue;
                    mpq_class cnt = d * subsets;
                    cnt.canonicalize();
                    if (cnt.get_den() != 1) throw std::logic_error("blowup: non-integral subgraph count");
                    rowsM[g].emplace_back(toReduced[it->second], A * cnt.get_num());
                }
            }
            // w(Q') - sum_Q d(Q', Q) w(Q) = 0, d(Q', Q) = (vertex deletions of Q giving Q') / |Q|
            for (int r = n - 1; levels && r >= 1; --r) {
                std::map<int, std::map<int, long>> del;    // row(Q') -> row(Q) -> count
                for (const auto& [Q, q] : levelIndex[r + 1])
                    for (int x = 0; x <= r; ++x) {
                        std::vector<int> keep;
                        for (int v = 0; v <= r; ++v)
                            if (v != x) keep.push_back(v);
                        auto it = levelIndex[r].find(T.canon(restrict(Q, keep, 0)));
                        if (it == levelIndex[r].end()) throw std::logic_error("blowup: deletion not in its level");
                        ++del[it->second][q];
                    }
                std::vector<int> rows;
                for (const auto& [Qp, row] : levelIndex[r]) rows.push_back(row);
                std::sort(rows.begin(), rows.end());
                for (int row : rows) {
                    std::vector<std::pair<int, mpq_class>> vec{{row, mpq_class(1)}};
                    for (const auto& [q, c] : del[row]) vec.emplace_back(q, mpq_class(-c, r + 1));
                    std::sort(vec.begin(), vec.end(), [](auto& a, auto& b) { return a.first < b.first; });
                    for (auto& e : vec) e.second.canonicalize();
                    levelEqs.push_back(std::move(vec));
                }
            }
            colsM.resize(reducedGraphs.size());
            for (std::size_t g = 0; g < graphs.size(); ++g)
                for (const auto& [q, c] : rowsM[g]) colsM[q].emplace_back(static_cast<int>(g), c);
            std::cerr << "blowup N=" << N << (levels ? " (levels)" : "") << ": " << reducedGraphs.size()
                      << " rows, " << levelEqs.size() << " level equations (" << seconds(t0) << "s)\n";
        }
        // c . y = (c M) . z, with M = (D M) / D
        auto reduceVec = [&](const std::vector<std::pair<int, mpq_class>>& vec) {
            std::map<int, mpq_class> acc;
            for (const auto& [g, v] : vec)
                for (const auto& [q, c] : rowsM[g]) acc[q] += v * c;
            std::vector<std::pair<int, mpq_class>> res;
            for (auto& [q, v] : acc)
                if (v != 0) {
                    mpq_class x = v / D;
                    x.canonicalize();
                    res.emplace_back(q, x);
                }
            return res;
        };
        const auto& outGraphs = reduced ? reducedGraphs : graphs;

        std::vector<std::pair<int, mpq_class>> objective;
        if (spec.objective) objective = toVector(T, evaluate(T, *spec.objective), n, gi, "objective");
        if (reduced) objective = reduceVec(objective);

        auto blocks = makeBlocks(T, n, spec.emptyType);
        std::cerr << "blocks: " << blocks.size() << " (" << seconds(t0) << "s)\n";

        // Localising blocks: entry (i, j) = lift([[ g F_i F_j ]], n), exact rationals.
        struct LocBlock {
            Flag sigma;
            int m;
            std::vector<Flag> flags;
            std::map<std::tuple<int, int, int>, mpq_class> entries;  // (graph, i, j)
        };
        std::vector<LocBlock> locBlocks;
        for (const auto& L : spec.localizing) {
            Combo g = evaluate(T, *L.g);
            if (g.empty()) throw std::invalid_argument("psd: multiplier is zero");
            const Flag& f0 = g.begin()->first;
            int maxk = 0;
            for (const auto& [f, v] : g) {
                if (f.s != f0.s || f.type() != f0.type())
                    throw std::invalid_argument("psd: multiplier mixes labelled types");
                maxk = std::max<int>(maxk, f.k);
            }
            LocBlock B;
            B.sigma = f0.type();
            B.m = L.m;
            const int s = f0.s;
            if (L.m < s || maxk + 2 * (L.m - s) > n)
                throw std::invalid_argument("psd: [[g F F]] does not fit in n vertices");
            B.flags = T.flags(B.sigma, L.m);
            for (std::size_t i = 0; i < B.flags.size(); ++i) {
                Combo gF = T.multiply(g, Combo{{B.flags[i], 1}});
                for (std::size_t j = i; j < B.flags.size(); ++j) {
                    Combo e = T.lift(T.average(T.multiply(gF, Combo{{B.flags[j], 1}})), n);
                    for (const auto& [G, v] : e) {
                        auto it = gi.find(G);
                        if (it == gi.end()) throw std::logic_error("psd: lifted graph not in the graph list");
                        B.entries[{it->second, static_cast<int>(i), static_cast<int>(j)}] = v;
                    }
                }
            }
            if (reduced) {
                std::map<std::tuple<int, int, int>, mpq_class> red;
                for (const auto& [key, v] : B.entries)
                    for (const auto& [q, c] : rowsM[std::get<0>(key)])
                        red[{q, std::get<1>(key), std::get<2>(key)}] += v * c;
                B.entries.clear();
                for (auto& [key, v] : red)
                    if (v != 0) {
                        mpq_class x = v / D;
                        x.canonicalize();
                        B.entries[key] = x;
                    }
            }
            std::cerr << "localising block: " << B.flags.size() << " flags (" << seconds(t0) << "s)\n";
            locBlocks.push_back(std::move(B));
        }

        std::ofstream out(argv[2]);
        out << "FLAGSDP 1\nn " << n << "\ncolors " << spec.colors << "\n";
        out << "graphs " << outGraphs.size() << "\n";
        for (const auto& g : outGraphs) out << g.str() << "\n";
        out << "blocks " << blocks.size() + locBlocks.size() << "\n";
        for (std::size_t b = 0; b < blocks.size(); ++b) {
            const auto& B = blocks[b];
            out << "block " << b << " " << B.m << " " << B.flags.size() << " "
                << (reduced ? mpz_class(mpz_class(static_cast<unsigned long>(B.denominator)) * D).get_str()
                            : std::to_string(B.denominator))
                << " plain "
                << B.sigma.str() << "\n";
            for (const auto& f : B.flags) out << f.str() << "\n";
        }
        for (std::size_t b = 0; b < locBlocks.size(); ++b) {
            const auto& B = locBlocks[b];
            out << "block " << blocks.size() + b << " " << B.m << " " << B.flags.size() << " 1 localizing "
                << B.sigma.str() << "\n";
            for (const auto& f : B.flags) out << f.str() << "\n";
        }
        out << "objective " << objective.size() << "\n";
        for (auto& [g, v] : objective) out << g << " " << v.get_str() << "\n";

        // Optional CSDP input, in the layout of python/flagcheck/sdp.py write_sdpa (same
        // entries; the order of entry lines differs, which the SDPA format allows):
        // LP block = lambda_G (m), d+, d-, (y+_j, y-_j) per equality, u_i per inequality, [nu].
        std::ofstream sdpa;
        const std::size_t m = outGraphs.size();
        std::size_t numEq = 0, numLe = 0;
        for (const auto& C : spec.constraints) ++(C.kind == Constraint::Eq ? numEq : numLe);
        numEq += levelEqs.size();
        const std::size_t lpSize = m + 2 + 2 * numEq + numLe + (normalize ? 1 : 0);
        const std::size_t lpBlock = blocks.size() + locBlocks.size() + 1;
        auto blockDim = [&](std::size_t b) {
            return b < blocks.size() ? blocks[b].flags.size() : locBlocks[b - blocks.size()].flags.size();
        };
        if (!sdpaPath.empty()) {
            sdpa.open(sdpaPath);
            sdpa.precision(17);
            sdpa << m + (normalize ? 1 : 0) << "\n" << lpBlock << "\n";
            for (std::size_t b = 0; b + 1 < lpBlock; ++b) sdpa << blockDim(b) << " ";
            sdpa << "-" << lpSize << "\n";
            std::vector<double> f(m, 0.0);
            for (auto& [g, v] : objective) f[g] = v.get_d();
            for (std::size_t g = 0; g < m; ++g) sdpa << f[g] << " ";
            if (normalize) sdpa << "1.0";
            sdpa << "\n";
            sdpa << "0 " << lpBlock << " " << m + 1 << " " << m + 1 << " 1.0\n";
            sdpa << "0 " << lpBlock << " " << m + 2 << " " << m + 2 << " -1.0\n";
            for (std::size_t g = 0; g < m; ++g) {
                const std::size_t row = g + 1;
                sdpa << row << " " << lpBlock << " " << g + 1 << " " << g + 1 << " 1.0\n";
                sdpa << row << " " << lpBlock << " " << m + 1 << " " << m + 1 << " 1.0\n";
                sdpa << row << " " << lpBlock << " " << m + 2 << " " << m + 2 << " -1.0\n";
            }
        }

        // Constraints, one at a time: each is lifted to n vertices, written to both files and
        // dropped, so memory does not grow with the number of constraints.
        std::size_t eqSeen = 0, leSeen = 0, lpEntries = 0;
        const std::size_t numConstraints = spec.constraints.size() + levelEqs.size();
        const Constraint levelEq{Constraint::Eq, mpq_class(0), nullptr};
        for (std::size_t c = 0; c < numConstraints; ++c) {
            const bool isLevel = c >= spec.constraints.size();
            const auto& C = isLevel ? levelEq : spec.constraints[c];
            const bool eq = C.kind == Constraint::Eq;
            std::vector<std::pair<int, mpq_class>> vec;
            if (isLevel) {
                vec = std::move(levelEqs[c - spec.constraints.size()]);
            } else {
                vec = toVector(T, evaluate(T, *C.expr), n, gi, "constraint");
                if (reduced) vec = reduceVec(vec);
            }
            out << "constraint " << (eq ? "eq" : "le") << " " << C.rhs.get_str() << " " << vec.size() << "\n";
            for (auto& [g, v] : vec) out << g << " " << v.get_str() << "\n";
            if (sdpa.is_open()) {
                // coefficient of graph g in this LP column: value - rhs (eq) or rhs - value (le)
                const std::size_t idx = eq ? m + 3 + 2 * eqSeen : m + 3 + 2 * numEq + leSeen;
                auto emit = [&](std::size_t g, const mpq_class& value) {
                    const double x = mpq_class(eq ? value - C.rhs : C.rhs - value).get_d();
                    if (x == 0) return;
                    sdpa << g + 1 << " " << lpBlock << " " << idx << " " << idx << " " << x << "\n";
                    if (eq) sdpa << g + 1 << " " << lpBlock << " " << idx + 1 << " " << idx + 1 << " " << -x << "\n";
                    lpEntries += eq ? 2 : 1;
                };
                if (C.rhs == 0) {
                    for (auto& [g, v] : vec) emit(g, v);
                } else {  // every graph has a nonzero coefficient unless its value equals rhs
                    const mpq_class zero(0);
                    std::size_t p = 0;
                    for (std::size_t g = 0; g < m; ++g) {
                        if (p < vec.size() && static_cast<std::size_t>(vec[p].first) == g) emit(g, vec[p++].second);
                        else emit(g, zero);
                    }
                }
            }
            ++(eq ? eqSeen : leSeen);
            if ((c + 1) % 1000 == 0 || c + 1 == numConstraints)
                std::cerr << "constraints: " << c + 1 << "/" << numConstraints << ", LP entries "
                          << lpEntries << " (" << seconds(t0) << "s)\n";
        }

        if (sdpa.is_open()) {
            if (normalize) {
                const std::size_t row = m + 1;
                for (std::size_t b = 0; b + 1 < lpBlock; ++b)
                    for (std::size_t i = 1; i <= blockDim(b); ++i) sdpa << row << " " << b + 1 << " " << i << " " << i << " 1.0\n";
                for (std::size_t idx = m + 3; idx <= lpSize; ++idx)
                    sdpa << row << " " << lpBlock << " " << idx << " " << idx << " 1.0\n";
            }
        }

        // Y, computed in parallel chunks and written in graph order. The body goes to a
        // temporary file (its length is only known at the end) instead of memory.
        const std::string yPath = std::string(argv[2]) + ".ybody";
        std::ofstream yBody(yPath);
        if (!yBody) throw std::runtime_error("cannot write " + yPath);
        std::size_t numEntries = 0;
        const std::size_t chunk = 256;
        std::vector<std::vector<YEntry>> allY(reduced && !levels ? graphs.size() : 0);  // full mode: kept
        for (std::size_t start = 0; start < (levels ? 0 : graphs.size()); start += chunk) {
            const std::size_t end = std::min(graphs.size(), start + chunk);
            std::vector<std::vector<YEntry>> results(end - start);
            std::string error;
#pragma omp parallel for schedule(dynamic, 1)
            for (std::size_t g = start; g < end; ++g) {
                try {
                    results[g - start] = computeY(T, blocks, graphs[g]);
                } catch (const std::exception& ex) {
#pragma omp critical
                    error = ex.what();
                }
            }
            if (!error.empty()) throw std::runtime_error(error);
            if (reduced) {
                for (std::size_t g = start; g < end; ++g) allY[g] = std::move(results[g - start]);
                continue;
            }
            for (std::size_t g = start; g < end; ++g)
                for (const auto& e : results[g - start]) {
                    yBody << g << " " << e.block << " " << e.i << " " << e.j << " " << e.count << "\n";
                    ++numEntries;
                    if (sdpa.is_open())
                        sdpa << g + 1 << " " << e.block + 1 << " " << e.i + 1 << " " << e.j + 1 << " "
                             << static_cast<double>(e.count) / static_cast<double>(blocks[e.block].denominator) << "\n";
                }
            std::cerr << "Y: " << end << "/" << graphs.size() << " graphs (" << seconds(t0) << "s)\n";
        }
        if (reduced) {
            // count'_Q(b, i, j) = sum_G (D M[G][Q]) count_G(b, i, j), over the block denominator
            // den_b D: still integer counts, so the file format and the checker are unchanged.
            std::vector<mpz_class> blockDen(blocks.size());
            for (std::size_t b = 0; b < blocks.size(); ++b)
                blockDen[b] = mpz_class(static_cast<unsigned long>(blocks[b].denominator)) * D;
            for (std::size_t start = 0; start < reducedGraphs.size(); start += chunk) {
                const std::size_t end = std::min(reducedGraphs.size(), start + chunk);
                std::vector<std::vector<std::pair<std::uint64_t, mpz_class>>> sums(end - start);
#pragma omp parallel for schedule(dynamic, 1)
                for (std::size_t q = start; q < end; ++q) {
                    std::unordered_map<std::uint64_t, mpz_class> acc;
                    for (const auto& [g, c] : colsM[q]) {
                        // levels: each G is in exactly one column, so computing Y_G here costs nothing extra
                        const std::vector<YEntry> own = levels ? computeY(T, blocks, graphs[g]) : std::vector<YEntry>();
                        for (const auto& e : levels ? own : allY[g]) {
                            const std::uint64_t key = (static_cast<std::uint64_t>(e.block) << 40) |
                                                      (static_cast<std::uint64_t>(e.i) << 20) | e.j;
                            mpz_addmul_ui(acc[key].get_mpz_t(), c.get_mpz_t(), static_cast<unsigned long>(e.count));
                        }
                    }
                    auto& v = sums[q - start];
                    v.assign(acc.begin(), acc.end());
                    std::sort(v.begin(), v.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
                }
                for (std::size_t q = start; q < end; ++q)
                    for (const auto& [key, cnt] : sums[q - start]) {
                        const std::uint64_t b = key >> 40, i = (key >> 20) & 0xFFFFF, j = key & 0xFFFFF;
                        yBody << q << " " << b << " " << i << " " << j << " " << cnt.get_str() << "\n";
                        ++numEntries;
                        if (sdpa.is_open())
                            sdpa << q + 1 << " " << b + 1 << " " << i + 1 << " " << j + 1 << " "
                                 << mpq_class(cnt, blockDen[b]).get_d() << "\n";
                    }
                std::cerr << "reduced Y: " << end << "/" << reducedGraphs.size() << " (" << seconds(t0) << "s)\n";
            }
        }
        for (std::size_t b = 0; b < locBlocks.size(); ++b)
            for (const auto& [key, v] : locBlocks[b].entries) {
                yBody << std::get<0>(key) << " " << blocks.size() + b << " " << std::get<1>(key) << " "
                      << std::get<2>(key) << " " << v.get_str() << "\n";
                ++numEntries;
                if (sdpa.is_open())
                    sdpa << std::get<0>(key) + 1 << " " << blocks.size() + b + 1 << " " << std::get<1>(key) + 1 << " "
                         << std::get<2>(key) + 1 << " " << v.get_d() << "\n";
            }
        if (sdpa.is_open() && !sdpa) throw std::runtime_error("SDPA write failed");
        yBody.close();
        if (!yBody) throw std::runtime_error("write of " + yPath + " failed");
        out << "Y " << numEntries << "\n";
        if (numEntries > 0) {
            std::ifstream yIn(yPath);
            out << yIn.rdbuf();
        }
        std::remove(yPath.c_str());
        out << "end\n";
        if (!out) throw std::runtime_error("write failed");
        std::cerr << "done (" << seconds(t0) << "s)\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "flagsdp: error: " << ex.what() << "\n";
        return 1;
    }
}
