// Test harness for the legacy (thesis-era) plain flag algebra code.
//
// Builds standard problems with the legacy Graph/Equation API and calls the legacy
// plainFlagAlgebra, which writes sdp.cbf and Duals.txt to the working directory and
// prints "Divide answers by <mult>". The Python side (flagcheck/legacy.py) converts the
// CBF file for CSDP and compares against the independent reference.
//
// Usage: harness <problem> <n> [args...]
//   goodman n | mantel n | ind_p3 n | ind_p4 n | ind_c5 n
//   ramsey_obj k l n
//   ramsey_feas k l N n          (green K_j = 1/N^(j-1) for j = 2..n; objective f = 1)
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <queue>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <omp.h>
#include <boost/functional/hash.hpp>

#include "general.h"
#include "fractions.h"
#include "flags_legacy.h"

static Graph clique(int k, int color, int numColors) {
    std::vector<Edge> edges;
    for (int i = 0; i < k; ++i)
        for (int j = i + 1; j < k; ++j) edges.push_back({i, j, color});
    return Graph(edges, k, numColors);
}

static std::vector<Graph> ramseyZeros(int k, int l) {
    const int c = 3;
    std::vector<Graph> zeros;
    zeros.push_back(Graph({{0, 1, 1}}, 3, c));
    zeros.push_back(Graph({{0, 1, 2}}, 3, c));
    zeros.push_back(Graph({{0, 1, 1}, {1, 2, 2}}, 3, c));
    zeros.push_back(clique(k, 1, c));
    zeros.push_back(clique(l, 2, c));
    return zeros;
}

int main(int argc, char** argv) {
    omp_set_num_threads(1);
    if (argc < 3) {
        std::cerr << "usage: harness <problem> <n> ...\n";
        return 2;
    }
    std::string prob = argv[1];
    std::vector<Graph> f, zeros;
    std::vector<Equation> known;
    bool maximize = false;
    int n = 0;

    if (prob == "goodman") {
        n = std::atoi(argv[2]);
        f = {clique(3, 1, 2), clique(3, 0, 2)};
    } else if (prob == "mantel") {
        n = std::atoi(argv[2]);
        zeros = {clique(3, 1, 2)};
        f = {Graph({{0, 1, 1}}, 2, 2)};
        maximize = true;
    } else if (prob == "ind_p3" || prob == "ind_p4" || prob == "ind_c5") {
        n = std::atoi(argv[2]);
        if (prob == "ind_p3") f = {Graph({{0, 1, 1}, {1, 2, 1}}, 3, 2)};
        if (prob == "ind_p4") f = {Graph({{0, 1, 1}, {1, 2, 1}, {2, 3, 1}}, 4, 2)};
        if (prob == "ind_c5") f = {Graph({{0, 1, 1}, {1, 2, 1}, {2, 3, 1}, {3, 4, 1}, {4, 0, 1}}, 5, 2)};
        maximize = true;
    } else if (prob == "ramsey_obj") {
        int k = std::atoi(argv[2]), l = std::atoi(argv[3]);
        n = std::atoi(argv[4]);
        zeros = ramseyZeros(k, l);
        f = {Graph({}, 2, 3)};
    } else if (prob == "ramsey_feas") {
        int k = std::atoi(argv[2]), l = std::atoi(argv[3]), N = std::atoi(argv[4]);
        n = std::atoi(argv[5]);
        zeros = ramseyZeros(k, l);
        // Objective f = 1 (all three coloured K2s). The legacy resize rejects an empty
        // objective. With the normalisation row added by the Python converter, a feasible
        // program gives exactly d* = 1 and an infeasible one gives d* > 1.
        f = {Graph({}, 2, 3), Graph({{0, 1, 1}}, 2, 3), Graph({{0, 1, 2}}, 2, 3)};
        for (int j = 2; j <= n; ++j) {
            known.push_back(Equation({Graph({}, j, 3)}, zeros, Frac(1, myPow(N, j - 1)), 0));
        }
    } else {
        std::cerr << "unknown problem " << prob << "\n";
        return 2;
    }

    plainFlagAlgebra(f, n, zeros, known, maximize);
    std::cout << "HARNESS maximize=" << (maximize ? 1 : 0) << std::endl;
    return 0;
}
