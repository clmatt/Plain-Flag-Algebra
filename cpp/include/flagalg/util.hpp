#pragma once

#include <cstdint>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace flagalg {

inline int omp_threads() {
#ifdef _OPENMP
    return omp_get_num_threads();
#else
    return 1;
#endif
}

inline int omp_thread() {
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}

// Advance a sorted k-subset of {0..n-1} to the next one in lexicographic order.
// Returns false (leaving `sub` unspecified) after the last subset. The empty subset
// has exactly one subset, so it returns false immediately.
inline bool nextSubset(std::vector<int>& sub, int n) {
    const int k = static_cast<int>(sub.size());
    int i = k - 1;
    while (i >= 0 && sub[i] == n - k + i) --i;
    if (i < 0) return false;
    ++sub[i];
    for (int j = i + 1; j < k; ++j) sub[j] = sub[j - 1] + 1;
    return true;
}

inline std::uint64_t binom(int n, int k) {
    if (k < 0 || k > n) return 0;
    std::uint64_t r = 1;
    for (int i = 1; i <= k; ++i) r = r * (n - k + i) / i;  // exact at every step
    return r;
}

inline std::uint64_t fallingFactorial(int n, int k) {
    std::uint64_t r = 1;
    for (int i = 0; i < k; ++i) r *= static_cast<std::uint64_t>(n - i);
    return r;
}

}  // namespace flagalg
