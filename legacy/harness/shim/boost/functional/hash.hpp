// Minimal stand-in for the one piece of Boost the legacy flags.h uses:
// boost::hash<std::pair<int,int>>. Only used to build the legacy code for testing.
#pragma once
#include <cstddef>
#include <functional>
#include <utility>

namespace boost {
template <class T>
struct hash : std::hash<T> {};

template <class A, class B>
struct hash<std::pair<A, B>> {
    std::size_t operator()(const std::pair<A, B>& p) const {
        std::size_t seed = std::hash<A>()(p.first);
        seed ^= std::hash<B>()(p.second) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};
}  // namespace boost
