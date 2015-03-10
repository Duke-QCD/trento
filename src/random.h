// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef RANDOM_H
#define RANDOM_H

#include <limits>
#include <random>

namespace trento { namespace random {

/// Mersenne Twister engine, 64-bit preset.
using Engine = std::mt19937_64;

/// Global variable defined in \c random.cxx.
extern Engine engine;

/// Helper function to easily generate random numbers in [0, 1).
template <typename RealType = double>
inline RealType canonical() {
  return std::generate_canonical
           <RealType, std::numeric_limits<RealType>::digits>
           (engine);
}

}}  // namespace trento::random

#endif  // RANDOM_H
