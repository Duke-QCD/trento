// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef RANDOM_H
#define RANDOM_H

#include <random>

namespace trento {

// Mersenne Twister engine, 64-bit preset.
// Global variable defined in random.cxx.
extern std::mt19937_64 random_engine;

}  // namespace trento::random

#endif  // RANDOM_H
