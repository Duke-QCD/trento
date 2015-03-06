// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "random.h"

// Seed random number generator from hardware device.
std::mt19937_64 trento::random_engine{std::random_device{}()};
