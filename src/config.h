// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef CONFIG_H
#define CONFIG_H

#include "fwd_decl.h"

namespace trento {

/// Set runtime static member variables from user input.
void set_static_vars(const po::variables_map&);

}  // namespace trento

#endif  // CONFIG_H
