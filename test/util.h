// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef UTIL_H
#define UTIL_H

// testing utilities

#include <boost/program_options/variables_map.hpp>

#include "../src/fwd_decl.h"

// Factory function; create a dummy boost::program_options::variables_map.
VarMap make_var_map(std::map<std::string, boost::any>&& args);

#endif  // UTIL_H
