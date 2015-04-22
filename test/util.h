// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef UTIL_H
#define UTIL_H

#include <iostream>
#include <sstream>

#include <boost/program_options/variables_map.hpp>

#include "../src/fwd_decl.h"

// testing utilities

// Factory function; create a dummy boost::program_options::variables_map.
VarMap make_var_map(std::map<std::string, boost::any>&& args);

// redirect stdout to a stringstream and safely restore upon destruction
struct capture_stdout {
  capture_stdout() {
    // replace stdout buffer
    cout_orig = std::cout.rdbuf(stream.rdbuf());
  }

  ~capture_stdout() {
    // restore stdout to working state
    std::cout.rdbuf(cout_orig);
  }

  std::streambuf* cout_orig;
  std::stringstream stream;
};

#endif  // UTIL_H
