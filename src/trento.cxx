// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include <iostream>

// CMake sets this definition.
// Fall back to a sane default.
#ifndef TRENTO_VERSION_STRING
#define TRENTO_VERSION_STRING "dev"
#endif

int main(int argc, char *argv[]) {
  if (argc == 2 && std::string(argv[1]) == "--version")
    std::cout << "TRENTO " << TRENTO_VERSION_STRING << '\n';

  return 0;
}
