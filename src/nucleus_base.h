// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef NUCLEUS_BASE_H
#define NUCLEUS_BASE_H

// Forward declaration of the NucleusBase interface class only.

#include <memory>
#include <string>

#include "fwd_decl.h"

namespace trento {

///
using NucleusPtr = std::unique_ptr<NucleusBase>;

///
class NucleusBase {
 public:
  ///
  static NucleusPtr create(const std::string& species);

  ///
  virtual void sample_nucleons(double offset) = 0;

  ///
  virtual ~NucleusBase() = default;

 protected:
  ///
  NucleusBase() = default;
};

}  // namespace trento

#endif  // NUCLEUS_BASE_H
