// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef COLLIDER_H
#define COLLIDER_H

#include <memory>

#include "fwd_decl.h"
#include "event.h"
#include "nucleon.h"
#include "output.h"

namespace trento {

///
class Collider {
 public:
  ///
  explicit Collider(const VarMap& var_map);

  /// Declare a destructor to properly handle the \c std::unique_ptr<Nucleus>
  /// members.  At this point in the code, \c Nucleus is an incomplete type so
  /// the compiler does not know how to delete it.  Therefore the destructor is
  /// defined (as default) in the implementation file, at which point \c Nucleus
  /// is fully defined.  See e.g. Item 22 of Effective Modern C++ by Scott
  /// Meyers and this stackoverflow answer: http://stackoverflow.com/a/6089065.
  ~Collider();

  ///
  void run_events();

 private:
  ///
  double sample_impact_param();

  ///
  std::unique_ptr<Nucleus> nucleusA_, nucleusB_;

  ///
  NucleonProfile nucleon_profile_;

  ///
  const int nevents_;

  ///
  const double bmin_, bmax_;

  ///
  const double asymmetry_;

  ///
  Event event_;

  ///
  Output output_;
};

}  // namespace trento

#endif  // COLLIDER_H
