// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef COLLIDER_H
#define COLLIDER_H

#include <memory>

#include "fwd_decl.h"
#include "event.h"
#include "nucleon.h"

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
  void new_event();

  ///
  void compute_nuclear_thickness_and_npart(const Nucleus& nucleus,
                                           Event::Grid& TX);

  ///
  std::function<void(void)> compute_reduced_thickness;

  ///
  void compute_observables();

  ///
  std::unique_ptr<Nucleus> nucleusA_, nucleusB_;

  ///
  Nucleon nucleon_;

  ///
  OutputFunctionVector output_functions_;

  ///
  const int num_events_;

  ///
  const double b_min_, b_max_;

  ///
  const double asymmetry_;

  ///
  const int grid_steps_;

  ///
  const double grid_half_width_, grid_delta_;

  ///
  Event event_;
};

}  // namespace trento

#endif  // COLLIDER_H
