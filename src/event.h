// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef EVENT_H
#define EVENT_H

#include <functional>
#include <map>

#ifdef NDEBUG
#define BOOST_DISABLE_ASSERTS
#endif
#include <boost/multi_array.hpp>

#include "fwd_decl.h"

namespace trento {

class NucleonProfile;

///
class Event {
 public:
  ///
  explicit Event(const VarMap& var_map);

  ///
  void compute(const Nucleus& nucleusA, const Nucleus& nucleusB,
               const NucleonProfile& profile);

  /// Alias for a two-dimensional thickness grid.
  using Grid = boost::multi_array<double, 2>;

  ///
  const int& npart() const
  { return npart_; }

  ///
  const double& multiplicity() const
  { return multiplicity_; }

  ///
  const std::map<int, double>& eccentricity() const
  { return eccentricity_; }

  ///
  const Grid& reduced_thickness_grid() const
  { return TR_; }

 private:
  ///
  void compute_nuclear_thickness(
      const Nucleus& nucleus, const NucleonProfile& profile, Grid& TX);

  ///
  template <typename GenMean>
  void compute_reduced_thickness(GenMean gen_mean);

  ///
  std::function<void()> compute_reduced_thickness_;

  ///
  void compute_observables();

  /// Normalization factor.
  const double norm_;

  /// Number of grid steps.
  const int nsteps_;

  /// Grid width and step size.
  const double width_, dxy_;

  /// Nuclear thickness grids TA, TB and reduced thickness grid TR.
  Grid TA_, TB_, TR_;

  /// Center of mass coordinates.
  double xcm_, ycm_;

  /// Number of participants.
  int npart_;

  /// Multiplicity (total entropy).
  double multiplicity_;

  /// Eccentricity harmonics.
  std::map<int, double> eccentricity_;
};

}  // namespace trento

#endif  // EVENT_H
