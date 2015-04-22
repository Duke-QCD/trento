// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef NUCLEON_H
#define NUCLEON_H

#include "fast_exp.h"
#include "fwd_decl.h"
#include "random.h"

namespace trento {

class Nucleon;

/// \rst
/// Encapsulates properties shared by all nucleons: transverse thickness
/// profile, cross section, fluctuations.  Responsible for sampling
/// nucleon-nucleon participation with given `\sigma_{NN}`.
/// \endrst
class NucleonProfile {
 public:
  /// Instantiate from the configuration.
  explicit NucleonProfile(const VarMap& var_map);

  /// The radius at which the nucleon profile is truncated.
  double radius() const;

  /// Compute the thickness function at a (squared) distance from the profile
  /// center.
  double thickness(double distance_squared) const;

  /// Randomly determine if a pair of nucleons participates.
  bool participate(Nucleon& A, Nucleon& B) const;

  /// Sample a random nucleon fluctuation.
  double sample_fluct() const;

 private:
  /// Truncate the Gaussian at this number of widths.
  // TODO: is this the optimal value?
  static constexpr double trunc_widths_ = 5.;

  /// Width of Gaussian thickness function.
  const double width_squared_;

  /// Truncate the Gaussian at this radius.
  const double trunc_radius_squared_;

  /// Fast exponential for calculating the thickness profile.
  FastExp<double> fast_exp_;

  /// Fluctuation distribution.
  /// mutable so that sample_fluct() can be const.
  mutable std::gamma_distribution<double> fluct_dist_;

  /// Dimensionless parameter set to reproduce the inelastic nucleon-nucleon
  /// cross section \sigma_{NN}.  Calculated in constructor.
  double cross_sec_param_;
};

/// \rst
/// Represents a single nucleon.  Stores its transverse position and whether or
/// not it's a participant.  These properties are globally readable, but can
/// only be set through ``Nucleus`` and ``NucleonProfile``.
/// \endrst
class Nucleon {
 public:
  /// Only a default constructor is necessary\---the class is designed to be
  /// contsructed once and repeatedly updated.
  Nucleon() = default;

  /// The transverse \em x position.
  double x() const;

  /// The transverse \em y position.
  double y() const;

  /// Whether or not this nucleon is a participant.
  bool is_participant() const;

 private:
  /// A Nucleus must be able to set its Nucleon positions.
  friend class Nucleus;

  /// The NucleonProfile samples participants so must be able to set
  /// participation status.
  friend bool NucleonProfile::participate(Nucleon&, Nucleon&) const;

  /// Set the transverse position and reset participant status to false.
  void set_position(double x, double y);

  /// Mark as a participant.
  void set_participant();

  /// Internal storage of the transverse position.
  double x_, y_;

  /// Internal storage of participant status.
  bool participant_;
};

// Nucleon inline member functions

inline double Nucleon::x() const {
  return x_;
}

inline double Nucleon::y() const {
  return y_;
}

inline bool Nucleon::is_participant() const {
  return participant_;
}

inline void Nucleon::set_position(double x, double y) {
  x_ = x;
  y_ = y;
  participant_ = false;
}

inline void Nucleon::set_participant() {
  participant_ = true;
}

// NucleonProfile inline member functions

inline double NucleonProfile::radius() const {
  return std::sqrt(trunc_radius_squared_);
}

inline double NucleonProfile::thickness(double distance_squared) const {
  if (distance_squared > trunc_radius_squared_)
    return 0.;
  return fast_exp_(-.5*distance_squared/width_squared_) / width_squared_;
}

inline bool NucleonProfile::participate(Nucleon& A, Nucleon& B) const {
  // If both nucleons are already participants, there's nothing to do.
  if (A.is_participant() && B.is_participant())
    return true;

  double dx = A.x() - B.x();
  double dy = A.y() - B.y();
  double distance_squared = dx*dx + dy*dy;

  // If distance > 2R, the nucleons are out of range and the participation
  // probability vanishes.
  if (distance_squared > 4.*trunc_radius_squared_)
    return false;

  // The probability is
  //   P = 1 - exp(...)
  // which can be calculated using std::expm1
  //   expm1() = exp() - 1 = -(1 - exp())
  auto prob = -std::expm1(
      -std::exp(cross_sec_param_ - .25*distance_squared/width_squared_));

  // Sample one random number and decide if this pair participates.
  if (prob > random::canonical<>()) {
    A.set_participant();
    B.set_participant();
    return true;
  }

  return false;
}

inline double NucleonProfile::sample_fluct() const {
  return fluct_dist_(random::engine);
}

}  // namespace trento

#endif  // NUCLEON_H
