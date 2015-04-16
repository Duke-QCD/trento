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

///
class NucleonProfile {
 public:
  /// \brief Initialize nucleon profile parameters, in particular determine
  /// \c cross_sec_param_.
  ///
  /// \param var_map Must contain keys
  /// \c "nucleon-width",
  /// \c "fluctation",
  /// \c "beam-energy",
  /// \c "cross-section".
  explicit NucleonProfile(const VarMap& var_map);

  /// The radius at which the nucleon profile is truncated.
  double radius() const;

  /// Calculate the thickness function.
  double thickness(double distance_squared) const;

  /// Randomly determine nucleon-nucleon participation.
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
  /// \c mutable so that \c sample_fluct() can be \c const.
  mutable std::gamma_distribution<double> fluct_dist_;

  /// Dimensionless parameter set to reproduce the inelastic nucleon-nucleon
  /// cross section \f$\sigma_{NN}\f$.  Calculated in constructor.
  double cross_sec_param_;
};

///
class Nucleon {
 public:
  ///
  Nucleon() = default;

  ///
  double x() const;

  ///
  double y() const;

  ///
  bool is_participant() const;

 private:
  ///
  friend class Nucleus;
  friend bool NucleonProfile::participate(Nucleon&, Nucleon&) const;

  ///
  void set_position(double x, double y);

  ///
  void set_participant();

  ///
  double x_, y_;

  ///
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
  // Nothing to do.
  if (A.is_participant() && B.is_participant())
    return true;

  double dx = A.x() - B.x();
  double dy = A.y() - B.y();
  double distance_squared = dx*dx + dy*dy;

  // Out of range.
  if (distance_squared > 4.*trunc_radius_squared_)
    return false;

  // The probability is
  //   P = 1 - exp(...)
  // which can be calculated using std::expm1
  //   expm1() = exp() - 1 = -(1 - exp())
  auto prob = -std::expm1(
      -std::exp(cross_sec_param_ - .25*distance_squared/width_squared_));

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
