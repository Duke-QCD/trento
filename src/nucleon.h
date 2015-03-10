// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef NUCLEON_H
#define NUCLEON_H

#include "fwd_decl.h"
#include "random.h"

namespace trento {

/// Nucleon class.
class Nucleon {
 public:
  /// \brief Initialize nucleon profile parameters, in particular determine
  /// \c cross_sec_param_.
  ///
  /// \param var_map Must contain keys
  /// \c "nucleon-width",
  /// \c "fluctation",
  /// \c "beam-energy",
  /// \c "cross-section".
  explicit Nucleon(const VarMap& var_map);

  /// The radius at which the nucleon profile is truncated.
  double radius() const;

  /// Sample a random nucleon fluctuation.
  double fluctuate() const;

  /// Calculate the thickness function.
  double thickness(double distance_squared) const;

  /// Randomly determine nucleon-nucleon participation.
  bool participate(double b_squared) const;

 private:
  /// Truncate the Gaussian at this number of widths.
  static constexpr double trunc_widths_ = 3.;

  /// Width of Gaussian thickness function.
  const double width_squared_;

  /// Truncate the Gaussian at this radius.
  const double trunc_radius_squared_;

  /// Dimensionless parameter set to reproduce the inelastic nucleon-nucleon
  /// cross section \f$\sigma_{NN}\f$.  Calculated in constructor.
  double cross_sec_param_;

  /// Fluctuation distribution.
  // mutable so that fluctuate() can be const
  mutable std::gamma_distribution<double> fluct_dist_;
};

inline double Nucleon::radius() const {
  return std::sqrt(trunc_radius_squared_);
}

inline double Nucleon::fluctuate() const {
  return fluct_dist_(random::engine);
}

inline double Nucleon::thickness(double distance_squared) const {
  if (distance_squared > trunc_radius_squared_)
    return 0.;

  return 1./width_squared_ * std::exp(-.5*distance_squared/width_squared_);
}

inline bool Nucleon::participate(double b_squared) const {
  if (b_squared > 4.*trunc_radius_squared_)
    return false;

  // The probability is
  //   P = 1 - exp(...)
  // which can be calculated using std::expm1
  //   expm1() = exp() - 1 = -(1 - exp())
  auto prob = -std::expm1(
      -std::exp(cross_sec_param_ - .25*b_squared/width_squared_));

  return prob > random::canonical<>();
}

}  // namespace trento

#endif  // NUCLEON_H
