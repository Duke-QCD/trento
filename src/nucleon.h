// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef NUCLEON_H
#define NUCLEON_H

#include <array>
#include <vector>

#include <boost/math/constants/constants.hpp>

#include "fast_exp.h"
#include "fwd_decl.h"
#include "random.h"

namespace trento {

// TODO explain this class structure

/// \rst
/// Represents a single nucleon.  Stores its transverse position and whether or
/// not it's a participant.  These properties are globally readable, but can
/// only be set through ``Nucleus`` and ``NucleonProfile``.
/// \endrst
class NucleonData {
 public:
  /// Only a default constructor is necessary\---the class is designed to be
  /// constructed once and repeatedly updated.
  NucleonData() = default;

  /// Whether or not this nucleon is a participant.
  bool is_participant() const;

  ///
  bool partons_sampled() const;

  /// The transverse \em x position.
  double x() const;

  /// The transverse \em y position.
  double y() const;

  /// The longitudinal \em z position.
  double z() const;

  ///
  double fluctuation() const;

 private:
  ///
  friend class NucleonCommon;

  /// A Nucleus must be able to set its Nucleon positions.
  friend class Nucleus;

  /// Set the transverse position and reset participant status to false.
  void set_position(double x, double y, double z);

  /// Internal storage of the transverse position.
  double x_, y_, z_;

  ///
  double fluctuation_;

  ///
  bool partons_exist = false;

  ///
  struct Parton {
    double x, y;
  };

  ///
  std::vector<Parton> partons_;
};

/// \rst
/// Encapsulates properties shared by all nucleons: transverse thickness
/// profile, cross section, fluctuations.  Responsible for sampling
/// nucleon-nucleon participation with given `\sigma_{NN}`.
/// \endrst
class NucleonCommon {
 public:
  /// Instantiate from the configuration.
  explicit NucleonCommon(const VarMap& var_map);

  /// The maximum impact parameter for participation.
  double max_impact() const;

  ///
  std::array<double, 4> boundary(const NucleonData& nucleon) const;

  ///
  double thickness(const NucleonData& nucleon, double x, double y) const;

  /// Randomly determine if a pair of nucleons participates.
  bool participate(NucleonData& A, NucleonData& B) const;

 private:
  ///
  void sample_parton_positions(NucleonData& nucleon) const;

  ///
  void set_participant(NucleonData& nucleon) const;

  /// Fast exponential for calculating the thickness profile.
  const FastExp<double> fast_exp_;

  /// Maximum impact parameter for participants.
  const double max_impact_sq_;

  ///
  const double parton_width_sq_, parton_radius_sq_;

  ///
  const int npartons_;

  ///
  const double sigma_qq_;

  /// Thickness function prefactor = 1/(npartons*2*pi*w^2) XXX
  const double prefactor_;

  /// Gamma distribution for nucleon fluctuations.
  mutable std::gamma_distribution<double> nucleon_fluctuation_dist_;

  ///
  mutable std::normal_distribution<double> parton_position_dist_;
};

// These functions are short, called very often, and account for a large
// fraction of the total computation time, so request inlining.

// Trivial helper function.
template <typename T>
inline constexpr T sqr(const T& value) {
  return value * value;
}

// NucleonData inline member functions

inline double NucleonData::x() const {
  return x_;
}

inline double NucleonData::y() const {
  return y_;
}

inline double NucleonData::z() const {
  return z_;
}

inline void NucleonData::set_position(double x, double y, double z) {
  x_ = x;
  y_ = y;
  z_ = z;
  fluctuation_ = -1.;
}

inline bool NucleonData::is_participant() const {
  return fluctuation_ > 0.;
}

inline bool NucleonData::partons_sampled() const {
  return partons_exist; 
}

// NucleonCommon inline member functions

inline double NucleonCommon::max_impact() const {
  return std::sqrt(max_impact_sq_);
}

inline std::array<double, 4>
NucleonCommon::boundary(const NucleonData& nucleon) const {
  auto parton = nucleon.partons_.begin();

  // Initialize the boundary with the position of the first parton.
  auto xmin = parton->x, xmax = parton->x;
  auto ymin = parton->y, ymax = parton->y;

  // Check the remaining partons and update the boundary accordingly.
  // Using this instead of something from std::algorithm because it finds all
  // four quantities {xmin, xmax, ymin, ymax} in a single pass over the partons.
  for (std::advance(parton, 1); parton != nucleon.partons_.end(); ++parton) {
    auto x = parton->x;
    if (x < xmin)
      xmin = x;
    else if (x > xmax)
      xmax = x;

    auto y = parton->y;
    if (y < ymin)
      ymin = y;
    else if (y > ymax)
      ymax = y;
  }

  // Remember to add and subtract the parton radius.
  auto r = std::sqrt(parton_radius_sq_);

  return {xmin - r, xmax + r, ymin - r, ymax + r};
}

inline double NucleonCommon::thickness(
    const NucleonData& nucleon, double x, double y) const {
  auto t = 0.;

  for (const auto& parton : nucleon.partons_) {
    auto distance_sq = std::pow(x - parton.x, 2) + std::pow(y - parton.y, 2);
    if (distance_sq < parton_radius_sq_)
      t += fast_exp_(-.5*distance_sq/parton_width_sq_);
  }

  return nucleon.fluctuation_ * prefactor_ * t;
}

inline bool NucleonCommon::participate(NucleonData& A, NucleonData& B) const {
  // If both nucleons are already participants, there's nothing to do.
  if (A.is_participant() && B.is_participant())
    return true;

  auto distance_sq = sqr(A.x() - B.x()) + sqr(A.y() - B.y());

  // Check if nucleons are out of range.
  if (distance_sq > max_impact_sq_)
    return false;

  sample_parton_positions(A);
  sample_parton_positions(B);

  auto overlap = 0.;
  for (const auto& qA : A.partons_) {
    for (const auto& qB : B.partons_) {
      auto distance_sq = sqr(qA.x - qB.x) + sqr(qA.y - qB.y);
      overlap += std::exp(-.25*distance_sq/parton_width_sq_);
    }
  }

  // The probability is
  //   P = 1 - exp(...)
  // which we could sample as
  //   P > U
  // where U is a standard uniform (0, 1) random number.  We can also compute
  //   1 - P = exp(...)
  // and then sample
  //   (1 - P) > (1 - U)
  // or equivalently
  //   (1 - P) < U
  auto one_minus_prob = std::exp(-sigma_qq_ * prefactor_/(2*npartons_) * overlap);

  // Sample one random number and decide if this pair participates.
  if (one_minus_prob < random::canonical<double>()) {
    set_participant(A);
    set_participant(B);
    return true;
  }

  return false;
}

inline void NucleonCommon::sample_parton_positions(NucleonData& nucleon) const {
  if (nucleon.partons_sampled())
    return;

  nucleon.partons_.resize(static_cast<std::size_t>(npartons_));

  for (auto&& parton : nucleon.partons_) {
    parton.x = nucleon.x() + parton_position_dist_(random::engine);
    parton.y = nucleon.y() + parton_position_dist_(random::engine);
  }
}

inline void NucleonCommon::set_participant(NucleonData& nucleon) const {
  if (nucleon.is_participant())
    return;

  nucleon.fluctuation_ = nucleon_fluctuation_dist_(random::engine);
}

}  // namespace trento

#endif  // NUCLEON_H
