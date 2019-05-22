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

  /// Whether or not its constituents have been sampled.
  bool constituents_exist() const;

  /// The transverse \em x position.
  double x() const;

  /// The transverse \em y position.
  double y() const;

  /// The longitudinal \em z position.
  double z() const;

 private:
  ///
  friend class NucleonCommon;

  /// A Nucleus must be able to set its Nucleon positions.
  friend class Nucleus;

  /// Set the transverse position and reset participant status to false.
  void set_position(double x, double y, double z);

  /// Internal storage of the transverse position.
  double x_, y_, z_;

  /// Whether or not nucleon is a participant.
  bool is_participant_ = false;

  /// Whether or not nucleon's constituents are sampled.
  bool constituents_exist_ = false;

  /// Constituent transverse position and fluctuation prefactor.
  struct Constituent {
    double x, y;
    double fluctuation;
  };

  /// Vector of constituent positions and fluctuation prefactors.
  std::vector<Constituent> constituents_;
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

  /// Corners of the tile enclosing the nucleon thickness.
  std::array<double, 4> boundary(const NucleonData& nucleon) const;

  /// Nucleon thickness as a function of transverse position.
  double thickness(const NucleonData& nucleon, double x, double y) const;

  /// Randomly determine if a pair of nucleons participates.
  bool participate(NucleonData& A, NucleonData& B) const;

 private:
  /// Sample constituent positions inside the nucleon.
  void sample_constituent_positions(NucleonData& nucleon) const;

  /// Flag the nucleon as a participant.
  void set_participant(NucleonData& nucleon) const;

  /// Fast exponential for calculating the thickness profile.
  const FastExp<double> fast_exp_;

  /// Gaussian nucleon width.
  const double nucleon_width_;

  /// Gaussian constituent width.
  const double constituent_width_;

  /// Number of constituents inside the nucleon.
  const std::size_t constituent_number_;

  /// Gaussian width of constituent position sampling distribution.
  const double sampling_width_;

  /// Maximum impact parameter for nucleon-nucleon participation.
  const double max_impact_sq_;

  /// Constituent width squared.
  const double constituent_width_sq_;

  /// Calculate thickness out to this distance from constituent center.
  const double constituent_radius_sq_;

  /// Nuclear opacity parameter
  double sigma_partonic_;

  /// Thickness function prefactor = 1/(constituent_number*2*pi*w^2) XXX
  const double prefactor_;

  /// Tracks binary collisions if true
  const bool calc_ncoll_;

  /// Gamma random variables used to weight each nucleon (or constituent)
  /// contribution. Controlled by the fluctuation parameter.
  mutable std::gamma_distribution<double> participant_fluctuation_dist_;

  /// Gaussian distribution for sampling constituent positions.
  mutable std::normal_distribution<double> constituent_position_dist_;

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
  is_participant_ = false;
  constituents_exist_ = false;
}

inline bool NucleonData::is_participant() const {
  return is_participant_;
}

inline bool NucleonData::constituents_exist() const {
  return constituents_exist_;
}

// NucleonCommon inline member functions

inline double NucleonCommon::max_impact() const {
  return std::sqrt(max_impact_sq_);
}

inline std::array<double, 4>
NucleonCommon::boundary(const NucleonData& nucleon) const {
  auto constituent = nucleon.constituents_.begin();

  // Initialize the boundary with the position of the first constituent.
  auto xmin = constituent->x, xmax = constituent->x;
  auto ymin = constituent->y, ymax = constituent->y;

  // Check the remaining constituents and update the boundary accordingly.
  // Using this instead of something from std::algorithm because it finds all
  // four quantities {xmin, xmax, ymin, ymax} in a single pass over the constituents.
  for (std::advance(constituent, 1); constituent != nucleon.constituents_.end(); ++constituent) {
    auto x = constituent->x;
    if (x < xmin)
      xmin = x;
    else if (x > xmax)
      xmax = x;

    auto y = constituent->y;
    if (y < ymin)
      ymin = y;
    else if (y > ymax)
      ymax = y;
  }

  // Remember to add and subtract the constituent radius.
  auto r = std::sqrt(constituent_radius_sq_);

  return {xmin - r, xmax + r, ymin - r, ymax + r};
}

inline double NucleonCommon::thickness(
    const NucleonData& nucleon, double x, double y) const {
  auto t = 0.;

  for (const auto& constituent : nucleon.constituents_) {
    auto fluct = constituent.fluctuation;
    auto distance_sq = sqr(x - constituent.x) + sqr(y - constituent.y);
    if (distance_sq < constituent_radius_sq_)
      t += fluct * fast_exp_(-.5*distance_sq/constituent_width_sq_);
  }

  return prefactor_ * t;
}

inline bool NucleonCommon::participate(NucleonData& A, NucleonData& B) const {
  // If both nucleons are already participants, there's nothing to do.
  if (A.is_participant() && B.is_participant() && !calc_ncoll_)
    return true;

  auto distance_sq = sqr(A.x() - B.x()) + sqr(A.y() - B.y());

  // Check if nucleons are out of range.
  if (distance_sq > max_impact_sq_)
    return false;

  sample_constituent_positions(A);
  sample_constituent_positions(B);

  auto overlap = 0.;
  for (const auto& qA : A.constituents_) {
    for (const auto& qB : B.constituents_) {
      auto distance_sq = sqr(qA.x - qB.x) + sqr(qA.y - qB.y);
      overlap += std::exp(-.25*distance_sq/constituent_width_sq_);
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
  auto one_minus_prob = std::exp(
      -sigma_partonic_ * prefactor_/(2.*constituent_number_) * overlap
      );

  // Sample one random number and decide if this pair participates.
  if (one_minus_prob < random::canonical<double>()) {
    set_participant(A);
    set_participant(B);
    return true;
  }

  return false;
}

inline void NucleonCommon::sample_constituent_positions(NucleonData& nucleon) const {
  if (nucleon.constituents_exist())
    return;

  nucleon.constituents_.resize(static_cast<std::size_t>(constituent_number_));

  double xcom = 0.0;
  double ycom = 0.0;

  // Sample nucleon constituent positions
  for (auto&& constituent : nucleon.constituents_) {
    auto xloc = constituent_position_dist_(random::engine);
    auto yloc = constituent_position_dist_(random::engine);

    constituent.x = xloc;
    constituent.y = yloc;
    constituent.fluctuation = participant_fluctuation_dist_(random::engine);

    xcom += xloc;
    ycom += yloc;
  }

  xcom /= constituent_number_;
  ycom /= constituent_number_;

  // Place nucleon at the desired position
  for (auto&& constituent : nucleon.constituents_) {
    constituent.x += (nucleon.x() - xcom);
    constituent.y += (nucleon.y() - ycom);
  }

  nucleon.constituents_exist_ = true;
}

inline void NucleonCommon::set_participant(NucleonData& nucleon) const {
  if (nucleon.is_participant())
    return;

  nucleon.is_participant_ = true;
}

class MonteCarloCrossSection {
 public:
  MonteCarloCrossSection(const VarMap& var_map);

  double operator() (const double sigma_partonic) const;

 private:
  const double nucleon_width_;
  const double constituent_width_;
  const std::size_t constituent_number_;
  const double sampling_width_;
  const double max_impact_;
  const double constituent_width_sq_;
  const double prefactor_;

  const std::size_t n_max = 10000000;
  const std::size_t cache_size = 1000000;
  const std::size_t n_loops = 10;
  const int n_pass = 10000;
  const double tolerance = 0.001;
};

}  // namespace trento

#endif  // NUCLEON_H
