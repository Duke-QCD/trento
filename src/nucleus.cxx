// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "nucleus.h"

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <boost/math/constants/constants.hpp>

#include "random.h"

namespace trento {

NucleusPtr Nucleus::create(const std::string& species) {
  // W-S params ref. in header
  // XXX: remember to add new species to the help output in main()
  if (species == "p")
    return NucleusPtr{new Proton{}};
  else if (species == "Cu")
    return NucleusPtr{new WoodsSaxonNucleus{62, 4.2, 0.596}};
  else if (species == "Cu2")
    return NucleusPtr{new DeformedWoodsSaxonNucleus{62, 4.2, 0.596, 0.162, -0.006}};
  else if (species == "Au")
    return NucleusPtr{new WoodsSaxonNucleus{197, 6.38, 0.535}};
  else if (species == "Au2")
    return NucleusPtr{new DeformedWoodsSaxonNucleus{197, 6.38, 0.535, -0.131, -0.031}};
  else if (species == "Pb")
    return NucleusPtr{new WoodsSaxonNucleus{208, 6.62, 0.546}};
  else if (species == "U")
    return NucleusPtr{new DeformedWoodsSaxonNucleus{238, 6.67, 0.44, 0.280, 0.093}};
  else
    throw std::invalid_argument{"unknown projectile species: " + species};
}

Nucleus::Nucleus(std::size_t A) : nucleons_(A) {}

template <typename... Args>
void Nucleus::set_nucleon_position(Nucleon& nucleon, Args&&... args) {
  nucleon.set_position(std::forward<Args>(args)...);
}

Proton::Proton() : Nucleus(1) {}

/// Always zero.
double Proton::radius() const {
  return 0.;
}

/// Always place the nucleon at (x, y) = (offset, 0).
void Proton::sample_nucleons(double offset) {
  set_nucleon_position(*begin(), offset,  0.);
}

// Extend the W-S dist out to R + 10a; for typical values of (R, a), the
// probability of sampling a nucleon beyond this radius is O(10^-5).
WoodsSaxonNucleus::WoodsSaxonNucleus(std::size_t A, double R, double a)
    : Nucleus(A),
      R_(R),
      a_(a),
      woods_saxon_dist_(1000, 0., R + 10.*a,
        [R, a](double r) { return r*r/(1.+std::exp((r-R)/a)); })
{}

/// Return something a bit smaller than the true maximum radius.  The
/// Woods-Saxon distribution falls off very rapidly (exponentially), and since
/// this radius determines the impact parameter range, the true maximum radius
/// would cause far too many events with zero participants.
double WoodsSaxonNucleus::radius() const {
  return R_ + 3.*a_;
}

/// Sample uncorrelated Woods-Saxon nucleon positions.
void WoodsSaxonNucleus::sample_nucleons(double offset) {
  for (auto&& nucleon : *this) {
    // Sample spherical radius from Woods-Saxon distribution.
    double r = woods_saxon_dist_(random::engine);

    // Sample isotropic spherical angles.
    double cos_theta = 2. * random::canonical<>() - 1.;
    double phi = math::double_constants::two_pi * random::canonical<>();

    double r_sin_theta = r * std::sqrt(1. - cos_theta*cos_theta);

    // Convert to transverse Cartesian coordinates
    double x = r_sin_theta * std::cos(phi) + offset;
    double y = r_sin_theta * std::sin(phi);

    set_nucleon_position(nucleon, x, y);
  }
  // XXX: re-center nucleon positions?
}

// Set rmax like the non-deformed case (R + 10a), but for the maximum
// "effective" radius.  The numerical coefficients for beta2 and beta4 are the
// approximate values of Y20 and Y40 at theta = 0.
DeformedWoodsSaxonNucleus::DeformedWoodsSaxonNucleus(
    std::size_t A, double R, double a, double beta2, double beta4)
    : Nucleus(A),
      R_(R),
      a_(a),
      beta2_(beta2),
      beta4_(beta4),
      rmax_(R*(1. + .63*std::fabs(beta2) + .85*std::fabs(beta4)) + 10.*a)
{}

/// Return something a bit smaller than the true maximum radius.  The
/// Woods-Saxon distribution falls off very rapidly (exponentially), and since
/// this radius determines the impact parameter range, the true maximum radius
/// would cause far too many events with zero participants.
double DeformedWoodsSaxonNucleus::radius() const {
  return rmax_ - 7.*a_;
}

double DeformedWoodsSaxonNucleus::deformed_woods_saxon_dist(
    double r, double cos_theta) const {
  auto cos_theta_sq = cos_theta*cos_theta;

  // spherical harmonics
  using math::double_constants::one_div_root_pi;
  auto Y20 = std::sqrt(5)/4. * one_div_root_pi * (3.*cos_theta_sq - 1.);
  auto Y40 = 3./16. * one_div_root_pi *
             (35.*cos_theta_sq*cos_theta_sq - 30.*cos_theta_sq + 3.);

  // "effective" radius
  auto Reff = R_ * (1. + beta2_*Y20 + beta4_*Y40);

  return 1. / (1. + std::exp((r - Reff) / a_));
}

/// Sample uncorrelated deformed Woods-Saxon nucleon positions.
void DeformedWoodsSaxonNucleus::sample_nucleons(double offset) {
  // The deformed W-S distribution is defined so the symmetry axis is aligned
  // with the Z axis, so e.g. the long axis of uranium coincides with Z.
  //
  // After sampling positions, they must be randomly rotated.  In general this
  // requires three Euler rotations, but in this case we only need two
  // because there is no use in rotating about the nuclear symmetry axis.
  //
  // The two rotations are:
  //  - a polar "tilt", i.e. rotation about the X axis
  //  - an azimuthal "spin", i.e. rotation about the original Z axis

  // "tilt" angle
  const double cos_a = 2. * random::canonical<>() - 1.;
  const double sin_a = std::sqrt(1. - cos_a*cos_a);

  // "spin" angle
  const double angle_b = math::double_constants::two_pi * random::canonical<>();
  const double cos_b = std::cos(angle_b);
  const double sin_b = std::sin(angle_b);

  for (auto&& nucleon : *this) {
    // Sample (r, theta) using a standard rejection method.
    // Remember to include the phase-space factors.
    double r, cos_theta;
    do {
      r = rmax_ * std::cbrt(random::canonical<>());
      cos_theta = 2. * random::canonical<>() - 1.;
    } while (random::canonical<>() > deformed_woods_saxon_dist(r, cos_theta));

    // Sample azimuthal angle.
    double phi = math::double_constants::two_pi * random::canonical<>();

    // Convert to Cartesian coordinates.
    double r_sin_theta = r * std::sqrt(1. - cos_theta*cos_theta);
    double x = r_sin_theta * std::cos(phi);
    double y = r_sin_theta * std::sin(phi);
    double z = r * cos_theta;

    // Rotate.
    // The rotation formula was derived by composing the "tilt" and "spin"
    // rotations described above.
    double x_rot = x*cos_b - y*cos_a*sin_b + z*sin_a*sin_b + offset;
    double y_rot = x*sin_b + y*cos_a*cos_b - z*sin_a*cos_b;

    set_nucleon_position(nucleon, x_rot, y_rot);
  }
}

}  // namespace trento
