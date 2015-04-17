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
  // W-S params ref. given in WoodsSaxonNucleus class def.
  // XXX: remember to add new species to the help output in main()
  if (species == "p")
    return NucleusPtr{new Proton{}};
  else if (species == "Au")
    return NucleusPtr{new WoodsSaxonNucleus{197, 6.38, 0.535}};
  else if (species == "Pb")
    return NucleusPtr{new WoodsSaxonNucleus{208, 6.62, 0.546}};
  else
    throw std::invalid_argument{"unknown projectile species: " + species};
}

Nucleus::Nucleus(std::size_t A) : nucleons_(A) {}

template <typename... Args>
void Nucleus::set_nucleon_position(Nucleon& nucleon, Args&&... args) {
  nucleon.set_position(std::forward<Args>(args)...);
}

Proton::Proton() : Nucleus(1) {}

double Proton::radius() const {
  return 0.;
}

void Proton::sample_nucleons(double offset) {
  // set (x, y) = (offset, 0)
  set_nucleon_position(*begin(), offset,  0.);
}

/// The Woods-Saxon distribution is
/// \f[
///   f(r) \propto \frac{1}{1 + \exp(\frac{r-R}{a})}.
/// \f]
/// Extend it out to \f$R + 10a\f$; for typical values of \f$(R, a)\f$, the
/// probability of sampling a nucleon beyond this radius is
/// \f$\mathcal O(10^{-5})\f$.
WoodsSaxonNucleus::WoodsSaxonNucleus(std::size_t A, double R, double a)
    : Nucleus(A),
      R_(R),
      a_(a),
      woods_saxon_dist_(1000, 0., R + 10.*a,
        [R, a](double r) { return r*r/(1.+std::exp((r-R)/a)); })
{}

double WoodsSaxonNucleus::radius() const {
  return R_ + 3.*a_;
}

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

}  // namespace trento
