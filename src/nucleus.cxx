// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "nucleus.h"

#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <boost/math/constants/constants.hpp>

#include "random.h"

namespace trento {

namespace {

using Coord = NucleusBase::TransverseCoord;

// Nucleon sampling policies.

/// Trivial sampler that always returns the same position.
/// For protons.
class FixedSampler {
 public:
  ///
  constexpr FixedSampler() = default;

  ///
  constexpr double radius() const noexcept { return 0.; };

  ///
  constexpr Coord sample() const noexcept { return {0., 0.}; }
};

/// Samples nucleons from a spherically symmetric Woods-Saxon distribution.
/// For non-deformed heavy nuclei such as lead.
class WoodsSaxonSampler {
 public:
  ///
  WoodsSaxonSampler(double R, double a);

  ///
  double radius() const { return R_ + 3.*a_; };

  ///
  Coord sample() const;

 private:
  ///
  const double R_, a_;

  ///
  mutable std::piecewise_linear_distribution<double> woods_saxon_dist_;
};

/// Extend Woods-Saxon distribution out to R + 10a.
/// For typical values of (R, a), the probability of sampling a nucleon beyond
/// this radius is O(10^-5).
WoodsSaxonSampler::WoodsSaxonSampler(double R, double a)
    : R_(R),
      a_(a),
      woods_saxon_dist_(1000, 0., R + 10.*a,
        [&R, &a](double r) { return r*r/(1.+std::exp((r-R)/a)); })
  {}

Coord WoodsSaxonSampler::sample() const {
  auto r = woods_saxon_dist_(random::engine);
  auto cos_theta = 2. * random::canonical<>() - 1.;
  auto phi = math::double_constants::two_pi * random::canonical<>();

  auto r_sin_theta = r * std::sqrt(1. - cos_theta*cos_theta);

  return {r_sin_theta*std::cos(phi), r_sin_theta*std::sin(phi)};
}

}  // unnamed namespace

NucleusPtr NucleusBase::create(const std::string& species) {
  if (species == "p")
    return NucleusPtr{new Nucleus<1, FixedSampler>{}};
  else if (species == "Pb")
    return NucleusPtr{new Nucleus<208, WoodsSaxonSampler>{6.67, 0.44}};
  else
    throw std::invalid_argument{"unknown projectile species: " + species};
}

template <std::size_t A, typename NucleonSampler>
template <typename... SamplerArgs>
Nucleus<A, NucleonSampler>::Nucleus(SamplerArgs&&... sampler_args)
    : nucleon_sampler_(std::forward<SamplerArgs>(sampler_args)...) {}

template <std::size_t A, typename NucleonSampler>
void Nucleus<A, NucleonSampler>::sample_nucleons(double offset) {
  for (auto& nucleon : nucleons_) {
    auto coord = nucleon_sampler_.sample();
    std::get<0>(coord) -= offset;
    nucleon = {std::move(coord), false};
  }
}

}  // namespace trento
