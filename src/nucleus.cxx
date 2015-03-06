// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "nucleus.h"

#include <array>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include "nucleon.h"

namespace trento {

namespace {

using TransverseCoord = std::array<double, 2>;

// Nucleon sampling policies.

/// Trivial sampler that always returns the same position.
/// For protons.
class FixedSampler {
 public:
  ///
  FixedSampler() : x_(0.), y_(0.) {}

  ///
  TransverseCoord sample() const {
    return {x_, y_};
  }

 private:
  ///
  const double x_, y_;
};

/// Samples nucleons from a spherically symmetric Woods-Saxon distribution.
/// For non-deformed heavy nuclei such as lead.
class WoodsSaxonSampler {
 public:
  ///
  WoodsSaxonSampler(double R, double a) : R_(R), a_(a) {}

  ///
  TransverseCoord sample() const {
    return {.5*R_, .5*R_};
  }

 private:
  ///
  double R_, a_;
};

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
    nucleon.set_position(std::get<0>(coord), std::get<1>(coord));
  }
}

// NucleonSampler policies

// FixedSampler::FixedSampler() : x_(0.), y_(0.) {}

// WoodsSaxonSampler::WoodsSaxonSampler(double R, double a) : R_(R), a_(a) {}

}  // namespace trento
