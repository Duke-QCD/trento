// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef NUCLEUS_H
#define NUCLEUS_H

#include <array>

#include "nucleus_base.h"

namespace trento {

using TransverseCoord = std::array<double, 2>;

///
template <std::size_t A, typename NucleonSampler>
class Nucleus : public NucleusBase {
  ///
  friend class NucleusBase;

 private:
  /// Private constructor, only accessible by base class.
  /// Perfect-forwards arguments to the NucleonSampler policy.
  template <typename... SamplerArgs>
  Nucleus(SamplerArgs&&... sampler_args);

  ///
  virtual void sample_nucleons(double offset) override;

  ///
  std::array<TransverseCoord, A> nucleons_;

  ///
  NucleonSampler nucleon_sampler_;
};

}  // namespace trento

#endif  // NUCLEUS_H
