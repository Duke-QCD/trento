// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef NUCLEUS_H
#define NUCLEUS_H

#include <array>
#include <memory>
#include <string>
#include <utility>

#include "fwd_decl.h"

namespace trento {

///
using NucleusPtr = std::unique_ptr<NucleusBase>;

///
class NucleusBase {
 public:
  ///
  static NucleusPtr create(const std::string& species);

  ///
  virtual void sample_nucleons(double offset) = 0;

  ///
  virtual double radius() const = 0;

  using TransverseCoord = std::array<double, 2>;
  using NucleonAttr = std::pair<TransverseCoord, bool>;

  template <std::size_t A>
  using NucleonArray = std::array<NucleonAttr, A>;

  using iterator = NucleonArray<1>::iterator;
  using const_iterator = NucleonArray<1>::const_iterator;

  virtual iterator begin() noexcept = 0;
  virtual const_iterator cbegin() const noexcept = 0;
  virtual iterator end() noexcept = 0;
  virtual const_iterator cend() const noexcept = 0;

  ///
  virtual ~NucleusBase() = default;

 protected:
  ///
  NucleusBase() = default;
};

///
template <std::size_t A, typename NucleonSampler>
class Nucleus : public NucleusBase {
 public:
  ///
  virtual void sample_nucleons(double offset) override;

  virtual double radius() const override
  { return nucleon_sampler_.radius(); }

  virtual iterator begin() noexcept override
  { return nucleons_.begin(); }

  virtual const_iterator cbegin() const noexcept override
  { return nucleons_.cbegin(); }

  virtual iterator end() noexcept override
  { return nucleons_.end(); }

  virtual const_iterator cend() const noexcept override
  { return nucleons_.cend(); }

 private:
  ///
  friend NucleusPtr NucleusBase::create(const std::string&);

  /// Private constructor, only accessible by \c NucleusBase::create.
  /// Perfect-forwards arguments to the NucleonSampler policy.
  template <typename... SamplerArgs>
  Nucleus(SamplerArgs&&... sampler_args);

  ///
  NucleonArray<A> nucleons_;

  ///
  NucleonSampler nucleon_sampler_;
};

}  // namespace trento

#endif  // NUCLEUS_H
