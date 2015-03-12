// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef NUCLEUS_H
#define NUCLEUS_H

#include <array>
#include <memory>
#include <string>

#include "fwd_decl.h"

namespace trento {

///
class NucleonData {
 public:
  ///
  NucleonData() = default;

  /// Disable all copy and move operations.
  NucleonData(const NucleonData&) = delete;
  NucleonData& operator=(const NucleonData&) = delete;
  NucleonData(NucleonData&&) = delete;
  NucleonData& operator=(NucleonData&&) = delete;

  ///
  const std::array<double, 2>& position() const { return position_; }

  ///
  void set_position(const std::array<double, 2>& new_position) {
    position_ = new_position;
    participant_ = false;
  }

  ///
  bool is_participant() const { return participant_; }

  ///
  void set_participant() { participant_ = true; }

 private:
  ///
  std::array<double, 2> position_;

  ///
  bool participant_ = false;
};

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

  ///
  template <std::size_t A>
  using NucleonDataArray = std::array<NucleonData, A>;

  ///
  using iterator = NucleonDataArray<1>::iterator;
  using const_iterator = NucleonDataArray<1>::const_iterator;

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
  { return nucleon_data_array_.begin(); }

  virtual const_iterator cbegin() const noexcept override
  { return nucleon_data_array_.cbegin(); }

  virtual iterator end() noexcept override
  { return nucleon_data_array_.end(); }

  virtual const_iterator cend() const noexcept override
  { return nucleon_data_array_.cend(); }

 private:
  ///
  friend NucleusPtr NucleusBase::create(const std::string&);

  /// Private constructor, only accessible by \c NucleusBase::create.
  /// Perfect-forwards arguments to the NucleonSampler policy.
  template <typename... SamplerArgs>
  Nucleus(SamplerArgs&&... sampler_args);

  ///
  NucleonDataArray<A> nucleon_data_array_;

  ///
  NucleonSampler nucleon_sampler_;
};

}  // namespace trento

#endif  // NUCLEUS_H
