// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef NUCLEUS_H
#define NUCLEUS_H

#include <memory>
#include <random>
#include <string>
#include <vector>

#include "fwd_decl.h"
#include "nucleon.h"

namespace trento {

///
using NucleusPtr = std::unique_ptr<Nucleus>;

///
class Nucleus {
 public:
  ///
  static NucleusPtr create(const std::string& species);

  ///
  virtual ~Nucleus() = default;

  ///
  virtual double radius() const = 0;

  ///
  virtual void sample_nucleons(double offset) = 0;

  ///
  using iterator = std::vector<Nucleon>::iterator;
  using const_iterator = std::vector<Nucleon>::const_iterator;

  ///
  iterator begin() noexcept
  { return nucleons_.begin(); }
  iterator end() noexcept
  { return nucleons_.end(); }

  ///
  const_iterator begin() const noexcept
  { return nucleons_.begin(); }
  const_iterator end() const noexcept
  { return nucleons_.end(); }

  ///
  const_iterator cbegin() const noexcept
  { return nucleons_.cbegin(); }
  const_iterator cend() const noexcept
  { return nucleons_.cend(); }

 protected:
  ///
  explicit Nucleus(std::size_t A);

  ///
  template <typename... Args>
  void set_nucleon_position(Nucleon& nucleon, Args&&... args);

 private:
  std::vector<Nucleon> nucleons_;
};

///
class Proton : public Nucleus {
 public:
  ///
  Proton();

  ///
  virtual double radius() const override;

  ///
  virtual void sample_nucleons(double offset) override;
};

/// Samples nucleons from a spherically symmetric Woods-Saxon distribution.
/// For non-deformed heavy nuclei such as lead.
class WoodsSaxonNucleus : public Nucleus {
 public:
  ///
  WoodsSaxonNucleus(std::size_t A, double R, double a);

  ///
  virtual double radius() const override;

  ///
  virtual void sample_nucleons(double offset) override;

 private:
  ///
  const double R_, a_;

  ///
  mutable std::piecewise_linear_distribution<double> woods_saxon_dist_;
};

}  // namespace trento

#endif  // NUCLEUS_H
