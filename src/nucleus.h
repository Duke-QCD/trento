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
  class NucleonData {
    friend class Nucleus;

   public:
    ///
    NucleonData() = default;

    /// Disable all copy and move operations.
    NucleonData(const NucleonData&) = delete;
    NucleonData& operator=(const NucleonData&) = delete;
    NucleonData(NucleonData&&) = delete;
    NucleonData& operator=(NucleonData&&) = delete;

    ///
    double x() const
    { return x_; }

    ///
    double y() const
    { return y_; }

    ///
    bool is_participant() const
    { return participant_; }

    ///
    void set_participant()
    { participant_ = true; }

   private:
    ///
    void set_position(double x, double y);

    ///
    double x_, y_;

    ///
    bool participant_ = false;
  };

  ///
  using iterator = std::vector<NucleonData>::iterator;
  using const_iterator = std::vector<NucleonData>::const_iterator;

  ///
  iterator begin() noexcept
  { return nucleon_data_vector_.begin(); }
  iterator end() noexcept
  { return nucleon_data_vector_.end(); }

  ///
  const_iterator begin() const noexcept
  { return nucleon_data_vector_.begin(); }
  const_iterator end() const noexcept
  { return nucleon_data_vector_.end(); }

  ///
  const_iterator cbegin() const noexcept
  { return nucleon_data_vector_.cbegin(); }
  const_iterator cend() const noexcept
  { return nucleon_data_vector_.cend(); }

 protected:
  ///
  explicit Nucleus(std::size_t A);

  ///
  template <typename... Args>
  void set_nucleon_position(NucleonData& nucleon_data, Args&&... args);

 private:
  std::vector<NucleonData> nucleon_data_vector_;
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
