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

// Alias for a smart pointer to a Nucleus.
using NucleusPtr = std::unique_ptr<Nucleus>;

/// \rst
/// Interface class to all nucleus types.  Stores an ensemble of nucleons
/// and randomly samples their transverse positions.  Implements a standard
/// iterator interface through ``begin()`` and ``end()`` functions.  Iterating over
/// a ``Nucleus`` means iterating over its ``Nucleon`` members.
/// \endrst
class Nucleus {
 public:
  /// \rst
  /// The canonical way to create a ``Nucleus``.  Constructs the appropriate
  /// subclass for the given species.  Sets the Woods-Saxon parameters for Au,
  /// Pb, U, etc; parameters copied from the PHOBOS Glauber model:
  ///
  /// - http://inspirehep.net/record/786828
  /// - http://inspirehep.net/record/1310629
  ///
  /// Example::
  ///
  ///   std::unique_ptr<Nucleus> lead_nucleus = Nucleus::create("Pb");
  ///   double radius = lead_nucleus->radius();
  ///   lead_nucleus->sample_nucleons(0);
  ///   for (const auto& nucleon : *lead_nucleus)
  ///     do_something(nucleon)
  ///
  /// \endrst
  ///
  /// \param species standard symbol, e.g. "p" for proton or "Pb" for lead-208
  ///
  /// \return a smart pointer \c std::unique_ptr<Nucleus>
  ///
  /// \throw std::invalid_argument for unknown species
  static NucleusPtr create(const std::string& species);

  /// Default virtual destructor for abstract base class.
  virtual ~Nucleus() = default;

  /// The "radius", i.e. the maximum distance at which a nucleon could be
  /// placed.
  virtual double radius() const = 0;

  /// Sample a new ensemble of nucleon positions.
  /// \param offset shift for each \em x position
  virtual void sample_nucleons(double offset) = 0;

  using iterator = std::vector<Nucleon>::iterator;
  using const_iterator = std::vector<Nucleon>::const_iterator;

  // non-const overload
  iterator begin() noexcept
  { return nucleons_.begin(); }
  iterator end() noexcept
  { return nucleons_.end(); }

  // const overload
  const_iterator begin() const noexcept
  { return nucleons_.begin(); }
  const_iterator end() const noexcept
  { return nucleons_.end(); }

  // forced const
  const_iterator cbegin() const noexcept
  { return nucleons_.cbegin(); }
  const_iterator cend() const noexcept
  { return nucleons_.cend(); }

 protected:
  /// Constructor only accessible by derived classes.
  /// \param A number of nucleons
  explicit Nucleus(std::size_t A);

  /// \rst
  /// Set a ``Nucleon`` position.  Perfect-forwards arguments to
  /// ``Nucleon::set_position``.  Implemented this way because while ``Nucleus``
  /// is a friend of ``Nucleon``, its derived classes are not.  This grants
  /// limited access to the derived classes.
  /// \endrst
  template <typename... Args>
  void set_nucleon_position(Nucleon& nucleon, Args&&... args);

 private:
  /// Internal storage of Nucleon objects.
  std::vector<Nucleon> nucleons_;
};

/// Trivial nucleus with a single nucleon.
class Proton : public Nucleus {
 public:
  /// Default constructor.
  Proton();

  virtual double radius() const override;

  virtual void sample_nucleons(double offset) override;
};

/// \rst
/// Samples nucleons from a spherically symmetric Woods-Saxon distribution
///
/// .. math::
///
///   f(r) \propto \frac{1}{1 + \exp(\frac{r-R}{a})}.
///
/// For non-deformed heavy nuclei such as lead.
///
/// \endrst
class WoodsSaxonNucleus : public Nucleus {
 public:
  /// ``Nucleus::create()`` sets these parameters for a given species.
  /// \param A number of nucleons
  /// \param R Woods-Saxon radius
  /// \param a Woods-Saxon surface thickness
  WoodsSaxonNucleus(std::size_t A, double R, double a);

  virtual double radius() const override;

  virtual void sample_nucleons(double offset) override;

 private:
  /// W-S parameters.
  const double R_, a_;

  /// W-S distribution object.  Since the dist does not have an analytic inverse
  /// CDF, approximate it as a piecewise linear dist.  For a large number of
  /// steps this is very accurate.
  mutable std::piecewise_linear_distribution<double> woods_saxon_dist_;
};

/// \rst
/// Samples nucleons from a deformed spheroidal Woods-Saxon distribution
///
/// .. math::
///
///   f(r, \theta) \propto
///   \frac{1}{1 + \exp(\frac{r-R(1+\beta_2Y_{20}+\beta_4Y_{40})}{a})}.
///
/// For deformed heavy nuclei such as uranium.
///
/// \endrst
class DeformedWoodsSaxonNucleus : public Nucleus {
 public:
  /// ``Nucleus::create()`` sets these parameters for a given species.
  /// \param A number of nucleons
  /// \param R Woods-Saxon radius
  /// \param a Woods-Saxon surface thickness
  /// \param beta2 Woods-Saxon deformation parameter
  /// \param beta4 Woods-Saxon deformation parameter
  DeformedWoodsSaxonNucleus(std::size_t A, double R, double a,
                            double beta2, double beta4);

  virtual double radius() const override;

  virtual void sample_nucleons(double offset) override;

 private:
  /// Evaluate the deformed Woods-Saxon distribution.
  double deformed_woods_saxon_dist(double r, double cos_theta) const;

  /// W-S parameters.
  const double R_, a_, beta2_, beta4_;

  /// Maximum radius.
  const double rmax_;
};

}  // namespace trento

#endif  // NUCLEUS_H
