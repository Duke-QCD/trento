// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "event.h"

#include <algorithm>
#include <cmath>

#include <boost/program_options/variables_map.hpp>

#include "nucleus.h"

namespace trento {

namespace {

constexpr double TINY = 1e-12;

// Generalized mean for p > 0.
// M_p(a, b) = (1/2*(a^p + b^p))^(1/p)
inline double positive_pmean(double p, double a, double b) {
  return std::pow(.5*(std::pow(a, p) + std::pow(b, p)), 1./p);
}

// Generalized mean for p < 0.
// Same as the positive version, except prevents division by zero.
inline double negative_pmean(double p, double a, double b) {
  if (a < TINY || b < TINY)
    return 0.;
  return positive_pmean(p, a, b);
}

// Generalized mean for p == 0.
inline double geometric_mean(double a, double b) {
  return std::sqrt(a*b);
}

}  // unnamed namespace

Event::Event(const VarMap& var_map)
    : norm_(var_map["normalization"].as<double>()),
      nsteps_(var_map["grid-steps"].as<int>()),
      width_(var_map["grid-width"].as<double>()),
      dxy_(width_/nsteps_),
      TA_(boost::extents[nsteps_][nsteps_]),
      TB_(boost::extents[nsteps_][nsteps_]),
      TR_(boost::extents[nsteps_][nsteps_]) {
  // Set reduced thickness function based on configuraton.
  auto p = var_map["reduced-thickness"].as<double>();

  // TODO: explain
  if (std::fabs(p) < TINY) {
    compute_reduced_thickness_ = [this]() {
      compute_reduced_thickness(geometric_mean);
    };
  } else if (p > 0.) {
    compute_reduced_thickness_ = [this, p]() {
      compute_reduced_thickness(
        [p](double a, double b) { return positive_pmean(p, a, b); });
    };
  } else {
    compute_reduced_thickness_ = [this, p]() {
      compute_reduced_thickness(
        [p](double a, double b) { return negative_pmean(p, a, b); });
    };
  }
}

void Event::compute(const Nucleus& nucleusA, const Nucleus& nucleusB,
                    const NucleonProfile& profile) {
  npart_ = 0;
  compute_nuclear_thickness(nucleusA, profile, TA_);
  compute_nuclear_thickness(nucleusB, profile, TB_);
  compute_reduced_thickness_();
  compute_observables();
}

namespace {

// Limit a value to a range.
// Used below to constrain grid indices.
template <typename T>
inline const T& clip(const T& value, const T& min, const T& max) {
  if (value < min)
    return min;
  if (value > max)
    return max;
  return value;
}

}  // unnamed namespace

void Event::compute_nuclear_thickness(
    const Nucleus& nucleus, const NucleonProfile& profile, Grid& TX) {
  // Wipe grid with zeros.
  std::fill(TX.origin(), TX.origin() + TX.num_elements(), 0.);

  const double r = profile.radius();

  // Deposit each participant onto the grid.
  for (const auto& nucleon : nucleus) {
    if (!nucleon.is_participant())
      continue;

    ++npart_;

    // Work in coordinates relative to (-width/2, -width/2).
    double x = nucleon.x() + .5*width_;
    double y = nucleon.y() + .5*width_;

    // Determine min & max indices of nucleon subgrid.
    int ixmin = clip(static_cast<int>((x-r)/dxy_), 0, nsteps_-1);
    int iymin = clip(static_cast<int>((y-r)/dxy_), 0, nsteps_-1);
    int ixmax = clip(static_cast<int>((x+r)/dxy_), 0, nsteps_-1);
    int iymax = clip(static_cast<int>((y+r)/dxy_), 0, nsteps_-1);

    double fluct = profile.sample_fluct();

    // Add profile to grid.
    for (auto iy = iymin; iy <= iymax; ++iy) {
      double dysq = std::pow(y - (static_cast<double>(iy)+.5)*dxy_, 2);
      for (auto ix = ixmin; ix <= ixmax; ++ix) {
        double dxsq = std::pow(x - (static_cast<double>(ix)+.5)*dxy_, 2);
        TX[iy][ix] += fluct * profile.thickness(dxsq + dysq);
      }
    }
  }
}

template <typename GenMean>
void Event::compute_reduced_thickness(GenMean gen_mean) {
  double sum = 0.;
  double xcm = 0.;
  double ycm = 0.;

  for (int iy = 0; iy < nsteps_; ++iy) {
    for (int ix = 0; ix < nsteps_; ++ix) {
      auto t = norm_ * gen_mean(TA_[iy][ix], TB_[iy][ix]);
      TR_[iy][ix] = t;
      sum += t;
      xcm += t * static_cast<double>(ix);
      ycm += t * static_cast<double>(iy);
    }
  }

  multiplicity_ = dxy_ * dxy_ * sum;
  xcm_ = xcm / sum;
  ycm_ = ycm / sum;
}

void Event::compute_observables() {
  // Compute eccentricity.

  // Simple helper class for use in the following loop.
  struct EccentricityAccumulator {
    double re = 0.;  // real part
    double im = 0.;  // imaginary part
    double wt = 0.;  // weight
    double finish() const  // compute final eccentricity
    { return std::sqrt(re*re + im*im) / std::fmax(wt, TINY); }
  } e2, e3, e4, e5;

  for (int iy = 0; iy < nsteps_; ++iy) {
    for (int ix = 0; ix < nsteps_; ++ix) {
      const auto& t = TR_[iy][ix];
      if (t < TINY)
        continue;

      auto x = static_cast<double>(ix) - xcm_;
      auto x2 = x*x;
      auto x3 = x2*x;
      auto x4 = x2*x2;

      auto y = static_cast<double>(iy) - ycm_;
      auto y2 = y*y;
      auto y3 = y2*y;
      auto y4 = y2*y2;

      auto r2 = x2 + y2;
      auto r = std::sqrt(r2);
      auto r4 = r2*r2;

      auto xy = x*y;
      auto x2y2 = x2*y2;

      // TODO: explain what is happening here
      e2.re += t * (y2 - x2);
      e2.im += t * 2.*xy;
      e2.wt += t * r2;

      e3.re += t * (y3 - 3.*y*x2);
      e3.im += t * (3.*x*y2 - x3);
      e3.wt += t * r2*r;

      e4.re += t * (x4 + y4 - 6.*x2y2);
      e4.im += t * 4.*xy*(y2 - x2);
      e4.wt += t * r4;

      e5.re += t * y*(5.*x4 - 10.*x2y2 + y4);
      e5.im += t * x*(x4 - 10.*x2y2 + 5.*y4);
      e5.wt += t * r4*r;
    }
  }

  eccentricity_[2] = e2.finish();
  eccentricity_[3] = e3.finish();
  eccentricity_[4] = e4.finish();
  eccentricity_[5] = e5.finish();
}

}  // namespace trento
