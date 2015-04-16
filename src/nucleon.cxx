// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "nucleon.h"

#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/program_options/variables_map.hpp>

#include "fwd_decl.h"

namespace trento {

namespace {

// Create ctor parameters for unit mean std::gamma_distribution.
//   mean = alpha*beta == 1  ->  beta = 1/alpha
// Used below in NucleonProfile ctor initializer list.

template <typename RealType> using param_type =
  typename std::gamma_distribution<RealType>::param_type;

template <typename RealType>
param_type<RealType> gamma_param_unit_mean(RealType alpha = 1.) {
  return param_type<RealType>{alpha, 1./alpha};
}

}  // unnamed namespace

/// \b TODO Derive the cross section parameter.
NucleonProfile::NucleonProfile(const VarMap& var_map)
    : width_squared_(std::pow(var_map["nucleon-width"].as<double>(), 2)),
      trunc_radius_squared_(std::pow(trunc_widths_, 2) * width_squared_),
      fast_exp_(-.5*trunc_widths_*trunc_widths_, 0., 1000),
      fluct_dist_(gamma_param_unit_mean(var_map["fluctuation"].as<double>())) {
  // Determine cross section.
  auto sigma_nn = var_map["cross-section"].as<double>();
  if (sigma_nn < 0.) {
    // TODO: automatically set from beam energy
    // auto sqrt_s = var_map["beam-energy"].as<double>();
    sigma_nn = 6.4;
  }

  // Initialize arguments for boost root finding function.

  // Bracket min and max.
  auto a = -10.;
  auto b = 20.;

  // Tolerance function.
  // Require 3/4 of double precision.
  math::tools::eps_tolerance<double> tol{
    (std::numeric_limits<double>::digits * 3) / 4};

  // Maximum iterations.
  // This is overkill -- in testing only 10-20 iterations were required
  // (but no harm in overestimating).
  boost::uintmax_t max_iter = 1000;

  // Cache some quantities.
  auto snn_div_four_pi_wsq = .5 * math::constants::one_div_two_pi<double>() *
                             sigma_nn / width_squared_;
  auto t = trunc_widths_ * trunc_widths_;  // for lack of a better name

  try {
    auto result = math::tools::toms748_solve(
      [&snn_div_four_pi_wsq, &t](double x) {
        using std::exp;
        using math::expint;
        return t - expint(-exp(x)) + expint(-exp(x-t)) - snn_div_four_pi_wsq;
      },
      a, b, tol, max_iter);

    cross_sec_param_ = .5*(result.first + result.second);
  }
  catch (const std::domain_error&) {
    // Root finding fails for very small nucleon widths, w^2/sigma_nn < ~0.01.
    throw std::domain_error{
      "unable to fit cross section -- nucleon width too small?"};
  }
}

}  // namespace trento
