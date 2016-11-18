// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "nucleon.h"

#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/expint.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/filesystem.hpp>


#include "fwd_decl.h"

namespace trento {

namespace {

// These constants define distances in terms of the width of the nucleon profile
// Gaussian thickness function.

// Truncation radius of the thickness function.
constexpr double radius_widths = 5.;

// Maximum impact parameter for participation.
constexpr double max_impact_widths = 6.;
// Create ctor parameters for unit mean std::gamma_distribution.
//   mean = alpha*beta == 1  ->  beta = 1/alpha
// Used below in NucleonProfile ctor initializer list.

template <typename RealType> using param_type =
  typename std::gamma_distribution<RealType>::param_type;

template <typename RealType>
param_type<RealType> gamma_param_unit_mean(RealType alpha = 1.) {
  return param_type<RealType>{alpha, 1./alpha};
}

// Return appropriate path to temporary file directory.
// Used to store cross section parameter \sigma_qq.
fs::path get_data_home() {
  const auto data_path = std::getenv("xdg_data_home");
  if(data_path == nullptr)
    return fs::path{std::getenv("HOME")} / ".local/share";
  return data_path;
}

// Test approximate cross section parameter equality
bool almost_equal(double a, double b) {
  return fabs(a - b) < 1e-6;
}

// Determine the effective parton-parton cross section for sampling participants.
// TODO: derive this
double analytic_sigma_qq(const VarMap& var_map) {
  // Read parameters from the configuration.
  auto sigma_nn = var_map["cross-section"].as<double>();
  auto width = var_map["nucleon-width"].as<double>();

  // TODO: automatically set from beam energy
  // if (sigma_nn < 0.) {
  //   auto sqrt_s = var_map["beam-energy"].as<double>();
  // }

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

  // The right-hand side of the equation.
  auto rhs = sigma_nn / (4 * math::double_constants::pi * sqr(width));

  // This quantity appears a couple times in the equation.
  auto c = sqr(max_impact_widths) / 4;

  try {
    auto result = math::tools::toms748_solve(
      [&rhs, &c](double x) {
        using std::exp;
        using math::expint;
        return c - expint(-exp(x)) + expint(-exp(x-c)) - rhs;
      },
      a, b, tol, max_iter);

    return .5*(result.first + result.second);
  }
  catch (const std::domain_error&) {
    // Root finding fails for very small nucleon widths, w^2/sigma_nn < ~0.01.
    throw std::domain_error{
      "unable to fit cross section -- nucleon width too small?"};
  }
}

//double monte_carlo_sigma_qq(const VarMap& var_map) {
//  // Read parameters from the configuration.
//  auto sigma_nn = var_map["cross-section"].as<double>();
//  auto nucleon_width = var_map["nucleon-width"].as<double>();
//  auto parton_width = var_map["parton-width"].as<double>();
//  auto nparton = var_map["parton-number"].as<int>();
//
//  return 1.;
//}

// Determine cross section parameter (partonic cross section), given
// nucleon width w, parton width v, cross section x, and parton number m.
double compute_sigma_qq(const VarMap& var_map) {
  // Read parameters from the configuration.
  auto sigma_nn = var_map["cross-section"].as<double>();
  auto nucleon_width = var_map["nucleon-width"].as<double>();
  auto parton_width = var_map["parton-width"].as<double>();
  auto nparton = var_map["parton-number"].as<int>();

  // create trento cache directory
  auto cache_dir = get_data_home() / "trento";
  fs::create_directory(cache_dir);

  // create dummy cross section file
  auto cache_path = cache_dir / "cross_section.dat";
  std::fstream cache_file(cache_path.string(),
      std::ios_base::out | std::ios_base::in | std::ios_base::app);

  double sigma_nn_;
  double nucleon_width_;
  double parton_width_;
  int nparton_;
  double xsec_param;

  while(!cache_file.eof()) {
    cache_file >> nparton_ >> nucleon_width_ >>
      parton_width_ >> sigma_nn_ >> xsec_param;
    if (almost_equal(nucleon_width, nucleon_width_) &&
        almost_equal(parton_width, parton_width_) &&
        almost_equal(sigma_nn, sigma_nn_) &&
        (nparton == nparton_)) {
      return xsec_param;
    }
  }

  xsec_param = analytic_sigma_qq(var_map);

  using std::fixed;
  using std::setprecision;
  using std::setw;
  using std::scientific;

  cache_file.clear();
  cache_file << setprecision(6) << std::left
             << setw(8) << fixed << nparton
             << setw(10) << fixed << nucleon_width
             << setw(10) << fixed << parton_width 
             << setw(10) << fixed << sigma_nn
             << setw(14) << scientific << xsec_param
             << std::endl;

  return xsec_param;
}

// Determine the width of the normal distribution for sampling parton positions.
double compute_parton_sampling_width(const VarMap& var_map) {
  // Read parameters from the configuration.
  auto nucleon_width = var_map["nucleon-width"].as<double>();
  auto parton_width = var_map["parton-width"].as<double>();

  // Compute Gaussian sampling width by deconvolving the parton profile from the
  // nucleon profile.  The convolution of two Gaussians is a third Gaussian
  // with the variances simply added together.
  auto sampling_width_sq = sqr(nucleon_width) - sqr(parton_width);

  if (sampling_width_sq < 0.)
    throw std::out_of_range{
      "the parton width cannot be larger than the nucleon width"};

  return std::sqrt(sampling_width_sq);
}

}  // unnamed namespace

NucleonCommon::NucleonCommon(const VarMap& var_map)
    : fast_exp_(-.5*sqr(radius_widths), 0., 1000),
      max_impact_sq_(sqr(max_impact_widths*var_map["nucleon-width"].as<double>())),
      parton_width_sq_(sqr(var_map["parton-width"].as<double>())),
      parton_radius_sq_(sqr(radius_widths)*parton_width_sq_),
      npartons_(var_map["parton-number"].as<int>()),
      sigma_qq_(compute_sigma_qq(var_map)),
      prefactor_(math::double_constants::one_div_two_pi/parton_width_sq_/npartons_),
      nucleon_fluctuation_dist_(gamma_param_unit_mean(var_map["fluctuation"].as<double>())),
      parton_position_dist_(0, compute_parton_sampling_width(var_map))
{}

}  // namespace trento
