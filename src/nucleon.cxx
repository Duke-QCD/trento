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
#include <boost/math/special_functions/erf.hpp>
#include <boost/filesystem.hpp>


#include "fwd_decl.h"

namespace trento {

namespace {

// These constants define distances in terms of the width of the nucleon profile
// Gaussian thickness function.

// Truncation radius of the thickness function.
constexpr double max_radius_widths = 5.;

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
// Used to store cross section parameter \sigma_partonic.
fs::path get_data_home() {
  const auto data_path = std::getenv("XDG_DATA_HOME");
  if(data_path == nullptr)
    return fs::path{std::getenv("HOME")} / ".local/share";
  return data_path;
}

// Test approximate cross section parameter equality.
bool almost_equal(double a, double b) {
  return fabs(a - b) < 1e-6;
}

// Calculate constituent position sampling width from CL arguments
double calc_sampling_width(const VarMap& var_map) {
  auto nucleon_width = var_map["nucleon-width"].as<double>();
  auto constituent_width = var_map["constit-width"].as<double>();
  auto constituent_number = var_map["constit-number"].as<int>();

  auto one_constituent = (constituent_number == 1);
  auto same_size = (nucleon_width - constituent_width < 1e-6);

  if (one_constituent || same_size)
    return 0.0;
  else {
    auto num = sqr(nucleon_width) - sqr(constituent_width);
    auto denom = 1. - 1./constituent_number;
    return std::sqrt(num / denom);
  }
}

// Determine cross section parameter for sampling nucleon participants.
// This semi-analytic method is only valid for one constituent.
double analytic_partonic_cross_section(const VarMap& var_map) {
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

    auto cross_section_param = .5*(result.first + result.second);
    return std::exp(cross_section_param) * 4 * math::double_constants::pi * sqr(width);
  }
  catch (const std::domain_error&) {
    // Root finding fails for very small nucleon widths, w^2/sigma_nn < ~0.01.
    throw std::domain_error{
      "unable to fit cross section -- nucleon width too small?"};
  }
}

// Determine cross section parameter for sampling nucleon participants.
// This Monte Carlo numeric method is used for more than one constituent.
double numeric_partonic_cross_section(const VarMap& var_map) {
  MonteCarloCrossSection mc_cross_section(var_map);
  auto sigma_nn = var_map["cross-section"].as<double>();
  auto width = var_map["constit-width"].as<double>();

  // Bracket min and max.
  auto a = -10.;
  auto b = 20.;

  // Tolerance function.
  math::tools::eps_tolerance<double> tol{8};

  // Maximum iterations.
  boost::uintmax_t max_iter = 20;

  // Dimensional constant [fm^-1] used to rescale sigma_partonic
  auto c = 4 * math::double_constants::pi * sqr(width);

  try {
    auto result = math::tools::toms748_solve(
      [&sigma_nn, &mc_cross_section, &c](double cross_section_param) {
        auto x = std::exp(cross_section_param) * c;
        return mc_cross_section(x) - sigma_nn;
      },
      a, b, tol, max_iter);

    auto cross_section_param = .5*(result.first + result.second);
    return std::exp(cross_section_param) * c;
  }
  catch (const std::domain_error&) {
    // Root finding fails for very small nucleon widths, w^2/sigma_nn < ~0.01.
    throw std::domain_error{
      "unable to fit cross section -- nucleon and/or constituent width too small?"};
  }
}

// Determine the cross section parameter given the nucleon_width,
// constituent_width, constituent number, and inelastic cross section.
double partonic_cross_section(const VarMap& var_map) {
  // Read parameters from the configuration.
  auto nucleon_width = var_map["nucleon-width"].as<double>();
  auto constituent_width = var_map["constit-width"].as<double>();
  auto constituent_number = var_map["constit-number"].as<int>();
  auto sigma_nn = var_map["cross-section"].as<double>();

  // Cross section parameters
  double nucleon_width_;
  double constituent_width_;
  int constituent_number_;
  double sigma_nn_;
  double sigma_partonic;

  // Establish trento cache path
  auto cache_dir = get_data_home() / "trento";
  auto cache_path = cache_dir / "cross_section.dat";
  fs::create_directory(cache_dir);

  // Check if cache exists
  if (fs::exists(cache_path.string())){

    // Open cache read only
    std::fstream cache_read(cache_path.string(), std::ios_base::in);

    // Check cache for previously tabulated cross section parameters
    while(!cache_read.eof()) {
      cache_read >> constituent_number_ >> nucleon_width_ >>
        constituent_width_ >> sigma_nn_ >> sigma_partonic;
      if (almost_equal(nucleon_width, nucleon_width_) &&
          almost_equal(constituent_width, constituent_width_) &&
          almost_equal(sigma_nn, sigma_nn_) &&
          (constituent_number == constituent_number_)) {
        return sigma_partonic;
      }
    }

    // Close read-only file stream
    cache_read.close();
  }

  // Open cache with write access
  std::fstream cache_write(cache_path.string(),
      std::ios_base::out | std::ios_base::app);

  // Use numeric method to determine sigma_partonic if there is nucleon
  // substructure, otherwise use semi-analytic method.
  if (constituent_number > 1)
    sigma_partonic = numeric_partonic_cross_section(var_map);
  else
    sigma_partonic = analytic_partonic_cross_section(var_map);

  // Save new cross section parameters to cache
  using std::fixed;
  using std::setprecision;
  using std::setw;
  using std::scientific;
  using std::endl;

  cache_write << setprecision(6)
              << std::left << setw(4) << fixed << constituent_number
              << std::right << setw(10) << fixed << nucleon_width
              << std::right << setw(10) << fixed << constituent_width
              << std::right << setw(10) << fixed << sigma_nn
              << std::right << setw(14) << scientific << sigma_partonic
              << endl;

  return sigma_partonic;
}

}  // unnamed namespace

NucleonCommon::NucleonCommon(const VarMap& var_map)
    : fast_exp_(-.5*sqr(max_radius_widths), 0., 1000),
      nucleon_width_(var_map["nucleon-width"].as<double>()),
      constituent_width_(var_map["constit-width"].as<double>()),
      constituent_number_(std::size_t(var_map["constit-number"].as<int>())),
      sampling_width_(calc_sampling_width(var_map)),
      max_impact_sq_(sqr(max_impact_widths*nucleon_width_)),
      constituent_width_sq_(sqr(constituent_width_)),
      constituent_radius_sq_(sqr(max_radius_widths*constituent_width_)),
      sigma_partonic_(partonic_cross_section(var_map)),
      prefactor_(math::double_constants::one_div_two_pi/constituent_width_sq_/constituent_number_),
      calc_ncoll_(var_map["ncoll"].as<bool>()),
      participant_fluctuation_dist_(gamma_param_unit_mean(var_map["fluctuation"].as<double>())),
      constituent_position_dist_(0, sampling_width_)
{}

MonteCarloCrossSection::MonteCarloCrossSection(const VarMap& var_map)
  : nucleon_width_(var_map["nucleon-width"].as<double>()),
    constituent_width_(var_map["constit-width"].as<double>()),
    constituent_number_(std::size_t(var_map["constit-number"].as<int>())),
    sampling_width_(calc_sampling_width(var_map)),
    max_impact_(max_impact_widths*nucleon_width_),
    constituent_width_sq_(sqr(constituent_width_)),
    prefactor_(math::double_constants::one_div_two_pi/constituent_width_sq_/constituent_number_)
{}

double MonteCarloCrossSection::operator() (const double sigma_partonic) const {

  random::CyclicNormal<> cyclic_normal{
    0., sampling_width_, cache_size, n_loops
  };

  struct Constituent {
    double x, y;
  };

  std::vector<Constituent> nucleonA(constituent_number_);
  std::vector<Constituent> nucleonB(constituent_number_);

  auto max_impact_sq_ = sqr(max_impact_);
  double ref_cross_section = 0.;
  double prob_miss = 0.;
  int pass_tolerance = 0;

  for (std::size_t n = 0; n < n_max; ++n) {
    // Sample b from P(b)db = 2*pi*b.
    auto b = max_impact_ * std::sqrt(random::canonical<double>());

    for (auto&& q : nucleonA) {
      q.x = cyclic_normal(random::engine);
      q.y = cyclic_normal(random::engine);
    }

    for (auto&& q : nucleonB) {
      q.x = cyclic_normal(random::engine);
      q.y = cyclic_normal(random::engine);
    }

    auto overlap = 0.;
    for (auto&& qA : nucleonA) {
      for (auto&& qB : nucleonB) {
        auto distance_sq = sqr(qA.x - qB.x + b) + sqr(qA.y - qB.y);
        overlap += std::exp(-.25*distance_sq/constituent_width_sq_);
      }
    }

    prob_miss +=
      std::exp(-sigma_partonic * prefactor_/(2.*constituent_number_) * overlap);

    auto prob_hit = 1. - (prob_miss/n);
    auto cross_section = M_PI*max_impact_sq_*prob_hit;

    auto update_difference = std::abs(cross_section - ref_cross_section);

    if (update_difference < tolerance) {
      ++pass_tolerance;
    } else {
      pass_tolerance = 0;
      ref_cross_section = cross_section;
    }

    if (pass_tolerance > n_pass){
      return cross_section;
    }
  }

  throw std::out_of_range{
    "Partonic cross section failed to converge \
      -- check nucleon width, constituent width, and constituent number"};
}

}  // namespace trento
