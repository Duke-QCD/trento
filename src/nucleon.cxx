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
// Used to store cross section parameter \sigma_partonic.
fs::path get_data_home() {
  const auto data_path = std::getenv("XDG_DATA_HOME");
  if(data_path == nullptr)
    return fs::path{std::getenv("HOME")} / ".local/share";
  return data_path;
}

// Test approximate cross section parameter equality
bool almost_equal(double a, double b) {
  return fabs(a - b) < 1e-6;
}

// Default parton width to nucleon width.
double infer_parton_width(const VarMap& var_map) {
  auto nucleon_width = var_map["nucleon-width"].as<double>();
  auto parton_width = var_map["parton-width"].as<double>();
  if (parton_width < 0)
    return nucleon_width;
  else
    return parton_width;
}

// Determine the effective parton-parton cross section for sampling participants.
// TODO: derive this
double analytic_partonic_cross_section(const VarMap& var_map) {
  // Read parameters from the configuration.

  // Use manual inelastic nucleon-nucleon cross section if specified.
  // Otherwise default to extrapolated cross section.
  auto sigma_nn = var_map["cross-section"].as<double>();
  if (sigma_nn < 0) {
    sigma_nn = cross_sec_from_energy(var_map["beam-energy"].as<double>());
  }
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
      "unable to fit cross section -- nucleon and/or parton width too small?"};
  }
}

// Use Monte Carlo methods to determine the partonic cross section
// \sigma_partonic. Solves MonteCarloCrossSection(sigma_partonic) = sigma_nn.
double numeric_partonic_cross_section(const VarMap& var_map) {
  MonteCarloCrossSection mc_cross_section(var_map);
  auto sigma_nn = var_map["cross-section"].as<double>();
  auto parton_width = infer_parton_width(var_map);

  // Bracket min and max.
  auto a = -10.;
  auto b = 20.;

  // Tolerance function.
  math::tools::eps_tolerance<double> tol{8};

  // Maximum iterations.
  boost::uintmax_t max_iter = 20;

  // Dimensional constant [fm^-1] used to rescale sigma_partonic
  auto c = 4 * math::double_constants::pi * sqr(parton_width);

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
      "unable to fit cross section -- nucleon and/or parton width too small?"};
  }
}

// Determine the cross section parameter sigma_partonic, given the
// nucleon width, parton width, nucleon-nucleon cross section, and parton number.
double partonic_cross_section(const VarMap& var_map) {
  // Read parameters from the configuration.
  auto parton_width = infer_parton_width(var_map);
  auto nucleon_width = var_map["nucleon-width"].as<double>();
  auto sigma_nn = var_map["cross-section"].as<double>();
  auto nparton = var_map["parton-number"].as<int>();

  // Create trento cache directory
  auto cache_dir = get_data_home() / "trento";
  fs::create_directory(cache_dir);

  // Create dummy cross section file
  auto cache_path = cache_dir / "cross_section.dat";
  std::fstream cache_file(cache_path.string(),
      std::ios_base::out | std::ios_base::in | std::ios_base::app);

  int nparton_;
  double nucleon_width_;
  double parton_width_;
  double sigma_nn_;
  double sigma_partonic;

  // Check parameter configuration cache
  while(!cache_file.eof()) {
    cache_file >> nparton_ >> nucleon_width_ >>
      parton_width_ >> sigma_nn_ >> sigma_partonic;
    if (almost_equal(nucleon_width, nucleon_width_) &&
        almost_equal(parton_width, parton_width_) &&
        almost_equal(sigma_nn, sigma_nn_) &&
        (nparton == nparton_)) {
      return sigma_partonic;
    }
  }

  // Use a semi-analytic method to determine sigma_partonic
  // if there is no subnucleonic structure---otherwise use MC method.
  if (nucleon_width == parton_width) 
    sigma_partonic = analytic_partonic_cross_section(var_map);
  else 
    sigma_partonic = numeric_partonic_cross_section(var_map);

  // Save new cross section parameters to cache
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
             << setw(14) << scientific << sigma_partonic
             << std::endl;

  return sigma_partonic;
}

// Determine the width of the normal distribution for sampling parton positions.
double compute_parton_sampling_width(const VarMap& var_map) {
  // Read parameters from the configuration.
  auto nucleon_width = var_map["nucleon-width"].as<double>();
  auto parton_width = var_map["parton-width"].as<double>();

  // Default parton width to nucleon width.
  if (parton_width < 0) return 0.;

  // Compute Gaussian sampling width by deconvolving the parton profile from the
  // nucleon profile.  The convolution of two Gaussians is a third Gaussian
  // with the variances simply added together.
  auto sampling_width_sq = sqr(nucleon_width) - sqr(parton_width);

  if (sampling_width_sq < 0.) {
    throw std::out_of_range{
      "the parton width cannot be larger than the nucleon width"};
  }

  return std::sqrt(sampling_width_sq);
}

}  // unnamed namespace

NucleonCommon::NucleonCommon(const VarMap& var_map)
    : fast_exp_(-.5*sqr(radius_widths), 0., 1000),
      max_impact_sq_(sqr(max_impact_widths*var_map["nucleon-width"].as<double>())),
      parton_width_sq_(sqr(infer_parton_width(var_map))),
      parton_radius_sq_(sqr(radius_widths)*parton_width_sq_),
      npartons_(std::size_t(var_map["parton-number"].as<int>())),
      sigma_partonic_(partonic_cross_section(var_map)),
      prefactor_(math::double_constants::one_div_two_pi/parton_width_sq_/npartons_),
      participant_fluctuation_dist_(gamma_param_unit_mean(var_map["fluctuation"].as<double>())),
      parton_position_dist_(0, compute_parton_sampling_width(var_map))
{}

MonteCarloCrossSection::MonteCarloCrossSection(const VarMap& var_map)
  : npartons(std::size_t(var_map["parton-number"].as<int>())),
    nucleon_width(var_map["nucleon-width"].as<double>()),
    parton_width(infer_parton_width(var_map)),
    parton_sampling_width(compute_parton_sampling_width(var_map)),
    parton_width_sq(sqr(parton_width)),
    prefactor(math::double_constants::one_div_two_pi/parton_width_sq/npartons)
{}

double MonteCarloCrossSection::operator() (const double sigma_partonic) const {

  random::CyclicNormal<> cyclic_normal{
    0., parton_sampling_width, cache_size, n_loops
  };

  struct Parton {
    double x, y;
  };

  std::vector<Parton> nucleonA(npartons);
  std::vector<Parton> nucleonB(npartons);

  auto bmax = max_impact_widths * nucleon_width;
  auto arg_max = 0.25*sqr(bmax)/parton_width_sq; 
  const FastExp<double> fast_exp(-arg_max, 0., 1000);

  double ref_cross_section = 0.;
  double prob_miss = 0.;
  int pass_tolerance = 0;

  for (std::size_t n = 0; n < n_max; ++n) {
    // Sample b from P(b)db = 2*pi*b.
    auto b = bmax * std::sqrt(random::canonical<double>());

    for (auto&& p : nucleonA) {
      p.x = cyclic_normal(random::engine);
      p.y = cyclic_normal(random::engine);
    }

    for (auto&& p : nucleonB) {
      p.x = cyclic_normal(random::engine);
      p.y = cyclic_normal(random::engine);
    }

    auto overlap = 0.;
    for (auto&& pa : nucleonA) {
      for (auto&& pb : nucleonB) {
        auto distance_sq = sqr(pa.x - pb.x + b) + sqr(pa.y - pb.y);
        auto arg = .25*distance_sq/parton_width_sq;
        if (arg < arg_max) overlap += fast_exp(-arg);
      }
    }

    prob_miss +=
      std::exp(-sigma_partonic * prefactor/(2.*npartons) * overlap);

    auto fraction_hit = 1. - (prob_miss/n);
    auto sample_area = math::double_constants::pi * sqr(bmax);
    auto cross_section = fraction_hit * sample_area;

    auto update_difference =
      std::abs((cross_section - ref_cross_section)/cross_section);

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
      -- check nucleon-width, parton-width and npartons"};

  // Code throws an error before it ever gets here.
  return -1;
}

}  // namespace trento
