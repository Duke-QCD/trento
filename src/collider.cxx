// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "collider.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

#include <boost/program_options/variables_map.hpp>

#include "fwd_decl.h"
#include "nucleus.h"

namespace trento {

namespace {

// Create one nucleus from the configuration.
// Helper function for the Collider ctor.
NucleusPtr create_nucleus(const VarMap& var_map, std::size_t index) {
  const auto& species = var_map["projectile"]
                        .as<std::vector<std::string>>().at(index);
  return Nucleus::create(species);
}

// Non-member functions to get the first and last element of a boost::multi_array.
// For iterating over the entire flattened grid.
// Provide const and non-const overloads.

template <typename T, std::size_t D, typename A>
inline T* flat_begin(boost::multi_array<T, D, A>& a) {
  return a.origin();
}

template <typename T, std::size_t D, typename A>
inline const T* flat_begin(const boost::multi_array<T, D, A>& a) {
  return a.origin();
}

template <typename T, std::size_t D, typename A>
inline const T* flat_cbegin(const boost::multi_array<T, D, A>& a) {
  return a.origin();
}

template <typename T, std::size_t D, typename A>
inline T* flat_end(boost::multi_array<T, D, A>& a) {
  return flat_begin(a) + a.num_elements();
}

template <typename T, std::size_t D, typename A>
inline const T* flat_end(const boost::multi_array<T, D, A>& a) {
  return flat_begin(a) + a.num_elements();
}

template <typename T, std::size_t D, typename A>
inline const T* flat_cend(const boost::multi_array<T, D, A>& a) {
  return flat_cbegin(a) + a.num_elements();
}

}  // unnamed namespace

// Lots of members to initialize...
Collider::Collider(const VarMap& var_map)
    : nucleusA_(create_nucleus(var_map, 0)),  // create_nucleus defined above
      nucleusB_(create_nucleus(var_map, 1)),
      nucleon_(var_map),
      output_functions_(create_output_functions(var_map)),
      num_events_(var_map["number-events"].as<int>()),
      b_min_(var_map["b-min"].as<double>()),
      b_max_([&var_map, this]() {  // create and immediately call a lambda
        auto bmax = var_map["b-max"].as<double>();
        if (bmax < 0.)
          bmax = nucleusA_->radius() +
                 nucleusB_->radius() +
                 2.*nucleon_.radius();
        return bmax;
      }()),
      asymmetry_([](double a, double b) {  // create and immediately call a lambda
        auto sum = a + b;
        if (sum < 0.1)
          return 0.5;
        else
          return a/sum;
      }(nucleusA_->radius(), nucleusB_->radius())),
      grid_steps_(var_map["grid-steps"].as<int>()),
      grid_half_width_(.5*var_map["grid-width"].as<double>()),
      grid_delta_(var_map["grid-width"].as<double>()/(grid_steps_ - 1)),
      event_(grid_steps_) {
  // Constructor body begins here.
  // Set random seed if requested.
  auto seed = var_map["random-seed"].as<int64_t>();
  if (seed > 0)
    random::engine.seed(static_cast<random::Engine::result_type>(seed));

  // Set reduced thickness function based on configuraton.
  auto norm = var_map["normalization"].as<double>();
  auto p = var_map["reduced-thickness"].as<double>();

  if (std::abs(p) < 1e-12) {
    compute_reduced_thickness = [this, norm]() {
      std::transform(flat_cbegin(event_.TA), flat_cend(event_.TA),
                     flat_cbegin(event_.TB),
                     flat_begin(event_.TR),
                     [norm](double a, double b) {
                       return norm * std::sqrt(a*b);
                     });
    };
  } else {
    compute_reduced_thickness = [this, norm, p]() {
      std::transform(flat_cbegin(event_.TA), flat_cend(event_.TA),
                     flat_cbegin(event_.TB),
                     flat_begin(event_.TR),
                     [norm, p](double a, double b) {
                       return norm *
                              std::pow(.5*(std::pow(a, p) + std::pow(b, p)),
                                       1./p);
                     });
    };
  }
}

Collider::~Collider() = default;

void Collider::run_events() {
  for (event_.num = 0; event_.num < num_events_; ++event_.num) {
    new_event();

    event_.npart = 0;
    compute_nuclear_thickness_and_npart(*nucleusA_, event_.TA);
    compute_nuclear_thickness_and_npart(*nucleusB_, event_.TB);

    compute_reduced_thickness();

    compute_observables();

    for (const auto& write : output_functions_)
      write(event_);
  }
}

void Collider::new_event() {
  for (bool collision = false; !collision; ) {
    event_.impact_param = b_min_ +
                          (b_max_ - b_min_) * std::sqrt(random::canonical<>());

    nucleusA_->sample_nucleons(asymmetry_ * event_.impact_param);
    nucleusB_->sample_nucleons((asymmetry_ - 1.) * event_.impact_param);

    for (auto&& A : *nucleusA_) {
      for (auto&& B : *nucleusB_) {
        if (!A.is_participant() || !B.is_participant()) {
          double dx = A.x() - B.x();
          double dy = A.y() - B.y();
          double distance_squared = dx*dx + dy*dy;
          if (nucleon_.participate(distance_squared)) {
            collision = true;
            A.set_participant();
            B.set_participant();
          }
        }
      }
    }
  }
}

void Collider::compute_nuclear_thickness_and_npart(const Nucleus& nucleus,
                                                   Event::Grid& TX) {
  // Wipe grid with zeros.
  std::fill(flat_begin(TX), flat_end(TX), 0.);

  const double nucleon_radius = nucleon_.radius();
  const double nucleon_xymax = grid_half_width_ + nucleon_radius;

  // Deposit each participant onto the grid.
  for (const auto& n : nucleus)
    if (n.is_participant()) {
      ++event_.npart;
      double x = n.x();
      double y = n.y();

      // Nucleon is completely off the grid (not even the edge touches).
      if ((std::abs(x) > nucleon_xymax) || (std::abs(y) > nucleon_xymax))
        continue;

      // Determine min & max indices of nucleon subgrid.
      int ixmin = std::max(static_cast<int>(std::floor(
                    (x - nucleon_radius + grid_half_width_)/grid_delta_
                  )), 0);
      int iymin = std::max(static_cast<int>(std::floor(
                    (y - nucleon_radius + grid_half_width_)/grid_delta_
                  )), 0);
      int ixmax = std::min(static_cast<int>(std::ceil(
                    (x + nucleon_radius + grid_half_width_)/grid_delta_
                  )), grid_steps_ - 1);
      int iymax = std::min(static_cast<int>(std::ceil(
                    (y + nucleon_radius + grid_half_width_)/grid_delta_
                  )), grid_steps_ - 1);

      double fluct = nucleon_.fluctuate();

      // Add profile to grid.
      for (auto iy = iymin; iy <= iymax; ++iy) {
        double dysq = std::pow(y - (iy*grid_delta_) + grid_half_width_, 2);
        for (auto ix = ixmin; ix <= ixmax; ++ix) {
          double dxsq = std::pow(x - (ix*grid_delta_) + grid_half_width_, 2);
          TX[iy][ix] += fluct * nucleon_.thickness(dxsq + dysq);
        }
      }
    }
}

void Collider::compute_observables() {
  double sum = 0.;
  double xcm = 0.;
  double ycm = 0.;

  // Compute multiplicity and center of mass coordinates.
  // TODO: include in initial TR loop
  for (int iy = 0; iy < grid_steps_; ++iy) {
    for (int ix = 0; ix < grid_steps_; ++ix) {
      const auto& t = event_.TR[iy][ix];
      sum += t;
      xcm += t * ix;
      ycm += t * iy;
    }
  }

  event_.multiplicity = grid_delta_ * grid_delta_ * sum;
  xcm /= sum;
  ycm /= sum;

  // Compute eccentricity.

  // Create map of (n, (real, imag, weight)).
  std::map<int, std::array<double, 3>> ecc;

  for (int iy = 0; iy < grid_steps_; ++iy) {
    for (int ix = 0; ix < grid_steps_; ++ix) {
      const auto& t = event_.TR[iy][ix];
      if (t < 1e-12)
        continue;

      auto x = static_cast<double>(ix) - xcm;
      auto x2 = x*x;
      auto x3 = x2*x;
      auto x4 = x2*x2;

      auto y = static_cast<double>(iy) - ycm;
      auto y2 = y*y;
      auto y3 = y2*y;
      auto y4 = y2*y2;

      auto r2 = x2 + y2;
      auto r = std::sqrt(r2);
      auto r4 = r2*r2;

      auto xy = x*y;
      auto x2y2 = x2*y2;

      // TODO: explain what is happening here
      ecc[2][0] += t * (y2 - x2);
      ecc[2][1] += t * 2.*xy;
      ecc[2][2] += t * r2;

      ecc[3][0] += t * (y3 - 3.*y*x2);
      ecc[3][1] += t * (3.*x*y2 - x3);
      ecc[3][2] += t * r2*r;

      ecc[4][0] += t * (x4 + y4 - 6.*x2y2);
      ecc[4][1] += t * 4.*xy*(y2 - x2);
      ecc[4][2] += t * r4;

      ecc[5][0] += t * y*(5.*x4 - 10.*x2y2 + y4);
      ecc[5][1] += t * x*(x4 - 10.*x2y2 + 5.*y4);
      ecc[5][2] += t * r4*r;
    }
  }

  for (const auto& e : ecc) {
    const auto& re = e.second[0];
    const auto& im = e.second[1];
    const auto w = std::fmax(e.second[2], 1e-12);
    event_.eccentricity[e.first] = std::sqrt(re*re + im*im) / w;
  }
}

}  // namespace trento
