// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "collider.h"

#include <algorithm>
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
      grid_delta_(var_map["grid-width"].as<double>()/grid_steps_),
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

  const double r = nucleon_.radius();

  // Deposit each participant onto the grid.
  for (const auto& n : nucleus) {
    if (!n.is_participant())
      continue;

    ++event_.npart;

    // Work in coordinates relative to (-width/2, -width/2).
    double x = n.x() + grid_half_width_;
    double y = n.y() + grid_half_width_;

    // Determine min & max indices of nucleon subgrid.
    int ixmin = std::max(static_cast<int>((x-r)/grid_delta_), 0);
    int iymin = std::max(static_cast<int>((y-r)/grid_delta_), 0);
    int ixmax = std::min(static_cast<int>((x+r)/grid_delta_), grid_steps_ - 1);
    int iymax = std::min(static_cast<int>((y+r)/grid_delta_), grid_steps_ - 1);

    double fluct = nucleon_.fluctuate();

    // Add profile to grid.
    for (auto iy = iymin; iy <= iymax; ++iy) {
      double dysq = std::pow(y - (iy+.5)*grid_delta_, 2);
      for (auto ix = ixmin; ix <= ixmax; ++ix) {
        double dxsq = std::pow(x - (ix+.5)*grid_delta_, 2);
        TX[iy][ix] += fluct * nucleon_.thickness(dxsq + dysq);
      }
    }
  }
}

namespace {

// Helper for eccentricity calculations.
// Used below in Collider::compute_observables();
struct EccentricityAccumulator {
  double re = 0.;  // real part
  double im = 0.;  // imaginary part
  double wt = 0.;  // weight
};

double compute_final_ecc(const EccentricityAccumulator& e) {
  return std::sqrt(e.re*e.re + e.im*e.im) / std::fmax(e.wt, 1e-12);
}

}  // unnamed namespace

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
  EccentricityAccumulator e2, e3, e4, e5;

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

  event_.eccentricity[2] = compute_final_ecc(e2);
  event_.eccentricity[3] = compute_final_ecc(e3);
  event_.eccentricity[4] = compute_final_ecc(e4);
  event_.eccentricity[5] = compute_final_ecc(e5);
}

}  // namespace trento
