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
// Useful for looping over the entire flattened grid.
// Named first() and last() to distinguish from begin() and end().
// Provide const and non-const overloads.

template <typename T, std::size_t D, typename A>
inline T* first(boost::multi_array<T, D, A>& a) {
  return a.origin();
}

template <typename T, std::size_t D, typename A>
inline const T* first(const boost::multi_array<T, D, A>& a) {
  return a.origin();
}

template <typename T, std::size_t D, typename A>
inline T* last(boost::multi_array<T, D, A>& a) {
  return a.origin() + a.num_elements();
}

template <typename T, std::size_t D, typename A>
inline const T* last(const boost::multi_array<T, D, A>& a) {
  return a.origin() + a.num_elements();
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
      grid_width_(var_map["grid-width"].as<double>()),
      grid_steps_(var_map["grid-steps"].as<int>()),
      grid_delta_(grid_width_/(grid_steps_ - 1)),
      event_(grid_steps_) {
  // Constructor body begins here.
  // Set random seed if requested.
  auto seed = var_map["random-seed"].as<int64_t>();
  if (seed > 0)
    random::engine.seed(static_cast<random::Engine::result_type>(seed));

  // Set reduced thickness function based on configuraton.
  auto norm = var_map["normalization"].as<double>();
  // auto p = var_map["reduced_thickness"].as<double>(); TODO

  compute_reduced_thickness = [this, norm]() {
    std::transform(first(event_.TA), last(event_.TA),
                   first(event_.TB),
                   first(event_.TR),
                   [norm](double a, double b) {
                     return norm * std::sqrt(a*b);
                   });
  };
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
                                                   Event::Grid& thickness) {
  // Wipe grid with zeros.
  std::fill(first(thickness), last(thickness), 0.);

  // Deposit each participant onto the grid.
  for (const auto& n : nucleus)
    if (n.is_participant()) {
      ++event_.npart;
      // TODO
    }
}

void Collider::compute_observables() {
    event_.multiplicity = grid_delta_ * grid_delta_ *
                          std::accumulate(first(event_.TR), last(event_.TR), 0.);

    for (int n = 2; n <= 5; ++n)
      event_.eccentricity[n] = random::canonical<>();
}

}  // namespace trento
