// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "collider.h"

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

}  // unnamed namespace

// Lots of members to initialize...
Collider::Collider(const VarMap& var_map)
    : nucleusA_(create_nucleus(var_map, 0)),  // create_nucleus defined above
      nucleusB_(create_nucleus(var_map, 1)),
      nucleon_profile_(var_map),
      nevents_(var_map["number-events"].as<int>()),
      bmin_(var_map["b-min"].as<double>()),
      bmax_([&var_map, this]() {  // create and immediately call a lambda
        auto bmax = var_map["b-max"].as<double>();
        if (bmax < 0.)
          bmax = nucleusA_->radius() +
                 nucleusB_->radius() +
                 2.*nucleon_profile_.radius();
        return bmax;
      }()),
      asymmetry_([](double a, double b) {  // create and immediately call a lambda
        auto sum = a + b;
        if (sum < 0.1)
          return 0.5;
        else
          return a/sum;
      }(nucleusA_->radius(), nucleusB_->radius())),
      event_(var_map),
      output_(var_map) {
  // Constructor body begins here.
  // Set random seed if requested.
  auto seed = var_map["random-seed"].as<int64_t>();
  if (seed > 0)
    random::engine.seed(static_cast<random::Engine::result_type>(seed));
}

Collider::~Collider() = default;

void Collider::run_events() {
  for (int n = 0; n < nevents_; ++n) {
    double b = sample_impact_param();
    event_.compute(*nucleusA_, *nucleusB_, nucleon_profile_);
    output_(n, b, event_);
  }
}

double Collider::sample_impact_param() {
  double b;
  bool collision = false;
  while (!collision) {
    b = bmin_ + (bmax_ - bmin_) * std::sqrt(random::canonical<>());

    nucleusA_->sample_nucleons(asymmetry_ * b);
    nucleusB_->sample_nucleons((asymmetry_ - 1.) * b);

    for (auto&& A : *nucleusA_) {
      for (auto&& B : *nucleusB_) {
        if (nucleon_profile_.participate(A, B))
          collision = true;
      }
    }
  }

  return b;
}

}  // namespace trento
