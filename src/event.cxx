// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "event.h"

#include <string>
#include <vector>

#include <boost/program_options/variables_map.hpp>

#include "fwd_decl.h"
#include "nucleus.h"

namespace trento {

namespace {

// Helper function to create NucleusPtr members of Event.
NucleusPtr make_nucleus(const VarMap& var_map, std::size_t index) {
  const auto& species = var_map["projectile"]
                        .as<std::vector<std::string>>().at(index);
  return NucleusBase::create(species);
}

}  // unnamed namespace

Event::Event(const VarMap& var_map)
  : nucA_(make_nucleus(var_map, 0)),  // make_nucleus defined above
    nucB_(make_nucleus(var_map, 1)),
    b_min_(var_map["b-min"].as<double>()),
    b_max_(var_map["b-max"].as<double>()) {}

// Explicitly-defined default destructor; see header for explanation.
Event::~Event() = default;

void Event::collide() {
  nucA_->sample_nucleons(-1.);
  nucB_->sample_nucleons(+1.);
}


}  // namespace trento
