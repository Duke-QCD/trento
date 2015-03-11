// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "../src/nucleus.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <map>
#include <sstream>

#include "catch.hpp"

using namespace trento;

TEST_CASE( "nucleus" ) {

  SECTION( "proton" ) {

    auto nucleus = NucleusBase::create("p");

    CHECK( std::distance(nucleus->begin(), nucleus->end()) == 1 );
    CHECK( std::distance(nucleus->cbegin(), nucleus->cend()) == 1 );

    CHECK( nucleus->radius() == 0. );

    nucleus->sample_nucleons(0.);

    auto proton = *(nucleus->begin());
    auto coord = std::get<0>(proton);
    auto part = std::get<1>(proton);
    auto x = std::get<0>(coord);
    auto y = std::get<1>(coord);

    CHECK( x == 0. );
    CHECK( y == 0. );
    CHECK( part == false );

  }

  SECTION( "lead" ) {

    auto nucleus = NucleusBase::create("Pb");

    CHECK( std::distance(nucleus->begin(), nucleus->end()) == 208 );
    CHECK( std::distance(nucleus->cbegin(), nucleus->cend()) == 208 );

    CHECK( nucleus->radius() > 6. );

    nucleus->sample_nucleons(0.);

    auto initial_participants = std::any_of(
        nucleus->cbegin(), nucleus->cend(),
        [](const NucleusBase::NucleonAttr& n) { return std::get<1>(n); });
    CHECK( initial_participants == false );

    constexpr auto n_events = 1000;
    constexpr auto n_samples = n_events*208;
    constexpr auto dr = 1.;
    std::map<int, int> hist{};

    for (auto i = 0; i < n_events; ++i) {
      nucleus->sample_nucleons(0.);
      for (const auto& nucleon : *nucleus) {
        auto coord = std::get<0>(nucleon);
        auto x = std::get<0>(coord);
        auto y = std::get<1>(coord);
        auto r = std::sqrt(x*x + y*y);
        ++hist[r/dr];
      }
    }

    std::ostringstream hist_output{};
    for (const auto& bin : hist) {
      hist_output << dr*bin.first << ": "
                  << std::string(300*bin.second/n_samples/dr, '*') << '\n';
    }
    INFO( hist_output.str() );
    CHECK( true );

  }

  SECTION( "unknown species" ) {
    CHECK_THROWS_AS( NucleusBase::create("hello"), std::invalid_argument );
  }

}
