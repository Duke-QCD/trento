// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "../src/nucleon.h"

#include "catch.hpp"

using namespace trento;


TEST_CASE( "Nucleons" ) {
  Nucleon::set_width(0.5);
  Nucleon::set_cross_sec(6.5);
  Nucleon::set_fluct_shape(1.);

  Nucleon n{};
  CHECK( 1 == 1 );
}
