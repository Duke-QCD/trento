// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "../src/output.h"

#include "catch.hpp"
#include "util.h"

#include <iostream>
#include <type_traits>

#include "../src/event.h"
#include "../src/nucleus.h"

using namespace trento;

// redirect stdout to a stringstream and safely restore upon destruction
struct capture_stdout {
  capture_stdout() {
    // replace stdout buffer
    cout_orig = std::cout.rdbuf(sstream.rdbuf());
  }

  ~capture_stdout() {
    // restore stdout to working state
    std::cout.rdbuf(cout_orig);
  }

  bool is_empty() {
    return sstream.peek() == std::stringstream::traits_type::eof();
  }

  std::streambuf* cout_orig;
  std::stringstream sstream;
};

TEST_CASE( "output" ) {
  auto var_map = make_var_map({
      {"normalization", 1.},
      {"reduced-thickness", 0.},
      {"grid-width", 18.},
      {"grid-steps", 90},
      {"fluctuation",   1.},
      {"cross-section", 6.4},
      {"nucleon-width", 0.5}
  });

  // create a test event
  Event event{var_map};
  NucleonProfile profile{var_map};

  auto nucleusA = Nucleus::create("Pb");
  auto nucleusB = Nucleus::create("Pb");

  auto b = 4.*std::sqrt(random::canonical<>());
  nucleusA->sample_nucleons(+.5*b);
  nucleusB->sample_nucleons(-.5*b);

  for (auto&& A : *nucleusA)
    for (auto&& B : *nucleusB)
      profile.participate(A, B);

  event.compute(*nucleusA, *nucleusB, profile);

  capture_stdout capture{};

  SECTION( "no output" ) {
    auto output_var_map = make_var_map({
        {"quiet", true},
        {"number-events", 1}
    });

    Output output{output_var_map};
    output(1, b, event);

    CHECK( capture.is_empty() );
  }

  SECTION( "stdout only" ) {
    auto nev = 11;
    auto output_var_map = make_var_map({
        {"quiet", false},
        {"number-events", nev}
    });

    Output output{output_var_map};
    for (int n = 0; n < nev; ++n)
      output(n, b, event);

    CHECK( capture.sstream.str().substr(0, 2) == " 0" );

    int num, npart;
    double impact, mult, e2, e3, e4, e5;
    capture.sstream >> num >> impact >> npart >> mult >> e2 >> e3 >> e4 >> e5;
    CHECK( num == 0 );
    CHECK( impact == Approx(b) );
    CHECK( npart == event.npart() );
    CHECK( mult == Approx(event.multiplicity()) );
    CHECK( e2 == Approx(event.eccentricity().at(2)) );
    CHECK( e3 == Approx(event.eccentricity().at(3)) );
    CHECK( e4 == Approx(event.eccentricity().at(4)) );
    CHECK( e5 == Approx(event.eccentricity().at(5)) );

    capture.sstream >> num;
    CHECK( num == 1 );

    std::string line;
    for (int k = 1; k < nev; ++k)
      std::getline(capture.sstream, line);
    CHECK( line.substr(0, 2) == "10" );

    CHECK( capture.is_empty() );
  }

  // SECTION( "text and stdout" ) {
  //   TODO
  // }
}
