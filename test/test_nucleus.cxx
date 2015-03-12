// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "../src/nucleus.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <map>

#include "catch.hpp"

#include "../src/random.h"

using namespace trento;

TEST_CASE( "nucleon data" ) {

  NucleonData nucleon_data{};

  // should not be a participant initially
  CHECK( !nucleon_data.is_participant() );

  // set and retrieve position
  nucleon_data.set_position({1, 2});
  std::array<double, 2> correct_position{1, 2};
  CHECK( nucleon_data.position() == correct_position );

  // should still not be a participant
  CHECK( !nucleon_data.is_participant() );

  // set and retrieve participant
  nucleon_data.set_participant();
  CHECK( nucleon_data.is_participant() );

  // setting new position resets participant
  nucleon_data.set_position({3, 4});
  CHECK( !nucleon_data.is_participant() );

}

TEST_CASE( "proton" ) {

  auto nucleus = NucleusBase::create("p");

  // contains one nucleon
  CHECK( std::distance(nucleus->begin(), nucleus->end()) == 1 );
  CHECK( std::distance(nucleus->cbegin(), nucleus->cend()) == 1 );

  CHECK( nucleus->radius() == 0. );

  double offset = random::canonical<>();
  nucleus->sample_nucleons(offset);

  auto&& proton = *(nucleus->begin());
  std::array<double, 2> correct_position{-offset, 0.};

  CHECK( proton.position() == correct_position );
  CHECK( !proton.is_participant() );

}

TEST_CASE( "lead nucleus" ) {

  auto nucleus = NucleusBase::create("Pb");

  constexpr auto A = 208;

  CHECK( std::distance(nucleus->begin(), nucleus->end()) == A );
  CHECK( std::distance(nucleus->cbegin(), nucleus->cend()) == A );

  CHECK( nucleus->radius() > 6. );

  double offset = nucleus->radius() * random::canonical<>();
  nucleus->sample_nucleons(offset);

  // average nucleon position
  double x = 0., y = 0.;
  for (const auto& nucleon : *nucleus) {
    auto coord = nucleon.position();
    x += std::get<0>(coord);
    y += std::get<1>(coord);
  }
  x /= A;
  y /= A;
  auto tolerance = .6;
  CHECK( std::abs(x + offset) < tolerance );
  CHECK( std::abs(y) < tolerance );

  // no initial participants
  bool initial_participants = std::any_of(
      nucleus->cbegin(), nucleus->cend(),
      [](const NucleonData& n) { return n.is_participant(); });
  CHECK( !initial_participants );

  // sample a bunch of nuclei and bin all the nucleon positions
  // using "p" for cylindrical radius to distinguish from spherical "r"
  constexpr auto dp = .2;
  std::map<int, int> hist{};
  for (auto i = 0; i < 5000; ++i) {
    nucleus->sample_nucleons(0.);
    for (const auto& nucleon : *nucleus) {
      auto coord = nucleon.position();
      auto x = std::get<0>(coord);
      auto y = std::get<1>(coord);
      auto p = std::sqrt(x*x + y*y);
      ++hist[p/dp];
    }
  }

  // calculate the W-S distribution at radius p by numerically integrating
  // through the z-direction
  auto woods_saxon = [](double p) -> double {
    double R = 6.67, a = 0.44;
    double zmax = R + 10.*a;
    int nstep = 1000;
    double result = 0.;
    for (auto i = 0; i < nstep; ++i) {
      auto z = (i+.5)*zmax/nstep;
      result += p/(1.+std::exp((std::sqrt(p*p + z*z) - R)/a));
    }
    return result;
  };

  // normalize the histogram and the numerical W-S dist by their values at a
  // point near their peaks
  double norm_point = 5.;
  double hist_norm = hist[norm_point/dp];
  double ws_norm = woods_saxon(norm_point);

  // check all histogram bins agree with the numerical dist within stat error
  auto all_bins_correct = true;
  std::ostringstream bin_output{};
  bin_output << std::fixed << std::boolalpha
             << "p        samples  dist     pass\n";
  for (const auto& bin : hist) {
    double p = (bin.first + .5)*dp;
    double hist_result = bin.second/hist_norm;
    double ws_result = woods_saxon(p)/ws_norm;
    double tolerance = 3/std::sqrt(hist_result);
    bool within_tol = std::abs(hist_result - ws_result) < tolerance;
    all_bins_correct = within_tol;
    bin_output << p << ' '
               << hist_result << ' '
               << ws_result << ' '
               << within_tol << '\n';
  }
  INFO( bin_output.str() );
  CHECK( all_bins_correct );
}

TEST_CASE( "unknown nucleus species" ) {
  CHECK_THROWS_AS( NucleusBase::create("hello"), std::invalid_argument );
}
