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

  // Test Woods-Saxon sampling.
  // This is honestly not a great test; while it does prove that the code
  // basically works as intended, it does not rigorously show that the generated
  // numbers are actually Woods-Saxon distributed.  Personally I am much more
  // convinced by plotting a histogram and visually comparing it to the smooth
  // curve.  The script 'plot-woods-saxon.py' in the 'woods-saxon' subdirectory
  // does this.

  // sample a bunch of nuclei and bin all the nucleon positions
  // using "p" for cylindrical radius to distinguish from spherical "r"
  constexpr auto dp = .5;
  constexpr auto nevents = 1000;
  constexpr auto nsamples = nevents * A;
  std::map<int, int> hist{};
  for (auto i = 0; i < nevents; ++i) {
    nucleus->sample_nucleons(0.);
    for (const auto& nucleon : *nucleus) {
      auto coord = nucleon.position();
      auto x = std::get<0>(coord);
      auto y = std::get<1>(coord);
      auto p = std::sqrt(x*x + y*y);
      ++hist[p/dp];
    }
  }

  // integrate the W-S dist from pmin to pmax and over all z
  constexpr double R = 6.67, a = 0.44;
  constexpr double rmax = R + 10.*a;
  auto integrate_woods_saxon = [R, a, rmax](double pmin, double pmax) -> double {
    int nzbins = 1000;
    auto npbins = static_cast<int>(nzbins*(pmax-pmin)/rmax);
    double result = 0.;
    for (auto ip = 0; ip <= npbins; ++ip) {
      auto p = pmin + (ip*(pmax-pmin))/(npbins - 1);
      for (auto iz = 0; iz < nzbins; ++iz) {
        auto z = (iz*rmax)/(nzbins-1);
        result += p/(1.+std::exp((std::sqrt(p*p + z*z) - R)/a));
      }
    }
    return result;
  };

  double ws_norm = integrate_woods_saxon(0, rmax);

  // check all histogram bins agree with numerical integration
  auto all_bins_correct = true;
  std::ostringstream bin_output{};
  bin_output << std::fixed << std::boolalpha
             << "pmin     pmax     prob     cprob    ratio    pass\n";
  for (const auto& bin : hist) {
    auto pmin = bin.first * dp;
    auto pmax = pmin + dp;
    auto prob = static_cast<double>(bin.second) / nsamples;
    auto correct_prob = integrate_woods_saxon(pmin, pmax) / ws_norm;
    bool within_tol = prob == Approx(correct_prob).epsilon(.1);
    if (!within_tol)
      all_bins_correct = false;
    bin_output << pmin << ' '
               << pmax << ' '
               << prob << ' '
               << correct_prob << ' '
               << prob/correct_prob << ' '
               << within_tol << '\n';
  }
  INFO( bin_output.str() );
  CHECK( all_bins_correct );
}

TEST_CASE( "unknown nucleus species" ) {
  CHECK_THROWS_AS( NucleusBase::create("hello"), std::invalid_argument );
}
