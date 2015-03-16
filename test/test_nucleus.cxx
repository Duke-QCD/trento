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

  Nucleus::NucleonData nucleon_data{};

  // should not be a participant initially
  CHECK( !nucleon_data.is_participant() );

  // set and retrieve participant
  nucleon_data.set_participant();
  CHECK( nucleon_data.is_participant() );

}

TEST_CASE( "proton" ) {

  auto nucleus = Nucleus::create("p");

  // proton contains one nucleon
  CHECK( std::distance(nucleus->begin(), nucleus->end()) == 1 );
  CHECK( std::distance(nucleus->cbegin(), nucleus->cend()) == 1 );

  // and has zero radius
  CHECK( nucleus->radius() == 0. );

  // sample position with random offset
  double offset = random::canonical<>();
  nucleus->sample_nucleons(offset);
  auto&& proton = *(nucleus->begin());

  // check correct position
  CHECK( proton.x() == offset );
  CHECK( proton.y() == 0. );

  // not a participant initially
  CHECK( !proton.is_participant() );

  // set as participant
  proton.set_participant();
  CHECK( proton.is_participant() );

  // resampling nucleons resets participant state
  nucleus->sample_nucleons(0.);
  CHECK( !proton.is_participant() );
}

TEST_CASE( "lead nucleus" ) {

  auto nucleus = Nucleus::create("Pb");

  int A = 208;
  CHECK( std::distance(nucleus->begin(), nucleus->end()) == A );
  CHECK( std::distance(nucleus->cbegin(), nucleus->cend()) == A );

  CHECK( nucleus->radius() > 6. );

  double offset = nucleus->radius() * random::canonical<>();
  nucleus->sample_nucleons(offset);

  // average nucleon position
  double x = 0., y = 0.;
  for (const auto& nucleon : *nucleus) {
    x += nucleon.x();
    y += nucleon.y();
  }
  x /= A;
  y /= A;
  auto tolerance = .6;
  CHECK( std::abs(x - offset) < tolerance );
  CHECK( std::abs(y) < tolerance );

  // no initial participants
  bool initial_participants = std::any_of(
      nucleus->cbegin(), nucleus->cend(),
      [](decltype(*nucleus->cbegin())& n) {
        return n.is_participant();
      });
  CHECK( !initial_participants );

}

TEST_CASE( "woods-saxon sampling" ) {

  int A = 200;
  double R = 6., a = .5;
  NucleusPtr nucleus{new WoodsSaxonNucleus{static_cast<std::size_t>(A), R, a}};

  // Test Woods-Saxon sampling.
  // This is honestly not a great test; while it does prove that the code
  // basically works as intended, it does not rigorously show that the generated
  // numbers are actually Woods-Saxon distributed.  Personally I am much more
  // convinced by plotting a histogram and visually comparing it to the smooth
  // curve.  The script 'plot-woods-saxon.py' in the 'woods-saxon' subdirectory
  // does this.

  // sample a bunch of nuclei and bin all the nucleon positions
  // using "p" for cylindrical radius to distinguish from spherical "r"
  auto dp = .5;
  auto nevents = 1000;
  auto nsamples = nevents * A;
  std::map<int, int> hist{};
  for (auto i = 0; i < nevents; ++i) {
    nucleus->sample_nucleons(0.);
    for (const auto& nucleon : *nucleus) {
      auto x = nucleon.x();
      auto y = nucleon.y();
      auto p = std::sqrt(x*x + y*y);
      ++hist[p/dp];
    }
  }

  // integrate the W-S dist from pmin to pmax and over all z
  double rmax = R + 10.*a;
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
  CHECK_THROWS_AS( Nucleus::create("hello"), std::invalid_argument );
}
