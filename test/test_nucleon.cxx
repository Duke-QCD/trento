// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "../src/nucleon.h"

#include <algorithm>
#include <cmath>

#include "catch.hpp"
#include "util.h"

#include "../src/random.h"

using namespace trento;

TEST_CASE( "nucleon properties",
           "[nucleon][thickness][cross-section][fluctuations]" ) {
  auto fluct = 1. + .5*random::canonical<>();
  auto xsec = 4. + 3.*random::canonical<>();
  auto width = .5 + .2*random::canonical<>();
  auto wsq = width*width;

  auto var_map = make_var_map({
      {"fluctuation",   fluct},
      {"cross-section", xsec},
      {"nucleon-width", width},
  });

  Nucleon nucleon{var_map};

  // truncation radius
  auto R = nucleon.radius();
  CHECK( R == Approx(3*width) );

  // thickness function
  CHECK( nucleon.thickness(0.) == Approx(1/wsq) );
  CHECK( nucleon.thickness(wsq) == Approx(1/wsq*std::exp(-.5)) );

  // random point inside radius
  auto dsq = std::pow(R*random::canonical<>(), 2);
  CHECK( nucleon.thickness(dsq) == Approx(1/wsq*std::exp(-.5*dsq/wsq)) );

  // random point outside radius
  dsq = std::pow(R*(1+random::canonical<>()), 2);
  CHECK( nucleon.thickness(dsq) == 0. );

  // cross section
  // min-bias impact params
  auto bsqmax = std::pow(1.5*2*R, 2);
  auto nev = 1e6;
  auto count = 0;
  for (auto i = 0; i < static_cast<int>(nev); ++i)
    if (nucleon.participate(bsqmax*random::canonical<>()))
      ++count;

  auto xsec_mc = M_PI*bsqmax * static_cast<double>(count)/nev;

  // precision is better than 1%, but let's be conservative
  CHECK( xsec_mc == Approx(xsec).epsilon(.01) );

  // fluctuations
  // just check they have unit mean -- the rest is handled by the C++ impl.
  auto total = 0.;
  auto n = 1e6;
  for (auto i = 0; i < static_cast<int>(n); ++i)
    total += nucleon.fluctuate();

  auto mean = total/n;
  CHECK( mean == Approx(1.).epsilon(.001) );
}

TEST_CASE( "exception when nucleon width too small",
           "[nucleon][exception]" ) {
  auto var_map = make_var_map({
      {"fluctuation",   1.},
      {"cross-section", 5.},
      {"nucleon-width", .1},
  });
  CHECK_THROWS_AS( Nucleon nucleon{var_map}, std::domain_error );
}
