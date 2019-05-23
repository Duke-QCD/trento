// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "../src/nucleon.h"

#include <cmath>

#include "catch.hpp"
#include "util.h"

#include "../src/nucleus.h"
#include "../src/random.h"
#include <boost/math/constants/constants.hpp>
#include <boost/multi_array.hpp>

#include <iostream>

using namespace trento;

TEST_CASE( "nucleon" ) {
  auto fluct = 1. + .5*random::canonical<>();
  auto xsec = 4. + 3.*random::canonical<>();
  auto width = .5 + .2*random::canonical<>();
  auto wsq = width*width;

  auto var_map = make_var_map({
      {"ncoll", false},
      {"fluctuation",   fluct},
      {"cross-section", xsec},
      {"nucleon-width", width},
      {"constit-width", width},
      {"constit-number", 1},
  });

  NucleonCommon nc{var_map};

  Proton proton{};
  proton.sample_nucleons(0.);
  auto& nucleon = *proton.begin();
  while (!nc.participate(nucleon, nucleon)) {}

  // truncation radius
  auto R = 5*width;

  // random angle
  auto phi = math::double_constants::two_pi * random::canonical<>();
  auto cs = std::cos(phi);
  auto sn = std::sin(phi);

  // thickness function
  // check relative to zero
  auto tzero = nc.thickness(nucleon, 0., 0.);
  CHECK( nc.thickness(nucleon, width*cs, width*sn) == Approx(tzero*std::exp(-.5)) );

  // random point inside radius
  auto x = R * random::canonical<>();
  auto y = R * random::canonical<>();
  auto dsq = x*x + y*y;
  CHECK( nc.thickness(nucleon, x, y) == Approx(tzero*std::exp(-.5*dsq/wsq)).epsilon(1e-5).margin(1e-5) );

  // random point outside radius
  auto d = R * (1 + random::canonical<>());
  CHECK( nc.thickness(nucleon, d*cs, d*sn) == 0. );

  // fluctuations
  // just check they have unit mean -- the rest is handled by the C++ impl.
  auto total = 0.;
  auto n = 1e6;
  for (auto i = 0; i < static_cast<int>(n); ++i) {
    proton.sample_nucleons(0.);
    nucleon = *proton.begin();
    while (!nc.participate(nucleon, nucleon)) {}
    total += nc.thickness(nucleon, 0., 0.) * (2*M_PI*wsq);
  }

  auto mean = total/n;
  CHECK( mean == Approx(1.).epsilon(.003) );

  // must use a Nucleus to set Nucleon position
  // Proton conveniently sets a deterministic position
  // a mock class would be better but this works fine
  Proton A{}, B{};
  A.sample_nucleons(0.);
  B.sample_nucleons(0.);
  auto& nA = *A.begin();
  auto& nB = *B.begin();
  CHECK( nA.x() == 0. );
  CHECK( nA.y() == 0. );
  CHECK( nA.z() == 0. );
  CHECK( !nA.is_participant() );

  // wait until the nucleons participate
  while (!nc.participate(nA, nB)) {}
  CHECK( nA.is_participant() );
  CHECK( nB.is_participant() );

  // resampling nucleons resets participant state
  A.sample_nucleons(0.);
  CHECK( !nA.is_participant() );

  // test cross section
  // min-bias impact params
  auto bmax = nc.max_impact();
  CHECK( bmax == Approx(6*width) );

  auto nev = 1e6;
  auto count = 0;
  for (auto i = 0; i < static_cast<int>(nev); ++i) {
    auto b = bmax * std::sqrt(random::canonical<>());
    A.sample_nucleons(.5*b);
    B.sample_nucleons(-.5*b);
    if (nc.participate(nA, nB))
      ++count;
  }

  auto xsec_mc = M_PI*bmax*bmax * static_cast<double>(count)/nev;

  // precision is better than this, but let's be conservative
  CHECK( xsec_mc == Approx(xsec).epsilon(.02) );

  // impact larger than max should never participate
  auto b = bmax + random::canonical<>();
  A.sample_nucleons(.5*b);
  B.sample_nucleons(-.5*b);
  CHECK( !nc.participate(nA, nB) );

  // very large fluctuation parameters mean no fluctuations
  auto no_fluct_var_map = make_var_map({
      {"ncoll", false},
      {"fluctuation",   1e12},
      {"cross-section", xsec},
      {"nucleon-width", width},
      {"constit-width", width},
      {"constit-number", 1},
  });

  NucleonCommon no_fluct_nc{no_fluct_var_map};

  proton.sample_nucleons(0.);
  nucleon = *proton.begin();
  while (!no_fluct_nc.participate(nucleon, nucleon)) {}
  CHECK( no_fluct_nc.thickness(nucleon, 0, 0) == Approx(1/(2*M_PI*wsq)) );

  CHECK_THROWS_AS([]() {
    // nucleon width too small
    auto bad_var_map = make_var_map({
        {"ncoll", false},
        {"fluctuation",   1.},
        {"cross-section", 5.},
        {"nucleon-width", .1},
        {"constit-width", .1},
        {"constit-number", 1},
    });
    NucleonCommon bad_profile{bad_var_map};
  }(),
  std::domain_error);

  // test nucleon-nucleon cross section with random constituent substructure
  auto nucleon_width = 0.5;
  auto constituent_number = std::uniform_int_distribution<>{2, 6}(random::engine);
  auto constituent_width = .3 + .2*random::canonical<>();

  // Coarse-ish p+p grid.
  auto grid_max = 3.;
  auto grid_step = 0.1;
  auto grid_nsteps = 60;

  // constituent var map
  var_map = make_var_map({
      {"ncoll", false},
      {"grid-max", grid_max},
      {"grid-step", grid_step},
      {"fluctuation",   1e12},
      {"cross-section", xsec},
      {"nucleon-width", nucleon_width},
      {"constit-width", constituent_width},
      {"constit-number", constituent_number},
  });

  NucleonCommon nc_constituent{var_map};
  bmax = nc_constituent.max_impact();

  // inelastic collision counter
  count = 0;

  // measure the cross section with min bias p+p collisions
  for (auto i = 0; i < static_cast<int>(nev); ++i) {
    b = bmax * std::sqrt(random::canonical<>());
    A.sample_nucleons(.5*b);
    B.sample_nucleons(-.5*b);
    if (nc_constituent.participate(*A.begin(), *B.begin()))
      ++count;
  }

  // calculate Monte Carlo cross section
  xsec_mc = M_PI*bmax*bmax * static_cast<double>(count)/nev;

  // assert a rather loose tolerance for the numerical cross section finder
  CHECK( xsec_mc/xsec == Approx(1).epsilon(.02) );

  // verify root-mean-square nucleon width for nucleons with substructure
  double mean_square_width_num = 0.0;
  double mean_square_width_denom = 0.0;

  for (auto i = 0; i < 1000; ++i) {
    A.sample_nucleons(0.);
    auto nucleon = *A.begin();

    while (!nc_constituent.participate(nucleon, nucleon)) {}

    for (auto iy = 0; iy < grid_nsteps; ++iy) {
      for (auto ix = 0; ix < grid_nsteps; ++ix) {
        auto x = (ix + .5) * 2 * grid_max/grid_nsteps - grid_max;
        auto y = (iy + .5) * 2 * grid_max/grid_nsteps - grid_max;
        auto thick = nc_constituent.thickness(nucleon, x, y);
        mean_square_width_num += (x*x + y*y)*thick;
        mean_square_width_denom += 2*thick;
      }
    }
  }

  auto rms_width = std::sqrt(mean_square_width_num / mean_square_width_denom);

  CHECK( rms_width == Approx(nucleon_width).epsilon(.01).margin(.01) );

}
