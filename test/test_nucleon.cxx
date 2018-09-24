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
      {"fluctuation",   fluct},
      {"cross-section", xsec},
      {"nucleon-width", width},
      {"constit-width", width},
      {"constit-number", 1},
      {"constit-position-radius", 0.0}
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
      {"fluctuation",   1e12},
      {"cross-section", xsec},
      {"nucleon-width", width},
      {"constit-width", width},
      {"constit-number", 1},
      {"constit-position-radius", 0.0}
  });

  NucleonCommon no_fluct_nc{no_fluct_var_map};

  proton.sample_nucleons(0.);
  nucleon = *proton.begin();
  while (!no_fluct_nc.participate(nucleon, nucleon)) {}
  CHECK( no_fluct_nc.thickness(nucleon, 0, 0) == Approx(1/(2*M_PI*wsq)) );

  CHECK_THROWS_AS([]() {
    // nucleon width too small
    auto bad_var_map = make_var_map({
        {"fluctuation",   1.},
        {"cross-section", 5.},
        {"nucleon-width", .1},
        {"constit-width", .1},
        {"constit-number", 1},
        {"constit-position-radius", 0.0}
    });
    NucleonCommon bad_profile{bad_var_map};
  }(),
  std::domain_error);

  // test nucleon-nucleon cross section with random constituent substructure
  auto constituent_number = std::uniform_int_distribution<>{1, 10}(random::engine);
  auto constituent_width = .3 + .2*random::canonical<>();
  auto constituent_position_radius = .5*random::canonical<>();
  auto envelope_width_sq = sqr(constituent_position_radius) + sqr(constituent_width);

  // Coarse-ish p+p grid.
  auto grid_max = 3.;
  auto grid_step = 0.1;
  auto grid_nsteps = 60;
  auto dxdy = std::pow(2*grid_max/grid_nsteps, 2);

  // constituent var map
  var_map = make_var_map({
      {"grid-max", grid_max},
      {"grid-step", grid_step},
      {"fluctuation",   1e12},
      {"cross-section", xsec},
      {"nucleon-width", 0.5},
      {"constit-width", constituent_width},
      {"constit-number", constituent_number},
      {"constit-position-radius", constituent_position_radius}
  });

  NucleonCommon nc_constituent{var_map};

  // inelastic collision counter
  auto ncoll = 0;

  // measure the cross section with min bias p+p collisions
  for (auto i = 0; i < static_cast<int>(nev); ++i) {
    b = bmax * std::sqrt(random::canonical<>());
    A.sample_nucleons(.5*b);
    B.sample_nucleons(-.5*b);
    if (nc_constituent.participate(*A.begin(), *B.begin()))
      ++ncoll;
  }

  // calculate Monte Carlo cross section
  auto fraction = double(ncoll)/nev;
  auto area = math::double_constants::pi * bmax * bmax;
  auto mc_xsec = fraction * area;

  // assert a rather loose tolerance for the numerical cross section finder
  CHECK( mc_xsec/xsec == Approx(1).epsilon(.02) );

  // Compute TR grid the slow way -- switch the order of grid and nucleon loops.
  boost::multi_array<double, 2> TR{boost::extents[grid_nsteps][grid_nsteps]};

  for (auto i = 0; i < 1000; ++i) {
    A.sample_nucleons(0.);
    auto nucleon = *A.begin();

    while (!nc_constituent.participate(nucleon, nucleon)) {}

    for (auto iy = 0; iy < grid_nsteps; ++iy) {
      for (auto ix = 0; ix < grid_nsteps; ++ix) {
        auto x = (ix + .5) * 2 * grid_max/grid_nsteps - grid_max;
        auto y = (iy + .5) * 2 * grid_max/grid_nsteps - grid_max;
        auto dsq = x*x + y*y;
        auto norm = math::double_constants::one_div_two_pi / envelope_width_sq;
        auto nucleon_thickness = norm * std::exp(-0.5*dsq/envelope_width_sq);

        TR[iy][ix] +=
          (nucleon_thickness - nc_constituent.thickness(nucleon, x, y)) * dxdy / nev;
      }
    }
  }

  auto gaussian_error = std::accumulate(TR.origin(), TR.origin() + TR.num_elements(), 0.);

  CHECK( gaussian_error == Approx(0.).epsilon(.01).margin(.01) );


}
