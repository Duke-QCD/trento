// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "nucleon.h"

#include <cmath>

#include "random.h"

namespace trento {

Nucleon::Nucleon() : x_(0.), y_(0.), fluct_(0.), participant_(false) {}

void Nucleon::set_width(double width) {
  width_squared_ = width*width;
  trunc_distance_squared_ = std::pow(trunc_widths_ * width, 2);
}

void Nucleon::set_cross_sec(double cross_sec) {
  sigma_gg_ = .5*cross_sec/width_squared_ - 0.577216;
}

void Nucleon::set_fluct_shape(double fluct_shape) {
  fluct_dist_.param(
    std::gamma_distribution<double>::param_type{fluct_shape, 1./fluct_shape});
}

double Nucleon::radius() {
  return trunc_widths_ * std::sqrt(width_squared_);
}

void Nucleon::set_position(double x, double y) {
  x_ = x;
  y_ = y;
  participant_ = false;
}

double Nucleon::thickness(double x, double y) const {
  x -= x_;
  y -= y_;

  auto distance_squared = x*x + y*y;
  if (distance_squared > trunc_distance_squared_)
    return 0.;

  return fluct_ / width_squared_ *
         std::exp(-.5*distance_squared/width_squared_);
}

bool Nucleon::is_participant() const {
  return participant_;
}

void Nucleon::set_participant() {
  participant_ = true;
  fluct_ = fluct_dist_(random_engine);
}

void Nucleon::participate_with(Nucleon& other) {
  auto dx = x_ - other.x_;
  auto dy = y_ - other.y_;
  auto distance_squared = dx*dx + dy*dy;

  auto prob = participation_prob(distance_squared);

  // If the probability is zero (less than epsilon), the first condition will
  // fail and no time will be wasted generating a random number.
  if (prob > 1e-12 && prob > unif_dist_(random_engine)) {
    // Participation is commutative.
    set_participant();
    other.set_participant();
  }
}

/// \f$ P(b) = 1 - \exp(-\sigma_{gg} \int dx \, dy \, T_A T_B) \f$
double Nucleon::participation_prob(double distance_squared) {
  if (distance_squared > 4.*trunc_distance_squared_) {
    return 0.;
  }
  else {
    return -std::expm1(
        -std::exp(sigma_gg_ - .25*distance_squared/width_squared_));
  }
}

double Nucleon::sigma_nn() {
  auto nsteps = 1000;
  auto db = 2.*std::sqrt(trunc_distance_squared_)/(nsteps - 1);

  auto sigma_nn = 0.;
  for (auto i = 0; i < nsteps; ++i) {
    auto b = (i+.5)*db;
    sigma_nn += b*participation_prob(b*b);
  }
  sigma_nn *= 2.*M_PI*db;

  return sigma_nn;
}

double Nucleon::width_squared_ = 1.;
double Nucleon::trunc_distance_squared_ = 1.;
double Nucleon::sigma_gg_ = 1.;
std::gamma_distribution<double> Nucleon::fluct_dist_{1., 1.};
std::uniform_real_distribution<double> Nucleon::unif_dist_{0., 1.};

}  // namespace trento
