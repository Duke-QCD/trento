#ifndef ETA_H
#define ETA_H

#include <cmath>
#include <vector>
#include <boost/math/constants/constants.hpp>

#include "fwd_decl.h"
#include <iostream>


constexpr double TINY = 1e-12;

///----------------------------------------------------------------------------------------------///
/// feel free to try you own parametrization of y-mean. y-std and y-skew as function of Ta and Tb///
///----------------------------------------------------------------------------------------------///
double inline mean_function(double ta, double tb, double exp_ybeam){
  if (ta < TINY && tb < TINY) return 0;
  return 0.5 * std::log((ta*exp_ybeam + tb/exp_ybeam) / (ta/exp_ybeam + tb*exp_ybeam));
}

double inline std_function(double ta, double tb){
  (void)ta;
  (void)tb;
  return 1.;
}

double inline skew_function(double ta, double tb){
  if (ta < TINY && tb < TINY) return 0.;
  return (ta - tb)/(ta + tb);
}
///----------------------------------------------------------------------------------------------///


double inline skew_normal_function(double eta, double mean, double width, double skew) {
  using math::double_constants::half_pi;
  using math::double_constants::one_div_root_two_pi;
  using math::double_constants::root_two;

  if (skew < 0) {
    skew *= -1.;
    mean *= -1.;
    eta *= -1.;
  }

  skew = std::fmin(skew, 0.99);

  auto zeta = std::pow(skew/(2. - half_pi), 2./3.);
  auto delta2 = half_pi*zeta/(1. + zeta); 
  auto omega = width/std::sqrt(1. - delta2/half_pi);
  auto xi = mean - omega*std::sqrt(delta2/half_pi);
  auto x = (eta - xi) / omega;
  auto alpha = 1./std::sqrt(1./delta2 - 1.);
  return one_div_root_two_pi/omega * std::exp(-x*x/2.)
    * (1. + std::erf(alpha * x / root_two));
}

double inline exp_gaussian(double x, double mean, double width, double skew) {
  if (skew < 0) {
    skew *= -1.;
    mean *= -1.;
    x *= -1.;
  }
  skew = std::fmin(skew, 1.9);

  auto g = pow(skew/2., 1./3.);
  auto mu = mean - width*g;
  auto sigma2 = width*width*(1. - g*g);
  auto lambda = 1./(width*g);
  return lambda/2. * std::exp(lambda/2. * (2.*mu + lambda*sigma2 - 2.*x))
    * std::erfc((mu + lambda*sigma2 - x)/std::sqrt(2.*sigma2));
}

/// A fast pseudorapidity to rapidity transformer using pretabulated values
/// Output both y(eta) and dy/deta(eta)
class fast_eta2y {
 private:
  double etamax_;
  double deta_;
  std::size_t neta_;
  std::vector<double> y_;
  std::vector<double> dydeta_;
 public:
  fast_eta2y(double J, double etamax, double deta)
      : etamax_(etamax),
        deta_(deta),
        neta_(std::ceil(2.*etamax_/deta_)),
        y_(neta_, 0.),
        dydeta_(neta_, 0.) {

    for (std::size_t ieta = 0; ieta < neta_; ++ieta) {
      double eta = -etamax_ + ieta*deta_;
      double Jsh = J*std::sinh(eta);
      double sq = std::sqrt(1. + Jsh*Jsh);
      y_[ieta] = std::log(sq + Jsh);
      dydeta_[ieta] = J*std::cosh(eta)/sq;
    }
  }
  
  double rapidity(double eta){
    double steps = (eta + etamax_)/deta_;
    double xi = std::fmod(steps, 1.);
    std::size_t index = std::floor(steps);
    return y_[index]*(1. - xi) + y_[index+1]*xi;
  }

  double Jacobian(double eta){
    double steps = (eta + etamax_)/deta_;
    double xi = std::fmod(steps, 1.);
    std::size_t index = std::floor(steps);
    return dydeta_[index]*(1. - xi) + dydeta_[index+1]*xi;
  }

};

#endif
