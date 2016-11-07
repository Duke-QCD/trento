#ifndef ETA_H
#define ETA_H

#include <cmath>
#include <vector>

double const sqrt2 = std::sqrt(2);
constexpr double TINY = 1e-12;

///----------------------------------------------------------------------------------------------///
/// feel free to try you own parametrization of y-mean. y-std and y-skew as function of Ta and Tb///
///----------------------------------------------------------------------------------------------///
double inline mean_function(double ta, double tb, double exp_ybeam){
  if (ta < TINY && tb < TINY) return 0.;
  return 0.5 * std::log((ta*exp_ybeam + tb/exp_ybeam) / std::max(ta/exp_ybeam + tb*exp_ybeam, TINY));
}

double inline std_function(double ta, double tb){
  (void)ta;
  (void)tb;
  return 1.;
}

double inline skew_function(double ta, double tb){
  return (ta - tb)/std::max(ta + tb, TINY);
}
///----------------------------------------------------------------------------------------------///


double inline skew_normal_function(double eta, double xi, double omega, double alpha){
  double x = (eta - xi) / omega;
  return std::exp(-x*x/2.) * (1. + std::erf(alpha * x / sqrt2));
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
