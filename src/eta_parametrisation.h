#ifndef ETA_H
#define ETA_H

#include <math.h>
#include <random>
#include <iostream>
#include <vector>

double const sqrt2 = sqrt(2);
///----------------------------------------------------------------------------------------------///
/// feel free to try you own parametrization of y-mean. y-std and y-skew as function of Ta and Tb///
///----------------------------------------------------------------------------------------------///
double inline mean_function(double ta, double tb){
  if (ta == 0)
    return -10.0;
  if (tb == 0)
    return 10.0;
  return 0.5*std::log(ta/tb);
}

double inline std_function(double ta, double tb){
  return 1.0;
}

double inline skew_function(double ta, double tb){
  return (ta-tb)/(ta+tb);
}
///----------------------------------------------------------------------------------------------///


double inline skew_normal_function(double eta, double xi, double omega, double alpha){
  double x = (eta-xi)/omega;
  return exp(-0.5*x*x)*(1.+erf(alpha*x/sqrt2));
}

/// A fast pseudorapidity to rapidity transformer using pretabulated values
/// Output both y(eta) and dy/deta(eta)
class fast_eta2y{
private:
  double J_, etamax_, deta_;
  int neta_;
  std::vector<double> y_;
  std::vector<double> dydeta_;
public:
  fast_eta2y(double J, double etamax, double deta){
    J_ = J;
    etamax_ = etamax;
    deta_ = deta;
    neta_ = std::ceil(2.*etamax_/deta_);
    y_.resize(neta_);
    dydeta_.resize(neta_);
    for (int i=0; i<neta_; i++){
      double eta = -etamax_ + i*deta_;
      double Jsh = J_*std::sinh(eta);
      double sq = std::sqrt(1.+Jsh*Jsh);
      y_[i] = std::log(sq+Jsh);
      dydeta_[i] = J*std::cosh(eta)/sq;
    }
  };
  
  double rapidity(double eta){
    double x = (eta + etamax_)/deta_;
    int i = std::floor(x);
    if (i<0 || i>=neta_-2)
      return eta;
    double ri = x - i;
    return y_[i]*(1.-ri) + y_[i+1]*ri; 
  };
  double Jacobi(double eta){
    double x = (eta + etamax_)/deta_;
    int i = std::floor(x);
    if (i<0 || i>=neta_-2)
      return 1.0;
    double ri = x - i;
    return dydeta_[i]*(1.-ri) + dydeta_[i+1]*ri;
  };

};

#endif
