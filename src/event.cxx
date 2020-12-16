// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#include "event.h"

#include <algorithm>
#include <cmath>

#include <boost/program_options/variables_map.hpp>
#include "nucleus.h"
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/special_functions/hypergeometric_1F1.hpp>
#include <boost/math/special_functions/erf.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
using boost::math::quadrature::trapezoidal;

namespace trento {

namespace {
constexpr double TINY = 1e-12;
constexpr double Mproton = 0.938272;
constexpr double Mneutron = 0.939565;
// Generalized mean for p > 0.
// M_p(a, b) = (1/2*(a^p + b^p))^(1/p)
double pmean(double p, double a, double b) {
  if (std::abs(p) < 1e-6) return std::sqrt(a*b);
  else{
    if (p>0.)
      return std::pow(.5*(std::pow(a,p)+std::pow(b,p)), 1./p);
    else{
      if (a<TINY || b<TINY) return 0.;
      else return std::pow(.5*(std::pow(a,p)+std::pow(b,p)), 1./p);
    }
  }
}



inline double EoverM_to_absrap(double p0_over_mT){
    if (p0_over_mT<=1.) return 0.;
    else return std::log( p0_over_mT + std::sqrt(std::pow(p0_over_mT,2)-1.) );
}

inline double EandP_to_absrap(double E, double P){
    if (E<(1.+TINY)*P) return 0.;
    else return .5*std::log((E+P)/(E-P));
}

inline double PandM_to_rap(double pz_over_mT){
    return std::log(std::sqrt(std::pow(pz_over_mT,2)+1.) + pz_over_mT);
}

}  // unnamed namespace

// Determine the grid parameters like so:
//   1. Read and set step size from the configuration.
//   2. Read grid max from the config, then set the number of steps as
//      nsteps = ceil(2*max/step).
//   3. Set the actual grid max as max = nsteps*step/2.  Hence if the step size
//      does not evenly divide the config max, the actual max will be marginally
//      larger (by at most one step size).
Event::Event(const VarMap& var_map)
    : mid_power_(var_map["mid-power"].as<double>()),
      mid_norm_(var_map["mid-norm"].as<double>()),
      shape_alpha_(var_map["shape-alpha"].as<double>()),
      shape_beta_(var_map["shape-beta"].as<double>()),
      kT_min_(var_map["kT-min"].as<double>()),
      sqrts_(var_map["sqrts"].as<double>()),
      nucleon_pabs_(std::sqrt(.25*sqrts_*sqrts_ - Mproton*Mproton)),
      eta_max_( std::asinh(.5*sqrts_/kT_min_) ),
      nsteps_etas_(var_map["nsteps-etas"].as<int>()),
      detas_(2.*eta_max_/nsteps_etas_),
      dxy_(var_map["grid-step"].as<double>()),
      nsteps_(std::ceil(2.*var_map["grid-max"].as<double>()/dxy_)),
      xymax_(.5*nsteps_*dxy_),
      TA_(boost::extents[nsteps_][nsteps_]),
      TB_(boost::extents[nsteps_][nsteps_]),
      Density_(boost::extents[nsteps_etas_][nsteps_][nsteps_]){
  Pplus_ = (sqrts_/2. + nucleon_pabs_)/2.;
  Pminus_ = (sqrts_/2. - nucleon_pabs_)/2.;
  // BeamRap width
  L_ = std::log(sqrts_/kT_min_);
  //Landau width
  LL_ = std::sqrt(L_);
  double sqrt2LL = std::sqrt(2.)*LL_;
  // Midrapidity norm
  norm_trento_ = mid_norm_ * Mproton * std::pow(sqrts_/Mproton, mid_power_);
  // Fragmentation region normalization
  auto f1 = [this](double x){
      return std::exp(-x*x/2./this->L_)*std::cosh(x)*std::pow(1.-std::pow(x/this->L_,4),2 );
  };
  double F1 = norm_trento_ * trapezoidal(f1, -L_, L_);
  // Energy fraction  needed to be depositied from the fragmentation region
  xloss_ = F1/sqrts_;
  if (xloss_>1) {
      std::cout << "central fireball too large!" << std::endl;
      exit(-1);
  }
  // Fragmentation region normalization
  auto f2 = [this](double x){
      return std::pow(eta_max_-x, this->shape_alpha_)
            *std::exp(-(eta_max_-x)*(this->shape_beta_+2.));
  };
  Nab_ = trapezoidal(f2, 0.0, eta_max_);
  // Cut fragmentation region using xloss_
  auto low = 0.0, high=100.;
  boost::uintmax_t max_iter = 1000;
  math::tools::eps_tolerance<double> tol{
     (std::numeric_limits<double>::digits * 3) / 4
  };
  try {
      auto result = math::tools::toms748_solve(
        [this](double x) {
            auto f3 = [this, x](double y){
                 return std::pow(this->eta_max_-y, this->shape_alpha_)
                       *std::exp(-(this->eta_max_-y)*(this->shape_beta_+2.)
                       -x*std::pow(this->eta_max_-y,2.) 
                       );
                };
        return 1.-trapezoidal(f3, 0.0, eta_max_)/this->Nab_ - this->xloss_;
      },
      low, high, tol, max_iter
      );
      shape_gamma_ = .5*(result.first + result.second);
  }
  catch (const std::domain_error&) {
     std::cout << "gamma = " << shape_gamma_
               << ", sqrts = " << sqrts_<<std::endl;
    throw std::domain_error{"unable to fit shape-gamma?"};
  }
  //std::cout << xloss_ << " " << shape_gamma_ << std::endl;
}

void Event::compute(const Nucleus& nucleusA, const Nucleus& nucleusB,
                    const NucleonCommon& nucleon_common) {
  // Reset npart; compute_nuclear_thickness() increments it.
  npart_ = 0;
  npartA_ = 0;
  npartB_ = 0;
  ixcm_.resize(nsteps_etas_);
  iycm_.resize(nsteps_etas_);
  dET_detas_.resize(nsteps_etas_);
  for (auto & it : ixcm_) it = 0.;
  for (auto & it : iycm_) it = 0.;
  for (auto & it : dET_detas_) it = 0.;

  for (const auto& nucleon : nucleusA) {
    if (nucleon.is_participant()) {
        npartA_++;
        npart_++;
    }
  }
  for (const auto& nucleon : nucleusB) {
    if (nucleon.is_participant()) {
        npartB_++;
        npart_++;
    }
  }
  compute_nuclear_thickness(nucleusA, nucleon_common, TA_);
  compute_nuclear_thickness(nucleusB, nucleon_common, TB_);
  compute_density();
  compute_observables();
}

namespace {

// Limit a value to a range.
// Used below to constrain grid indices.
template <typename T>
inline const T& clip(const T& value, const T& min, const T& max) {
  if (value < min)
    return min;
  if (value > max)
    return max;
  return value;
}

}  // unnamed namespace


void Event::compute_nuclear_thickness(
    const Nucleus& nucleus, const NucleonCommon& nucleon_common, Grid& TX) {
  // Construct the thickness grid by looping over participants and adding each
  // to a small subgrid within its radius.  Compared to the other possibility
  // (grid cells as the outer loop and participants as the inner loop), this
  // reduces the number of required distance-squared calculations by a factor of
  // ~20 (depending on the nucleon size).  The Event unit test verifies that the
  // two methods agree.

  // Wipe grid with zeros.
  std::fill(TX.origin(), TX.origin() + TX.num_elements(), 0.);

  // Deposit each participant onto the grid.
  for (const auto& nucleon : nucleus) {
    if (!nucleon.is_participant())
      continue;

    // Get nucleon subgrid boundary {xmin, xmax, ymin, ymax}.
    const auto boundary = nucleon_common.boundary(nucleon);

    // Determine min & max indices of nucleon subgrid.
    int ixmin = clip(static_cast<int>((boundary[0]+xymax_)/dxy_), 0, nsteps_-1);
    int ixmax = clip(static_cast<int>((boundary[1]+xymax_)/dxy_), 0, nsteps_-1);
    int iymin = clip(static_cast<int>((boundary[2]+xymax_)/dxy_), 0, nsteps_-1);
    int iymax = clip(static_cast<int>((boundary[3]+xymax_)/dxy_), 0, nsteps_-1);

    // Add profile to grid.
    for (auto iy = iymin; iy <= iymax; ++iy) {
      for (auto ix = ixmin; ix <= ixmax; ++ix) {
        TX[iy][ix] += nucleon_common.thickness(
          nucleon, (ix+.5)*dxy_ - xymax_, (iy+.5)*dxy_ - xymax_
        );
      }
    }
  }
}



void Event::compute_density() {
  std::fill(Density_.origin(), Density_.origin()+Density_.num_elements(), 0.);  

  for (int iy = 0; iy < nsteps_; ++iy) {
    for (int ix = 0; ix < nsteps_; ++ix) {
      double ta = TA_[iy][ix];
      double tb = TB_[iy][ix];
      if (ta<TINY || tb<TINY) continue;
      double TR = pmean(0., ta, tb);
      double etacm = .5*std::log((ta*Pplus_+tb*Pminus_)
                                /(tb*Pplus_+ta*Pminus_));
      for (int iz = 0; iz < nsteps_etas_; ++iz) {
        double etas = -eta_max_ + (iz+.5) * detas_; 
        double absetas = std::abs(etas);
        double TaTb = (ta*(etas>0)+tb*(etas<=0));
        double mid_profile = 
                  std::exp(-std::pow(etas-etacm,2)/2./L_)
                  * std::pow(1.-std::pow((etas-etacm)/L_, 4), 2) 
                  * (std::abs(etas-etacm)<L_);
        double x = std::abs(2*kT_min_*std::sinh(etas)/sqrts_);
        double fb_profile = std::pow(eta_max_-absetas, shape_alpha_)
                        *std::exp(
                            - (shape_beta_+1.) * (eta_max_-absetas)
                            - shape_gamma_*std::pow(eta_max_-absetas,2)
                         )   
                        /Nab_
         *(1.-std::exp(-std::pow(shape_gamma_*std::abs(etas)/eta_max_,2)));
         Density_[iz][iy][ix] = TR*norm_trento_*mid_profile
                             + TaTb*kT_min_*fb_profile;
      }
    }
  }
  // compute multiplicity and com
  for (int iz = 0; iz < nsteps_etas_; ++iz) { 
    double sum = 0.;
    for (int iy = 0; iy < nsteps_; ++iy) {
      for (int ix = 0; ix < nsteps_; ++ix) {
        double dE = Density_[iz][iy][ix];
        if (dE<TINY) continue;
        ixcm_[iz] += dxy_*ix*dE;
        iycm_[iz] += dxy_*iy*dE;
        sum += dE;
      }
    }
    ixcm_[iz] /= sum;
    iycm_[iz] /= sum;
    dET_detas_[iz] = sum * dxy_ * dxy_;
  }
}

void Event::compute_observables() {
  // Compute eccentricity.

  // Simple helper class for use in the following loop.
  struct EccentricityAccumulator {
    std::vector<double> re; // real part
    std::vector<double> im; // imaginary part
    std::vector<double> wt; // weight
    void init(int N){
         re.resize(N);
         im.resize(N);
         wt.resize(N);
         std::fill(re.begin(), re.end(), 0.);
         std::fill(im.begin(), im.end(), 0.);
         std::fill(wt.begin(), wt.end(), 0.);
    }
    std::vector<double> magnitudes() { 
         std::vector<double> enabs;
         enabs.resize(re.size());
         for (size_t i=0; i<enabs.size(); i++)
             enabs[i] = std::sqrt(re[i]*re[i] + im[i]*im[i]) 
                      / std::fmax(wt[i], TINY);
         return enabs;
    }
    std::vector<double> angles() { 
         std::vector<double> enphi;
         enphi.resize(re.size());
         for (size_t i=0; i<enphi.size(); i++)
             enphi[i] = std::atan2(im[i], re[i]);
         return enphi;
    }
  } e2, e3, e4, e5;
  e2.init(nsteps_etas_);
  e3.init(nsteps_etas_);
  e4.init(nsteps_etas_);
  e5.init(nsteps_etas_);

  for (size_t iz = 0; iz < nsteps_etas_; ++iz){
    for (size_t iy = 0; iy < nsteps_; ++iy) {
      for (size_t ix = 0; ix < nsteps_; ++ix) {
        const auto& t = Density_[iz][iy][ix];
        if (t < TINY)
          continue;

        // Compute (x, y) relative to the CM and cache powers of x, y, r.
        auto x = ix*dxy_ - ixcm_[iz];
        auto x2 = x*x;
        auto x3 = x2*x;
        auto x4 = x2*x2;

        auto y = iy*dxy_ - iycm_[iz];
        auto y2 = y*y;
        auto y3 = y2*y;
        auto y4 = y2*y2;

        auto r2 = x2 + y2;
        auto r = std::sqrt(r2);
        auto r4 = r2*r2;

        auto xy = x*y;
        auto x2y2 = x2*y2;

        // The eccentricity harmonics are weighted averages of r^n*exp(i*n*phi)
        // over the entropy profile (reduced thickness). 
        // The naive way to compute
        // exp(i*n*phi) at a given (x, y) point is essentially:
        //
        //   phi = arctan2(y, x)
        //   real = cos(n*phi)
        //   imag = sin(n*phi)
        //
        // However this implementation uses three unnecessary trig functions; a
        // much faster method is to express 
        // the cos and sin directly in terms of x
        // and y.  For example, it is trivial to show 
        // (by drawing a triangle and
        // using rudimentary trig) that
        //
        //   cos(arctan2(y, x)) = x/r = x/sqrt(x^2 + y^2)
        //   sin(arctan2(y, x)) = y/r = x/sqrt(x^2 + y^2)
        //
        // This is easily generalized to cos and sin of (n*phi) by invoking the
        // multiple angle formula, e.g. sin(2x) = 2sin(x)cos(x), and hence
        //
        //   sin(2*arctan2(y, x)) = 2*sin(arctan2(y, x))*cos(arctan2(y, x))
        //                        = 2*x*y / r^2
        //
        // Which not only eliminates the trig functions, but also naturally
        // cancels the r^2 weight.  This cancellation occurs for all n.
        //
        // The Event unit test verifies that the two methods agree.
        e2.re[iz] += t * (y2 - x2);
        e2.im[iz] += t * 2.*xy;
        e2.wt[iz] += t * r2;

        e3.re[iz] += t * (y3 - 3.*y*x2);
        e3.im[iz] += t * (3.*x*y2 - x3);
        e3.wt[iz] += t * r2*r;

        e4.re[iz] += t * (x4 + y4 - 6.*x2y2);
        e4.im[iz] += t * 4.*xy*(y2 - x2);
        e4.wt[iz] += t * r4;

        e5.re[iz] += t * y*(5.*x4 - 10.*x2y2 + y4);
        e5.im[iz] += t * x*(x4 - 10.*x2y2 + 5.*y4);
        e5.wt[iz] += t * r4*r;
      }
    }
  }

  ecc_mag_[2] = e2.magnitudes();
  ecc_ang_[2] = e2.angles();

  ecc_mag_[3] = e3.magnitudes();
  ecc_ang_[3] = e3.angles();

  ecc_mag_[4] = e4.magnitudes();
  ecc_ang_[4] = e4.angles();

  ecc_mag_[5] = e5.magnitudes();
  ecc_ang_[5] = e5.angles();
}

}  // namespace trento
