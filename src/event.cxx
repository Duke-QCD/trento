// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License


#include <algorithm>
#include <cmath>
#include <boost/program_options/variables_map.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include "event.h"
#include "nucleus.h"

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
      flatness_(var_map["flatness"].as<double>()),
      shape_alpha_(var_map["shape-alpha"].as<double>()),
      shape_beta_(var_map["shape-beta"].as<double>()),
      kT_min_(var_map["kT-min"].as<double>()),
      sqrts_(var_map["sqrts"].as<double>()),
      nucleon_pabs_(std::sqrt(.25*sqrts_*sqrts_ - Mproton*Mproton)),
      // eta_max_(sqrts, kTmin) is the physically kinetiac range considered
      eta_max_( std::asinh(.5*sqrts_/kT_min_) ),  
      // eta_grid_max_(sqrts) is the range of the computing grid
      // it also has a physical meaning though
      eta_grid_max_( std::asinh(.5*sqrts_/0.2) ),
      nsteps_etas_(var_map["nsteps-etas"].as<int>()),
      detas_(2.*eta_grid_max_/nsteps_etas_),
      dxy_(var_map["grid-step"].as<double>()),
      nsteps_(std::ceil(2.*var_map["grid-max"].as<double>()/dxy_)),
      xymax_(.5*nsteps_*dxy_),
      // Thickness function to deposite into the middle fireball
      TA_(boost::extents[nsteps_][nsteps_]),
      TB_(boost::extents[nsteps_][nsteps_]),
      // Thickness function fragments at large x
      FA_(boost::extents[nsteps_][nsteps_]),
      FB_(boost::extents[nsteps_][nsteps_]),
      // 3D initial energy density at tau=0+
      Density_(boost::extents[nsteps_etas_][nsteps_][nsteps_])
  {
  // light-cone moemtnum compoents of the projectile (per nucleon)
  // P+ = (E+|Pz|)/2 = (sqrts/2.+pabs)/2. = Mproton * exp(ybeam)
  Pplus_ = (sqrts_/2. + nucleon_pabs_)/2.;
  // P- = (E-|Pz|)/2 = (sqrts/2.-pabs)/2. = Mproton * exp(-ybeam)
  Pminus_ = (sqrts_/2. - nucleon_pabs_)/2.;
  // Midrapidity norm
  norm_trento_ = mid_norm_ * Mproton * std::pow(sqrts_/Mproton, mid_power_);
  // Compute central fireball normalization Ncentral
  auto f1 = [this](double eta){
      return this->central_profile(eta)*std::cosh(eta);
  };
  double Ncentral = norm_trento_ * trapezoidal(f1, -eta_max_, eta_max_);
  // Energy fraction of the projectile/target to be deposited in 
  // the central fireball, cannot be larger than one!
  xloss_ = Ncentral/sqrts_; 
  if (xloss_>1.-TINY) {
      std::cout << "central fireball too large!" << std::endl;
      exit(-1);
  }
  // Normalziation for the fragmentation region, Nfrag_
  auto f2 = [this](double x){
      return std::pow(eta_max_-x, this->shape_alpha_)
            *std::exp(-(eta_max_-x)*(this->shape_beta_+2.));
  };
  Nfrag_ = trapezoidal(f2, -eta_max_, eta_max_);
}

void Event::compute(const Nucleus& nucleusA, 
                    const Nucleus& nucleusB,
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
  // compute and decompose thickness function into two parts
  // Toriginal = T_ + F_
  // T_: the fraction to be deposited into central fireball
  // F_: the fraction remained at large x limiting fragmentation region
  compute_nuclear_thickness(nucleusA, nucleon_common, TA_, FA_);
  compute_nuclear_thickness(nucleusB, nucleon_common, TB_, FB_);
  // compute 3D energy density using TA, TB, FA, FB
  compute_density();
  // compute eccentricity 
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
    const Nucleus& nucleus, const NucleonCommon& nucleon_common, 
    Grid& TX, Grid& FX) {
  // Construct the thickness grid by looping over participants and adding each
  // to a small subgrid within its radius.  Compared to the other possibility
  // (grid cells as the outer loop and participants as the inner loop), this
  // reduces the number of required distance-squared calculations by a factor of
  // ~20 (depending on the nucleon size).  The Event unit test verifies that the
  // two methods agree.

  // Wipe grid with zeros.
  std::fill(TX.origin(), TX.origin() + TX.num_elements(), 0.);
  std::fill(FX.origin(), FX.origin() + FX.num_elements(), 0.);


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
        FX[iy][ix] += nucleon_common.fragmentation(
          nucleon, (ix+.5)*dxy_ - xymax_, (iy+.5)*dxy_ - xymax_
        );
      }
    }
  }
}



void Event::compute_density() {
  std::fill(Density_.origin(), 
            Density_.origin()+Density_.num_elements(), 0.);  

  for (int iy = 0; iy < nsteps_; ++iy) {
    for (int ix = 0; ix < nsteps_; ++ix) {
      double ta = TA_[iy][ix], fa = FA_[iy][ix];
      double tb = TB_[iy][ix], fb = FB_[iy][ix];
      if (ta<TINY && tb<TINY && fa < TINY && fb < TINY) continue;
      double TR = pmean(0., ta, tb)/xloss_;
      double etacm = center_of_mass_eta(ta, tb);
      for (int iz = 0; iz < nsteps_etas_; ++iz) {
        double etas = - eta_grid_max_ + (iz+.5) * detas_; 
        double e_central = (std::abs(etas-etacm)<eta_grid_max_)? 
                (norm_trento_*TR*central_profile(etas-etacm)) : 0.;
        double xA = std::min(std::max(std::exp(-eta_max_+etas), TINY),1.);
        double xB = std::min(std::max(std::exp(-eta_max_-etas), TINY),1.);
        double e_fb = kT_min_*(ta*frag_profile(xA) + tb*frag_profile(xB));
        Density_[iz][iy][ix] = e_central + e_fb;
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
        ixcm_[iz] += ix*dE;
        iycm_[iz] += iy*dE;
        sum += dE;
      }
    }
    ixcm_[iz] *= dxy_/sum;
    iycm_[iz] *= dxy_/sum;
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
    int order;
    void init(int N, int order_){
         re.resize(N);
         im.resize(N);
         wt.resize(N);
         std::fill(re.begin(), re.end(), 0.);
         std::fill(im.begin(), im.end(), 0.);
         std::fill(wt.begin(), wt.end(), 0.);
         order = order_;
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
             enphi[i] = std::atan2(im[i], re[i])/order;
         return enphi;
    }
  } e2, e3, e4, e5;
  e2.init(nsteps_etas_,2);
  e3.init(nsteps_etas_,3);
  e4.init(nsteps_etas_,4);
  e5.init(nsteps_etas_,5);

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
