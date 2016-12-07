// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef RANDOM_H
#define RANDOM_H

#include <limits>
#include <random>
#include <vector>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/erf.hpp>

#include "fwd_decl.h"
#include <iostream>

namespace trento { namespace random {

/// Mersenne Twister engine, 64-bit preset.
using Engine = std::mt19937_64;

/// Global variable defined in \c random.cxx.
extern Engine engine;

/// Helper function to easily generate random numbers in [0, 1).
template <typename RealType = double>
inline RealType canonical() {
  return std::generate_canonical
           <RealType, std::numeric_limits<RealType>::digits>
           (engine);
}

/// Sample a spherical polar angle cos(theta).
template <typename RealType = double>
inline double cos_theta() {
  return 2 * canonical<RealType>() - 1;
}

/// Sample a spherical azimuthal angle phi.
template <typename RealType = double>
inline double phi() {
  return math::constants::two_pi<RealType>() * canonical<RealType>();
}

/// Generates and cycles through uniform samples from a normal distribution
template <typename RealType = double>
class CyclicNormal {
 public:
  CyclicNormal(
      RealType mean = 0.0, RealType stddev = 1.0,
      std::size_t cache_size = 1000, std::size_t n_loops = 10
      );

  const RealType& operator()() {
    if (loop_count > n_loops_) {
      std::shuffle(cache_.begin(), cache_.end(), random::engine); 
      loop_count = 0;
    }
    if (iter_ == cache_.end()) {
      iter_ = cache_.begin();
      ++loop_count;
    }

    return *iter_++;
  }

  template <class Generator>
    const RealType& operator()(Generator&) {
      return operator()();
    }

 private:
  std::vector<RealType> cache_;
  std::size_t n_loops_;
  std::size_t loop_count = 0;
  typename std::vector<RealType>::const_iterator iter_;
};

template <typename RealType>
CyclicNormal<RealType>::CyclicNormal(
      RealType mean, RealType stddev,
      std::size_t cache_size, std::size_t n_loops
    )
  : cache_(cache_size),
    n_loops_(n_loops),
    iter_(cache_.begin()) {

  for (std::size_t i = 1; i < cache_size + 1; ++i) {
    auto z = -1. + 2.*double(i)/(cache_size + 1);
    auto normal = math::double_constants::root_two * boost::math::erf_inv(z);
    cache_[i-1] = mean + stddev * normal;
  }

  std::shuffle(cache_.begin(), cache_.end(), random::engine); 
}

}}  // namespace trento::random


#endif  // RANDOM_H
