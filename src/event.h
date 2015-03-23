// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef EVENT_H
#define EVENT_H

#include <functional>
#include <map>
#include <vector>

#ifdef NDEBUG
#define BOOST_DISABLE_ASSERTS
#endif
#include <boost/multi_array.hpp>

#include "fwd_decl.h"

namespace trento {

///
struct Event {
  ///
  Event(int grid_size);

  /// The event number within the batch.
  int num;

  /// Impact parameter.
  double impact_param;

  /// Number of participants.
  int npart;

  /// Multiplicity (total entropy).
  double multiplicity;

  /// Eccentricity harmonics.
  std::map<int, double> eccentricity;

  /// Alias for a two-dimensional thickness grid.
  using Grid = boost::multi_array<double, 2>;

  /// Nuclear thickness grids TA, TB and reduced thickness grid TR.
  Grid TA, TB, TR;
};

///
using OutputFunctionVector = std::vector<std::function<void(const Event&)>>;

///
OutputFunctionVector create_output_functions(const VarMap& var_map);

}  // namespace trento

#endif  // EVENT_H
