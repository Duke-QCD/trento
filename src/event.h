// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef EVENT_H
#define EVENT_H

#include <memory>

#include "fwd_decl.h"

namespace trento {

///
class Event {
 public:
  ///
  explicit Event(const VarMap& var_map);

  ///
  ~Event();

  ///
  void collide();

 private:
  ///
  std::unique_ptr<NucleusBase> nucA_, nucB_;

  ///
  double b_min_, b_max_;
};

}  // namespace trento

#endif  // EVENT_H
