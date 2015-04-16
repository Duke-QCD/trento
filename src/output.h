// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef OUTPUT_H
#define OUTPUT_H

#include <functional>
#include <utility>
#include <vector>

#include "fwd_decl.h"

namespace trento {

class Event;

///
class Output {
 public:
  ///
  Output(const VarMap& var_map);

  ///
  template <typename... Args>
  void operator()(Args&&... args) const;

 private:
  ///
  std::vector<std::function<void(int, double, const Event&)>> writers_;
};

template <typename... Args>
void Output::operator()(Args&&... args) const {
  for (const auto& write : writers_)
    write(std::forward<Args>(args)...);
}

}  // namespace trento

#endif  // OUTPUT_H
