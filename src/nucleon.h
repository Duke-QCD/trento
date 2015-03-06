// TRENTO: Reduced Thickness Event-by-event Nuclear Topology
// Copyright 2015 Jonah E. Bernhard, J. Scott Moreland
// MIT License

#ifndef NUCLEON_H
#define NUCLEON_H

#include <random>

namespace trento {

/// Nucleon class.
class Nucleon {
 public:
  /// Default constructor.
  Nucleon();

  ///
  static void set_width(double width);
  static void set_cross_sec(double cross_sec);
  static void set_fluct_shape(double fluct_shape);

  static double radius();

  ///
  void set_position(double x, double y);

  ///
  void participate_with(Nucleon& other);

  ///
  double thickness(double x, double y) const;

  ///
  bool is_participant() const;

 private:
  ///
  static double participation_prob(double distance_squared);

  ///
  static double sigma_nn();

  /// Truncate the thickness function at this number of nucleon widths.
  static constexpr double trunc_widths_ = 3.;

  /// Width of Gaussian thickness function.
  static double width_squared_;

  ///
  static double trunc_distance_squared_;

  ///
  static double sigma_gg_;

  ///
  static std::gamma_distribution<double> fluct_dist_;

  ///
  static std::uniform_real_distribution<double> unif_dist_;

  ///
  void set_participant();

  ///
  double x_, y_;

  ///
  double fluct_;

  ///
  bool participant_;
};

}  // namespace trento

#endif  // NUCLEON_H
