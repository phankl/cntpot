#ifndef bicubicSpline_header
#define bicubicSpline_header

#include <vector>

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

#include "cubicSpline.h"

class BicubicSpline {
private:
  std::vector<CubicSpline> xSplines;
  double xSpacing;
  double xStart;
  int points;

public:
  BicubicSpline();
  BicubicSpline(const std::vector<CubicSpline> &, double, double);

  double operator()(double, double) const;
};

#endif
