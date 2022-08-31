#ifndef cubicSpline_header
#define cubicSpline_header

#include <vector>

#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>

class CubicSpline {
private:
  boost::math::interpolators::cardinal_cubic_b_spline<double> spline;

public:
  CubicSpline();
  CubicSpline(const std::vector<double> &, double, double);

  double operator()(double) const;
};

#endif
