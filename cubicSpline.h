#ifndef cubicSpline_header
#define cubicSpline_header

#include <cmath>
#include <vector>

using namespace std;

class CubicSpline {
private:
  int size;
  double xstart, dx;
  vector<vector<double>> coeff;

  double spline(double) const;

  vector<vector<double>> splineCoeff(const vector<double> &);

public:
  CubicSpline();
  CubicSpline(const vector<double> &, double, double);

  double operator()(double) const;
};

#endif
