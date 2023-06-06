#ifndef bicubicSpline_header
#define bicubicSpline_header

#include <cmath>
#include <vector>

using namespace std;

class BicubicSpline {
private:
  int size;
  double xstart, ystart, dx, dy;
  vector<vector<vector<vector<double>>>> coeff;

  double spline(double, double) const;

  vector<vector<vector<vector<double>>>>
  splineCoeff(const vector<vector<double>> &);

public:
  BicubicSpline();
  BicubicSpline(const vector<vector<double>> &, double, double, double, double);

  double operator()(double, double) const;
};

#endif
