#include "cubicSpline.h"

using namespace std;

CubicSpline::CubicSpline() {}

CubicSpline::CubicSpline(const vector<double> &data, double xstart, double dx)
    : xstart(xstart), dx(dx), size(data.size()), coeff(splineCoeff(data)) {}

vector<vector<double>> CubicSpline::splineCoeff(const vector<double> &data) {
  vector<vector<double>> g(size, vector<double>(4));

  vector<vector<double>> b(size, vector<double>(size));
  vector<double> p(size);
  vector<double> bprime(size);
  vector<double> dprime(size);

  double dx_inv = 1.0 / dx;
  double dxsq_inv = dx_inv * dx_inv;
  double dxcb_inv = dx_inv * dxsq_inv;

  double ax[4][4] = {{1, 0, 0, 0},
                     {0, 1, 0, 0},
                     {-3 * dxsq_inv, -2 * dx_inv, 3 * dxsq_inv, -dx_inv},
                     {2 * dxcb_inv, dxsq_inv, -2 * dxcb_inv, dxsq_inv}};

  // compute finite difference derivatives at boundaries

  p[0] = (data[1] - data[0]) * dx_inv;
  p[size - 1] = (data[size - 1] - data[size - 2]) * dx_inv;

  // compute derivatives inside domain

  for (int i = 1; i < size - 1; i++) {
    if (i > 1)
      b[i][i - 1] = dx;
    b[i][i] = 4 * dx;
    if (i < size - 2)
      b[i][i + 1] = dx;
  }
  bprime[1] = b[1][1];
  for (int i = 2; i < size - 1; i++)
    bprime[i] = b[i][i] - b[i][i - 1] * b[i - 1][i] / bprime[i - 1];

  for (int i = 1; i < size - 1; i++) {
    double d = 3 * (data[i + 1] - data[i - 1]);
    if (i == 1)
      d -= dx * p[i - 1];
    if (i == size - 2)
      d -= dx * p[i + 1];
    dprime[i] = d;
    if (i != 1)
      dprime[i] -= b[i][i - 1] * dprime[i - 1] / bprime[i - 1];
  }

  p[size - 2] = dprime[size - 2] / bprime[size - 2];
  for (int i = size - 3; i > 0; i--)
    p[i] = (dprime[i] - b[i][i + 1] * p[i + 1]) / bprime[i];

  // compute spline coefficients

  for (int i = 1; i < size; i++) {
    for (int j = 0; j < 4; j++)
      g[i][j] = 0;

    double k[4] = {data[i - 1], p[i - 1], data[i], p[i]};

    for (int j = 0; j < 4; j++)
      for (int l = 0; l < 4; l++)
        g[i][j] += ax[j][l] * k[l];
  }

  return g;
}

double CubicSpline::spline(double x) const {
  int i = ceil((x - xstart) / dx);

  // linear extrapolation

  if (i < 1) {
    return coeff[1][0] + coeff[1][1] * (x - xstart);

    // constant extrapolation

  } else if (i > size - 1) {
    i = size - 1;
    x = xstart + (size - 1) * dx;
  }

  // cubic interpolation

  double xlo = xstart + (i - 1) * dx;
  double xbar = x - xlo;

  return coeff[i][0] +
         xbar * (coeff[i][1] + xbar * (coeff[i][2] + xbar * coeff[i][3]));
}

double CubicSpline::operator()(double x) const { return spline(x); }
