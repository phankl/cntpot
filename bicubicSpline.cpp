#include "bicubicSpline.h"

using namespace std;

BicubicSpline::BicubicSpline() {}

BicubicSpline::BicubicSpline(const vector<vector<double>> &data, double xstart,
                             double ystart, double dx, double dy)
    : xstart(xstart), ystart(ystart), dx(dx), dy(dy), size(data.size()),
      coeff(splineCoeff(data)) {}

double BicubicSpline::spline(double x, double y) const {
  int i = ceil((x - xstart) / dx);
  int j = ceil((y - ystart) / dy);

  // constant extrapolation

  if (i < 1) {
    i = 1;
    x = xstart;
  } else if (i > size - 1) {
    i = size - 1;
    x = xstart + (size - 1) * dx;
  }

  if (j < 1) {
    j = 1;
    y = ystart;
  } else if (j > size - 1) {
    j = size - 1;
    y = ystart + (size - 1) * dy;
  }

  // cubic interpolation

  double xlo = xstart + (i - 1) * dx;
  double ylo = ystart + (j - 1) * dy;
  double xbar = x - xlo;
  double ybar = y - ylo;

  double y0 = coeff[i][j][0][0] +
              ybar * (coeff[i][j][0][1] +
                      ybar * (coeff[i][j][0][2] + ybar * (coeff[i][j][0][3])));
  double y1 = coeff[i][j][1][0] +
              ybar * (coeff[i][j][1][1] +
                      ybar * (coeff[i][j][1][2] + ybar * (coeff[i][j][1][3])));
  double y2 = coeff[i][j][2][0] +
              ybar * (coeff[i][j][2][1] +
                      ybar * (coeff[i][j][2][2] + ybar * (coeff[i][j][2][3])));
  double y3 = coeff[i][j][3][0] +
              ybar * (coeff[i][j][3][1] +
                      ybar * (coeff[i][j][3][2] + ybar * (coeff[i][j][3][3])));

  return y0 + xbar * (y1 + xbar * (y2 + xbar * y3));
}

vector<vector<vector<vector<double>>>>
BicubicSpline::splineCoeff(const vector<vector<double>> &data) {
  vector<vector<vector<vector<double>>>> g(
      size, vector<vector<vector<double>>>(
                size, vector<vector<double>>(4, vector<double>(4))));

  double d;
  vector<double> bprime(size);
  vector<double> dprime(size);
  vector<vector<double>> p(size, vector<double>(size));
  vector<vector<double>> q(size, vector<double>(size));
  vector<vector<double>> s(size, vector<double>(size));
  vector<vector<double>> b(size, vector<double>(size));

  double dx_inv = 1.0 / dx;
  double dy_inv = 1.0 / dy;
  double dxsq_inv = dx_inv * dx_inv;
  double dysq_inv = dy_inv * dy_inv;
  double dxcb_inv = dx_inv * dxsq_inv;
  double dycb_inv = dy_inv * dysq_inv;

  double ax[4][4] = {{1, 0, 0, 0},
                     {0, 1, 0, 0},
                     {-3 * dxsq_inv, -2 * dx_inv, 3 * dxsq_inv, -dx_inv},
                     {2 * dxcb_inv, dxsq_inv, -2 * dxcb_inv, dxsq_inv}};
  double ay[4][4] = {{1, 0, 0, 0},
                     {0, 1, 0, 0},
                     {-3 * dysq_inv, -2 * dy_inv, 3 * dysq_inv, -dy_inv},
                     {2 * dycb_inv, dysq_inv, -2 * dycb_inv, dysq_inv}};

  // compute finite difference derivatives at boundaries

  for (int i = 0; i < size; i++) {
    p[0][i] = (data[1][i] - data[0][i]) * dx_inv;
    p[size - 1][i] = (data[size - 1][i] - data[size - 2][i]) * dx_inv;
  }

  for (int i = 0; i < size; i++) {
    q[i][0] = (data[i][1] - data[i][0]) * dy_inv;
    q[i][size - 1] = (data[i][size - 1] - data[i][size - 2]) * dy_inv;
  }

  s[0][0] = (p[0][1] - p[0][0]) * dy_inv;
  s[0][size - 1] = (p[0][size - 1] - p[0][size - 2]) * dy_inv;
  s[size - 1][0] = (p[size - 1][1] - p[size - 1][0]) * dy_inv;
  s[size - 1][size - 1] =
      (p[size - 1][size - 1] - p[size - 1][size - 2]) * dy_inv;

  // compute derivatives inside domain

  // sweep in x

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

  // compute p

  for (int j = 0; j < size; j++) {
    for (int i = 1; i < size - 1; i++) {
      d = 3 * (data[i + 1][j] - data[i - 1][j]);
      if (i == 1)
        d -= dx * p[i - 1][j];
      if (i == size - 2)
        d -= dx * p[i + 1][j];
      dprime[i] = d;
      if (i != 1)
        dprime[i] -= b[i][i - 1] * dprime[i - 1] / bprime[i - 1];
    }

    p[size - 2][j] = dprime[size - 2] / bprime[size - 2];
    for (int i = size - 3; i > 0; i--)
      p[i][j] = (dprime[i] - b[i][i + 1] * p[i + 1][j]) / bprime[i];
  }

  // compute s

  for (int j = 0; j < size; j += size - 1) {
    for (int i = 1; i < size - 1; i++) {
      d = 3 * (q[i + 1][j] - q[i - 1][j]);
      if (i == 1)
        d -= dx * s[i - 1][j];
      if (i == size - 2)
        d -= dx * s[i + 1][j];
      dprime[i] = d;
      if (i != 1)
        dprime[i] -= b[i][i - 1] * dprime[i - 1] / bprime[i - 1];
    }

    s[size - 2][j] = dprime[size - 2] / bprime[size - 2];
    for (int i = size - 3; i > 0; i--)
      s[i][j] = (dprime[i] - b[i][i + 1] * s[i + 1][j]) / bprime[i];
  }

  // sweep in y

  for (int i = 1; i < size - 1; i++) {
    if (i > 1)
      b[i][i - 1] = dy;
    b[i][i] = 4 * dy;
    if (i < size - 2)
      b[i][i + 1] = dy;
  }
  bprime[1] = b[1][1];
  for (int i = 2; i < size - 1; i++)
    bprime[i] = b[i][i] - b[i][i - 1] * b[i - 1][i] / bprime[i - 1];

  // compute q

  for (int i = 0; i < size; i++) {
    for (int j = 1; j < size - 1; j++) {
      d = 3 * (data[i][j + 1] - data[i][j - 1]);
      if (j == 1)
        d -= dy * q[i][j - 1];
      if (j == size - 2)
        d -= dy * q[i][j + 1];
      dprime[j] = d;
      if (j != 1)
        dprime[j] -= b[j][j - 1] * dprime[j - 1] / bprime[j - 1];
    }

    q[i][size - 2] = dprime[size - 2] / bprime[size - 2];
    for (int j = size - 3; j > 0; j--)
      q[i][j] = (dprime[j] - b[j][j + 1] * q[i][j + 1]) / bprime[j];
  }

  // compute s

  for (int i = 0; i < size; i++) {
    for (int j = 1; j < size - 1; j++) {
      d = 3 * (p[i][j + 1] - p[i][j - 1]);
      if (j == 1)
        d -= dy * s[i][j - 1];
      if (j == size - 2)
        d -= dy * s[i][j + 1];
      dprime[j] = d;
      if (j != 1)
        dprime[j] -= b[j][j - 1] * dprime[j - 1] / bprime[j - 1];
    }

    s[i][size - 2] = dprime[size - 2] / bprime[size - 2];
    for (int j = size - 3; j > 0; j--)
      s[i][j] = (dprime[j] - b[j][j + 1] * s[i][j + 1]) / bprime[j];
  }

  for (int i = 1; i < size; i++)
    for (int j = 1; j < size; j++) {
      for (int l = 0; l < 4; l++)
        for (int m = 0; m < 4; m++)
          g[i][j][l][m] = 0;

      double k[4][4] = {
          {data[i - 1][j - 1], q[i - 1][j - 1], data[i - 1][j], q[i - 1][j]},
          {p[i - 1][j - 1], s[i - 1][j - 1], p[i - 1][j], s[i - 1][j]},
          {data[i][j - 1], q[i][j - 1], data[i][j], q[i][j]},
          {p[i][j - 1], s[i][j - 1], p[i][j], s[i][j]}};

      for (int l = 0; l < 4; l++)
        for (int m = 0; m < 4; m++)
          for (int n = 0; n < 4; n++)
            for (int o = 0; o < 4; o++)
              g[i][j][l][m] += ax[l][n] * k[n][o] * ay[m][n];
    }

  return g;
}

double BicubicSpline::operator()(double x, double y) const {
  return spline(x, y);
}
