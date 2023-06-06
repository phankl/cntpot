#ifndef utility_header
#define utility_header

#include <cmath>
#include <iostream>
#include <utility>

#include "constants.h"
#include "cubicSpline.h"

using namespace std;

// General functions
double s5(double);
int heaviside(double);
int sgn(double);

// Specific functions for potential
double gammaFunction(double, double, const CubicSpline &);
double omegaFunction(double);
double thetaFunction(double);

double ljPair(double);
double surfaceElementDistance(double, double, double, double, double, double);
double smoothCutoff(double);
double zetaMin(double);

// Minimisation with gradient descent with line search
template <typename F>
std::pair<double, double> gradientDescent(F f, double x0, double dx,
                                          double eps) {
  double x = x0;
  double prevX = x0;
  double fX = f(x);
  double prevFX = fX;
  double df = (f(x + dx) - fX) / dx;
  double alpha = 1;

  // Perform line search
  do {
    if (std::abs(df) < eps)
      break;

    double fAlphaX = f(x - alpha * df);
    double prevFAlphaX = fAlphaX;

    while (fAlphaX <= prevFAlphaX && fAlphaX < fX) {
      prevFAlphaX = fAlphaX;
      alpha *= 2.0;
      fAlphaX = f(x - alpha * df);
    }
    fAlphaX = prevFAlphaX;
    alpha *= 0.5;
    while (fAlphaX >= fX) {
      alpha *= 0.5;
      fAlphaX = f(x - alpha * df);
    }

    prevX = x;
    prevFX = fX;
    x -= alpha * df;
    fX = fAlphaX;
    df = (f(x + dx) - fX) / dx;
  } while (std::abs((fX - prevFX) / prevFX) > eps);

  if (fX < prevFX) {
    return {x, fX};
  } else {
    return {prevX, prevFX};
  }
}

#endif
