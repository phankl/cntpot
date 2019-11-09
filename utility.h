#ifndef utility_header
#define utility_header

#include <cmath>
#include <boost/math/interpolators/cubic_b_spline.hpp>

#include "constants.h"

//General functions
double s5(double);
int heaviside(double);
int sgn(double);

//Specific functions for potential
double gamma(double, double, const boost::math::cubic_b_spline<double>&);
double omega(double);
double theta(double);

double ljPair(double);
double surfaceElementDistance(double, double, double, double, double, double);
double smoothCutoff(double);
double zetaMin(double);

#endif
