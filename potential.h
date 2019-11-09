#ifndef potential_header
#define potential_header

#include <vector>
#include <cmath>
#include <fstream>

#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <boost/math/quadrature/gauss.hpp>

#include "constants.h"
#include "utility.h"

//File generation functions
boost::math::cubic_b_spline<double> uInfGeneration();
boost::math::cubic_b_spline<double> gammaOrthGeneration(const boost::math::cubic_b_spline<double>&);
void phiGeneration(const boost::math::cubic_b_spline<double>&);
void uSemiGeneration();

//Exact potential line density functions
double uExactInfDensity(double, double, double);
double uExactSemiDensity(double, double, double);

//Exact potential functions
double uExactInf(double, double, double, double);
double uExactSemi();

//Approximate potential functions
double uApproximateInf();
double uApproximateSemi();

#endif
