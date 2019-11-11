#ifndef potential_header
#define potential_header

#include <vector>
#include <cmath>
#include <fstream>
#include <utility>

#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/tools/minima.hpp>

#include "constants.h"
#include "utility.h"
#include "io.h"
#include "bicubicSpline.h"

//File generation functions
boost::math::cubic_b_spline<double> uInfGeneration();
boost::math::cubic_b_spline<double> gammaOrthGeneration(const boost::math::cubic_b_spline<double>&);
BicubicSpline phiGeneration(const boost::math::cubic_b_spline<double>&);
BicubicSpline uSemiGeneration();

//Exact potential line density functions
double uExactInfDensity(double, double, double);
double uExactSemiDensity(double, double, double, double);

//Exact potential functions
double uExactInf(double, double, double, double);
double uExactSemi(double, double, double, double, double);

//Approximate potential functions
double uApproximateInf(double, double, double, double, const boost::math::cubic_b_spline<double>&, BicubicSpline);
double uApproximateSemi();

#endif
