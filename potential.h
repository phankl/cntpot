#ifndef potential_header
#define potential_header

#include <cmath>
#include <fstream>
#include <iostream>
#include <utility>
#include <vector>

#include <boost/math/quadrature/gauss.hpp>
#include <boost/math/tools/minima.hpp>

#include "bicubicSpline.h"
#include "constants.h"
#include "cubicSpline.h"
#include "io.h"
#include "utility.h"

// File generation functions
CubicSpline uInfGeneration();
CubicSpline gammaOrthGeneration(const CubicSpline &);
BicubicSpline phiGeneration(const CubicSpline &);
BicubicSpline uSemiGeneration();

// Exact potential line density functions
double uExactInfDensity(double, double, double);
double uExactSemiDensity(double, double, double, double);

// Exact potential functions
double uExactInf(double, double, double, double);
double uExactSemi(double, double, double, double, double);

// Approximate potential functions
double uApproximateInf(double, double, double, double, const CubicSpline &,
                       const BicubicSpline &);
double uApproximateSemi(double, double, double, double, double,
                        const CubicSpline &, const BicubicSpline &);

#endif
