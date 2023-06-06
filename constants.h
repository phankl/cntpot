#ifndef constants_header
#define constants_header

#include <cmath>
#include <fstream>
#include <string>

#include "integrator.h"

// Global non-constant variables
// CNT chiral vector
extern int N;
extern int M;

// Derived CNT variables
extern double R_CNT;
extern double R_C;
extern double R_C0;
extern double C_OMEGA;
extern double C_THETA;

// Filestream
extern std::ofstream POTFILE;

// True constants
// Mathematical constants
extern const double PI;

// Potential and CNT constants
extern const double SIGMA;
extern const double EPSILON;
extern const double N_SIGMA;
extern const double A_C;

// Numerical parameters
extern const int TRAPEZOIDAL_POINTS;
extern const int OUTPUT_PRECISION;

extern const double DELTA0;
extern const double DELTA1;
extern const double DELTA2;
extern const double DELTA3;

// Potential file point numbers
extern const int GAMMA_POINTS;
extern const int UINF_POINTS;
extern const int PHI_POINTS;
extern const int USEMI_POINTS;
extern const int MINIMUM_PRECISION_BITS;

// Potential file strings
extern const std::string FILE_PREFIX;
extern const std::string FILE_SUFFIX;
extern const std::string CONTRIBUTOR;
extern const std::string EMAIL;
extern const std::string CITATION;

// Gauss-Legendre integrators
extern const Integrator ANGLE_INTEGRATOR;
extern const Integrator LENGTH_INTEGRATOR;

// Gradient descent
extern const double GD_DX;
extern const double GD_EPS;

#endif
