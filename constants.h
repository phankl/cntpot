#ifndef constants_header
#define constants_header

#include <cmath>
#include <string>

#include <boost/math/quadrature/gauss.hpp>

//Chiral vector components
extern const int N;
extern const int M;

//Mathematical constants
extern const double PI;

//Potential and CNT constants
extern const double SIGMA;
extern const double EPSILON;
extern const double N_SIGMA;
extern const double A_C;
extern const double R_CNT;
extern const double R_C;
extern const double R_C0;
extern const double C_OMEGA;
extern const double C_THETA;

//Numerical spacing parameters
extern const int MINIMUM_PRECISION_BITS;

extern const double DELTA0;
extern const double DELTA1;
extern const double DELTA2;
extern const double DELTA3;
extern const double XI_LOWER;

//Potential file point numbers
extern const int GAMMA_POINTS;
extern const int UINF_POINTS;
extern const int PHI_POINTS;
extern const int USEMI_POINTS;
extern const int MINIMUM_PRECISION_BITS;

//Potential file names
extern const std::string GAMMA_FILE_NAME;
extern const std::string UINF_FILE_NAME;
extern const std::string PHI_FILE_NAME;
extern const std::string USEMI_FILE_NAME;
extern const std::string FILE_SUFFIX;

//Gauss-Legendre integrators
extern const boost::math::quadrature::gauss<double, 128> ANGLE_INTEGRATOR;
extern const boost::math::quadrature::gauss<double, 129> LENGTH_INTEGRATOR;

#endif
