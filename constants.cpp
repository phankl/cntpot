#include "constants.h"

// Global non-constant variables
// CNT chiral vector
int N;
int M;

// Derived CNT variables
double R_CNT;
double R_C;
double R_C0;
double C_OMEGA;
double C_THETA;

// Mathematical constants
const double PI = 3.14159265359;

// Filestream
std::ofstream POTFILE;

// Potential and CNT constants
const double SIGMA = 3.40;
const double EPSILON = 2.84e-3;
const double N_SIGMA = 0.381;
const double A_C = 1.421;

// Numerical parameters
const int TRAPEZOIDAL_POINTS = 100;
const int OUTPUT_PRECISION = 8;

const double DELTA0 = 0.3;
const double DELTA1 = 0.3;
const double DELTA2 = 2.0;
const double DELTA3 = 1.0;

// Potential file point numbers
const int GAMMA_POINTS = 1001;
const int UINF_POINTS = 1001;
const int PHI_POINTS = 1001;
const int USEMI_POINTS = 1001;

// Potential file strings
const std::string FILE_PREFIX = "C_";
const std::string FILE_SUFFIX = ".mesocnt";
const std::string CONTRIBUTOR = "Philipp Kloza";
const std::string EMAIL = "pak37@cam.ac.uk";
const std::string CITATION =
    "Volkov and Zhigilei, J Phys Chem C, 114, 5513 (2010)";

// Gauss-Legendre integrators
extern const int ANGLE_POINTS = 128;
extern const int LENGTH_POINTS = 129;
extern const int BISECTION_STEPS = 1000000;
extern const double BISECTION_EPS = 1.0e-15;
extern const Integrator ANGLE_INTEGRATOR(ANGLE_POINTS, BISECTION_STEPS,
                                         BISECTION_EPS);
extern const Integrator LENGTH_INTEGRATOR(LENGTH_POINTS, BISECTION_STEPS,
                                          BISECTION_EPS);

// Gradient descent
extern const double GD_DX = 1.0e-8;
extern const double GD_EPS = 1.0e-6;
