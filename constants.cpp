#include "constants.h"

//Chiral vector components
const int N = 10;
const int M = 10;

//Mathematical constants
const double PI = 3.14159265359;

//Potential and CNT constants
const double SIGMA = 3.40;
const double EPSILON = 2.84e-3;
const double N_SIGMA = 0.381;
const double A_C = 1.421;
const double R_CNT = 0.5 / PI * sqrt(3 * (N*N + M*M + N*M)) * A_C;
const double R_C = 3.0 * SIGMA;
const double R_C0 = 2.16 * SIGMA;
const double C_OMEGA = 0.275 * (1.0 - 1.0/(1.0 + 0.59*R_CNT));
const double C_THETA = 0.35 + 0.0226*(R_CNT - 6.785);

//Numerical parameters
const int MINIMUM_PRECISION_BITS = 26;

const double DELTA0 = 0.3;
const double DELTA1 = 0.1;
const double DELTA2 = 0.2;
const double DELTA3 = 0.3;
const double XI_LOWER = -25.0;

//Potential file point numbers
const int GAMMA_POINTS = 26;
const int UINF_POINTS = 1001;
const int PHI_POINTS = 1001;
const int USEMI_POINTS = 1001;

//Potential file names
const std::string GAMMA_FILE_NAME = "gamma";
const std::string UINF_FILE_NAME = "uInf";
const std::string PHI_FILE_NAME = "phi";
const std::string USEMI_FILE_NAME = "uSemi";
const std::string FILE_SUFFIX = ".dat";
