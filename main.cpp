#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>

#include <boost/math/interpolators/cubic_b_spline.hpp>

#include "potential.h"
#include "constants.h"
#include "cubicSpline.h"
#include "bicubicSpline.h"

int main(int argc, char* argv[]){
  if(argc != 3){ 
    std::cout << "Wrong use of cntpot command. Please provide "
      "CNT chiral vector (M,N) with the following syntax: cntpot M N" << std::endl;
    return 1;
  }

  //Initialise variables
  std::stringstream arg1(argv[1]);
  std::stringstream arg2(argv[2]);
  arg1 >> M;
  arg2 >> N;
  std::cout << "Computing potential parameters for mesoscopic LJ potential for "
    "CNTs with chiral vector (" << M << "," << N << ")" << std::endl;

  R_CNT = 0.5 / PI * sqrt(3 * (N*N + M*M + N*M)) * A_C;
  R_C = 3.0 * SIGMA;
  R_C0 = 2.16 * SIGMA;
  C_OMEGA = 0.275 * (1.0 - 1.0/(1.0 + 0.59*R_CNT));
  C_THETA = 0.35 + 0.0226*(R_CNT - 6.785);

  std::string fileName = std::to_string(M) + "_" + std::to_string(N) + FILE_SUFFIX;
  POTFILE.open(fileName);

  //Generate file
  io::printHeader();
	CubicSpline uInfParaSpline = uInfGeneration();
	std::cout << "uInfPara computed!" << std::endl;
	CubicSpline gammaOrthSpline = gammaOrthGeneration(uInfParaSpline);
	std::cout << "gammaOrth computed!" << std::endl;
	BicubicSpline phiSpline = phiGeneration(uInfParaSpline);
	std::cout << "phi computed!" << std::endl;
	BicubicSpline uSemiParaSpline = uSemiGeneration();
	std::cout << "uSemiPara compted!" << std::endl;

  POTFILE.close();

	return 0;
}
