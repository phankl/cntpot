#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include <boost/math/interpolators/cubic_b_spline.hpp>

#include "potential.h"
#include "constants.h"
#include "cubicSpline.h"
#include "bicubicSpline.h"

int main(int argc, char* argv[]){
	CubicSpline uInfParaSpline = uInfGeneration();
	std::cout << "uInfPara computed!" << std::endl;
	CubicSpline gammaOrthSpline = gammaOrthGeneration(uInfParaSpline);
	std::cout << "gammaOrth computed!" << std::endl;
	BicubicSpline phiSpline = phiGeneration(uInfParaSpline);
	std::cout << "phi computed!" << std::endl;
	BicubicSpline uSemiParaSpline = uSemiGeneration();
	std::cout << "uSemiPara compted!" << std::endl;

	return 0;
}
