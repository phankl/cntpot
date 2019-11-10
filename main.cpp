#include <mpi.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

#include <boost/math/interpolators/cubic_b_spline.hpp>

#include "potential.h"
#include "constants.h"

int main(int argc, char* argv[]){
	MPI_Init(&argc, &argv);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	boost::math::cubic_b_spline<double> uInfParaSpline = uInfGeneration();
	std::cout << "uInfPara computed!" << std::endl;
	boost::math::cubic_b_spline<double> gammaOrthSpline = gammaOrthGeneration(uInfParaSpline);
	std::cout << "gammaOrth computed!" << std::endl;
	BicubicSpline phiSpline = phiGeneration(uInfParaSpline);
	std::cout << "phi computed!" << std::endl;
	
	std::ofstream outFile("phiSpline.dat");

	int points = 1001;
	double h = 2*R_CNT + 3.15;
	double xi1 = -10.0;
	double xi2 = 10.0;
	double alphaSpacing = PI / (points - 1);
	for(int i = 0; i < points; i++){
		double alpha = i * alphaSpacing;
		double uInf = uApproximateInf(h, alpha, xi1, xi2, gammaOrthSpline, phiSpline);

		outFile << alpha * 180.0 / PI << " " << uInf << std::endl;
	}

	outFile.close();

	MPI_Finalize();

	return 0;
}
