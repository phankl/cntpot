#include "potential.h"

boost::math::cubic_b_spline<double> uInfGeneration(){
	double hStart = 2*R_CNT + DELTA0;
	double hEnd = 2*R_CNT + R_C;
	double hSpacing = (hEnd - hStart) / (UINF_POINTS - 1);

	std::vector<double> uInfData(UINF_POINTS, 0);
	std::vector<std::vector<double>> fileData(UINF_POINTS, std::vector<double>(2, 0));

	#pragma omp parallel for
	for(int i = 0; i < UINF_POINTS; i++){
		double h = hStart + i*hSpacing;
		double uDensity = uExactInfDensity(h, 0, 0);
		uInfData[i] = uDensity;
		fileData[i][0] = h;
		fileData[i][1] = uDensity;
	}

	std::string uInfFileName = UINF_FILE_NAME + std::to_string(N) + "x" + std::to_string(M) + FILE_SUFFIX;
	io::printDataFile(fileData, uInfFileName);

	boost::math::cubic_b_spline<double> spline(uInfData.begin(), uInfData.end(), hStart, hSpacing);
	return spline;
}

boost::math::cubic_b_spline<double> gammaOrthGeneration(const boost::math::cubic_b_spline<double>& uInfParaSpline){
	double hStart = 2*R_CNT + DELTA0;
	double hEnd = 2*R_CNT + R_C - DELTA0;
	double hSpacing = (hEnd - hStart) / (GAMMA_POINTS - 1);

	std::vector<double> gammaOrthData(GAMMA_POINTS, 0);
	std::vector<std::vector<double>> fileData(GAMMA_POINTS, std::vector<double>(2, 0));

	#pragma omp parallel for
	for(int i = 0; i < GAMMA_POINTS; i++){
		double h = hStart + i*hSpacing;

		//Define lambda functions for minimisation
		auto uInfLambda = [h](double xi)->double{
			return uExactInfDensity(h, 0.5*PI, xi);
		};
		auto uInfBarLambda = [h, uInfParaSpline](double xi)->double{
			double hBar = sqrt(h*h + xi*xi);
			if(hBar >= 2*R_CNT + R_C) return 0;
			else return uInfParaSpline(hBar);
		};

		std::pair<double, double> uInfMinimum = boost::math::tools::brent_find_minima(uInfLambda, XI_LOWER, 0.0, MINIMUM_PRECISION_BITS);
		std::pair<double, double> uInfBarMinimum = boost::math::tools::brent_find_minima(uInfBarLambda, XI_LOWER, 0.0, MINIMUM_PRECISION_BITS);

		double gamma = uInfMinimum.second / uInfBarMinimum.second;

		gammaOrthData[i] = gamma;
		fileData[i][0] = h;
		fileData[i][1] = gamma;
	}

	std::string gammaFileName = GAMMA_FILE_NAME + std::to_string(N) + "x" + std::to_string(M) + FILE_SUFFIX;
	io::printDataFile(fileData, gammaFileName);

	boost::math::cubic_b_spline<double> spline(gammaOrthData.begin(), gammaOrthData.end(), hStart, hSpacing);
	return spline;
}

BicubicSpline phiGeneration(const boost::math::cubic_b_spline<double>& uInfParaSpline){
	double hStart = 0.0;
	double hEnd = 2*R_CNT + R_C;
	double hSpacing = (hEnd - hStart) / (PHI_POINTS - 1);
	double psiSpacing = 1.0 / (PHI_POINTS - 1);

	std::vector<std::vector<double>> fileData(PHI_POINTS*PHI_POINTS, std::vector<double>(3, 0));
	std::vector<boost::math::cubic_b_spline<double>> hSplines(PHI_POINTS, boost::math::cubic_b_spline<double>());

	#pragma omp parallel for
	for(int i = 0; i < PHI_POINTS; i++){
		double h = hStart + i*hSpacing;

		auto zetaIntegrand = [h, uInfParaSpline](double zeta)->double{
			double hBar = sqrt(h*h + zeta*zeta);
			if(hBar >= 2*R_CNT + R_C) return 0;
			else return uInfParaSpline(hBar);
		};

		double hMax = 2*R_CNT + R_C;
		double zetaLower = zetaMin(h);
		double zetaUpper = sqrt(hMax*hMax - h*h);
		double zetaSpacing = (zetaUpper - zetaLower) / (PHI_POINTS - 1);

		std::vector<double> constantHData(PHI_POINTS, 0);

		for(int j = 0; j < PHI_POINTS; j++){
			double zeta = zetaLower + j*zetaSpacing;
			double psi = j * psiSpacing;
			double phi;
			if(zeta == zetaLower) phi = 0;
			else phi = LENGTH_INTEGRATOR.integrate(zetaIntegrand, zetaLower, zeta);

			constantHData[j] = phi;
			fileData[i*PHI_POINTS + j][0] = h;
			fileData[i*PHI_POINTS + j][1] = psi;
			fileData[i*PHI_POINTS + j][2] = phi;
		}

		hSplines[i] = boost::math::cubic_b_spline<double>(constantHData.begin(), constantHData.end(), 0.0, psiSpacing);
	}

	std::string phiFileName = PHI_FILE_NAME + std::to_string(N) + "x" + std::to_string(M) + FILE_SUFFIX;
	io::printDataFile(fileData, phiFileName);

	BicubicSpline phiSpline(hSplines, hStart, hSpacing);
	return phiSpline;
}

double uExactInfDensity(double h, double alpha, double xi){
	//3D integral with nested lambda functions
	auto phi1Integrand = 
		[h, alpha, xi]
		(double phi1)->double{
		auto etaIntegrand = 
			[h, alpha, xi, phi1]
			(double eta)->double{
			auto phi2Integrand =
				[h, alpha, xi, eta, phi1]
				(double phi2)->double{
				double distance = surfaceElementDistance(h, alpha, xi, eta, phi1, phi2);
				return ljPair(distance);
			};
			return ANGLE_INTEGRATOR.integrate(phi2Integrand, 0, 2.0*PI);
		};
		double sinAlpha = sin(alpha);
		double axisDistance = sqrt(h*h + xi*xi*sinAlpha*sinAlpha);
		double deltaH = axisDistance - 2*R_CNT;
		if(fabs(deltaH) >= R_C) return 0;
		
		//Compute integral range
		double etaMax; 
		if(alpha == 0) etaMax = sqrt(R_C*R_C - deltaH*deltaH);
		//Conversative choice for general case
		else etaMax = 2*R_CNT + R_C;

		return LENGTH_INTEGRATOR.integrate(etaIntegrand, xi-etaMax, xi+etaMax);
	};
	return R_CNT * R_CNT * N_SIGMA * N_SIGMA * ANGLE_INTEGRATOR.integrate(phi1Integrand, 0, 2.0*PI);
}

double uExactInf(double h, double alpha, double xi1, double xi2){
	auto xiIntegrand =
		[h, alpha]
		(double xi)->double{
			return uExactInfDensity(h, alpha, xi);
		};
		return LENGTH_INTEGRATOR.integrate(xiIntegrand, xi1, xi2);
}

double uApproximateInf(double h, double alpha, double xi1, double xi2, const boost::math::cubic_b_spline<double>& gammaOrthSpline, BicubicSpline phiSpline){
	double sinAlpha = sin(alpha);
	double gamma = gammaFunction(h, alpha, gammaOrthSpline);
	double omega = omegaFunction(alpha);
	double a = omega * sinAlpha;
	double hMax = 2*R_CNT + R_C;

	double zeta1 = xi1 * a;
	double zeta2 = xi2 * a;
	double zetaLower = zetaMin(h);
	double zetaUpper = sqrt(hMax*hMax - h*h);
	double psi1 = (fabs(zeta1) - zetaLower) / (zetaUpper - zetaLower);
	double psi2 = (fabs(zeta2) - zetaLower) / (zetaUpper - zetaLower);
	double phi1;
	double phi2;
	if(zeta1 < 0) phi1 = -phiSpline(h, psi1);
	else phi1 = phiSpline(h, psi1);
	if(zeta2 < 0) phi2 = -phiSpline(h, psi2);
	else phi2 = phiSpline(h, psi2);

	return gamma / a * (phi2 - phi1);
}
