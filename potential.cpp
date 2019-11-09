#include "potential.h"

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
		double etaMax = sqrt(R_C*R_C - deltaH*deltaH);

		return LENGTH_INTEGRATOR.integrate(etaIntegrand, xi-etaMax, xi+etaMax);
	};
	return R_CNT * R_CNT * N_SIGMA * N_SIGMA * ANGLE_INTEGRATOR.integrate(phi1Integrand, 0, 2*PI);
}
