#include "utility.h"

double s5(double x){
	return heaviside(-x) + heaviside(x)*heaviside(1-x)*(1 - x*x*x*(6*x*x - 15*x + 10));
}

int heaviside(double x){
	if(x > 0) return 1;
	else return 0;
}

int sgn(double x){
	if(x > 0) return 1;
	else if(x < 0) return -1;
	else return 0;
}

double gamma(double h, double alpha, const boost::math::cubic_b_spline<double>& gammaOrthSpline){
	double sinAlpha = sin(alpha);
	return 1.0 + sinAlpha*sinAlpha*(gammaOrthSpline(h) - 1.0);	
}

double omega(double alpha){
	double sinAlpha = sin(alpha);
	return 1.0 / (1.0 - C_OMEGA*sinAlpha*sinAlpha);
}

double theta(double alpha){
	double sinAlpha = sin(alpha);
	return 1.0 - C_THETA*sinAlpha*sinAlpha;
}

double ljPair(double r){
	double lj6 = pow(SIGMA/r, 6);
	double lj12 = lj6 * lj6;
	
	return 4.0 * EPSILON * (lj12 - lj6) * smoothCutoff(r);
}

double surfaceElementDistance(double h, double alpha, double xi, double eta, double phi1, double phi2){
	double x = h + R_CNT*(cos(phi2) - cos(phi1));
	double y = R_CNT*(sin(phi2)*cos(alpha) - sin(phi1)) - eta*sin(alpha);
	double z = R_CNT*sin(phi2)*sin(alpha) + eta*cos(alpha) - xi;

	return sqrt(x*x + y*y + z*z);
}

double smoothCutoff(double r){
	double tau = (r - R_C0) / (R_C - R_C0);
	return heaviside(-tau) + heaviside(tau)*heaviside(1-tau)*(1 - tau*tau*(3 - 2*tau));
}

double zetaMin(double h){
	double limit = 2.0*R_CNT + DELTA2;
	if(h > limit) return 0;
	else{
		double cutoff = s5((h - 2.0*R_CNT - DELTA1)/(DELTA2 - DELTA1));
		return cutoff * sqrt(limit*limit - h*h);
	}
}
