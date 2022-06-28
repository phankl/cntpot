#include "cubicSpline.h"

CubicSpline::CubicSpline(){
}

CubicSpline::CubicSpline(const std::vector<double>& inData, double inStart, double inSpacing):
	spline(boost::math::interpolators::cardinal_cubic_b_spline<double>(inData.begin(), inData.end(), inStart, inSpacing)){
}

double CubicSpline::operator ()(double x) const{
	return spline(x);
}
