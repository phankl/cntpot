#include "bicubicSpline.h"

BicubicSpline::BicubicSpline(const std::vector<boost::math::cubic_b_spline<double>>& inXSplines, double inXStart, double inXSpacing):
	xSplines(inXSplines),
	xStart(inXStart),
	xSpacing(inXSpacing),
	points(inXSplines.size())
{
}

double BicubicSpline::operator ()(double x, double y){
	std::vector<double> constantYData(points, 0);
	for(int i = 0; i < points; i++) constantYData[i] = xSplines[i](y);
	
	boost::math::cubic_b_spline<double> ySpline(constantYData.begin(), constantYData.end(), xStart, xSpacing);
	return ySpline(x);
}
