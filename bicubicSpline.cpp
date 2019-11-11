#include "bicubicSpline.h"

BicubicSpline::BicubicSpline(){
}

BicubicSpline::BicubicSpline(const std::vector<CubicSpline>& inXSplines, double inXStart, double inXSpacing):
	xSplines(inXSplines),
	xStart(inXStart),
	xSpacing(inXSpacing),
	points(inXSplines.size()){
}

double BicubicSpline::operator ()(double x, double y) const{
	std::vector<double> constantYData(points, 0);
	for(int i = 0; i < points; i++) constantYData[i] = xSplines[i](y);
	
	CubicSpline ySpline(constantYData, xStart, xSpacing);
	return ySpline(x);
}
