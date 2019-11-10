#ifndef bicubicSpline_header
#define bicubicSpline_header

#include <vector>

#include <boost/math/interpolators/cubic_b_spline.hpp>

class BicubicSpline{
	private:
		const std::vector<boost::math::cubic_b_spline<double>> xSplines;
		const double xSpacing;
		const double xStart;
		const int points;

	public:
		BicubicSpline(const std::vector<boost::math::cubic_b_spline<double>>&, double, double);

		double operator ()(double, double);
};

#endif
