#ifndef linalg_header
#define linalg_header

#include <cmath>
#include <vector>
#include <iostream>

#include "io.h"

namespace linalg{
	//Linear algebra functions
	void qrFactorisation(const std::vector<std::vector<double>>&, std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);
	void qrEigenvalue(const std::vector<std::vector<double>>&, std::vector<double>&, std::vector<std::vector<double>>&);
	void projection(const std::vector<double>&, const std::vector<double>&, std::vector<double>&);

	void transpose(const std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);
	void matrixProduct(const std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);

	void add(const std::vector<double>&, const std::vector<double>&, std::vector<double>&);
	void subtract(const std::vector<double>&, const std::vector<double>&, std::vector<double>&);
	void normalise(const std::vector<double>&, std::vector<double>&);
	void scale(double, const std::vector<double>&, std::vector<double>&);

	double dot(const std::vector<double>&, const std::vector<double>&);
	double norm(const std::vector<double>&);

	void subtractMatrix(const std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&, std::vector<std::vector<double>>&);

	double dotMatrix(const std::vector<std::vector<double>>&, const std::vector<std::vector<double>>&);
	double normMatrix(const std::vector<std::vector<double>>&);
}

#endif
