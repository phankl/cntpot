#include <mpi.h>
#include <vector>
#include <iostream>

#include "integrator.h"
#include "linalg.h"
#include "io.h"

int main(int argc, char* argv[]){
	MPI_Init(&argc, &argv);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int n = 3;

	std::vector<std::vector<double>> a{{2, 0, 0},{0, 3, 4}, {0, 4, 9}};
	std::vector<std::vector<double>> q(n, std::vector<double>(n, 0));
	std::vector<std::vector<double>> r(n, std::vector<double>(n, 0));

	std::vector<double> eigenvalues(n, 0);
	std::vector<std::vector<double>> eigenvectors(n, std::vector<double>(n, 0));

	linalg::qrFactorisation(a,q,r);
	
	std::cout << "QR Factorisation" << std::endl;
	io::printMatrix(a);
	io::printMatrix(q);
	io::printMatrix(r);

	std::cout << "Check correctness of QR Factorisation" << std::endl;
	linalg::matrixProduct(q,r,a);
	io::printMatrix(a);

	linalg::qrEigenvalue(a,eigenvalues,eigenvectors);
	std::cout << "Eigenvalues" << std::endl;
	io::printVector(eigenvalues);
	
	std::cout << "Eigenvectors" << std::endl;
	io::printMatrix(eigenvectors);

	std::cout << std::endl;
	Integrator integrator(3);

	MPI_Finalize();

	return 0;
}
