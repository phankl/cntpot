#include "integrator.h"

Integrator::Integrator(int n){
	std::vector<std::vector<double>> j(n, std::vector<double>(n, 0));
 
	for(int i = 0; i < n; i++){
		if(i != 0) j[i][i-1] = sqrt(i*i / (2.0*i + 1.0) / (2.0*i - 1.0));
		if(i != n-1) j[i][i+1] = sqrt((i + 1.0)*(i + 1.0) / (2.0*i + 3.0) / (2.0*i + 1.0));
	}

	io::printMatrix(j);
	std::cout << std::endl;

	std::vector<double> eigenvalues(n, 0);
	std::vector<std::vector<double>> eigenvectors(n, std::vector<double>(n, 0));
	std::vector<std::vector<double>> q(n, std::vector<double>(n, 0));
	std::vector<std::vector<double>> qTranspose(n, std::vector<double>(n, 0));
	std::vector<std::vector<double>> unity(n, std::vector<double>(n, 0));
	std::vector<std::vector<double>> r(n, std::vector<double>(n, 0));

	linalg::qrFactorisation(j, q, r);
	linalg::qrEigenvalue(j, eigenvalues, eigenvectors);
	
	std::cout << "Q of J:" << std::endl;
	io::printMatrix(q);

	std::cout << "R of J:" << std::endl;
	io::printMatrix(r);

	std::cout << "Orthogonality check of Q:" << std::endl;
	linalg::transpose(q,qTranspose);
	linalg::matrixProduct(q,qTranspose,unity);
	io::printMatrix(unity);
	linalg::matrixProduct(qTranspose,q,unity);
	io::printMatrix(unity);

	std::cout << "Factorisation check of J:" << std::endl;
	linalg::matrixProduct(q,r,j);
	io::printMatrix(j);

	std::cout << "Eigenvalues of J:" << std::endl;
	io::printVector(eigenvalues);

	std::cout << "Eigenvectors of J:" << std::endl;
	io::printMatrix(eigenvectors);

	std::cout << "Weights:" << std::endl;
	for(int i = 0; i < n; i++){
		std::cout << 2*eigenvectors[i][0]*eigenvectors[i][0] << " ";
	}
	std::cout << std::endl;
}
