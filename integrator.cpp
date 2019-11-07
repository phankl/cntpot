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
	
	linalg::qrEigenvalue(j, eigenvalues, eigenvectors);

	std::cout << "Eigenvalues of J:" << std::endl;
	io::printVector(eigenvalues);

	std::cout << "Eigenvectors of J:" << std::endl;
	io::printMatrix(eigenvectors);
}
