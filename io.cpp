#include "io.h"

namespace io{
	/****************************************************
	Print vector to standard output
	
	Input: vector a
	****************************************************/

	void printVector(const std::vector<double>& a){
		int n = a.size();
		for(int i = 0; i < n; i++){
			std::cout << a[i] << " ";
		}
		std::cout << std::endl;
	};


	/****************************************************
	Print matrix to standard output
	
	Input: matrix a
	****************************************************/

	void printMatrix(const std::vector<std::vector<double>>& a){
		int n = a.size();
		for(int i = 0; i < n; i++){
			int m = a[i].size();
			for(int j = 0; j < m; j++){
				std::cout << a[i][j] << " ";
			}
			std::cout << std::endl;
		}
	};
}
