#include "io.h"

namespace io{
	/****************************************************
	Write data to file
	
	Input: 	2D data vector fileData
		file name string fileName
	****************************************************/

	void printDataFile(const std::vector<std::vector<double>>& fileData, std::string fileName){
		std::ofstream outFile(fileName);
		
		int m = fileData.size();
		for(int i = 0; i < m; i++){
			int n = fileData[i].size();
			for(int j = 0; j < n; j++){
				outFile << std::setprecision(OUTPUT_PRECISION) << fileData[i][j] << " ";
			}
			outFile << std::endl;
		}

		outFile.close();
	};
}
