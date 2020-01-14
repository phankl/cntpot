#include "io.h"

namespace io{
	/****************************************************
	Write potential header to file
	****************************************************/

	void printHeader(){
    time_t now = std::time(0);
    tm *ltime = localtime(&now);

    POTFILE << "DATE: " 
      << 1900 + ltime->tm_year << "-";
    if(ltime->tm_mon < 10) POTFILE << "0";
    POTFILE << 1 + ltime->tm_mon << "-";
    if(ltime->tm_mday < 10) POTFILE << "0";
    POTFILE << ltime->tm_mday
      << " CONTRIBUTOR: " << CONTRIBUTOR
      << ", " << EMAIL
      << " CITATION: " << CITATION
      << " COMMENT: Chiral vector ("
      << M << "," << N << ")"
      << std::endl;
    POTFILE << UINF_POINTS << " "
      << GAMMA_POINTS << " "
      << PHI_POINTS << " "
      << USEMI_POINTS
      << std::endl;
    POTFILE << R_CNT << " "
      << SIGMA << " " 
      << DELTA1 << " "
      << DELTA2
      << std::endl;
	};

	/****************************************************
	Write data to file
	
	Input: Data vector fileData
	****************************************************/

	void printData(const std::vector<std::vector<double>>& fileData){
	  POTFILE << std::endl;
		int m = fileData.size();
		for(int i = 0; i < m; i++){
			int n = fileData[i].size();
			for(int j = 0; j < n; j++){
				POTFILE << std::setprecision(OUTPUT_PRECISION) << fileData[i][j] << " ";
			}
			POTFILE << std::endl;
		}
	};
}
