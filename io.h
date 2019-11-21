#ifndef io_header
#define io_header

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "constants.h"

namespace io{
	void printDataFile(const std::vector<std::vector<double>>&, std::string);
}

#endif
