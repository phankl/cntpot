#ifndef io_header
#define io_header

#include <ctime>
#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

#include "constants.h"

namespace io{
  void printHeader();
	void printData(const std::vector<std::vector<double>>&);
}

#endif
