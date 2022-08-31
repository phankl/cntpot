#ifndef io_header
#define io_header

#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "constants.h"

namespace io {
void printHeader();
void printData(const std::vector<std::vector<double>> &);
} // namespace io

#endif
