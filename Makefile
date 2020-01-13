#Makefile to compile the adaptive CNT aerogelation model for LAMMPS

BOOSTDIR = /usr/local/include/boost/
BOOSTLIB = /usr/local/lib/

INCLUDES = -I$(BOOSTDIR)
LFLAGS = -L$(BOOSTLIB)
LIBS =

CXX = g++
CXXFLAGS = -O3 -g -fopenmp -std=c++11 $(INCLUDES)

MAIN = cntpot

SRC = $(wildcard *.cpp)
OBJ = $(SRC:.cpp=.o)

$(MAIN): $(OBJ)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LFLAGS) $(LIBS)

.PHONY: clean
clean:
	rm -f $(OBJ) $(MAIN)
