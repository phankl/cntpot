CNTpot is a code for generating data files for the mesoscopic carbon nanotube potential [mesoCNT](https://docs.lammps.org/pair_mesocnt.html) 
for the MD software package LAMMPS. It is an OpenMP parallelised code which generates a single file for immediate use in LAMMPS.

The only input parameters it needs are the two components of the CNT chiral vector (M,N).
Optionally, the number of data points can also be specified.
The file generation takes about 1-4h on Intel Skylake CPUs.

Building the code is straightforward by simply running `make` in the `cntpot` directory.
The parameters of the potential other than the chiral vector can be modified by changing the values the `constants.cpp` file and recompiling the code.
