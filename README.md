CNTPOT is a code for generating data files for the mesoscopic carbon nanotube potential mesoCNT 
for the MD software package LAMMPS. It is an OpenMP parallelised code which generates four files for the immediate use 
in LAMMPS.

The only input parameters it needs are the two components of the CNT chiral vector (M,N).
Optionally, the number of data points can also be specified.
The file generation takes about 1-4h on modern Intel Skylake CPUs.
