# analysis_fibm15
Analysis codes in .kloe extension used before on fibm15 server, now can be used for prod2root files.

Compilators: g++ (v.11.4.0), gfortran (v.11.4.0)
Libraries: Lapack, BLAS, KERNLIB

Prerequisites:
1) GFortran: sudo apt-get install gfortran
2) Lapack and BLAS: sudo apt-get install libblas-dev liblapack-dev
3) KERNLIB: 

The analysis is used for:
1) Reconstruction of variables
2) Preselection of signal events

To compile:
1) Copy content of the repository to fibm15 src_kskl folder
2) source compile_klspm00.csh

How to pass simple variables, 1D arrays and 2D arrays from C++ to Fortran:
1) 

Warning: the hard coded paths have to be changed.
