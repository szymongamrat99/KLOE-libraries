g++ -c fortran_to_cpp_linking.cpp
gfortran -c analysis_library.f
g++ -o klspm00_analysis.exe analysis_library.o fortran_to_cpp_linking.o -lgfortran -llapack -lblas -L$CERN -lkernlib
