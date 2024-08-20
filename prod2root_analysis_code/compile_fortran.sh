g++ -c fortran_to_cpp_linking.cpp
gfortran -c analysis_library.f
gfortran -c analysismodule.f
g++ -o klspm00_analysis.exe analysismodule.o analysis_library.o fortran_to_cpp_linking.o -lgfortran -llapack -lblas -L$CERN -lkernlib
