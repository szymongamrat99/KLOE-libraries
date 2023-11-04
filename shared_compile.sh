#!/bin/bash

#Trilateration class
g++ -fPIC `root-config --cflags --glibs` -lm -c Codes/reconstructor.cpp -o Compiled/reconstructor.o

#IP reconstruction
g++ -fPIC `root-config --cflags --glibs` -lm -c Codes/plane_intersection.cpp -o Compiled/plane_intersection.o
g++ -fPIC `root-config --cflags --glibs` -lm -c Codes/closest_approach.cpp -o Compiled/closest_approach.o

#Interference function
g++ -fPIC `root-config --cflags --glibs` -lm -c Codes/interf_function.cpp -o Compiled/interf_function.o

#Uncertainties
g++ -fPIC `root-config --cflags --glibs` -lm -c Codes/uncertainties.cpp -o Compiled/uncertainties.o

#Charged momenta
g++ -fPIC `root-config --cflags --glibs` -lm -c Codes/charged_mom.cpp -o Compiled/charged_mom.o

#Neutral momenta
g++ -fPIC `root-config --cflags --glibs` -lm -c Codes/neutral_mom.cpp -o Compiled/neutral_mom.o

#Chi2 distribution
g++ -fPIC `root-config --cflags --glibs` -lm -c Codes/chi2_dist.cpp -o Compiled/chi2_dist.o

#Lorentz transf
g++ -fPIC `root-config --cflags --glibs` -lm -c Codes/lorentz_transf.cpp -o Compiled/lorentz_transf.o

#Constraints tri
g++ -fPIC `root-config --cflags --glibs` -lm -c Codes/constraints_tri.cpp -o Compiled/constraints_tri.o

#Shared library creation
g++ -shared -g Compiled/*.o -o librec.so
