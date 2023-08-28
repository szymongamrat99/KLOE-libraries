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

#Kinematic fits
gfortran -fPIC -llapack -lstdc++ -c Codes/kinematic_fits.f -o Compiled/kinematic_fits.o

#Shared library creation
g++ -shared -g Compiled/*.o -o librec.so
