CXX = g++
CXXFLAGS = -fPIC `root-config --cflags --glibs` -lm
LDFLAGS = -shared -g
PATH = Codes

kloe_library.so: kloe_class.o interference.o kinematic_fits.o
	$(CXX) -o kloe_library.so $(LDFLAGS) *.o

kloe_class.o: kloe_class.cpp kloe_class.h
	$(CXX) $(CXXFLAGS) -c kloe_class.cpp -o kloe_class.o

interference.o: interference.cpp interference.h
	$(CXX) $(CXXFLAGS) -c interference.cpp -o interference.o

kinematic_fits.o: kinematic_fits.cpp kinematic_fits.h
	$(CXX) $(CXXFLAGS) -c kinematic_fits.cpp -o kinematic_fits.o

clean:
	rm -f *.o
