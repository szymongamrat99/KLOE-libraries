CXX = g++
CXXFLAGS = -fPIC `root-config --cflags --glibs` -lm
LDFLAGS = -shared -g

kloe_library.so: kloe_class.o interference.o 
	$(CXX) -o kloe_library.so $(LDFLAGS) *.o

kloe_class.o: kloe_class.cpp kloe_class.h
	$(CXX) $(CXXFLAGS) -c kloe_class.cpp -o kloe_class.o

interference.o: interference.cpp interference.h
	$(CXX) $(CXXFLAGS) -c interference.cpp -o interference.o

clean:
	rm -f *.o
