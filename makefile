CXX = g++
CXXFLAGS = -fPIC `root-config --cflags --glibs` -lm
LDFLAGS = -shared -g

SRCPATH = Codes
OBJPATH = Compiled
INCPATH = Inc

SRC = $(wildcard $(SRCPATH)/*.cpp)
DICT = $(wildcard $(INCPATH)/*.h)
OBJ = 

SHARED = librec.so

.PHONY: all clean
all: $(OBJ) $(SHARED)

$(SHARED): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@

%.o: %.c %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean: 
	rm -f $(OBJ)


