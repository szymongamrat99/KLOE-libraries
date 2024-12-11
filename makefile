CXX = g++
CXXFLAGS = -fPIC `root-config --cflags --glibs` -lm -MD -MP
LDFLAGS = -shared -g

SRCPATH = Codes
OBJPATH = Compiled
INCPATH = Inc

SRC = $(wildcard $(SRCPATH)/*.cpp)
DICT = $(wildcard $(SRCPATH)/*.h)
OBJ = $(SRC:$(SRCPATH)/%.cpp=$(OBJPATH)/%.o)

SHARED = librec.so

.PHONY: all
all: $(OBJ) $(SHARED)

$(SHARED): $(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@

$(OBJPATH)/%.o: Codes/%.cpp Codes/%.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean: 
	rm -f $(OBJ)


