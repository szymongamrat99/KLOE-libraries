# Compiler
CXX = g++

# Directories
SRC_DIR = src
INC_DIR = inc
OBJ_DIR = obj
LIB_DIR = lib

# Output shared library
TARGET = $(LIB_DIR)/libmylibrary.so

# Root flags
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS = $(shell root-config --libs)

# Compiler flags
CXXFLAGS = -fPIC -I$(INC_DIR) $(ROOTCFLAGS)
LDFLAGS = -shared $(ROOTLIBS)

# Source and object files
SOURCES = $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SOURCES))

# Default target
all: $(TARGET)

# Link objects to create shared library
$(TARGET): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^

# Compile source files to object files
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -c -o $@ $<

# Clean up
clean:
	rm -rf $(OBJ_DIR) $(TARGET)

.PHONY: all clean