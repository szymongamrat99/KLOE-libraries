CFF = gfortran
LDFLAGS = -lgfortran -llapack -lblas -lkernlib

VPATH = src:obj
OBJDIR = obj
MODDIR = mod

all: klspm00_analysis.exe $(OBJDIR)/fortran_to_cpp_linking.o $(OBJDIR)/analysis_library.o $(OBJDIR)/analysismodule.o $(MODDIR)/analysismodule.mod

klspm00_analysis.exe: $(OBJDIR)/fortran_to_cpp_linking.o $(OBJDIR)/analysis_library.o
	$(CXX) $^ -o $@ $(LDFLAGS)

$(OBJDIR)/fortran_to_cpp_linking.o: fortran_to_cpp_linking.cpp
	$(CXX) -c $^ -o $@

$(OBJDIR)/analysismodule.o: analysismodule.f
	$(CFF) -c $< -o $@ -J$(MODDIR)

$(OBJDIR)/analysis_library.o: analysis_library.f $(OBJDIR)/analysismodule.o
	$(CFF) -c $< -o $@ -J$(MODDIR)

clean:
	rm -f klspm00_analysis.exe $(OBJDIR)/*.o $(MODDIR)/*.mod

