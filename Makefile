CXX      = g++

CXXFLAGS= $(shell root-config --cflags) -I/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/boost/1.63.0/include
LIBS    = $(shell root-config --libs) -L/cvmfs/cms.cern.ch/slc6_amd64_gcc630/external/boost/1.63.0/lib

SOURCES = WeightCalculatorFromHistogram.cc LeptonEfficiencyCorrector.cc RoccoR.cc MainEvent.cc HmmAnalyzer.cc
HEADERS = WeightCalculatorFromHistogram.h LeptonEfficiencyCorrector.h RoccoR.h MainEvent.h HmmAnalyzer.h

OBJECTS = $(SOURCES:.cc=.o)

EXECUTABLE = analyzeHmm

all: $(SOURCES) $(EXECUTABLE)

%.o: %.cc $(HEADERS)
	@echo Compiling $<...
	$(CXX) $(CXXFLAGS) -c -o $@ $< 

$(EXECUTABLE): $(OBJECTS)
	@echo "Linking $(EXECUTABLE) ..."
	@echo "@$(CXX) $(LIBS) $(OBJECTS) -o $@"
	@$(CXX) -o $@ $^ $(LIBS) 
	@echo "done"

# Specifying the object files as intermediates deletes them automatically after the build process.
.INTERMEDIATE:  $(OBJECTS)

# The default target, which gives instructions, can be called regardless of whether or not files need to be updated.
.PHONY : clean
clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

###
MainEvent.o: MainEvent.h
AnalyzeHmm.o:MainEvent.h HmmAnalyzer.h
