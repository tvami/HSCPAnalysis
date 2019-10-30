CXX = g++
CXXFLAGS = -g -Wno-unused

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

CXX_LIBS = -lm

HSCPAnalyzer.exe : HSCPAnalyzer.o pixelTree.o 
	$(CXX) $(CXXFLAGS) -o HSCPAnalyzer.exe HSCPAnalyzer.o \
	pixelTree.o \
	$(ROOTGLIBS) \
	$(CXX_LIBS)

HSCPAnalyzer2D.exe : HSCPAnalyzer2D.o pixelTree.o
	$(CXX) $(CXXFLAGS) -o HSCPAnalyzer2D.exe HSCPAnalyzer2D.o \
	pixelTree.o \
	$(ROOTGLIBS) \
	$(CXX_LIBS)

HSCPAnalyzer.o : HSCPAnalyzer.cc 
	$(CXX) -c $(CXXFLAGS) $(ROOTCFLAGS) -I.  $<

HSCPAnalyzer2D.o : HSCPAnalyzer2D.cc 
	$(CXX) -c $(CXXFLAGS) $(ROOTCFLAGS) -I.  $<

pixelTree.o : pixelTree.C pixelTree.h
	$(CXX) -c $(CXXFLAGS) $(ROOTCFLAGS) -I.  $<

clean	:
	rm *.o *.exe
	
all	:
	make HSCPAnalyzer.exe HSCPAnalyzer2D.exe