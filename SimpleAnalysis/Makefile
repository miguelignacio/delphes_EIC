CXX = g++
CXXFLAGS = -std=c++17 $(shell root-config --cflags --ldflags --libs)
CFILES   = $(wildcard *.cc)
INCLUDE  = -I../../delphes
LIBS     = -L../../delphes -lDelphes

SimpleAnalysis.exe: *.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBS) -I. -o $@ $(CFILES) 


clean:
	rm -f SimpleAnalysis.exe
