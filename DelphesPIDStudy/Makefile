#delphes_path = $(DELPHES_PATH)
CXX = g++
CXXFLAGS = -std=c++17 $(shell root-config --cflags --ldflags --libs) -lEG -lTMVA
CFILES   = $(wildcard *.cc)
INCLUDE  = -I$(DELPHES_PATH) 
LIBS     = -L$(DELPHES_PATH) -lDelphes

.PHONY: build check-env


build: check-env DelphesPIDStudy.exe

DelphesPIDStudy.exe: *.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBS) -I. -o $@ $(CFILES) 


clean:
	rm -f DelphesPIDStudy.exe


check-env:
ifndef DELPHES_PATH
	$(error DELPHES_PATH is undefined)
endif
