#ifndef MODULE_HH
#define MODULE_HH

#include <map>
#include <string>
#include <any>

#include "TTree.h"
#include "TClonesArray.h"

#include "external/ExRootAnalysis/ExRootTreeReader.h"

#include "AnalysisFunctions.cc"


class Module {

 public:

  Module(ExRootTreeReader* data);

  ~Module();

  virtual void initialize();
  virtual bool execute(std::map<std::string, std::any>* DataStore);
  virtual void finalize();

  ExRootTreeReader* getData() { return _data;};


  // Particle Objects
  void setJets(TClonesArray* jets) { _jets = jets; };
  void setElectrons(TClonesArray* electrons) { _electrons = electrons; };
  void setPhotons(TClonesArray* photons) { _photons = photons; };
  void setNeutralHadrons(TClonesArray* neutralhadrons) { _neutralhadrons = neutralhadrons; };
  void setMuons(TClonesArray* muons) { _muons = muons; };
  void setTracks(TClonesArray* tracks) { _tracks = tracks; };
  void setEFlowTracks(TClonesArray* tracks) { _eflowtracks = tracks; };
  void setMET(TClonesArray* met) { _met = met; };

  void setGenParticles(TClonesArray* particles) { _genparticles = particles; };
  void setGenJets(TClonesArray* genjets) { _genjets = genjets; };

  TClonesArray* getJets() { return _jets;};
  TClonesArray* getElectrons() { return _electrons;};
  TClonesArray* getPhotons() { return _photons;};
  TClonesArray* getNeutralHadrons() { return _neutralhadrons;};
  TClonesArray* getMuons() { return _muons;};
  TClonesArray* getEFlowTracks() { return _eflowtracks;};
  TClonesArray* getTracks() { return _tracks;};
  TClonesArray* getMET() { return _met;};
  TClonesArray* getGenParticles() { return _genparticles;};
  TClonesArray* getGenJets() { return _genjets;};

 private:

  ExRootTreeReader* _data = nullptr; 


  // Particle Object Array Pointers
  TClonesArray* _jets = nullptr;
  TClonesArray* _electrons = nullptr;
  TClonesArray* _photons = nullptr;
  TClonesArray* _neutralhadrons = nullptr;
  TClonesArray* _muons = nullptr;
  TClonesArray* _tracks = nullptr;
  TClonesArray* _eflowtracks = nullptr;

  TClonesArray* _met = nullptr;

  TClonesArray* _genparticles = nullptr;
  TClonesArray* _genjets   = nullptr;


};

#endif
