#ifndef TAGGINGMODULE_HH
#define TAGGINGMODULE_HH

#include "classes/DelphesClasses.h"

#include "Module.h"

class TaggingModule : public Module {

 public:

  TaggingModule(ExRootTreeReader* data);

  ~TaggingModule();

  void initialize() override;
  bool execute(std::map<std::string, std::any>* DataStore) override;
  void finalize() override;


 private:

  float sIP3D(Jet* jet, Track* track);

  bool Tagged_Kaon(Jet* jet, std::vector<Track*> kaons, float minSignif = 3.0, float minPT = 1.0, int minKaons = 1);
  bool Tagged_Electron(Jet* jet, std::vector<Track*> electrons, float minSignif = 3.0, float minPT = 1.0, int minElectrons = 1);
  bool Tagged_Muon(Jet* jet, std::vector<Track*> muons, float minSignif = 3.0, float minPT = 1.0, int minMuons = 1);
  bool Tagged_sIP3D(Jet* jet, float minSignif=2.2, float minPT = 1.0, int minTracks = 3);


  // Branch variables for storage to disk
  float _jet_pt;
  float _jet_eta;
  float _jet_flavor;
  float _jet_sip3dtagged;
  float _jet_ktagged;
  float _jet_etagged;
  float _jet_mutagged;
  float _jet_btag;

  float _jet_nk;



};

#endif
