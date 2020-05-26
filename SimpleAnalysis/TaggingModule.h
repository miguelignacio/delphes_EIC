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

  bool Tagged(Jet* jet, float minSignif=2.2, float minPT = 1.0, int minTracks = 3);

  // Branch variables for storage to disk
  float _jet_pt;
  float _jet_eta;
  float _jet_flavor;
  float _jet_tagged;
  float _jet_btag;



};

#endif
