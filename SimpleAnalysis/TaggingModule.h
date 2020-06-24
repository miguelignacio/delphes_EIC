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

  // Branch variables for storage to disk
  float _jet_pt;
  float _jet_eta;
  float _jet_flavor;
  float _jet_sip3dtagged;
  float _jet_ktagged;
  float _jet_etagged;
  float _jet_mutagged;
  float _jet_btag;

  float _bjorken_x;
  float _bjorken_Q2;
  float _JB_x;
  float _JB_Q2;

 private:
  // Methods internal to the class for flavor tagging applications


};

#endif
