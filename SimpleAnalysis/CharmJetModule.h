#ifndef CHARMJETMODULE_HH
#define CHARMJETMODULE_HH

#include "Module.h"

class CharmJetModule : public Module {

 public:

  CharmJetModule(ExRootTreeReader* data);

  ~CharmJetModule();

  void initialize() override {};
  bool execute(std::map<std::string, std::any>* DataStore) override;
  void finalize() override {};

};

#endif
