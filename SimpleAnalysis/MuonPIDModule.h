#ifndef MUONPIDMODULE_HH
#define MUONPIDMODULE_HH

#include "Module.h"
#include "classes/DelphesClasses.h"


class MuonPIDModule : public Module {

 public:

  MuonPIDModule(ExRootTreeReader* data);

  ~MuonPIDModule();

  void initialize() override {};
  bool execute(std::map<std::string, std::any>* DataStore) override;
  void finalize() override {};

 private:

  bool MuonPID(Track* track, float muIDprob = 0.95, float separation = 2.0);
};

#endif
