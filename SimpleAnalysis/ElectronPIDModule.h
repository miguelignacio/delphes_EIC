#ifndef ELECTRONPIDMODULE_HH
#define ELECTRONPIDMODULE_HH

#include "Module.h"
#include "classes/DelphesClasses.h"


class ElectronPIDModule : public Module {

 public:

  ElectronPIDModule(ExRootTreeReader* data);

  ~ElectronPIDModule();

  void initialize() override {};
  bool execute(std::map<std::string, std::any>* DataStore) override;
  void finalize() override {};

 private:

  bool ElectronPID(Track* track, float eIDprob = 0.90, float separation = 3.0);
};

#endif
