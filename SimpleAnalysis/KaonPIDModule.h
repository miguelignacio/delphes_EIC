#ifndef KAONPIDMODULE_HH
#define KAONPIDMODULE_HH

#include "Module.h"
#include "classes/DelphesClasses.h"


class KaonPIDModule : public Module {

 public:

  KaonPIDModule(ExRootTreeReader* data);

  ~KaonPIDModule();

  void initialize() override {};
  bool execute(std::map<std::string, std::any>* DataStore) override;
  void finalize() override {};

 private:

  bool KaonPID(Track* track, float kIDprob = 0.90, float separation = 3.0);
};

#endif
