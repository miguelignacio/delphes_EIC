#ifndef EVENTSELECTIONMODULE_HH
#define EVENTSELECTIONMODULE_HH

#include "Module.h"
#include "classes/DelphesClasses.h"
#include <vector>
#include <map>
#include <utility>


class EventSelectionModule : public Module {

 public:

  EventSelectionModule(ExRootTreeReader* data);

  ~EventSelectionModule();

  void initialize() override;
  bool execute(std::map<std::string, std::any>* DataStore) override;
  void finalize() override;

 private:

  // Branch variables for storage to disk
  unsigned int       _charmjet_n;
  std::vector<float> _charmjet_pt;
  std::vector<float> _charmjet_eta;

  float              _met_et;

  // DIS variables
  float              _bjorken_x;
  float              _jb_x;
  float              _jb_Q2;

  // Cut flow
  //std::vector<std::pair<std::string, int>> _cut_flow;
  std::map<std::string, int> _cut_flow;

};

#endif
