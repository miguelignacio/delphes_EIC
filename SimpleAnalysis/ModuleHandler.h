#ifndef MODULEHANDLER_HH
#define MODULEHANDLER_HH

#include <iostream>
#include <vector>

#include "TTree.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

#include "Module.h"
#include "CharmJetModule.h"
#include "KaonPIDModule.h"
#include "ElectronPIDModule.h"
#include "MuonPIDModule.h"
#include "TaggingModule.h"
#include "TaggingStudyModule.h"
#include "EventSelectionModule.h"

using namespace std;

class ModuleHandler {
  static ModuleHandler *instance;
  std::vector<Module*> module_sequence;
  
  // Private constructor so that no objects can be created.
  ModuleHandler(ExRootTreeReader* data) {
    module_sequence = std::vector<Module*>();
    _data = data;
  }
  
 public:

  static ModuleHandler *getInstance(ExRootTreeReader* data) {
    if (!instance)
      instance = new ModuleHandler(data);
    return instance;
  }

  std::vector<Module*> getModules() {
    return this -> module_sequence;
  }

  void addModule(std::string name) {
    Module* module = nullptr;

    if (name == "CharmJetModule") {
      module = new CharmJetModule(_data);
    }
    else if (name == "KaonPIDModule") {
      module = new KaonPIDModule(_data);
    }
    else if (name == "ElectronPIDModule") {
      module = new ElectronPIDModule(_data);
    }
    else if (name == "MuonPIDModule") {
      module = new MuonPIDModule(_data);
    }
    else if (name == "TaggingModule") {
      module = new TaggingModule(_data);
    }
    else if (name == "TaggingStudyModule") {
      module = new TaggingStudyModule(_data);
    }
    else if (name == "EventSelectionModule") {
      module = new EventSelectionModule(_data);
    } else {
      std::cout << "ModuleHandler(): The requested module, " << name << ", is unknown to the ModuleHandler!" << std::endl;
      assert(1==1);
    }

    if (module != nullptr)
      this -> module_sequence.push_back(module);
  }

 private:
  ExRootTreeReader* _data = nullptr;

};

#endif
