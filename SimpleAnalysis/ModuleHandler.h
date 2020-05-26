#ifndef MODULEHANDLER_HH
#define MODULEHANDLER_HH

#include <iostream>
#include <vector>

#include "TTree.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

#include "Module.h"
#include "CharmJetModule.h"
#include "TaggingModule.h"
#include "TaggingStudyModule.h"

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
    else if (name == "TaggingModule") {
      module = new TaggingModule(_data);
    }
    else if (name == "TaggingStudyModule") {
      module = new TaggingStudyModule(_data);
    }

    if (module != nullptr)
      this -> module_sequence.push_back(module);
  }

 private:
  ExRootTreeReader* _data = nullptr;

};

#endif
