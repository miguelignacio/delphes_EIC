#ifndef TREEHANDLER_HH
#define TREEHANDLER_HH

#include <iostream>
#include <vector>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

using namespace std;

class TreeHandler {
  static TreeHandler *instance;
  
  TFile* _file = nullptr;
  TTree* _tree = nullptr;

  // Private constructor so that no objects can be created.
  TreeHandler(std::string filename, std::string treename) {
    _filename = filename;
    _treename = treename;
  }
  
 public:

  static TreeHandler *getInstance(std::string filename="", std::string treename="") {
    if (!instance)
      instance = new TreeHandler(filename, treename);
    return instance;
  }

  TFile* getFile() {
    return this -> _file;
  }

  TTree* getTree() {
    return this -> _tree;
  }

  void initialize() {
    _file = new TFile(_filename.c_str(), "RECREATE");
    _file->cd();
    _tree = new TTree(_treename.c_str(), "");
  }

  void execute() {
    if (_tree == nullptr)
      return;
    _tree->Fill();
  }

  void finalize() {
    if (_tree == nullptr)
      return;
    if (_file == nullptr)
      return;
    _file->cd();
    _tree->Write();
    _file->Close();
  }

 private:
  std::string _filename;
  std::string _treename;

};

#endif
