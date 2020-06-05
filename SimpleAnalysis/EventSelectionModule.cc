#include "EventSelectionModule.h"

#include "TClonesArray.h"
#include "TRandom3.h"
#include "Math/PdfFuncMathCore.h"

#include "AnalysisFunctions.cc"
#include "TreeHandler.h"

#include <iostream>
#include <iomanip>  
#include <fstream>

EventSelectionModule::EventSelectionModule(ExRootTreeReader* data)
  : Module(data)
{

}

EventSelectionModule::~EventSelectionModule()
{

}

void EventSelectionModule::initialize()
{
  TreeHandler *tree_handler = tree_handler->getInstance();

  if (tree_handler->getTree() != nullptr) {
    _charmjet_n = 0;
    _charmjet_pt = std::vector<float>();
    _charmjet_eta = std::vector<float>();
    tree_handler->getTree()->Branch("charmjet_n", &_charmjet_n, "charmjet_n/I");
    tree_handler->getTree()->Branch("charmjet_pt", "std::vector<float>", &_charmjet_pt);
    tree_handler->getTree()->Branch("charmjet_eta", "std::vector<float>", &_charmjet_eta);
  }


  // Initialize the cut flow
  // _cut_flow.push_back(std::pair<std::string, int>("All events", 0));
  // _cut_flow.push_back(std::pair<std::string, int>("MET > 12 GeV", 0));
  // _cut_flow.push_back(std::pair<std::string, int>("Fiducial Jets >= 1", 0));
  // _cut_flow.push_back(std::pair<std::string, int>("Charm Jet == 1", 0));

  _cut_flow["1: All events"] = 0;
  _cut_flow["2: MET > 10 GeV"] = 0;
  _cut_flow["3: Fiducial Jets >= 1"] = 0;
  _cut_flow["4: Charm Jet == 1"] = 0;

  

}

void EventSelectionModule::finalize()
{
  ofstream csvfile;
  csvfile.open("cut_flow.csv");
  
  csvfile << "Cut,Yield" << std::endl;

  for (auto& [cut,yield] : _cut_flow) {
    csvfile << "\"" << cut << "\"," << int(yield) << std::endl;
  }

  csvfile.close();

}


bool EventSelectionModule::execute(std::map<std::string, std::any>* DataStore)
{
  auto data = getData();

  // Initialize output variables
  _charmjet_pt = std::vector<float>();
  _charmjet_eta = std::vector<float>();
  _charmjet_n = _charmjet_pt.size();

  // Cut flow
  _cut_flow["1: All events"] += 1;
  bool passed = true;

  // Get the MET object and use it
  MissingET* MET = nullptr;
  for (int imet = 0; imet < getMET()->GetEntries(); imet++) {
    MET = static_cast<MissingET*>(getMET()->At(imet));
  }
  
  if (MET == nullptr) {
    return false;
  }
  
  if (passed == true && MET->MET > 10.0) {
    _cut_flow["2: MET > 10 GeV"] += 1;
  } else {
    passed = false;
  }


  // If event contains at least 1 jet
  std::vector<Jet*> all_jets;
  for (int ijet = 0; ijet < getJets()->GetEntries(); ijet++) 
    {
      // Take first jet
      Jet *jet = (Jet*) getJets()->At(ijet);
      all_jets.push_back(jet);
    }

  std::vector<Jet*> fiducial_jets = SelectorFcn<Jet>(all_jets, [](Jet* j){ return (TMath::Abs(j->Eta) < 3.0 && j->PT > 5.0); });

  if (passed == true && fiducial_jets.size() > 0) {
    _cut_flow["3: Fiducial Jets >= 1"] += 1;
  } else {
    passed = false;
  }

  //std::vector<Jet*> charmJets = SelectorFcn<Jet>(fiducial_jets, [](Jet* j){ return (j->Flavor == 4); });
  
  std::vector<Jet*> charmJets;
  if (DataStore->find("CharmJets") != DataStore->end()) {
    charmJets = std::any_cast<std::vector<Jet*>>((*DataStore)["CharmJets"]);
  }

  if (passed == true && charmJets.size() > 0) {
    _cut_flow["4: Charm Jet == 1"] += 1;
  } else {
    passed = false;
  }


  // Store charm jet information
  if (passed) {
    for (auto jet : charmJets) {
      _charmjet_pt.push_back( jet->PT );
      _charmjet_eta.push_back( jet->Eta );
    }
    _charmjet_n = _charmjet_pt.size();
  }


  return true;
}




