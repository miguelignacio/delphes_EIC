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

    _jet_n = 0;
    _jet_pt = std::vector<float>();
    _jet_eta = std::vector<float>();
    _jet_flavor = std::vector<int>();
    _jet_sip3dtag = std::vector<int>();
    _jet_ktag = std::vector<int>();
    _jet_etag = std::vector<int>();
    _jet_mutag = std::vector<int>();

    tree_handler->getTree()->Branch("jet_n",        &_jet_n, "jet_n/I");
    tree_handler->getTree()->Branch("jet_pt",       "std::vector<float>", &_jet_pt);
    tree_handler->getTree()->Branch("jet_eta",      "std::vector<float>", &_jet_eta);
    tree_handler->getTree()->Branch("jet_flavor",   "std::vector<int>", &_jet_flavor);
    tree_handler->getTree()->Branch("jet_sip3dtag", "std::vector<int>", &_jet_sip3dtag);
    tree_handler->getTree()->Branch("jet_ktag",     "std::vector<int>", &_jet_ktag);
    tree_handler->getTree()->Branch("jet_etag",     "std::vector<int>", &_jet_etag);
    tree_handler->getTree()->Branch("jet_mutag",    "std::vector<int>", &_jet_mutag);


    _charmjet_n = 0;
    _charmjet_pt = std::vector<float>();
    _charmjet_eta = std::vector<float>();
    tree_handler->getTree()->Branch("charmjet_n", &_charmjet_n, "charmjet_n/I");
    tree_handler->getTree()->Branch("charmjet_pt", "std::vector<float>", &_charmjet_pt);
    tree_handler->getTree()->Branch("charmjet_eta", "std::vector<float>", &_charmjet_eta);

    _met_et = 0.0;
    tree_handler->getTree()->Branch("met_et", &_met_et, "met_et/F");

    _bjorken_x = 0.0;
    _bjorken_Q2 = 0.0;
    _jb_x = 0.0;
    _jb_Q2 = 0.0;
    tree_handler->getTree()->Branch("bjorken_x", &_bjorken_x, "bjorken_x/F");
    tree_handler->getTree()->Branch("bjorken_Q2", &_bjorken_Q2, "bjorken_Q2/F");
    tree_handler->getTree()->Branch("jb_x", &_jb_x, "jb_x/F");
    tree_handler->getTree()->Branch("jb_Q2", &_jb_Q2, "jb_Q2/F");

  }


  // Initialize the cut flow

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

  // Compute global DIS variables
  auto dis_variables = DISVariables(getParticles());
  _bjorken_x = dis_variables["x"];
  _bjorken_Q2 = dis_variables["Q2"];

  auto jb_variables = DISJacquetBlondel(getTracks(), getElectrons(), getPhotons(), getNeutralHadrons());
  _jb_x = jb_variables["x_JB"];
  _jb_Q2 = jb_variables["Q2_JB"];


  // Initialize output variables
  // _charmjet_pt = std::vector<float>();
  // _charmjet_eta = std::vector<float>();
  _charmjet_pt.clear();
  _charmjet_eta.clear();
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
    passed = false;
  }

  _met_et = MET->MET;
  
  if (passed == true && MET->MET > 10.0) {
    _cut_flow["2: MET > 10 GeV"] += 1;
  } else {
    passed = false;
  }


  // If event contains at least 1 jet


  _jet_n = 0;
  // _jet_pt = std::vector<float>();
  // _jet_eta = std::vector<float>();
  // _jet_flavor = std::vector<int>();
  // _jet_sip3dtag = std::vector<int>();
  // _jet_ktag = std::vector<int>();
  // _jet_etag = std::vector<int>();
  // _jet_mutag = std::vector<int>();
  _jet_pt.clear();
  _jet_eta.clear();
  _jet_flavor.clear();
  _jet_sip3dtag.clear();
  _jet_ktag.clear();
  _jet_etag.clear();
  _jet_mutag.clear();
  
  
  auto tracks = getTracks();


  bool use_kaons = false;
  if (DataStore->find("Kaons") != DataStore->end()) {
    // store the number of kaons in the jets
    use_kaons = true;
  }

  bool use_electrons = false;
  if (DataStore->find("Electrons") != DataStore->end()) {
    // store the number of electrons in the jets
    use_electrons = true;
  }

  bool use_muons = false;
  if (DataStore->find("Muons") != DataStore->end()) {
    // store the number of muons in the jets
    use_muons = true;
  }



  std::vector<Jet*> all_jets;
  for (int ijet = 0; ijet < getJets()->GetEntries(); ijet++) 
    {
      // Take first jet
      Jet *jet = (Jet*) getJets()->At(ijet);
      all_jets.push_back(jet);

      _jet_pt.push_back( jet->PT );
      _jet_eta.push_back( jet->Eta );
      _jet_flavor.push_back( jet->Flavor );
      _jet_sip3dtag.push_back( Tagged_sIP3D(jet, *tracks, 3.75, 0.75, 2.0) );
      if (use_electrons) {
	_jet_etag.push_back( Tagged_Electron(jet, std::any_cast<std::vector<Track*>>((*DataStore)["Electrons"]), 3.0, 1.0, 1) );
      } else { 
	_jet_etag.push_back(0.0);
      }


      if (use_muons) {
	_jet_mutag.push_back( Tagged_Muon(jet, std::any_cast<std::vector<Track*>>((*DataStore)["Muons"]), 3.0, 1.0, 1) );
      } else {
	_jet_mutag.push_back( 0.0 );
      }

      if (use_kaons) {
	_jet_ktag.push_back( Tagged_Kaon(jet, std::any_cast<std::vector<Track*>>((*DataStore)["Kaons"]), 3.0, 1.0, 1) );
      } else {
	_jet_ktag.push_back( 0.0 );
      }
    }

  _jet_n = _jet_pt.size();

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




