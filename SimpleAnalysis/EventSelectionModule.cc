#include "EventSelectionModule.h"

#include "TClonesArray.h"
#include "TRandom3.h"
#include "TDatabasePDG.h"
#include "Math/PdfFuncMathCore.h"

#include "AnalysisFunctions.cc"
#include "TreeHandler.h"
#include "TSystem.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

EventSelectionModule::EventSelectionModule(ExRootTreeReader *data)
  : Module(data)
{}

EventSelectionModule::~EventSelectionModule()
{}

void EventSelectionModule::initialize()
{
  // Get the pointer to the mRICH tracks
  _branch_mRICHTracks      = getData()->UseBranch("mRICHTrack");
  _branch_barrelDircTracks = getData()->UseBranch("barrelDIRCTrack");
  _branch_dualRICHagTracks    = getData()->UseBranch("dualRICHagTrack");
  _branch_dualRICHcfTracks    = getData()->UseBranch("dualRICHcfTrack");


  // Tree handler initialization
  TreeHandler *tree_handler = tree_handler->getInstance();

  if (tree_handler->getTree() != nullptr) {
    _jet_n             = 0;
    _jet_pt            = std::vector<float>();
    _jet_eta           = std::vector<float>();
    _jet_flavor        = std::vector<int>();
    _jet_nconstituents = std::vector<int>();
    _jet_sip3dtag      = std::vector<int>();

    _jet_ktag     = std::vector<int>();
    _jet_k1_pt    = std::vector<float>();
    _jet_k1_sIP3D = std::vector<float>();
    _jet_k1_pid   = std::vector<int>();
    _jet_k2_pt    = std::vector<float>();
    _jet_k2_sIP3D = std::vector<float>();
    _jet_k2_pid   = std::vector<int>();

    _jet_e1_pt    = std::vector<float>();
    _jet_e1_sIP3D = std::vector<float>();
    _jet_e2_pt    = std::vector<float>();
    _jet_e2_sIP3D = std::vector<float>();

    _jet_mu1_pt    = std::vector<float>();
    _jet_mu1_sIP3D = std::vector<float>();
    _jet_mu2_pt    = std::vector<float>();
    _jet_mu2_sIP3D = std::vector<float>();

    _jet_t1_pt    = std::vector<float>();
    _jet_t1_sIP3D = std::vector<float>();
    _jet_t2_pt    = std::vector<float>();
    _jet_t2_sIP3D = std::vector<float>();
    _jet_t3_pt    = std::vector<float>();
    _jet_t3_sIP3D = std::vector<float>();
    _jet_t4_pt    = std::vector<float>();
    _jet_t4_sIP3D = std::vector<float>();


    _jet_etag  = std::vector<int>();
    _jet_mutag = std::vector<int>();

    _jet_mlp_ktagger    = std::vector<float>();
    _jet_mlp_eltagger   = std::vector<float>();
    _jet_mlp_mutagger   = std::vector<float>();
    _jet_mlp_ip3dtagger = std::vector<float>();


    _jet_charge = std::vector<float>();


    tree_handler->getTree()->Branch("jet_n",        &_jet_n, "jet_n/I");
    tree_handler->getTree()->Branch("jet_pt",       "std::vector<float>", &_jet_pt);
    tree_handler->getTree()->Branch("jet_eta",      "std::vector<float>", &_jet_eta);
    tree_handler->getTree()->Branch("jet_flavor",   "std::vector<int>", &_jet_flavor);
    tree_handler->getTree()->Branch("jet_nconstituents",   "std::vector<int>", &_jet_nconstituents);
    tree_handler->getTree()->Branch("jet_sip3dtag", "std::vector<int>", &_jet_sip3dtag);

    tree_handler->getTree()->Branch("jet_ktag",     "std::vector<int>", &_jet_ktag);
    tree_handler->getTree()->Branch("jet_k1_pt",    "std::vector<float>", &_jet_k1_pt);
    tree_handler->getTree()->Branch("jet_k1_pid",    "std::vector<int>", &_jet_k1_pid);
    tree_handler->getTree()->Branch("jet_k1_sIP3D", "std::vector<float>", &_jet_k1_sIP3D);
    tree_handler->getTree()->Branch("jet_k2_pt",    "std::vector<float>", &_jet_k2_pt);
    tree_handler->getTree()->Branch("jet_k2_pid",    "std::vector<int>", &_jet_k2_pid);
    tree_handler->getTree()->Branch("jet_k2_sIP3D", "std::vector<float>", &_jet_k2_sIP3D);

    tree_handler->getTree()->Branch("jet_e1_pt",    "std::vector<float>", &_jet_e1_pt);
    tree_handler->getTree()->Branch("jet_e1_sIP3D", "std::vector<float>", &_jet_e1_sIP3D);
    tree_handler->getTree()->Branch("jet_e2_pt",    "std::vector<float>", &_jet_e2_pt);
    tree_handler->getTree()->Branch("jet_e2_sIP3D", "std::vector<float>", &_jet_e2_sIP3D);

    tree_handler->getTree()->Branch("jet_mu1_pt",    "std::vector<float>", &_jet_mu1_pt);
    tree_handler->getTree()->Branch("jet_mu1_sIP3D", "std::vector<float>", &_jet_mu1_sIP3D);
    tree_handler->getTree()->Branch("jet_mu2_pt",    "std::vector<float>", &_jet_mu2_pt);
    tree_handler->getTree()->Branch("jet_mu2_sIP3D", "std::vector<float>", &_jet_mu2_sIP3D);

    tree_handler->getTree()->Branch("jet_t1_sIP3D", "std::vector<float>", &_jet_t1_sIP3D);
    tree_handler->getTree()->Branch("jet_t2_sIP3D", "std::vector<float>", &_jet_t2_sIP3D);
    tree_handler->getTree()->Branch("jet_t3_sIP3D", "std::vector<float>", &_jet_t3_sIP3D);
    tree_handler->getTree()->Branch("jet_t4_sIP3D", "std::vector<float>", &_jet_t4_sIP3D);

    tree_handler->getTree()->Branch("jet_t1_pt", "std::vector<float>", &_jet_t1_pt);
    tree_handler->getTree()->Branch("jet_t2_pt", "std::vector<float>", &_jet_t2_pt);
    tree_handler->getTree()->Branch("jet_t3_pt", "std::vector<float>", &_jet_t3_pt);
    tree_handler->getTree()->Branch("jet_t4_pt", "std::vector<float>", &_jet_t4_pt);

    tree_handler->getTree()->Branch("jet_mlp_ip3dtagger", "std::vector<float>", &_jet_mlp_ip3dtagger);
    tree_handler->getTree()->Branch("jet_mlp_ktagger", "std::vector<float>",    &_jet_mlp_ktagger);
    tree_handler->getTree()->Branch("jet_mlp_eltagger", "std::vector<float>",   &_jet_mlp_eltagger);
    tree_handler->getTree()->Branch("jet_mlp_mutagger", "std::vector<float>",   &_jet_mlp_mutagger);
    tree_handler->getTree()->Branch("jet_mlp_globaltagger", "std::vector<float>",   &_jet_mlp_globaltagger);

    tree_handler->getTree()->Branch("jet_etag",     "std::vector<int>", &_jet_etag);
    tree_handler->getTree()->Branch("jet_mutag",    "std::vector<int>", &_jet_mutag);
    tree_handler->getTree()->Branch("jet_charge",    "std::vector<float>", &_jet_charge);

    _jet_Ks_mass            = std::vector<std::vector<float> >();
    _jet_Ks_p               = std::vector<std::vector<float> >();
    _jet_Ks_flightlength    = std::vector<std::vector<float> >();
    _jet_Ks_sumpt           = std::vector<float>();
    _jet_K_sumpt            = std::vector<float>();
    _jet_Ks_zhadron         = std::vector<std::vector<float> >();
    _jet_K_zhadron          = std::vector<std::vector<float> >();
    _jet_Ks_leading_zhadron = std::vector<float>();
    _jet_K_leading_zhadron  =  std::vector<float>();
    _jet_ehadoveremratio    = std::vector<float>();

    tree_handler->getTree()->Branch("jet_Ks_mass",          "std::vector<std::vector<float>>", &_jet_Ks_mass);
    tree_handler->getTree()->Branch("jet_Ks_p",             "std::vector<std::vector<float>>", &_jet_Ks_p);
    tree_handler->getTree()->Branch("jet_Ks_flightlength",  "std::vector<std::vector<float>>", &_jet_Ks_flightlength);
    tree_handler->getTree()->Branch("jet_Ks_sumpt",         "std::vector<float>", &_jet_Ks_sumpt);
    tree_handler->getTree()->Branch("jet_K_sumpt",          "std::vector<float>", &_jet_K_sumpt);
    tree_handler->getTree()->Branch("jet_Ks_zhadron",             "std::vector<std::vector<float>>", &_jet_Ks_zhadron);
    tree_handler->getTree()->Branch("jet_K_zhadron",             "std::vector<std::vector<float>>", &_jet_K_zhadron);
    tree_handler->getTree()->Branch("jet_Ks_leading_zhadron",             "std::vector<float>", &_jet_Ks_leading_zhadron);
    tree_handler->getTree()->Branch("jet_K_leading_zhadron",             "std::vector<float>", &_jet_K_leading_zhadron);
    tree_handler->getTree()->Branch("jet_ehadoveremratio",  "std::vector<float>", &_jet_ehadoveremratio);

    _charmjet_n   = 0;
    _charmjet_pt  = std::vector<float>();
    _charmjet_eta = std::vector<float>();
    tree_handler->getTree()->Branch("charmjet_n", &_charmjet_n, "charmjet_n/I");
    tree_handler->getTree()->Branch("charmjet_pt", "std::vector<float>", &_charmjet_pt);
    tree_handler->getTree()->Branch("charmjet_eta", "std::vector<float>", &_charmjet_eta);

    _met_et = 0.0;
    tree_handler->getTree()->Branch("met_et", &_met_et, "met_et/F");

    _bjorken_x  = 0.0;
    _bjorken_Q2 = 0.0;
    _bjorken_y  = 0.0;
    _jb_x       = 0.0;
    _jb_Q2      = 0.0;
    tree_handler->getTree()->Branch("bjorken_x", &_bjorken_x, "bjorken_x/F");
    tree_handler->getTree()->Branch("bjorken_Q2", &_bjorken_Q2, "bjorken_Q2/F");
    tree_handler->getTree()->Branch("bjorken_y", &_bjorken_y, "bjorken_y/F");
    tree_handler->getTree()->Branch("jb_x", &_jb_x, "jb_x/F");
    tree_handler->getTree()->Branch("jb_Q2", &_jb_Q2, "jb_Q2/F");

    // pid tracks
    _pid_track_n         = 0;
    _pid_track_jetmother = std::vector<int>();
    _pid_track_jet_pt    = std::vector<float>();
    _pid_track_jet_eta   = std::vector<float>();
    _pid_track_true_pid  = std::vector<int>();
    _pid_track_reco_pid  = std::vector<int>();
    _pid_track_eta       = std::vector<float>();
    _pid_track_pt        = std::vector<float>();

    tree_handler->getTree()->Branch("pid_track_n", &_pid_track_n, "pid_track_n/I");
    tree_handler->getTree()->Branch("pid_track_true_pid", "std::vector<int>", &_pid_track_true_pid);
    tree_handler->getTree()->Branch("pid_track_reco_pid", "std::vector<int>", &_pid_track_reco_pid);
    tree_handler->getTree()->Branch("pid_track_pt", "std::vector<float>", &_pid_track_pt);
    tree_handler->getTree()->Branch("pid_track_eta", "std::vector<float>", &_pid_track_eta);
    tree_handler->getTree()->Branch("pid_track_jetmother", "std::vector<int>", &_pid_track_jetmother);
    tree_handler->getTree()->Branch("pid_track_jet_pt", "std::vector<float>", &_pid_track_jet_pt);
    tree_handler->getTree()->Branch("pid_track_jet_eta", "std::vector<float>", &_pid_track_jet_eta);
  }


  // Initialize the cut flow

  _cut_flow["1: All events"]         = 0;
  _cut_flow["2: MET > 10 GeV"]       = 0;
  _cut_flow["3: Fiducial Jets >= 1"] = 0;
  _cut_flow["4: Charm Jet == 1"]     = 0;

  // Global variables
  _mpi = TDatabasePDG().GetParticle(211)->Mass();
  _mK  = TDatabasePDG().GetParticle(321)->Mass();


  // MVA Algorithms
  _mva_reader_ip3dtagger   = new TMVA::Reader("Silent");
  _mva_reader_ktagger      = new TMVA::Reader("Silent");
  _mva_reader_eltagger     = new TMVA::Reader("Silent");
  _mva_reader_mutagger     = new TMVA::Reader("Silent");
  _mva_reader_globaltagger = new TMVA::Reader("Silent");


  _mva_inputs_float             = std::map<TString, Float_t>();
  _mva_inputs_int               = std::map<TString, Int_t>();
  _mva_inputs_float["jet_pt"]   = 0.0;
  _mva_inputs_float["jet_eta"]  = 0.0;
  _mva_inputs_int["jet_flavor"] = 0;
  _mva_inputs_float["met_et"]   = 0.0;


  _mva_inputs_float["jet_e1_pt"]    = 0.0;
  _mva_inputs_float["jet_e1_sIP3D"] = 0.0;
  _mva_inputs_float["jet_e2_pt"]    = 0.0;
  _mva_inputs_float["jet_e2_sIP3D"] = 0.0;

  _mva_inputs_float["jet_t1_pt"]    = 0.0;
  _mva_inputs_float["jet_t1_sIP3D"] = 0.0;
  _mva_inputs_float["jet_t2_pt"]    = 0.0;
  _mva_inputs_float["jet_t2_sIP3D"] = 0.0;
  _mva_inputs_float["jet_t3_pt"]    = 0.0;
  _mva_inputs_float["jet_t3_sIP3D"] = 0.0;
  _mva_inputs_float["jet_t4_pt"]    = 0.0;
  _mva_inputs_float["jet_t4_sIP3D"] = 0.0;

  _mva_inputs_float["jet_k1_pt"]    = 0.0;
  _mva_inputs_float["jet_k1_sIP3D"] = 0.0;
  _mva_inputs_float["jet_k1_pid"]   = 0;
  _mva_inputs_float["jet_k2_pt"]    = 0.0;
  _mva_inputs_float["jet_k2_sIP3D"] = 0.0;
  _mva_inputs_float["jet_k2_pid"]   = 0;

  _mva_inputs_float["jet_mu1_pt"]    = 0.0;
  _mva_inputs_float["jet_mu1_sIP3D"] = 0.0;
  _mva_inputs_float["jet_mu2_pt"]    = 0.0;
  _mva_inputs_float["jet_mu2_sIP3D"] = 0.0;


  _mva_inputs_float["jet_mlp_ip3dtagger"] = 0.0;
  _mva_inputs_float["jet_mlp_ktagger"]    = 0.0;
  _mva_inputs_float["jet_mlp_eltagger"]   = 0.0;
  _mva_inputs_float["jet_mlp_mutagger"]   = 0.0;


  _mva_reader_ip3dtagger->AddSpectator("jet_pt",     &(_mva_inputs_float["jet_pt"]));
  _mva_reader_ip3dtagger->AddSpectator("jet_eta",    &(_mva_inputs_float["jet_eta"]));
  _mva_reader_ip3dtagger->AddSpectator("jet_flavor", &(_mva_inputs_int["jet_flavor"]));
  _mva_reader_ip3dtagger->AddSpectator("met_et",     &(_mva_inputs_float["met_et"]));

  _mva_reader_ktagger->AddSpectator("jet_pt",     &(_mva_inputs_float["jet_pt"]));
  _mva_reader_ktagger->AddSpectator("jet_eta",    &(_mva_inputs_float["jet_eta"]));
  _mva_reader_ktagger->AddSpectator("jet_flavor", &(_mva_inputs_int["jet_flavor"]));
  _mva_reader_ktagger->AddSpectator("met_et",     &(_mva_inputs_float["met_et"]));

  _mva_reader_eltagger->AddSpectator("jet_pt",     &(_mva_inputs_float["jet_pt"]));
  _mva_reader_eltagger->AddSpectator("jet_eta",    &(_mva_inputs_float["jet_eta"]));
  _mva_reader_eltagger->AddSpectator("jet_flavor", &(_mva_inputs_int["jet_flavor"]));
  _mva_reader_eltagger->AddSpectator("met_et",     &(_mva_inputs_float["met_et"]));

  _mva_reader_mutagger->AddSpectator("jet_pt",     &(_mva_inputs_float["jet_pt"]));
  _mva_reader_mutagger->AddSpectator("jet_eta",    &(_mva_inputs_float["jet_eta"]));
  _mva_reader_mutagger->AddSpectator("jet_flavor", &(_mva_inputs_int["jet_flavor"]));
  _mva_reader_mutagger->AddSpectator("met_et",     &(_mva_inputs_float["met_et"]));

  _mva_reader_globaltagger->AddSpectator("jet_pt",     &(_mva_inputs_float["jet_pt"]));
  _mva_reader_globaltagger->AddSpectator("jet_eta",    &(_mva_inputs_float["jet_eta"]));
  _mva_reader_globaltagger->AddSpectator("jet_flavor", &(_mva_inputs_int["jet_flavor"]));
  _mva_reader_globaltagger->AddSpectator("met_et",     &(_mva_inputs_float["met_et"]));


  _mva_reader_eltagger->AddVariable("jet_e1_pt",    &(_mva_inputs_float["jet_e1_pt"]));
  _mva_reader_eltagger->AddVariable("jet_e1_sIP3D", &(_mva_inputs_float["jet_e1_sIP3D"]));
  _mva_reader_eltagger->AddVariable("jet_e2_pt",    &(_mva_inputs_float["jet_e2_pt"]));
  _mva_reader_eltagger->AddVariable("jet_e2_sIP3D", &(_mva_inputs_float["jet_e2_sIP3D"]));

  _mva_reader_ip3dtagger->AddVariable("jet_t1_pt",    &(_mva_inputs_float["jet_t1_pt"]));
  _mva_reader_ip3dtagger->AddVariable("jet_t1_sIP3D", &(_mva_inputs_float["jet_t1_sIP3D"]));
  _mva_reader_ip3dtagger->AddVariable("jet_t2_pt",    &(_mva_inputs_float["jet_t2_pt"]));
  _mva_reader_ip3dtagger->AddVariable("jet_t2_sIP3D", &(_mva_inputs_float["jet_t2_sIP3D"]));
  _mva_reader_ip3dtagger->AddVariable("jet_t3_pt",    &(_mva_inputs_float["jet_t3_pt"]));
  _mva_reader_ip3dtagger->AddVariable("jet_t3_sIP3D", &(_mva_inputs_float["jet_t3_sIP3D"]));
  _mva_reader_ip3dtagger->AddVariable("jet_t4_pt",    &(_mva_inputs_float["jet_t4_pt"]));
  _mva_reader_ip3dtagger->AddVariable("jet_t4_sIP3D", &(_mva_inputs_float["jet_t4_sIP3D"]));

  _mva_reader_ktagger->AddVariable("jet_k1_pt",    &(_mva_inputs_float["jet_k1_pt"]));
  _mva_reader_ktagger->AddVariable("jet_k1_sIP3D", &(_mva_inputs_float["jet_k1_sIP3D"]));
  _mva_reader_ktagger->AddVariable("jet_k2_pt",    &(_mva_inputs_float["jet_k2_pt"]));
  _mva_reader_ktagger->AddVariable("jet_k2_sIP3D", &(_mva_inputs_float["jet_k2_sIP3D"]));

  _mva_reader_mutagger->AddVariable("jet_mu1_pt",    &(_mva_inputs_float["jet_mu1_pt"]));
  _mva_reader_mutagger->AddVariable("jet_mu1_sIP3D", &(_mva_inputs_float["jet_mu1_sIP3D"]));
  _mva_reader_mutagger->AddVariable("jet_mu2_pt",    &(_mva_inputs_float["jet_mu2_pt"]));
  _mva_reader_mutagger->AddVariable("jet_mu2_sIP3D", &(_mva_inputs_float["jet_mu2_sIP3D"]));

  _mva_reader_globaltagger->AddVariable("jet_mlp_ip3dtagger", &(_mva_inputs_float["jet_mlp_ip3dtagger"]));
  _mva_reader_globaltagger->AddVariable("jet_mlp_ktagger",    &(_mva_inputs_float["jet_mlp_ktagger"]));
  _mva_reader_globaltagger->AddVariable("jet_mlp_eltagger",   &(_mva_inputs_float["jet_mlp_eltagger"]));
  _mva_reader_globaltagger->AddVariable("jet_mlp_mutagger",   &(_mva_inputs_float["jet_mlp_mutagger"]));

  _mva_reader_ip3dtagger->BookMVA("CharmIP3DTagger", "mva_taggers/TMVAClassification_CharmIP3DTagger.weights.xml");
  _mva_reader_eltagger->BookMVA("CharmETagger", "mva_taggers/TMVAClassification_CharmETagger.weights.xml");
  _mva_reader_ktagger->BookMVA("CharmKTagger", "mva_taggers/TMVAClassification_CharmKTagger.weights.xml");
  _mva_reader_mutagger->BookMVA("CharmMuTagger", "mva_taggers/TMVAClassification_CharmMuTagger.weights.xml");
  _mva_reader_globaltagger->BookMVA("CharmGlobalTagger", "mva_taggers/TMVAClassification_CharmGlobalTagger.weights.xml");
}

void EventSelectionModule::finalize()
{
  ofstream csvfile;

  csvfile.open("cut_flow.csv");

  csvfile << "Cut,Yield" << std::endl;

  for (auto&[cut, yield] : _cut_flow) {
    csvfile << "\"" << cut << "\"," << int(yield) << std::endl;
  }

  csvfile.close();

  if (_mva_reader_ip3dtagger) {
    delete _mva_reader_ip3dtagger;
  }

  if (_mva_reader_ktagger) {
    delete _mva_reader_ktagger;
  }

  if (_mva_reader_eltagger) {
    delete _mva_reader_eltagger;
  }

  if (_mva_reader_mutagger) {
    delete _mva_reader_mutagger;
  }

  if (_mva_reader_globaltagger) {
    delete _mva_reader_globaltagger;
  }
}

bool EventSelectionModule::execute(std::map<std::string, std::any> *DataStore)
{
  auto data = getData();


  // Compute global DIS variables
  auto dis_variables = DISVariables(getGenParticles());

  _bjorken_x  = dis_variables["x"];
  _bjorken_Q2 = dis_variables["Q2"];
  _bjorken_y  = dis_variables["y"];

  auto jb_variables = DISJacquetBlondel(getEFlowTracks(), getElectrons(), getPhotons(), getNeutralHadrons());
  _jb_x  = jb_variables["x_JB"];
  _jb_Q2 = jb_variables["Q2_JB"];


  // Initialize output variables
  // _charmjet_pt = std::vector<float>();
  // _charmjet_eta = std::vector<float>();
  _charmjet_pt.clear();
  _charmjet_eta.clear();
  _charmjet_n = _charmjet_pt.size();

  _pid_track_jetmother.clear();
  _pid_track_jet_pt.clear();
  _pid_track_jet_eta.clear();
  _pid_track_true_pid.clear();
  _pid_track_reco_pid.clear();
  _pid_track_pt.clear();
  _pid_track_eta.clear();
  _pid_track_n = _pid_track_pt.size();


  // Cut flow
  _cut_flow["1: All events"] += 1;
  bool passed = true;

  // Get the MET object and use it
  MissingET *MET = nullptr;

  for (int imet = 0; imet < getMET()->GetEntries(); imet++) {
    MET = static_cast<MissingET *>(getMET()->At(imet));
  }

  if (MET == nullptr) {
    passed = false;
  }

  _met_et = MET->MET;

  if ((passed == true) && (MET->MET > 10.0)) {
    _cut_flow["2: MET > 10 GeV"] += 1;
  } else {
    passed = false;
  }


  // If event contains at least 1 jet


  _jet_n = 0;
  _jet_pt.clear();
  _jet_eta.clear();
  _jet_flavor.clear();
  _jet_nconstituents.clear();
  _jet_sip3dtag.clear();
  _jet_ktag.clear();
  _jet_k1_pt.clear();
  _jet_k1_pid.clear();
  _jet_k1_sIP3D.clear();
  _jet_k2_pt.clear();
  _jet_k2_pid.clear();
  _jet_k2_sIP3D.clear();

  _jet_e1_pt.clear();
  _jet_e1_sIP3D.clear();
  _jet_e2_pt.clear();
  _jet_e2_sIP3D.clear();

  _jet_mu1_pt.clear();
  _jet_mu1_sIP3D.clear();
  _jet_mu2_pt.clear();
  _jet_mu2_sIP3D.clear();

  _jet_t1_sIP3D.clear();
  _jet_t2_sIP3D.clear();
  _jet_t3_sIP3D.clear();
  _jet_t4_sIP3D.clear();

  _jet_t1_pt.clear();
  _jet_t2_pt.clear();
  _jet_t3_pt.clear();
  _jet_t4_pt.clear();


  _jet_mlp_ktagger.clear();
  _jet_mlp_eltagger.clear();
  _jet_mlp_mutagger.clear();
  _jet_mlp_ip3dtagger.clear();
  _jet_mlp_globaltagger.clear();


  _jet_etag.clear();
  _jet_mutag.clear();
  _jet_charge.clear();

  _jet_Ks_mass.clear();
  _jet_Ks_p.clear();
  _jet_Ks_flightlength.clear();
  _jet_Ks_sumpt.clear();
  _jet_K_sumpt.clear();
  _jet_Ks_zhadron.clear();
  _jet_K_zhadron.clear();
  _jet_Ks_leading_zhadron.clear();
  _jet_K_leading_zhadron.clear();
  _jet_ehadoveremratio.clear();

  auto tracks = getTracks();


  bool use_kaons = true; // Kaons, for now, are constructed from PID maps in this code.

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


  // Loop over tracks in the event; make opposite-sign pairs; compute mass and
  // call them Ks if within some window
  // of the Ks0 mass.

  // auto particles = getGenParticles();

  std::vector<Track *> all_tracks;

  // std::vector<GenParticle*> all_tracks;
  // std::vector<TLorentzVector> all_tracks;
  for (auto obj_track : *tracks)
  {
    auto track = static_cast<Track *>(obj_track);

    if (track->PT < 0.1) continue;
    all_tracks.push_back(track);


    // if (TMath::Abs(track->PID) == 321 && track->P4().Vect().Mag() > 3.5 && track->P4().Vect().Mag() < 4.50 && track->Eta < -1.0 ) {

    //   std::cout << "  True Kaon" << std::endl;

    //   if (_branch_mRICHTracks != nullptr) {
    // 	for (int itrk = 0; itrk < _branch_mRICHTracks->GetEntries(); itrk++) {
    // 	  Track *pid_track = (Track *)_branch_mRICHTracks->At(itrk);
    // 	  if (pid_track->Particle == track->Particle) {
    // 	    std::cout << "     ID'd as " << pid_track->PID << std::endl;
    // 	  }
    // 	}
    //   }
    // }



  }


  // Build all Ks candidates
  std::vector<Candidate> Ks_candidates;

  for (int i = 0; i < all_tracks.size(); i++)
  {
    for (int j = i + 1; j < all_tracks.size(); j++)
    {
      auto track1 = all_tracks[i];
      auto track2 = all_tracks[j];

      if (track1->Charge * track2->Charge != -1) continue;

      // Treat the tracks under the charged pion hypothesis
      auto track1P4 = TLorentzVector();
      track1P4.SetPtEtaPhiM(track1->PT, track1->Eta, track1->Phi, _mpi);
      auto track2P4 = TLorentzVector();
      track2P4.SetPtEtaPhiM(track2->PT, track2->Eta, track2->Phi, _mpi);

      // std::cout <<
      // "===========================================================" <<
      // std::endl;
      // track1->P4().Print();
      // track2->P4().Print();

      // TLorentzVector Ks_candidate = track1->P4() + track2->P4();
      TLorentzVector Ks_candidate = track1P4 + track2P4;

      if (TMath::Abs(Ks_candidate.M() - 0.497) > 0.250) continue;

      // Spatial coincidence requirement
      TVector3 track1_POCA(track1->X, track1->Y, track1->Z);
      TVector3 track2_POCA(track2->X, track2->Y, track2->Z);

      auto  intertrack_displacement = track1_POCA - track2_POCA;
      float intertrack_distance     = intertrack_displacement.Mag();

      float track1_d0err = track1->ErrorD0;
      float track1_z0err = track1->ErrorDZ;
      float track2_d0err = track2->ErrorD0;
      float track2_z0err = track2->ErrorDZ;

      float err_3D = TMath::Sqrt(TMath::Power(track1_d0err, 2.0) + TMath::Power(track1_z0err, 2.0) +
                                 TMath::Power(track2_d0err, 2.0) + TMath::Power(track2_z0err, 2.0));

      float intertrack_distance_signif = intertrack_distance / err_3D;

      // std::cout << " Intertrack distance, error, significance: " <<
      // intertrack_distance << " mm"
      //            << ", " << err_3D
      //            << ", " << intertrack_distance_signif << std::endl;

      // To be from a common decay, their displacement significance should be
      // small

      if (intertrack_distance_signif > 1.5) continue;

      // build a new Candidate
      Candidate Ks;
      Ks.PID      = 310;
      Ks.Mass     = Ks_candidate.M();
      Ks.Momentum = Ks_candidate;

      // midpoint between the pocas of the two tracks
      Ks.Position = TLorentzVector(track2_POCA + 0.5 * intertrack_displacement, 0.0);

      Ks_candidates.push_back(Ks);
    }
  }


  // Get Kaons from the EIC Detector PID Systems
  std::vector<Track *> kaon_candidates;
  std::vector<Track>   pid_candidates;

  // The mRICH and barrel DIRC are physicall distinct systems covering different
  // eta ranges
  // The dualRICH system has two components - an aerogel detector for low-momentum
  // performance
  // and a CF6 system for high-momentum performance.


  if (_branch_mRICHTracks != nullptr) {
    for (int itrk = 0; itrk < _branch_mRICHTracks->GetEntries(); itrk++) {
      Track *track = (Track *)_branch_mRICHTracks->At(itrk);

      if ((track->Eta  < -3.5) || (-1.0 < track->Eta)) continue;
      
      pid_candidates.push_back(*track);
      
      Int_t reco_pid = track->PID;

      if (TMath::Abs(reco_pid) == 321) {
	kaon_candidates.push_back(track);
      }
    }
  }

  if (_branch_barrelDircTracks != nullptr) {
    for (int itrk = 0; itrk < _branch_barrelDircTracks->GetEntries(); itrk++) {
      Track *track = (Track *)_branch_barrelDircTracks->At(itrk);
      
      if ((track->Eta  < -1.0) || (1.0 < track->Eta)) continue;
      
      pid_candidates.push_back(*track);
      
      Int_t reco_pid = track->PID;

      if (TMath::Abs(reco_pid) == 321) {
	kaon_candidates.push_back(track);
      }
    }
  }


  // Handle tracks in the forward direction (dualRICH)
  for (auto track : all_tracks) {
    if (1.0 <= track->Eta && track->Eta <= 3.5) {
      Int_t final_pid = 0;
      Double_t p_track = track->P4().Vect().Mag();
      Double_t ag_p_threshold = 12.0;

      if (p_track < ag_p_threshold) {
	// region of sensitivity for Aerogel
	if (_branch_dualRICHagTracks != nullptr) {
	  for (int itrk = 0; itrk < _branch_dualRICHagTracks->GetEntries(); itrk++) {
	    Track *track_ag = (Track *)_branch_dualRICHagTracks->At(itrk);
	    if (track_ag->Particle == track->Particle) {
	      final_pid = track_ag->PID;
	      break;
	    }
	  }
	}
      } else {
	if (_branch_dualRICHcfTracks != nullptr) {
	  for (int itrk = 0; itrk < _branch_dualRICHcfTracks->GetEntries(); itrk++) {
	    Track *track_cf = (Track *)_branch_dualRICHcfTracks->At(itrk);
	    if (track_cf->Particle == track->Particle) {
	      final_pid = track_cf->PID;
	      break;
	    }
	  }
	}
      }
      
      Track drich_track = *track;
      drich_track.PID = final_pid;
      pid_candidates.push_back(drich_track);
      
      if (TMath::Abs(final_pid) == 321) 
	kaon_candidates.push_back(track);

    }
  }

  _pid_track_n = pid_candidates.size();

  for (auto track : all_tracks) {
    // try to identify the jet that might contain this track and save its
    // mother's identity (e.g. 1, 4, 5, etc.)
    float jet_min_dR     = 1e9;
    int   jet_mother_id  = 0;
    float jet_mother_pt  = 0;
    float jet_mother_eta = 1e6;

    for (int ijet = 0; ijet < getJets()->GetEntries(); ijet++)
    {
      // Take first jet
      Jet  *jet     = (Jet *)getJets()->At(ijet);
      float this_dR = jet->P4().DeltaR(track->P4());

      if ((this_dR < 1.0) && (this_dR < jet_min_dR)) {
        jet_min_dR     = this_dR;
        jet_mother_id  = jet->Flavor;
        jet_mother_pt  = jet->PT;
        jet_mother_eta = jet->Eta;
      }
    }
    _pid_track_jetmother.push_back(jet_mother_id);
    _pid_track_jet_pt.push_back(jet_mother_pt);
    _pid_track_jet_eta.push_back(jet_mother_eta);

    Int_t pid_id = 0;
    for (auto pid_track : pid_candidates) {
      if ( pid_track.Particle == track->Particle ) {
	pid_id = pid_track.PID;
	break;
      }
    }

    // if (TMath::Abs(track->PID) == 321 && track->P4().Vect().Mag() > 3.5 && track->P4().Vect().Mag() < 4.50 && track->Eta < -1.0 ) {

    //   std::cout << "  > True Kaon" << std::endl;

    //   if (_branch_mRICHTracks != nullptr) {
    // 	for (int itrk = 0; itrk < _branch_mRICHTracks->GetEntries(); itrk++) {
    // 	  Track *pid_track = (Track *)_branch_mRICHTracks->At(itrk);
    // 	  if (pid_track->Particle == track->Particle) {
    // 	    std::cout << "       ID'd as      " << pid_track->PID << std::endl;
    // 	  }
    // 	}
    //   }
    //   std::cout << "       Alt. ID'd as " << pid_id << std::endl;
    // }

    _pid_track_true_pid.push_back(track->PID);
    _pid_track_reco_pid.push_back( pid_id );
    _pid_track_pt.push_back(track->PT);
    _pid_track_eta.push_back(track->Eta);
  }


  std::vector<Jet *> all_jets;

  for (int ijet = 0; ijet < getJets()->GetEntries(); ijet++)
  {
    // Take first jet
    Jet *jet = (Jet *)getJets()->At(ijet);
    all_jets.push_back(jet);

    _jet_pt.push_back(jet->PT);
    _jet_eta.push_back(jet->Eta);
    _jet_flavor.push_back(jet->Flavor);
    _jet_nconstituents.push_back(jet->Constituents.GetEntries());

    // pre-long-lived particle decay optimization
    // _jet_sip3dtag.push_back( Tagged_sIP3D(jet, *getEFlowTracks(), 3.75, 1.00,
    // 2.0) );
    // after including long-lived particle decay, optimization
    // _jet_sip3dtag.push_back( Tagged_sIP3D(jet, *getEFlowTracks(), 2.75, 0.25,
    // 2.0) );
    // _jet_sip3dtag.push_back( Tagged_sIP3D(jet, *tagging_tracks, 2.75, 0.25,
    // 2.0) );
    // _jet_sip3dtag.push_back( Tagged_sIP3D(jet, *tagging_tracks, 3.75, 1.00,
    // 2.0) );

    // After using the significance of the charm tag yield in 100/fb to optimize
    // the tagging,
    // and after fixing the PYTHIA ctau_max bug:
    _jet_sip3dtag.push_back(Tagged_sIP3D(jet, *getEFlowTracks(), 3.00, 0.25, 2.0));


    // Find all the good flavor-tagging tracks in this jet
    std::vector<Track *> jet_tracks;

    for (int itrk = 0; itrk < getEFlowTracks()->GetEntries(); itrk++) {
      auto track = (Track *)getEFlowTracks()->At(itrk);

      if ((track->P4().DeltaR(jet->P4()) < 0.5) && IsTaggingTrack(track)) {
        jet_tracks.push_back(track);
      }
    }


    // Rank by sIP3D
    std::sort(jet_tracks.begin(), jet_tracks.end(), [jet](auto& lhs, const auto& rhs)
    {
      return sIP3D(jet, lhs) > sIP3D(jet, rhs);
    });

    // store the HIP track information
    if (jet_tracks.size() > 0) {
      _jet_t1_sIP3D.push_back(sIP3D(jet, jet_tracks[0]));
      _jet_t1_pt.push_back(jet_tracks[0]->PT);
    } else {
      _jet_t1_sIP3D.push_back(-199.0);
      _jet_t1_pt.push_back(-1.0);
    }

    if (jet_tracks.size() > 1) {
      _jet_t2_sIP3D.push_back(sIP3D(jet, jet_tracks[1]));
      _jet_t2_pt.push_back(jet_tracks[1]->PT);
    } else {
      _jet_t2_sIP3D.push_back(-199.0);
      _jet_t2_pt.push_back(-1.0);
    }

    if (jet_tracks.size() > 2) {
      _jet_t3_sIP3D.push_back(sIP3D(jet, jet_tracks[2]));
      _jet_t3_pt.push_back(jet_tracks[2]->PT);
    } else {
      _jet_t3_sIP3D.push_back(-199.0);
      _jet_t3_pt.push_back(-1.0);
    }

    if (jet_tracks.size() > 3) {
      _jet_t4_sIP3D.push_back(sIP3D(jet, jet_tracks[3]));
      _jet_t4_pt.push_back(jet_tracks[3]->PT);
    } else {
      _jet_t4_sIP3D.push_back(-199.0);
      _jet_t4_pt.push_back(-1.0);
    }


    if (use_electrons) {
      _jet_etag.push_back(Tagged_Electron(jet, std::any_cast<std::vector<Track *> >((*DataStore)["Electrons"]), 3.0, 1.0, 1));
    } else {
      _jet_etag.push_back(0.0);
    }


    if (use_muons) {
      _jet_mutag.push_back(Tagged_Muon(jet, std::any_cast<std::vector<Track *> >((*DataStore)["Muons"]), 3.0, 1.0, 1));
    } else {
      _jet_mutag.push_back(0.0);
    }

    if (use_kaons) {
      _jet_ktag.push_back(Tagged_Kaon(jet, std::any_cast<std::vector<Track *> >(kaon_candidates), 3.0, 1.0, 1));
    } else {
      _jet_ktag.push_back(0.0);
    }

    // PID-system kaons

    // Sort by PT
    std::sort(kaon_candidates.begin(), kaon_candidates.end(), [](auto& lhs, const auto& rhs)
    {
      return lhs->PT > rhs->PT;
    });

    auto kaons_list     = kaon_candidates;
    auto electrons_list = std::any_cast<std::vector<Track *> >((*DataStore)["Electrons"]);
    auto muons_list     = std::any_cast<std::vector<Track *> >((*DataStore)["Muons"]);

    // if (kaon_candidates.size() > 0) {
    if (kaons_list.size() > 0) {
      auto k1 = kaons_list.at(0);
      _jet_k1_pt.push_back(k1->PT);
      _jet_k1_sIP3D.push_back(sIP3D(jet, k1));
      _jet_k1_pid.push_back(k1->PID);
    } else {
      _jet_k1_pt.push_back(-1.0);
      _jet_k1_sIP3D.push_back(-199.0);
      _jet_k1_pid.push_back(0);
    }

    // if (kaon_candidates.size() > 1) {
    if (kaons_list.size() > 1) {
      auto k2 = kaons_list.at(1);
      _jet_k2_pt.push_back(k2->PT);
      _jet_k2_sIP3D.push_back(sIP3D(jet, k2));
      _jet_k2_pid.push_back(k2->PID);
    } else {
      _jet_k2_pt.push_back(-1.0);
      _jet_k2_sIP3D.push_back(-199.0);
      _jet_k2_pid.push_back(0);
    }

    TVector3 Ks_sumpt;
    _jet_Ks_mass.push_back(std::vector<float>());
    _jet_Ks_p.push_back(std::vector<float>());
    _jet_Ks_flightlength.push_back(std::vector<float>());
    _jet_Ks_zhadron.push_back(std::vector<float>());

    for (auto Ks_candidate : Ks_candidates) {
      if (Ks_candidate.Position.Rho() < 5) // 5mm minimum displacement from IP
        continue;

      if (Ks_candidate.Momentum.DeltaR(jet->P4()) < 0.5) {
        Ks_sumpt += Ks_candidate.Momentum.Vect();
        _jet_Ks_mass[ijet].push_back(Ks_candidate.Mass);
        _jet_Ks_p[ijet].push_back(Ks_candidate.Momentum.Rho());
        _jet_Ks_flightlength[ijet].push_back(Ks_candidate.Position.Rho());

        _jet_Ks_zhadron[ijet].push_back(TMath::Abs((Ks_candidate.Momentum.Vect() * jet->P4().Vect()) / jet->P4().Vect().Mag2()));
      }
    }

    if (_jet_Ks_zhadron[ijet].size() == 0) {
      _jet_Ks_zhadron[ijet].push_back(-1.0); // TMVA hates empty vector
                                             // branches!
    }

    _jet_Ks_leading_zhadron.push_back(*max_element(_jet_Ks_zhadron[ijet].begin(), _jet_Ks_zhadron[ijet].end()));


    _jet_Ks_sumpt.push_back(Ks_sumpt.Perp());

    // handle charged kaons
    _jet_K_zhadron.push_back(std::vector<float>());
    TVector3 K_sumpt;

    if (use_kaons) {
      auto kaons = kaon_candidates;

      for (auto kaon : kaons) {
        if (kaon->P4().DeltaR(jet->P4()) < 0.5) {
          K_sumpt += kaon->P4().Vect();
          _jet_K_zhadron[ijet].push_back(TMath::Abs((kaon->P4().Vect() * jet->P4().Vect()) / jet->P4().Vect().Mag2()));
        }
      }
    }

    if (_jet_K_zhadron[ijet].size() == 0) {
      _jet_K_zhadron[ijet].push_back(-1.0); // TMVA hates empty vector branches!
    }

    _jet_K_leading_zhadron.push_back(*max_element(_jet_K_zhadron[ijet].begin(), _jet_K_zhadron[ijet].end()));

    _jet_K_sumpt.push_back(K_sumpt.Perp());

    // Electrons
    if (electrons_list.size() > 0) {
      auto e1 = electrons_list.at(0);
      _jet_e1_pt.push_back(e1->PT);
      _jet_e1_sIP3D.push_back(sIP3D(jet, e1));
    } else {
      _jet_e1_pt.push_back(-1.0);
      _jet_e1_sIP3D.push_back(-199.0);
    }

    if (electrons_list.size() > 1) {
      auto e2 = electrons_list.at(1);
      _jet_e2_pt.push_back(e2->PT);
      _jet_e2_sIP3D.push_back(sIP3D(jet, e2));
    } else {
      _jet_e2_pt.push_back(-1.0);
      _jet_e2_sIP3D.push_back(-199.0);
    }

    // Muons
    if (muons_list.size() > 0) {
      auto mu1 = muons_list.at(0);
      _jet_mu1_pt.push_back(mu1->PT);
      _jet_mu1_sIP3D.push_back(sIP3D(jet, mu1));
    } else {
      _jet_mu1_pt.push_back(-1.0);
      _jet_mu1_sIP3D.push_back(-199.0);
    }

    if (muons_list.size() > 1) {
      auto mu2 = muons_list.at(1);
      _jet_mu2_pt.push_back(mu2->PT);
      _jet_mu2_sIP3D.push_back(sIP3D(jet, mu2));
    } else {
      _jet_mu2_pt.push_back(-1.0);
      _jet_mu2_sIP3D.push_back(-199.0);
    }

    // MLP Tagger Decisions
    _mva_inputs_int["jet_flavor"] = static_cast<Int_t>(jet->Flavor);
    _mva_inputs_float["jet_pt"]   = jet->PT;
    _mva_inputs_float["jet_eta"]  = jet->Eta;
    _mva_inputs_float["met_et"]   = _met_et;


    _mva_inputs_float["jet_t1_pt"]    = _jet_t1_pt.back();
    _mva_inputs_float["jet_t1_sIP3D"] = _jet_t1_sIP3D.back();
    _mva_inputs_float["jet_t2_pt"]    = _jet_t2_pt.back();
    _mva_inputs_float["jet_t2_sIP3D"] = _jet_t2_sIP3D.back();
    _mva_inputs_float["jet_t3_pt"]    = _jet_t3_pt.back();
    _mva_inputs_float["jet_t3_sIP3D"] = _jet_t3_sIP3D.back();
    _mva_inputs_float["jet_t4_pt"]    = _jet_t4_pt.back();
    _mva_inputs_float["jet_t4_sIP3D"] = _jet_t4_sIP3D.back();

    _mva_inputs_float["jet_k1_pt"]    = _jet_k1_pt.back();
    _mva_inputs_float["jet_k2_pt"]    = _jet_k2_pt.back();
    _mva_inputs_float["jet_k1_sIP3D"] = _jet_k1_sIP3D.back();
    _mva_inputs_float["jet_k2_sIP3D"] = _jet_k2_sIP3D.back();

    _mva_inputs_float["jet_e1_pt"]    = _jet_e1_pt.back();
    _mva_inputs_float["jet_e2_pt"]    = _jet_e2_pt.back();
    _mva_inputs_float["jet_e1_sIP3D"] = _jet_e1_sIP3D.back();
    _mva_inputs_float["jet_e2_sIP3D"] = _jet_e2_sIP3D.back();

    _mva_inputs_float["jet_mu1_pt"]    = _jet_mu1_pt.back();
    _mva_inputs_float["jet_mu2_pt"]    = _jet_mu2_pt.back();
    _mva_inputs_float["jet_mu1_sIP3D"] = _jet_mu1_sIP3D.back();
    _mva_inputs_float["jet_mu2_sIP3D"] = _jet_mu2_sIP3D.back();


    _jet_mlp_ip3dtagger.push_back(_mva_reader_ip3dtagger->EvaluateMVA("CharmIP3DTagger"));
    _jet_mlp_ktagger.push_back(_mva_reader_ktagger->EvaluateMVA("CharmKTagger"));
    _jet_mlp_eltagger.push_back(_mva_reader_eltagger->EvaluateMVA("CharmETagger"));
    _jet_mlp_mutagger.push_back(_mva_reader_mutagger->EvaluateMVA("CharmMuTagger"));

    _mva_inputs_float["jet_mlp_ip3dtagger"] = _mva_reader_ip3dtagger->EvaluateMVA("CharmIP3DTagger");
    _mva_inputs_float["jet_mlp_ktagger"]    = _mva_reader_ktagger->EvaluateMVA("CharmKTagger");
    _mva_inputs_float["jet_mlp_eltagger"]   = _mva_reader_eltagger->EvaluateMVA("CharmETagger");
    _mva_inputs_float["jet_mlp_mutagger"]   = _mva_reader_mutagger->EvaluateMVA("CharmMuTagger");

    _jet_mlp_globaltagger.push_back(_mva_reader_globaltagger->EvaluateMVA("CharmGlobalTagger"));

    // Calorimeter energy ratios
    _jet_ehadoveremratio.push_back(jet->EhadOverEem);

    // Jet charge!
    float jet_Q = -10.0;
    float kappa = 0.5;

    for (auto obj_track : *getEFlowTracks()) {
      auto track = static_cast<Track *>(obj_track);

      if (track->P4().DeltaR(jet->P4()) < 0.5) {
        if (jet_Q == -10.0) {
          jet_Q = 0.0;
        }
        jet_Q += track->Charge * TMath::Power(track->PT, kappa);
      }
    }
    jet_Q /= TMath::Power(jet->PT, kappa);
    _jet_charge.push_back(jet_Q);
  }

  _jet_n = _jet_pt.size();

  std::vector<Jet *> fiducial_jets = SelectorFcn<Jet>(all_jets, [](Jet *j) {
    return TMath::Abs(j->Eta) < 3.0 && j->PT > 5.0;
  });

  if ((passed == true) && (fiducial_jets.size() > 0)) {
    _cut_flow["3: Fiducial Jets >= 1"] += 1;
  } else {
    passed = false;
  }

  // std::vector<Jet*> charmJets = SelectorFcn<Jet>(fiducial_jets, [](Jet* j){
  // return (j->Flavor == 4); });

  std::vector<Jet *> charmJets;

  if (DataStore->find("CharmJets") != DataStore->end()) {
    charmJets = std::any_cast<std::vector<Jet *> >((*DataStore)["CharmJets"]);
  }

  if ((passed == true) && (charmJets.size() > 0)) {
    _cut_flow["4: Charm Jet == 1"] += 1;
  } else {
    passed = false;
  }


  // Store charm jet information
  if (passed) {
    for (auto jet : charmJets) {
      _charmjet_pt.push_back(jet->PT);
      _charmjet_eta.push_back(jet->Eta);
    }
    _charmjet_n = _charmjet_pt.size();
  }


  return true;
}
