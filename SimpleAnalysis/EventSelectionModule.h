#ifndef EVENTSELECTIONMODULE_HH
#define EVENTSELECTIONMODULE_HH

#include "Module.h"
#include "classes/DelphesClasses.h"
#include <vector>
#include <map>
#include <utility>
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"


class EventSelectionModule : public Module {
public:

  EventSelectionModule(ExRootTreeReader *data);

  ~EventSelectionModule();

  void initialize() override;
  bool execute(std::map<std::string, std::any> *DataStore) override;
  void finalize() override;

private:

  // Branch variables for storage to disk
  unsigned int _charmjet_n;
  std::vector<float>_charmjet_pt;
  std::vector<float>_charmjet_eta;


  unsigned int _pid_track_n;
  std::vector<int>_pid_track_true_pid;
  std::vector<int>_pid_track_reco_pid;
  std::vector<float>_pid_track_pt;
  std::vector<float>_pid_track_eta;
  std::vector<int>_pid_track_jetmother;
  std::vector<float>_pid_track_jet_eta;
  std::vector<float>_pid_track_jet_pt;

  unsigned int _jet_n;
  std::vector<float>_jet_pt;
  std::vector<float>_jet_eta;
  std::vector<int>_jet_flavor;
  std::vector<int>_jet_nconstituents;
  std::vector<int>_jet_sip3dtag;

  // kaon information
  std::vector<int>_jet_ktag;

  std::vector<float>_jet_k1_pt;
  std::vector<float>_jet_k1_sIP3D;
  std::vector<int>_jet_k1_pid;
  std::vector<float>_jet_k2_pt;
  std::vector<float>_jet_k2_sIP3D;
  std::vector<int>_jet_k2_pid;

  std::vector<float>_jet_e1_pt;
  std::vector<float>_jet_e1_sIP3D;
  std::vector<float>_jet_e2_pt;
  std::vector<float>_jet_e2_sIP3D;

  std::vector<float>_jet_mu1_pt;
  std::vector<float>_jet_mu1_sIP3D;
  std::vector<float>_jet_mu2_pt;
  std::vector<float>_jet_mu2_sIP3D;

  std::vector<float>_jet_t1_sIP3D;
  std::vector<float>_jet_t2_sIP3D;
  std::vector<float>_jet_t3_sIP3D;
  std::vector<float>_jet_t4_sIP3D;

  std::vector<float>_jet_t1_pt;
  std::vector<float>_jet_t2_pt;
  std::vector<float>_jet_t3_pt;
  std::vector<float>_jet_t4_pt;


  std::vector<int>_jet_etag;
  std::vector<int>_jet_mutag;
  std::vector<float>_jet_charge;
  std::vector<std::vector<float> >_jet_Ks_mass;
  std::vector<std::vector<float> >_jet_Ks_p;
  std::vector<std::vector<float> >_jet_Ks_flightlength;
  std::vector<float>_jet_Ks_sumpt;
  std::vector<float>_jet_K_sumpt;
  std::vector<std::vector<float> >_jet_Ks_zhadron;
  std::vector<std::vector<float> >_jet_K_zhadron;
  std::vector<float>_jet_Ks_leading_zhadron;
  std::vector<float>_jet_K_leading_zhadron;
  std::vector<float>_jet_ehadoveremratio;

  std::vector<float>_jet_mlp_ktagger;
  std::vector<float>_jet_mlp_eltagger;
  std::vector<float>_jet_mlp_mutagger;
  std::vector<float>_jet_mlp_ip3dtagger;
  std::vector<float>_jet_mlp_globaltagger;


  float _met_et;

  // DIS variables
  float _bjorken_x;
  float _bjorken_Q2;
  float _bjorken_y;
  float _jb_x;
  float _jb_Q2;

  // Cut flow
  // std::vector<std::pair<std::string, int>> _cut_flow;
  std::map<std::string, int>_cut_flow;

  // Global variables
  float _mpi;
  float _mK;

  // Branch pointers - only initialize once per run
  TClonesArray *_branch_mRICHTracks;
  TClonesArray *_branch_barrelDircTracks;
  TClonesArray *_branch_dualRICHagTracks;
  TClonesArray *_branch_dualRICHcfTracks;

  // High-level multivariate tagging
  TMVA::Reader *_mva_reader_ip3dtagger;
  TMVA::Reader *_mva_reader_ktagger;
  TMVA::Reader *_mva_reader_eltagger;
  TMVA::Reader *_mva_reader_mutagger;
  TMVA::Reader *_mva_reader_globaltagger;

  std::map<TString, Float_t>_mva_inputs_float;
  std::map<TString, Int_t>_mva_inputs_int;
};

#endif // ifndef EVENTSELECTIONMODULE_HH
