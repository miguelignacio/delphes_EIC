#ifndef EVENTSELECTIONMODULE_HH
#define EVENTSELECTIONMODULE_HH

#include "Module.h"
#include "classes/DelphesClasses.h"
#include <vector>
#include <map>
#include <utility>


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
  std::vector<float>_jet_k2_pt;
  std::vector<float>_jet_k2_sIP3D;

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
};

#endif // ifndef EVENTSELECTIONMODULE_HH
