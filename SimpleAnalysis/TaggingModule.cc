#include "TaggingModule.h"

#include "TClonesArray.h"

#include "AnalysisFunctions.cc"

#include <iostream>
#include <iomanip>  
#include <fstream>

#include "TreeHandler.h"

TaggingModule::TaggingModule(ExRootTreeReader* data)
  : Module(data)
{

}

TaggingModule::~TaggingModule()
{

}

void TaggingModule::initialize()
{
  TreeHandler *tree_handler = tree_handler->getInstance();

  if (tree_handler->getTree() != nullptr) {
    tree_handler->getTree()->Branch("jet_pt", &_jet_pt, "jet_pt/F");
    tree_handler->getTree()->Branch("jet_eta", &_jet_eta, "jet_eta/F");
    tree_handler->getTree()->Branch("jet_flavor", &_jet_flavor, "jet_flavor/F");
    tree_handler->getTree()->Branch("jet_sip3dtagged", &_jet_sip3dtagged, "jet_sip3dtagged/F");
    tree_handler->getTree()->Branch("jet_ktagged", &_jet_ktagged, "jet_ktagged/F");
    tree_handler->getTree()->Branch("jet_etagged", &_jet_etagged, "jet_etagged/F");
    tree_handler->getTree()->Branch("jet_mutagged", &_jet_mutagged, "jet_mutagged/F");
    tree_handler->getTree()->Branch("jet_btag", &_jet_btag, "jet_btag/F");
    tree_handler->getTree()->Branch("bjorken_x", &_bjorken_x, "bjorken_x/F");
    tree_handler->getTree()->Branch("bjorken_Q2", &_bjorken_Q2, "bjorken_Q2/F");
    tree_handler->getTree()->Branch("JB_x", &_JB_x, "JB_x/F");
    tree_handler->getTree()->Branch("JB_Q2", &_JB_Q2, "JB_Q2/F");
  }

}

void TaggingModule::finalize()
{

}
bool TaggingModule::execute(std::map<std::string, std::any>* DataStore)
{
  auto data = getData();
  TreeHandler *tree_handler = tree_handler->getInstance();

  // Compute global DIS variables
  auto dis_variables = DISVariables(getParticles());
  _bjorken_x = dis_variables["x"];
  _bjorken_Q2 = dis_variables["Q2"];

  auto jb_variables = DISJacquetBlondel(getTracks(), getElectrons(), getPhotons(), getNeutralHadrons());
  _JB_x = jb_variables["x_JB"];
  _JB_Q2 = jb_variables["Q2_JB"];



  // If event contains at least 1 jet
  std::vector<Jet*> all_jets;
  for (int ijet = 0; ijet < getJets()->GetEntries(); ijet++) 
    {
      // Take first jet
      Jet *jet = (Jet*) getJets()->At(ijet);
      all_jets.push_back(jet);
    }

  std::vector<Jet*> fiducial_jets = SelectorFcn<Jet>(all_jets, [](Jet* j){ return (TMath::Abs(j->Eta) < 3.0 && j->PT > 5.0); });

  std::vector<Jet*> charmJets = SelectorFcn<Jet>(fiducial_jets, [](Jet* j){ return (j->Flavor == 4); });

  std::vector<Jet*> lightJets = SelectorFcn<Jet>(fiducial_jets, [](Jet* j){ return (j->Flavor < 4 || j->Flavor == 21); });

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


  
  for (auto charmjet : charmJets) {
    _jet_pt = charmjet->PT;
    _jet_eta = charmjet->Eta;
    _jet_flavor = charmjet->Flavor;
    _jet_sip3dtagged = Tagged_sIP3D(charmjet, 3.75, 0.75, 2.0);
    _jet_btag = charmjet->BTag;

    _jet_ktagged = 0;
    if (use_kaons == true) {
      _jet_ktagged = Tagged_Kaon(charmjet, std::any_cast<std::vector<Track*>>((*DataStore)["Kaons"]), 3.0, 1.0, 1);
    }

    _jet_etagged = 0;
    if (use_electrons == true) {
      _jet_etagged = Tagged_Electron(charmjet, std::any_cast<std::vector<Track*>>((*DataStore)["Electrons"]), 3.0, 1.0, 1);
    }

    _jet_mutagged = 0;
    if (use_muons == true) {
      _jet_mutagged = Tagged_Muon(charmjet, std::any_cast<std::vector<Track*>>((*DataStore)["Muons"]), 3.0, 1.0, 1);
    }

    tree_handler->getTree()->Fill();
  }

  for (auto lightjet : lightJets) {
    _jet_pt = lightjet->PT;
    _jet_eta = lightjet->Eta;
    _jet_flavor = lightjet->Flavor;
    _jet_sip3dtagged = Tagged_sIP3D(lightjet, 3.75, 0.75, 2.0);
    _jet_btag = lightjet->BTag;

    _jet_ktagged = 0;
    if (use_kaons == true) {
      _jet_ktagged = Tagged_Kaon(lightjet, std::any_cast<std::vector<Track*>>((*DataStore)["Kaons"]), 3.0, 1.0, 1);
    }

    _jet_etagged = 0;
    if (use_electrons == true) {
      _jet_etagged = Tagged_Electron(lightjet, std::any_cast<std::vector<Track*>>((*DataStore)["Electrons"]), 3.0, 1.0, 1);
    }

    _jet_mutagged = 0;
    if (use_muons == true) {
      _jet_mutagged = Tagged_Muon(lightjet, std::any_cast<std::vector<Track*>>((*DataStore)["Muons"]), 3.0, 1.0, 1);
    }

    tree_handler->getTree()->Fill();
  }


  return true;
}

float TaggingModule::sIP3D(Jet* jet, Track* track)
{
  const TLorentzVector &jetMomentum = jet->P4();
  float jpx = jetMomentum.Px();
  float jpy = jetMomentum.Py();
  float jpz = jetMomentum.Pz();

  float d0 = TMath::Abs(track->D0);

  float xd = track->Xd;
  float yd = track->Yd;
  float zd = track->Zd;
  float dd0 = TMath::Abs(track->ErrorD0);
  float dz = TMath::Abs(track->DZ);
  float ddz = TMath::Abs(track->ErrorDZ);
  
  int sign = (jpx * xd + jpy * yd + jpz * zd > 0.0) ? 1 : -1;
  //add transverse and longitudinal significances in quadrature
  float sip = sign * TMath::Sqrt(TMath::Power(d0 / dd0, 2) + TMath::Power(dz / ddz, 2));
  
  return sip;
}


bool TaggingModule::Tagged_Kaon(Jet* jet, std::vector<Track*> kaons, float minSignif, float minPT, int minKaons)
{
  bool tagged = false;
  int kaon_count = 0;
  for (auto kaon : kaons) {
    if (kaon->P4().DeltaR( jet->P4() ) > 0.5)
      continue;
    
    if (kaon->PT < minPT)
      continue;

    if (sIP3D(jet, kaon) < minSignif)
      continue;

    kaon_count++;

    if (kaon_count >= minKaons)
      break;
  }

  tagged = (kaon_count >= minKaons); 

  return tagged;
}

bool TaggingModule::Tagged_Electron(Jet* jet, std::vector<Track*> electrons, float minSignif, float minPT, int minElectrons)
{
  bool tagged = false;
  int electron_count = 0;
  for (auto electron : electrons) {
    if (electron->P4().DeltaR( jet->P4() ) > 0.5)
      continue;
    
    if (electron->PT < minPT)
      continue;

    if (sIP3D(jet, electron) < minSignif)
      continue;

    electron_count++;

    if (electron_count >= minElectrons)
      break;
  }

  tagged = (electron_count >= minElectrons); 

  return tagged;
}

bool TaggingModule::Tagged_Muon(Jet* jet, std::vector<Track*> muons, float minSignif, float minPT, int minMuons)
{
  bool tagged = false;
  int muon_count = 0;
  for (auto muon : muons) {
    if (muon->P4().DeltaR( jet->P4() ) > 0.5)
      continue;
    
    if (muon->PT < minPT)
      continue;
    
    if (sIP3D(jet, muon) < minSignif)
      continue;

    muon_count++;

    if (muon_count >= minMuons)
      break;
  }

  tagged = (muon_count >= minMuons); 

  return tagged;
}

bool TaggingModule::Tagged_sIP3D(Jet* jet, 
				float minSignif, float minPT, int minTracks)
{
  bool tagged = false;

  const TLorentzVector &jetMomentum = jet->P4();
  float jpx = jetMomentum.Px();
  float jpy = jetMomentum.Py();
  float jpz = jetMomentum.Pz();
  
  auto jet_constituents = *(getTracks());

  int N_sIPtrack = 0;

  for (int iconst = 0; iconst < jet_constituents.GetEntries(); iconst++) {


    if (N_sIPtrack >= minTracks) 
      break;

    auto constituent = jet_constituents.At(iconst);
    
    if(constituent == 0) continue;

    if (constituent->IsA() == Track::Class()) {
      auto track = static_cast<Track*>(constituent);
    
      const TLorentzVector &trkMomentum = track->P4();
      float tpt = trkMomentum.Pt();
      if(tpt < minPT) continue;

      if (trkMomentum.DeltaR(jetMomentum) > 0.5)
	continue;

      float d0 = TMath::Abs(track->D0);
      if (d0 > 3.0) 
	continue;

      float sip = sIP3D(jet, track);


      if(sip > minSignif) N_sIPtrack++;
    }

      


  }
  
  tagged = (N_sIPtrack >= minTracks);

  return tagged;
}
