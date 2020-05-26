#include "TaggingStudyModule.h"

#include "TClonesArray.h"

#include "AnalysisFunctions.cc"

#include <iostream>
#include <iomanip>  
#include <fstream>

#include "TreeHandler.h"

TaggingStudyModule::TaggingStudyModule(ExRootTreeReader* data)
  : Module(data)
{

}

TaggingStudyModule::~TaggingStudyModule()
{

}

void TaggingStudyModule::initialize()
{
  TreeHandler *tree_handler = tree_handler->getInstance();

  if (tree_handler->getTree() != nullptr) {
    tree_handler->getTree()->Branch("jet_pt", &_jet_pt, "jet_pt/F");
    tree_handler->getTree()->Branch("jet_eta", &_jet_eta, "jet_eta/F");
    tree_handler->getTree()->Branch("jet_flavor", &_jet_flavor, "jet_flavor/F");
    tree_handler->getTree()->Branch("jet_tagged", &_jet_tagged, "jet_tagged/F");
    tree_handler->getTree()->Branch("jet_btag", &_jet_btag, "jet_btag/F");
  }

  // create the study variations map
  study_variations = std::map<std::string, std::map<std::string, int>>();

  study_variations["all"]["charm"] = 0;
  study_variations["all"]["light"] = 0;
  for (auto minTrack : minTrackVar) {
    for (auto minTrkPt : minTrkPTVar) {
      for (auto minSignif : minSignifVar) {
	std::string varName = std::string(Form("MinTrk: %d;TrkPT: %.2f;MinSig: %.2f",minTrack,minTrkPt,minSignif));
	//if (study_variations.find(varName) == study_variations.end()
	study_variations[varName]["charm"] = 0;
	study_variations[varName]["light"] = 0;
      }
    }
  }
}

void TaggingStudyModule::finalize()
{
  std::vector<std::string> jets{"light", "charm"};
  // Report the yield for each variation as well as the Punzi Figure of Merit
  
  float allLight = float(study_variations["all"]["light"]);
  float allCharm = float(study_variations["all"]["charm"]);

  float maxFOM = 0.0;
  std::string bestVariation = "";

  std::cout << std::setw(40) << "VARIATION" << std::setw(6) << "LIGHT" 
	    << std::setw(6) << "CHARM" << std::setw(8) << "PUNZI" << std::endl;


  ofstream csvfile;
  csvfile.open("tagging_study.csv");
  
  csvfile << "Variation,Light,Charm" << std::endl;

  for (auto& [variation, jet] : study_variations) {
    float nLight = float(study_variations[variation]["light"]);
    float nCharm = float(study_variations[variation]["charm"]);
    float PunziFOM = (nCharm/allCharm)/( 1.5 + TMath::Sqrt(nLight));

    std::cout << std::setw(40) << variation << std::setw(6) << int(nLight)
	      << std::setw(6) << int(nCharm) << std::setw(8) << Form("%.4f", PunziFOM) << std::endl;

    csvfile << "\"" << variation << "\"," << int(nLight) << "," << nCharm << std::endl;

    if (PunziFOM > maxFOM) {
      maxFOM = PunziFOM;
      bestVariation = variation;
    }
  }

  std::cout << "======================================================================" << std::endl;
  std::cout << "Best Variation: " << bestVariation << " (FOM: " << maxFOM << ")" << std::endl;

  csvfile.close();
  

}
bool TaggingStudyModule::execute(std::map<std::string, std::any>* DataStore)
{
  auto data = getData();
  TreeHandler *tree_handler = tree_handler->getInstance();

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

  
  study_variations["all"]["charm"] += charmJets.size();
  for (auto charmjet : charmJets) {
    _jet_pt = charmjet->PT;
    _jet_eta = charmjet->Eta;
    _jet_flavor = charmjet->Flavor;
    _jet_tagged = Tagged(charmjet);
    _jet_btag = charmjet->BTag;
    tree_handler->getTree()->Fill();

    for (auto minTrack : minTrackVar) {
      for (auto minTrkPt : minTrkPTVar) {
	for (auto minSignif : minSignifVar) {
	  std::string varName = std::string(Form("MinTrk: %d;TrkPT: %.2f;MinSig: %.2f",minTrack,minTrkPt,minSignif));
	  if (Tagged(charmjet, minSignif, minTrkPt, minTrack))
	    study_variations[varName]["charm"] += 1;
	}
      }
    }
    
  }

  study_variations["all"]["light"] += lightJets.size();
  for (auto lightjet : lightJets) {
    _jet_pt = lightjet->PT;
    _jet_eta = lightjet->Eta;
    _jet_flavor = lightjet->Flavor;
    _jet_tagged = Tagged(lightjet);
    _jet_btag = lightjet->BTag;
    tree_handler->getTree()->Fill();

    for (auto minTrack : minTrackVar) {
      for (auto minTrkPt : minTrkPTVar) {
	for (auto minSignif : minSignifVar) {
	  std::string varName = std::string(Form("MinTrk: %d;TrkPT: %.2f;MinSig: %.2f",minTrack,minTrkPt,minSignif));
	  if (Tagged(lightjet, minSignif, minTrkPt, minTrack))
	    study_variations[varName]["light"] += 1;
	}
      }
    }
  }


  return true;
}


bool TaggingStudyModule::Tagged(Jet* jet, 
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

      float xd = track->Xd;
      float yd = track->Yd;
      float zd = track->Zd;
      float dd0 = TMath::Abs(track->ErrorD0);
      float dz = TMath::Abs(track->DZ);
      float ddz = TMath::Abs(track->ErrorDZ);

      int sign = (jpx * xd + jpy * yd + jpz * zd > 0.0) ? 1 : -1;
      //add transverse and longitudinal significances in quadrature
      float sip = sign * TMath::Sqrt(TMath::Power(d0 / dd0, 2) + TMath::Power(dz / ddz, 2));

      if(sip > minSignif) N_sIPtrack++;
    }

      


  }
  
  tagged = (N_sIPtrack >= minTracks);

  return tagged;
}
