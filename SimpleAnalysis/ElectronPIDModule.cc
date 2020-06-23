#include "ElectronPIDModule.h"

#include "TClonesArray.h"
#include "TRandom3.h"
#include "Math/PdfFuncMathCore.h"

#include "AnalysisFunctions.cc"

#include <iostream>
#include <algorithm>

ElectronPIDModule::ElectronPIDModule(ExRootTreeReader* data)
  : Module(data)
{

}

ElectronPIDModule::~ElectronPIDModule()
{

}



bool ElectronPIDModule::execute(std::map<std::string, std::any>* DataStore)
{
  auto data = getData();

  // If event contains at least 1 jet
  // std::vector<Jet*> all_jets;
  // for (int ijet = 0; ijet < getJets()->GetEntries(); ijet++) 
  //   {
  //     // Take first jet
  //     Jet *jet = (Jet*) getJets()->At(ijet);
  //     all_jets.push_back(jet);
  //   }

  // std::vector<Jet*> fiducial_jets = SelectorFcn<Jet>(all_jets, [](Jet* j){ return (TMath::Abs(j->Eta) < 3.0 && j->PT > 5.0); });
  // std::vector<Jet*> charmJets = SelectorFcn<Jet>(fiducial_jets, [](Jet* j){ return (j->Flavor == 4); });

  std::vector<Track*> all_electrons;

  // If the DataStore contains already a list of tracks to be used for PID assignment,
  // use that. If not, use the getEFlowTracks() method to get the general list of all tracks.
  // TracksForPID is special - it contains tracks NOT already used by another PID algorithm,
  // to avoid using the same track (pion) twice in two categories.
  
  float e_efficiency = 0.90;
  float epi_separation = 2.4; //sigma
  //float e_efficiency = 1.00;
  //float epi_separation = 100.0; //sigma

  if (DataStore->find("TracksForPID") != DataStore->end()) {
    for (auto track : std::any_cast<std::vector<Track*>>((*DataStore)["TracksForPID"])) {
      if (ElectronPID(track, e_efficiency, epi_separation))
	all_electrons.push_back(track);
    } 
    
    
  } else {      
    for (int itrk = 0; itrk < getEFlowTracks()->GetEntries(); itrk++)
      {
	Track* track = (Track*) getEFlowTracks()->At(itrk);
	if (ElectronPID(track, e_efficiency, epi_separation))
	  all_electrons.push_back(track);
      }

    std::vector<Track*> tracks_for_PID;
    for (int itrk = 0; itrk < getEFlowTracks()->GetEntries(); itrk++) 
      {
	Track* track = (Track*) getEFlowTracks()->At(itrk);
	tracks_for_PID.push_back(track);
      }
    (*DataStore)["TracksForPID"] = tracks_for_PID;
  }

  //auto reconstructed_electrons = all_electrons;
  std::vector<Track*> reconstructed_electrons = SelectorFcn<Track>(all_electrons, 
								   [](Track* t){ return (t->PT >= 0.1); });
  
  std::vector<Track*> tracks_for_PID = std::any_cast<std::vector<Track*>>((*DataStore)["TracksForPID"]);
  for (auto electron : reconstructed_electrons) {
    tracks_for_PID.erase(std::find(tracks_for_PID.begin(), tracks_for_PID.end(), electron));
  }
  (*DataStore)["TracksForPID"] = tracks_for_PID;
  
  
  // store the electrons
  (*DataStore)["Electrons"] = reconstructed_electrons;


  return true;
}


bool ElectronPIDModule::ElectronPID(Track* track, float eIDprob, float separation)
{
  bool TrackIsElectron = false;


  // Apply a basic parameterized electron PID to tracks. If the track is really a
  // electron from truth information, apply a flat ID probability. If it's a pion,
  // use the separation (in Gaussian sigma) to determine if it's mis-identified.

  if (track->PT < 0.1) 
    return TrackIsElectron;

  int track_truth = track->PID;

  if (TMath::Abs(track_truth) == 11) {
    // true charged electron
    if (gRandom->Uniform(0, 1) <= eIDprob) {
      TrackIsElectron = true;
    } else {
      TrackIsElectron = false;
    }

    

  } else if (TMath::Abs(track_truth) == 211) {
    // true charged pion
    if (gRandom->Uniform(0,1) <= ROOT::Math::gaussian_pdf(separation)) {
      TrackIsElectron = true;
    } else {
      TrackIsElectron = false;
    }

  } else {
    // ignore ALL other species for now
    TrackIsElectron = false;
  }


  return TrackIsElectron;
}




