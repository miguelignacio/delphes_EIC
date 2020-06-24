#include "MuonPIDModule.h"

#include "TClonesArray.h"
#include "TRandom3.h"
#include "Math/PdfFuncMathCore.h"

#include "AnalysisFunctions.cc"

#include <iostream>

MuonPIDModule::MuonPIDModule(ExRootTreeReader* data)
  : Module(data)
{

}

MuonPIDModule::~MuonPIDModule()
{

}



bool MuonPIDModule::execute(std::map<std::string, std::any>* DataStore)
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

  std::vector<Track*> all_muons;

  // If the DataStore contains already a list of tracks to be used for PID assignment,
  // use that. If not, use the getEFlowTracks() method to get the general list of all tracks.
  // TracksForPID is special - it contains tracks NOT already used by another PID algorithm,
  // to avoid using the same track (pion) twice in two categories.

  float mu_efficiency = 0.95;
  float mupi_separation = 2.0; //sigma
  // float mu_efficiency = 1.0;
  // float mupi_separation = 100.0; //sigma
  
  
  if (DataStore->find("TracksForPID") != DataStore->end()) {
    for (auto track : std::any_cast<std::vector<Track*>>((*DataStore)["TracksForPID"])) {
      if (MuonPID(track, mu_efficiency, mupi_separation))
	all_muons.push_back(track);
    } 
    
    
  } else {      
    for (int itrk = 0; itrk < getEFlowTracks()->GetEntries(); itrk++)
      {
	Track* track = (Track*) getEFlowTracks()->At(itrk);
	if (MuonPID(track, mu_efficiency, mupi_separation))
	  all_muons.push_back(track);
      }

    std::vector<Track*> tracks_for_PID;
    for (int itrk = 0; itrk < getEFlowTracks()->GetEntries(); itrk++) 
      {
	Track* track = (Track*) getEFlowTracks()->At(itrk);
	tracks_for_PID.push_back(track);
      }
    (*DataStore)["TracksForPID"] = tracks_for_PID;
  }

  //auto reconstructed_muons = all_muons;
  std::vector<Track*> reconstructed_muons = SelectorFcn<Track>(all_muons, [](Track* t){ return (t->PT >= 0.1); });
  
  std::vector<Track*> tracks_for_PID = std::any_cast<std::vector<Track*>>((*DataStore)["TracksForPID"]);
  for (auto muon : reconstructed_muons) {
    tracks_for_PID.erase(std::find(tracks_for_PID.begin(), tracks_for_PID.end(), muon));
  }
  (*DataStore)["TracksForPID"] = tracks_for_PID;
   
  (*DataStore)["Muons"] = reconstructed_muons;


  return true;
}


bool MuonPIDModule::MuonPID(Track* track, float muIDprob, float separation)
{
  bool TrackIsMuon = false;


  // Apply a basic parameterized muon PID to tracks. If the track is really a
  // muon from truth information, apply a flat ID probability. If it's a pion,
  // use the separation (in Gaussian sigma) to determine if it's mis-identified.

  if (track->PT < 0.1) 
    return TrackIsMuon;

  int track_truth = track->PID;

  if (TMath::Abs(track_truth) == 13) {
    // true charged muon
    if (gRandom->Uniform(0, 1) <= muIDprob) {
      TrackIsMuon = true;
    } else {
      TrackIsMuon = false;
    }

    

  } else if (TMath::Abs(track_truth) == 211) {
    // true charged pion
    if (gRandom->Uniform(0,1) <= ROOT::Math::gaussian_pdf(separation)) {
      TrackIsMuon = true;
    } else {
      TrackIsMuon = false;
    }

  } else {
    // ignore ALL other species for now
    TrackIsMuon = false;
  }


  return TrackIsMuon;
}




