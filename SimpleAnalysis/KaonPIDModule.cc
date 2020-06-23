#include "KaonPIDModule.h"

#include "TClonesArray.h"
#include "TRandom3.h"
#include "Math/PdfFuncMathCore.h"

#include "AnalysisFunctions.cc"

#include <iostream>

KaonPIDModule::KaonPIDModule(ExRootTreeReader* data)
  : Module(data)
{

}

KaonPIDModule::~KaonPIDModule()
{

}



bool KaonPIDModule::execute(std::map<std::string, std::any>* DataStore)
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

  std::vector<Track*> all_kaons;


  if (DataStore->find("TracksForPID") != DataStore->end()) {
    for (auto track : std::any_cast<std::vector<Track*>>((*DataStore)["TracksForPID"])) {
      if (KaonPID(track, 0.90, 3.0))
	all_kaons.push_back(track);
    } 
    
    
  } else {      
    for (int itrk = 0; itrk < getEFlowTracks()->GetEntries(); itrk++)
      {
	Track* track = (Track*) getEFlowTracks()->At(itrk);
	if (KaonPID(track, 0.90, 3.0))
	  all_kaons.push_back(track);
      }

    std::vector<Track*> tracks_for_PID;
    for (int itrk = 0; itrk < getEFlowTracks()->GetEntries(); itrk++) 
      {
	Track* track = (Track*) getEFlowTracks()->At(itrk);
	tracks_for_PID.push_back(track);
      }
    (*DataStore)["TracksForPID"] = tracks_for_PID;
  }

  //auto reconstructed_kaons = all_kaons;
   std::vector<Track*> reconstructed_kaons = SelectorFcn<Track>(all_kaons, [](Track* t){ return (t->PT >= 0.1); });

   std::vector<Track*> tracks_for_PID = std::any_cast<std::vector<Track*>>((*DataStore)["TracksForPID"]);
   for (auto kaon : reconstructed_kaons) {
     tracks_for_PID.erase(std::find(tracks_for_PID.begin(), tracks_for_PID.end(), kaon));
   }
   (*DataStore)["TracksForPID"] = tracks_for_PID;
   

  (*DataStore)["Kaons"] = reconstructed_kaons;


  return true;
}


bool KaonPIDModule::KaonPID(Track* track, float kIDprob, float separation)
{
  bool TrackIsKaon = false;


  // Apply a basic parameterized kaon PID to tracks. If the track is really a
  // kaon from truth information, apply a flat ID probability. If it's a pion,
  // use the separation (in Gaussian sigma) to determine if it's mis-identified.

  if (track->PT < 0.1) 
    return TrackIsKaon;

  int track_truth = track->PID;

  if (TMath::Abs(track_truth) == 321) {
    // true charged kaon
    if (gRandom->Uniform(0, 1) <= kIDprob) {
      TrackIsKaon = true;
    } else {
      TrackIsKaon = false;
    }

    

  } else if (TMath::Abs(track_truth) == 211) {
    // true charged pion
    if (gRandom->Uniform(0,1) <= ROOT::Math::gaussian_pdf(separation)) {
      TrackIsKaon = true;
    } else {
      TrackIsKaon = false;
    }

  } else {
    // ignore ALL other species for now
    TrackIsKaon = false;
  }


  return TrackIsKaon;
}




