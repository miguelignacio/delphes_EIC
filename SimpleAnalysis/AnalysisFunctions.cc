#ifndef ANALYSISFUNCTIONS
#define ANALYSISFUNCTIONS

#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TClonesArray.h"

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <math.h>

#include "classes/DelphesClasses.h"

// Template functions for applying cuts and selecting subsets of particle lists

template <class T> std::vector<T*> SelectorFcn(std::vector<T*> particles, bool selector(T*)) 
{
  std::vector<T*> output_list;
  for (T* particle : particles) {
    if (selector(particle)) {
      output_list.push_back(particle);
    }
  }

  return output_list;
}


// Computations of fundamental quantities, like Bjorken x, etc.


inline std::map<std::string, float> DISVariables(TClonesArray *branchParticle)
{
  // four-momenta of proton, electron, virtual photon/Z^0/W^+-.
  auto pProton      = static_cast<GenParticle *>(branchParticle->At(0))->P4(); //these numbers 0 , 3, 5 are hardcoded in Pythia8
  auto pleptonIn    = static_cast<GenParticle *>(branchParticle->At(3))->P4();
  auto pleptonOut   = static_cast<GenParticle *>(branchParticle->At(5))->P4();
  auto pPhoton      = pleptonIn - pleptonOut;
  
  // Q2, W2, Bjorken x, y, nu.
  float Q2 = -pPhoton.M2();
  float W2 = (pProton + pPhoton).M2();
  float x = Q2 / (2. * pProton.Dot(pPhoton));
  float y = (pProton.Dot(pPhoton)) / (pProton.Dot(pleptonIn));

  
  std::map<std::string, float> dis_variables;
  dis_variables["Q2"] = Q2;
  dis_variables["W2"] = W2;
  dis_variables["x"] = x;
  dis_variables["y"] = y;
  
  return dis_variables;
}

inline std::map<std::string, float> DISJacquetBlondel(TClonesArray* tracks, 
						      TClonesArray* electrons,
						      TClonesArray* photons,
						      TClonesArray* neutral_hadrons)
{
  // Jacquet-Blondel method:
  float delta_track = 0.0;
  auto temp_p = TVector3();
  for (int i = 0 ; i < tracks->GetEntries(); i++) {
    auto track_mom = static_cast<Track*>(tracks->At(i))->P4();
    if (isnan(track_mom.E()))
      continue;
    delta_track += (track_mom.E() - track_mom.Pz());
    temp_p = temp_p + track_mom.Vect();
  }

  float delta_track_noel = 0.0;
  if (electrons->GetEntries()>0) {
    auto e = static_cast<Track*>(electrons->At(0))->P4();
    if (!isnan(e.E())) {
      delta_track_noel = delta_track - (e.E() - e.Pz());
    }
  }
  
  float delta_photon = 0.0;
  for (int i = 0; i < photons->GetEntries(); i++) {
    auto pf_mom = static_cast<Photon*>(photons->At(i))->P4();
    if (isnan(pf_mom.E()))
      continue;
    delta_photon += (pf_mom.E() - pf_mom.Pz());
    temp_p = temp_p + pf_mom.Vect();
  }

  float delta_neutral = 0;
  float delta_neutral_noBarrel = 0;
  auto ptmiss_noBarrel = TVector3();
    
  for (int i = 0; i < neutral_hadrons->GetEntries(); i++) {
    auto pf_mom = static_cast<Tower*>(neutral_hadrons->At(i))->P4();
    if (isnan(pf_mom.E()))
      continue;
    delta_neutral += (pf_mom.E() - pf_mom.Pz());
    temp_p = temp_p+ pf_mom.Vect();
    if (TMath::Abs(pf_mom.Eta())>1.0) {
      delta_neutral_noBarrel += (pf_mom.E() - pf_mom.Pz());
      ptmiss_noBarrel = ptmiss_noBarrel + pf_mom.Vect();
    }
  }

  float delta = delta_track + delta_photon + delta_neutral;
  //delta_noel = delta_track_noel + delta_photon + delta_neutral
  //delta_noel_noBarrel = delta_track_noel + delta_photon + delta_neutral_noBarrel
  float delta_noBarrel = delta_track + delta_photon + delta_neutral_noBarrel;
    
  float y_JB   = delta/(2.0*10.0);
  auto ptmiss = temp_p.Perp();
  float Q2_JB  = (ptmiss*ptmiss)/(1.0-y_JB);
  float s     = 4.0*10.0*275.0;
  float x_JB  = Q2_JB/(s*y_JB);

  std::map<std::string, float> dis_variables;
  dis_variables["x_JB"] = x_JB;
  dis_variables["Q2_JB"] = Q2_JB;
  dis_variables["y_JB"] = y_JB;
  
  return dis_variables;
}

#endif
