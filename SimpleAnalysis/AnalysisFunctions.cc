#ifndef ANALYSISFUNCTIONS
#define ANALYSISFUNCTIONS

#include "TROOT.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TClonesArray.h"
#include "TMath.h"

#include <vector>
#include <string>
#include <map>
#include <iostream>
#include <math.h>

#include "classes/DelphesClasses.h"

// Template functions for applying cuts and selecting subsets of particle lists

template<class T>std::vector<T *>                                  SelectorFcn(std::vector<T *>particles,
                                                              bool                             selector(T *))
{
  std::vector<T *> output_list;

  for (T *particle : particles) {
    if (selector(particle)) {
      output_list.push_back(particle);
    }
  }

  return output_list;
}

// Computations of fundamental quantities, like Bjorken x, etc.


inline std::map<std::string, float>DISVariables(TClonesArray *branchParticle)
{
  // four-momenta of proton, electron, virtual photon/Z^0/W^+-.
  auto pProton = static_cast<GenParticle *>(branchParticle->At(0))->P4(); // these
                                                                          // numbers
                                                                          // 0
                                                                          // ,
                                                                          // 3,
                                                                          // 5
                                                                          // are
                                                                          // hardcoded
                                                                          // in
                                                                          // Pythia8
  auto pleptonIn  = static_cast<GenParticle *>(branchParticle->At(3))->P4();
  auto pleptonOut = static_cast<GenParticle *>(branchParticle->At(5))->P4();
  auto pPhoton    = pleptonIn - pleptonOut;

  // Q2, W2, Bjorken x, y, nu.
  float Q2 = -pPhoton.M2();
  float W2 = (pProton + pPhoton).M2();
  float x  = Q2 / (2. * pProton.Dot(pPhoton));
  float y  = (pProton.Dot(pPhoton)) / (pProton.Dot(pleptonIn));


  std::map<std::string, float> dis_variables;
  dis_variables["Q2"] = Q2;
  dis_variables["W2"] = W2;
  dis_variables["x"]  = x;
  dis_variables["y"]  = y;

  return dis_variables;
}

inline std::map<std::string, float>DISJacquetBlondel(TClonesArray *tracks,
                                                     TClonesArray *electrons,
                                                     TClonesArray *photons,
                                                     TClonesArray *neutral_hadrons)
{
  // Jacquet-Blondel method:
  float delta_track = 0.0;
  auto  temp_p      = TVector3();

  for (int i = 0; i < tracks->GetEntries(); i++) {
    auto track_mom = static_cast<Track *>(tracks->At(i))->P4();

    if (isnan(track_mom.E()))
      continue;
    delta_track += (track_mom.E() - track_mom.Pz());
    temp_p       = temp_p + track_mom.Vect();
  }

  float delta_track_noel = 0.0;

  if (electrons->GetEntries() > 0) {
    auto e = static_cast<Track *>(electrons->At(0))->P4();

    if (!isnan(e.E())) {
      delta_track_noel = delta_track - (e.E() - e.Pz());
    }
  }

  float delta_photon = 0.0;

  for (int i = 0; i < photons->GetEntries(); i++) {
    auto pf_mom = static_cast<Photon *>(photons->At(i))->P4();

    if (isnan(pf_mom.E()))
      continue;
    delta_photon += (pf_mom.E() - pf_mom.Pz());
    temp_p        = temp_p + pf_mom.Vect();
  }

  float delta_neutral          = 0;
  float delta_neutral_noBarrel = 0;
  auto  ptmiss_noBarrel        = TVector3();

  for (int i = 0; i < neutral_hadrons->GetEntries(); i++) {
    auto pf_mom = static_cast<Tower *>(neutral_hadrons->At(i))->P4();

    if (isnan(pf_mom.E()))
      continue;
    delta_neutral += (pf_mom.E() - pf_mom.Pz());
    temp_p         = temp_p + pf_mom.Vect();

    if (TMath::Abs(pf_mom.Eta()) > 1.0) {
      delta_neutral_noBarrel += (pf_mom.E() - pf_mom.Pz());
      ptmiss_noBarrel         = ptmiss_noBarrel + pf_mom.Vect();
    }
  }

  float delta = delta_track + delta_photon + delta_neutral;

  // delta_noel = delta_track_noel + delta_photon + delta_neutral
  // delta_noel_noBarrel = delta_track_noel + delta_photon +
  // delta_neutral_noBarrel
  float delta_noBarrel = delta_track + delta_photon + delta_neutral_noBarrel;

  float y_JB   = delta / (2.0 * 10.0);
  auto  ptmiss = temp_p.Perp();
  float Q2_JB  = (ptmiss * ptmiss) / (1.0 - y_JB);
  float s      = 4.0 * 10.0 * 275.0;
  float x_JB   = Q2_JB / (s * y_JB);

  std::map<std::string, float> dis_variables;
  dis_variables["x_JB"]  = x_JB;
  dis_variables["Q2_JB"] = Q2_JB;
  dis_variables["y_JB"]  = y_JB;

  return dis_variables;
}

//
// Tagging Methods
//

inline bool IsTaggingTrack(Track *track)
{
  bool good_track = true;


  // Avoid decays that are too far from the beamspot (e.g. Ks, etc.)
  float d0 = TMath::Abs(track->D0);
  float z0 = TMath::Abs(track->DZ);

  // In the original approach, we merely restricted d0 to be < 3mm
  // However, with Ks and Lambda, etc. long-lived light decays turned on,
  // more displaced tracks started contaminating this in light jets.
  // (e.g. with small d0 but large z0)
  // Use a 3D impact parameter maximum to constrain these tracks.

  float r0 = TMath::Sqrt(d0 * d0 + z0 * z0);

  if (r0 > 3.0) // mm
    good_track &= false;

  return good_track;
}

inline float sIP3D(Jet *jet, Track *track)
{
  const TLorentzVector& jetMomentum = jet->P4();
  float jpx                         = jetMomentum.Px();
  float jpy                         = jetMomentum.Py();
  float jpz                         = jetMomentum.Pz();

  float d0 = TMath::Abs(track->D0);

  float xd  = track->Xd;
  float yd  = track->Yd;
  float zd  = track->Zd;
  float dd0 = TMath::Abs(track->ErrorD0);
  float dz  = TMath::Abs(track->DZ);
  float ddz = TMath::Abs(track->ErrorDZ);

  int sign = (jpx * xd + jpy * yd + jpz * zd > 0.0) ? 1 : -1;

  // add transverse and longitudinal significances in quadrature
  float sip = sign * TMath::Sqrt(TMath::Power(d0 / dd0, 2) + TMath::Power(dz / ddz, 2));

  // if (d0 > 3) { // 3mm maximum in d0
  if (!IsTaggingTrack(track) || TMath::IsNaN(sip) || !TMath::Finite(sip)) {
    // Tracks far outside the beamspot region should be ignored (Ks, Lambda,
    // etc.); tracks with bad d0, z0 error values, etc. should be ignored.
    sip = -199.0;
  }

  return sip;
}

inline bool Tagged_Kaon(Jet *jet, std::vector<Track *>kaons, float minSignif, float minPT, int minKaons)
{
  bool tagged     = false;
  int  kaon_count = 0;

  for (auto kaon : kaons) {
    if (kaon->P4().DeltaR(jet->P4()) > 0.5)
      continue;

    if (kaon->PT < minPT)
      continue;

    if (!IsTaggingTrack(kaon))
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

inline bool Tagged_Electron(Jet *jet, std::vector<Track *>electrons, float minSignif, float minPT, int minElectrons)
{
  bool tagged         = false;
  int  electron_count = 0;

  for (auto electron : electrons) {
    if (electron->P4().DeltaR(jet->P4()) > 0.5)
      continue;

    if (electron->PT < minPT)
      continue;

    if (!IsTaggingTrack(electron))
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

inline bool Tagged_Muon(Jet *jet, std::vector<Track *>muons, float minSignif, float minPT, int minMuons)
{
  bool tagged     = false;
  int  muon_count = 0;

  for (auto muon : muons) {
    if (muon->P4().DeltaR(jet->P4()) > 0.5)
      continue;

    if (muon->PT < minPT)
      continue;

    if (!IsTaggingTrack(muon))
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

inline bool Tagged_sIP3D(Jet *jet, TClonesArray tracks,
                         float minSignif, float minPT, int minTracks)
{
  bool tagged = false;

  const TLorentzVector& jetMomentum = jet->P4();
  float jpx                         = jetMomentum.Px();
  float jpy                         = jetMomentum.Py();
  float jpz                         = jetMomentum.Pz();

  auto jet_constituents = tracks;

  int N_sIPtrack = 0;

  // std::cout << "------------ Processing Jet " << jet << " ------------" <<
  // std::endl;

  for (int iconst = 0; iconst < jet_constituents.GetEntries(); iconst++) {
    if (N_sIPtrack >= minTracks)
      break;

    auto constituent = jet_constituents.At(iconst);

    if (constituent == 0) continue;

    if (constituent->IsA() == Track::Class()) {
      auto track = static_cast<Track *>(constituent);

      const TLorentzVector& trkMomentum = track->P4();
      float tpt                         = trkMomentum.Pt();

      if (tpt < minPT) continue;

      if (trkMomentum.DeltaR(jetMomentum) > 0.5)
        continue;

      // Avoid decays that are too far from the beamspot (e.g. Ks, etc.)
      // float d0 = TMath::Abs(track->D0);
      // float z0 = TMath::Abs(track->DZ);
      if (!IsTaggingTrack(track))
        continue;

      float sip = sIP3D(jet, track);


      if (sip > minSignif) N_sIPtrack++;
    }
  }

  tagged = (N_sIPtrack >= minTracks);

  // std::cout << "------------ xxxxxxxxxxxxxx ------------" << std::endl;

  return tagged;
}

#endif // ifndef ANALYSISFUNCTIONS
