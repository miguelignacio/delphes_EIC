#!/usr/bin/env python

import sys

import ROOT

try:
  input = raw_input
except:
  pass

if len(sys.argv) < 2:
  print(" Usage: Example1.py input_file")
  sys.exit(1)

ROOT.gSystem.Load("libDelphes")

try:
  ROOT.gInterpreter.Declare('#include "classes/DelphesClasses.h"')
  ROOT.gInterpreter.Declare('#include "external/ExRootAnalysis/ExRootTreeReader.h"')
except:
  pass

inputFile = sys.argv[1]

# Create chain of root trees
chain = ROOT.TChain("Delphes")
chain.Add(inputFile)

# Create object of class ExRootTreeReader
treeReader = ROOT.ExRootTreeReader(chain)
numberOfEntries = treeReader.GetEntries()

# Get pointers to branches used in this analysis
branchParticle = treeReader.UseBranch("Particle")
branchJet = treeReader.UseBranch("Jet")
branchGenJet = treeReader.UseBranch("GenJet")
branchElectron = treeReader.UseBranch("Electron")
branchMet = treeReader.UseBranch("MissingET")
branchGenMet = treeReader.UseBranch("GenMissingET")

# Book histograms
histJetPT = ROOT.TH1F("jet_pt", "jet P_{T}", 100, 0.0, 100.0)
METMatrix = ROOT.TH2F("METMatrix", "Met Matrix", 100, 10.0, 40.0, 100, 10.0, 40.0)
JetMatrix = ROOT.TH2F("JetMatrix", "Jet Matrix", 100, 10.0, 40.0, 100, 10.0, 40.0)
ElectronMatrix = ROOT.TH2F("ElectronMatrix", "Electron Matrix", 100,10.0,40.0,100,10.0,40.0)




JetMatrix_profile = ROOT.TH2F("JetMatrix_profile", "Jet Matrix profile", 100, 10.0, 40.0, 100, -1.0, 1.0)
profile_Jet = ROOT.TProfile("profile_Jet", "profile_Jet", 30, 10, 40, -1.0,1.0,"s") 

METMatrix_profile = ROOT.TH2F("METMatrix_profile", "MET Matrix profile", 100, 10.0, 40.0, 100, -1.0, 1.0)
profile_MET       = ROOT.TProfile("profile_MET", "profile_MET", 30, 10, 40, -1.0,1.0,"s")     

ElectronMatrix_profile = ROOT.TH2F("ElectronMatrix_profile", "Electron Matrix profile", 100, 10.0, 40.0, 100, -.20, .20)
profile_Electron       = ROOT.TProfile("profile_Electron", "profile_Electron", 30, 10, 40, -.20,0.20,"s")



ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetPalette(112)

c = ROOT.TCanvas("c", "c", 800,600)

JetMatrix.SetTitle("; generated p_{T} [GeV]; reconstructed p_{T} [GeV]")
METMatrix.SetTitle("; generated MET [GeV]; reconstructed MET [GeV]")
ElectronMatrix.SetTitle("; generated p_{T} [GeV]; reconstructed p_{T} [GeV]")

JetMatrix_profile.SetTitle("Jet response matrix; p_{T}^{gen} [GeV]; (p_{T}^{reco}-p_{T}^{gen})/p_{T}^{gen}")
METMatrix_profile.SetTitle("Met response matrix; MET_{T}^{gen} [GeV]; (MET_{T}^{reco}-MET_{T}^{gen})/MET_{T}^{gen}")          
ElectronMatrix_profile.SetTitle("Electron response matrix; p_{T}^{gen} [GeV]; (p_{T}^{reco}-p_{T}^{gen})/p_{T}^{gen}")  


# Loop over all events
for entry in range(0, numberOfEntries):
  # Load selected branches with data from specified event
  treeReader.ReadEntry(entry)

  ## Electron response matrix:
  if branchElectron.GetEntries()>0:
      electron = branchElectron.At(0)
      if branchParticle.GetEntries()>0:
          for i in range(branchParticle.GetEntries()):
              particle = branchParticle.At(i)
              gen_mom = particle.P4()
              pid = particle.PID
              status = particle.Status
              if (pid==11 and status==1 ):
                  gen_electron = particle
                  #print 'electron pT', gen_electron.Pt()
                  ElectronMatrix.Fill(gen_electron.PT, electron.PT)
                  res = (electron.PT-gen_electron.PT)/gen_electron.PT
                  profile_Electron.Fill(gen_electron.PT, res)
                  ElectronMatrix_profile.Fill(gen_electron.PT,res)
                                                                                            
                                                                                            
  ##Jet response matrix
  if branchJet.GetEntries() > 0:
      # Take first jet
      jet = branchJet.At(0)
      genjet = branchGenJet.At(0)
      # Plot jet transverse momentum
      histJetPT.Fill(jet.PT)
      # Print jet transverse momentum
      if branchGenJet.GetEntries()>0:
          JetMatrix.Fill(genjet.PT, jet.PT)
          res = (jet.PT-genjet.PT)/genjet.PT
          profile_Jet.Fill(genjet.PT, res)
          JetMatrix_profile.Fill(genjet.PT,res)
  ##Met response matrix
  if branchMet.GetEntries() > 0:
      met = branchMet.At(0)
      genmet = branchGenMet.At(0)

      if branchGenMet.GetEntries()>0:
          METMatrix.Fill(genmet.MET, met.MET)
          res = (met.MET-genmet.MET)/genmet.MET
          profile_MET.Fill(genmet.MET,res)
          METMatrix_profile.Fill(genmet.MET, res)

###Colors
profile_Jet.SetLineColor(2)
profile_Jet.SetLineWidth(2)
profile_MET.SetLineColor(2)
profile_MET.SetLineWidth(2)
profile_Electron.SetLineColor(2)
profile_Electron.SetLineWidth(2)

# Show resulting histograms
histJetPT.Draw()
input("Press Enter to continue...")
c.Clear()
JetMatrix.Draw("colz")
c.SaveAs("JetResponse.png")

input("Press Enter to continue...")
c.Clear()
ElectronMatrix.Draw("colz")
c.SaveAs("ElectronResponse.png")


input("Press enter to continue...")
c.Clear()
JetMatrix_profile.Draw("colz")
profile_Jet.Draw("same")
c.SaveAs("profile_jet.png")

input("Press enter to continue...")
c.Clear()
METMatrix_profile.Draw("colz")
profile_MET.Draw("same")
c.SaveAs("profile_MET%s.png"%inputFile)

input("Press enter to continue...")
c.Clear()
ElectronMatrix_profile.Draw("colz")
profile_Electron.Draw("same")
c.SaveAs("profile_Electron%.png"%inputFile)


input("Press enter to continue")
c.Clear()
METMatrix.Draw("colz")
c.SaveAs("METResponse.png")

