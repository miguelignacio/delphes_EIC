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

isCC = False
isNC = False

inputFile = sys.argv[1]
if 'CC' in inputFile:
    print 'CC analysis'
    isCC = True
elif 'NC in inputFile':
    print 'NC analysis'
    isNC = True

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

##Tracks/Towers
branchTrack = treeReader.UseBranch("Track")
branchTower = treeReader.UseBranch("Tower")
##Particle flow objects
branchEFlowTrack =  treeReader.UseBranch("EFlowTrack")
branchEFlowPhoton = treeReader.UseBranch("EFlowPhoton")
branchEFlowNeutralHadron = treeReader.UseBranch("EFlowNeutralHadron")


# Book histograms
histJetPT = ROOT.TH1F("jet_pt", "jet P_{T}", 100, 0.0, 100.0)
METMatrix = ROOT.TH2F("METMatrix", "Met Matrix", 100, 10.0, 40.0, 100, 10.0, 40.0)
JetMatrix = ROOT.TH2F("JetMatrix", "Jet Matrix", 100, 10.0, 40.0, 100, 10.0, 40.0)
ElectronMatrix = ROOT.TH2F("ElectronMatrix", "Electron Matrix", 100,10.0,40.0,100,10.0,40.0)

y_Matrix = ROOT.TH2F("y_Matrix", "inelasticity response matrix, JB method", 100, 0.0,1.0, 100,0.0,1.0)
x_Matrix = ROOT.TH2F("x_Matrix", "Bjorken x response matrix, JB method", 100, 0.01,1.0, 100,0.01,1.0)
Q2_Matrix = ROOT.TH2F("Q2_Matrix", "Q2 response matrix, JB method", 100, 100,10000, 100,100,10000) 

JetMatrix_profile = ROOT.TH2F("JetMatrix_profile", "Jet Matrix profile", 100, 10.0, 40.0, 100, -1.0, 1.0)
profile_Jet = ROOT.TProfile("profile_Jet", "profile_Jet", 30, 10, 40, -1.0,1.0,"s") 

METMatrix_profile = ROOT.TH2F("METMatrix_profile", "MET Matrix profile", 100, 10.0, 40.0, 100, -1.0, 1.0)
profile_MET       = ROOT.TProfile("profile_MET", "profile_MET", 30, 10, 40, -1.0,1.0,"s")     

ElectronMatrix_profile = ROOT.TH2F("ElectronMatrix_profile", "Electron Matrix profile", 100, 10.0, 40.0, 100, -.20, .20)
profile_Electron       = ROOT.TProfile("profile_Electron", "profile_Electron", 30, 10, 40, -.20,0.20,"s")

dphi_reco_bin1 = ROOT.TH1F("dphi_reco_bin1", "dphi reco #1", 20, 2.8, ROOT.TMath.Pi())
dphi_reco_bin2 = ROOT.TH1F("dphi_reco_bin2", "dphi reco #2", 20, 2.8, ROOT.TMath.Pi())
dphi_reco_bin3 = ROOT.TH1F("dphi_reco_bin3", "dphi reco #3", 20, 2.8, ROOT.TMath.Pi())
dphi_reco_bin4 = ROOT.TH1F("dphi_reco_bin4", "dphi reco #4", 20, 2.8, ROOT.TMath.Pi())

dphi_truth_bin1 = ROOT.TH1F("dphi_truth_bin1", "dphi truth #1", 20, 2.8, ROOT.TMath.Pi())
dphi_truth_bin2 = ROOT.TH1F("dphi_truth_bin2", "dphi truth #2", 20, 2.8, ROOT.TMath.Pi())
dphi_truth_bin3 = ROOT.TH1F("dphi_truth_bin3", "dphi truth #3", 20, 2.8, ROOT.TMath.Pi())  
dphi_truth_bin4 = ROOT.TH1F("dphi_truth_bin4", "dphi truth #4", 20, 2.8, ROOT.TMath.Pi())

dphi_res_bin1 = ROOT.TH1F("dphi_res_bin1", "", 50, -0.3,0.3)
dphi_res_bin2 =  ROOT.TH1F("dphi_res_bin2", "", 50, -0.3,0.3)
dphi_res_bin3 = ROOT.TH1F("dphi_res_bin3", "", 50, -0.3,0.3)
dphi_res_bin4 = ROOT.TH1F("dphi_res_bin4", "", 50, -0.3,0.3)  

ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetPalette(112)

c = ROOT.TCanvas("c", "c", 800,600)

JetMatrix.SetTitle("; generated p_{T} [GeV]; reconstructed p_{T} [GeV]")
METMatrix.SetTitle("; generated MET [GeV]; reconstructed MET [GeV]")
ElectronMatrix.SetTitle("; generated p_{T} [GeV]; reconstructed p_{T} [GeV]")

JetMatrix_profile.SetTitle("Jet response matrix; p_{T}^{gen} [GeV]; (p_{T}^{reco}-p_{T}^{gen})/p_{T}^{gen}")
METMatrix_profile.SetTitle("Met response matrix; MET_{T}^{gen} [GeV]; (MET_{T}^{reco}-MET_{T}^{gen})/MET_{T}^{gen}")          
ElectronMatrix_profile.SetTitle("Electron response matrix; p_{T}^{gen} [GeV]; (p_{T}^{reco}-p_{T}^{gen})/p_{T}^{gen}")  

Q2_Matrix.SetTitle(" Q^{2} response matrix, JB method ; generated Q^{2} [GeV^{2}]; reconstructed Q^{2} [GeV^{2}]")
x_Matrix.SetTitle(" x response matrix, JB method; generated x; reconstructed x ")
y_Matrix.SetTitle(" y response matrix, JB method; generated y, reconstructed y")

# Loop over all events
for entry in range(0, numberOfEntries):
  if entry%10000==0:
      print 'event ' , entry
  # Load selected branches with data from specified event
  treeReader.ReadEntry(entry)

  # four-momenta of proton, electron, virtual photon/Z^0/W^+-.
  pProton      = branchParticle.At(0).P4(); #these numbers 0 , 3, 5 are hardcoded in Pythia8
  pleptonIn    = branchParticle.At(3).P4();
  pleptonOut   = branchParticle.At(5).P4();
  pPhoton      = pleptonIn - pleptonOut;
  
   #Q2, W2, Bjorken x, y, nu.
  Q2 = -pPhoton.M2()
  W2 = (pProton + pPhoton).M2()
  x = Q2 / (2. * pProton.Dot(pPhoton))
  y = (pProton.Dot(pPhoton)) / (pProton.Dot(pleptonIn))
  
  #Jacquet Blondet method:
  temp = 0
  temp_p = ROOT.TVector3()
  for i in range(branchEFlowTrack.GetEntries()):
      track_mom = branchEFlowTrack.At(i).P4()
      temp = temp + (track_mom.E() - track_mom.Pz())
      temp_p = temp_p + track_mom.Vect()
      
  for i in range(branchEFlowPhoton.GetEntries()):
      pf_mom = branchEFlowPhoton.At(i).P4()        
      temp = temp + (pf_mom.E() - pf_mom.Pz())
      temp_p = temp_p + pf_mom.Vect()

  for i in range(branchEFlowNeutralHadron.GetEntries()):
      pf_mom = branchEFlowNeutralHadron.At(i).P4()
      temp = temp + (pf_mom.E() - pf_mom.Pz())
      temp_p = temp_p + pf_mom.Vect()
      

  y_JB   = temp/(2.0*10.0)
  ptmiss = temp_p.Perp()
  Q2_JB  = (ptmiss*ptmiss)/(1-y_JB)
  s     = 4*10.0*275.0
  x_JB  = Q2_JB/(s*y_JB)

 # print 'Kinematic variables ', x_JB, ' x, ', Q2_JB, ' , Q2 , ', y_JB,' , y'
 # print 'Kinematic variables ', x, ' x, ', Q2, ' , Q2 , ', y,' , y'
 # print 'Missing momentum ' , ptmiss

  x_Matrix.Fill(x, x_JB)
  Q2_Matrix.Fill(Q2, Q2_JB)
  y_Matrix.Fill(y, y_JB)
 
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

      #print 'MET', met.MET
      #print 'GenMET', genmet.MET
      if branchGenMet.GetEntries()>0:
          METMatrix.Fill(genmet.MET, met.MET)
          res = (met.MET-genmet.MET)/genmet.MET
          profile_MET.Fill(genmet.MET,res)
          METMatrix_profile.Fill(genmet.MET, res)

  leptonOK = False
  if(isCC and branchMet.GetEntries()>0 and branchGenMet.GetEntries()>0):
      lepton = branchMet.At(0).P4()
      genlepton = branchGenMet.At(0).P4()
      lepton_pt = lepton.Pt()
      genlepton_pt = lepton.Pt()
      leptonOK = True
  elif(isNC and branchElectron.GetEntries()>0):
      lepton = branchElectron.At(0).P4()
      genlepton = branchElectron.At(0).P4()
      lepton_pt = lepton.Pt()
      genlepton_pt = lepton.Pt() 
      leptonOK = True
      
  if leptonOK and branchJet.GetEntries() > 0 and branchGenJet.GetEntries()>0:
      jet = branchJet.At(0)
      genjet = branchGenJet.At(0)
      genlepton_phi = genlepton.Phi()
      lepton_phi =    lepton.Phi()
      
      jet_phi =    jet.P4().Phi()
      genjet_phi=  genjet.P4().Phi()
      
      dphi_gen = abs(ROOT.TVector2.Phi_mpi_pi(genlepton_phi-genjet_phi))
      dphi_reco = abs(ROOT.TVector2.Phi_mpi_pi(lepton_phi-jet_phi))

      #print genjet.P4().Pt(), ' ' ,jet.P4().Pt()
      if(genlepton_pt>10 and genlepton_pt<15):
          dphi_reco_bin1.Fill(dphi_reco)
          dphi_truth_bin1.Fill(dphi_gen)
          dphi_res_bin1.Fill(dphi_reco-dphi_gen)     
      elif(genlepton_pt>15 and genlepton_pt<20):
          dphi_reco_bin2.Fill(dphi_reco)
          dphi_truth_bin2.Fill(dphi_gen)
          dphi_res_bin2.Fill(dphi_reco-dphi_gen)     
      elif(genlepton_pt>20 and genlepton_pt<30):
          dphi_reco_bin3.Fill(dphi_reco)
          dphi_truth_bin3.Fill(dphi_gen)
          dphi_res_bin3.Fill(dphi_reco-dphi_gen)     
      elif(genlepton_pt>30 and genlepton_pt<40):
          dphi_reco_bin4.Fill(dphi_reco)
          dphi_truth_bin4.Fill(dphi_gen)    
          dphi_res_bin4.Fill(dphi_reco-dphi_gen)
                  

###Colors
profile_Jet.SetLineColor(2)
profile_Jet.SetLineWidth(2)
profile_MET.SetLineColor(2)
profile_MET.SetLineWidth(2)
profile_Electron.SetLineColor(2)
profile_Electron.SetLineWidth(2)

y_Matrix.Draw("colz")
c.SaveAs("y_response_%s.png"%inputFile)
c.SaveAs("y_response_%s.pdf"%inputFile)

c.Clear()
x_Matrix.Draw("colz")
c.SaveAs("x_response_%s.png"%inputFile)
c.SaveAs("x_response_%s.pdf"%inputFile)

c.Clear()

Q2_Matrix.Draw("colz")
c.SetLogy(1)
c.SetLogx(1)
c.SaveAs("Q2_response_%s.png"%inputFile)
c.SaveAs("Q2_response_%s.pdf"%inputFile)
c.SetLogy(0)
c.SetLogx(0)


# Show resulting histograms
histJetPT.Draw()
#input("Press Enter to continue...")
c.Clear()
JetMatrix.Draw("colz")
c.SaveAs("JetResponse_%s.png"%inputFile)
c.SaveAs("JetResponse_%s.pdf"%inputFile)

#input("Press Enter to continue...")
c.Clear()
ElectronMatrix.Draw("colz")
c.SaveAs("ElectronResponse_%s.png"%inputFile)    
c.SaveAs("ElectronResponse_%s.pdf"%inputFile)

#input("Press enter to continue...")
c.Clear()
JetMatrix_profile.Draw("colz")
profile_Jet.Draw("same")
c.SaveAs("profile_jet_%s.png"%inputFile)    
c.SaveAs("profile_jet_%s.pdf"%inputFile)

#input("Press enter to continue...")
c.Clear()
METMatrix_profile.Draw("colz")
profile_MET.Draw("same")
c.SaveAs("profile_MET_%s.png"%inputFile)
c.SaveAs("profile_MET_%s.pdf"%inputFile)


#input("Press enter to continue...")
c.Clear()
ElectronMatrix_profile.Draw("colz")
profile_Electron.Draw("same")
c.SaveAs("profile_Electron_%s.png"%inputFile)
c.SaveAs("profile_Electron_%s.pdf"%inputFile)

c.Clear()
METMatrix.Draw("colz")
c.SaveAs("METResponse_%s.png"%inputFile)
c.SaveAs("METResponse_%s.pdf"%inputFile)



c.Clear()
dphi_truth_bin1.SetTitle("; #Delta#phi [rad]; 1/#sigma d#sigma/d#Delta#phi")
dphi_truth_bin1.SetLineColor(2)
dphi_truth_bin1.SetLineWidth(2)
dphi_truth_bin1.DrawNormalized()
dphi_reco_bin1.DrawNormalized("same")
c.SaveAs("Dphi_bin1_%s.png"%inputFile)
c.SaveAs("Dphi_bin1_%s.pdf"%inputFile)


c.Clear()
dphi_truth_bin2.SetTitle("; #Delta#phi [rad]; 1/#sigma d#sigma/d#Delta#phi")
dphi_truth_bin2.SetLineColor(2)
dphi_truth_bin2.SetLineWidth(2)
dphi_truth_bin2.DrawNormalized()
dphi_reco_bin2.DrawNormalized("same")
c.SaveAs("Dphi_bin2_%s.png"%inputFile)
c.SaveAs("Dphi_bin2_%s.pdf"%inputFile)

c.Clear()
dphi_truth_bin3.SetTitle("; #Delta#phi [rad]; 1/#sigma d#sigma/d#Delta#phi")
dphi_truth_bin3.SetLineColor(2)
dphi_truth_bin3.SetLineWidth(2)
dphi_truth_bin3.DrawNormalized()
dphi_reco_bin3.DrawNormalized("same")
c.SaveAs("Dphi_bin3_%s.png"%inputFile)
c.SaveAs("Dphi_bin3_%s.pdf"%inputFile)     

c.Clear()
dphi_truth_bin4.SetTitle("; #Delta#phi [rad]; 1/#sigma d#sigma/d#Delta#phi")
dphi_truth_bin4.SetLineColor(2)
dphi_truth_bin4.SetLineWidth(2)
dphi_truth_bin4.DrawNormalized()
dphi_reco_bin4.DrawNormalized("same")

c.SaveAs("Dphi_bin4_%s.png"%inputFile)
c.SaveAs("Dphi_bin4_%s.pdf"%inputFile)


c.Clear()

c.Divide(2,2)

dphi_res_bin1.SetTitle("; #Delta#phi_{reco} - #Delta#phi_{truth}; entries")
dphi_res_bin2.SetTitle("; #Delta#phi_{reco} - #Delta#phi_{truth}; entries")
dphi_res_bin3.SetTitle("; #Delta#phi_{reco} - #Delta#phi_{truth}; entries")
dphi_res_bin4.SetTitle("; #Delta#phi_{reco} - #Delta#phi_{truth}; entries") 

dphi_res_bin1.SetLineColor(1)
dphi_res_bin2.SetLineColor(2)
dphi_res_bin3.SetLineColor(4)
dphi_res_bin4.SetLineColor(8)

c.cd(1)
dphi_res_bin1.Fit("gaus")
dphi_res_bin1.Draw()
sigma = dphi_res_bin1.GetFunction("gaus").GetParameter(2)
print sigma
c.cd(2)
dphi_res_bin2.DrawNormalized("same")
c.cd(3)
dphi_res_bin3.Fit("gaus")
dphi_res_bin3.DrawNormalized("same")
c.cd(4)
dphi_res_bin4.DrawNormalized("same")
c.SaveAs("DphiResolution_%s.png"%inputFile)
c.SaveAs("DphiResolution_%s.pdf"%inputFile)
