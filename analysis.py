#!/usr/bin/env python

import sys

import ROOT
from AtlasCommonUtils import *
#from Legend import Legend

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


#from ROOT import *
#ROOT.gROOT.LoadMacro("AtlasUtils.C")    

def DrawText(x, y, text, color=1, size=0.05):
    l = ROOT.TLatex()
    # l.SetTextAlign(12)
    l.SetTextSize(size)
    l.SetNDC()
    l.SetTextColor(color)
    l.DrawLatex(x, y, text)

isCC = False
isNC = False
SetAtlasStyle()

inputFile = sys.argv[1]
if 'CC' in inputFile:
    print 'CC analysis'
    isCC = True
elif 'NC' in inputFile:
    print 'NC analysis'
    isNC = True

if(isNC or isCC):
    print ' Either CC or NC analysis'
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

neutral_E = {}
photon_E  = {}
track_E   = {}
for i in range(1,5):
    neutral_E['etabin%i'%i] = ROOT.TH1F("neutral_E_eta%i"%i, "", 100, 0.0, 100.0)
    photon_E['etabin%i'%i]  = ROOT.TH1F("photon_E_eta%i"%i , "" , 100, 0.0,100.0)
    track_E['etabin%i'%i]   = ROOT.TH1F("track_E_eta%i"%i, "", 100, 0.0,100.0)


#Kinematic
y_Matrix = ROOT.TH2F("y_Matrix", "inelasticity response matrix, JB method", 10, 0.0,1.0, 10,0.0,1.0)
x_Matrix = ROOT.TH2F("x_Matrix", "Bjorken x response matrix, JB method", 10, 0.0,1.0, 10,0.0,1.0)
Q2_Matrix = ROOT.TH2F("Q2_Matrix", "Q2 response matrix, JB method", 20, 10,100, 20,10,100) 


## met vs x
Ngen_met_x = ROOT.TH2F("Ngen_met_x", "", 6, 10.0, 40.0, 5, 0.0, 1.0)
Nout_met_x = ROOT.TH2F("Nout_met_x", "", 6, 10.0, 40.0, 5, 0.0, 1.0)
Nin_met_x  = ROOT.TH2F("Nin_met_x", "" , 6, 10.0, 40.0, 5, 0.0, 1.0)

## y vs x 
Ngen_y_x = ROOT.TH2F("Ngen_y_x", "", 5,0.0,1.0, 5, 0.0, 1.0)
Nout_y_x = ROOT.TH2F("Nout_y_x", "", 5,0.0,1.0, 5, 0.0, 1.0)
Nin_y_x  = ROOT.TH2F("Nin_y_x", "" , 5,0.0,1.0, 5, 0.0, 1.0)

## Q2 vs x

import numpy as np
from array import array
binsQ2 = np.logspace(1.69897,4.2,8)
binsx = np.logspace(-2,0,11)
print 'BINS Q2', binsQ2
print 'BINS x ', binsx

Ngen_Q2_x = ROOT.TH2F("Ngen_Q2_x", "", 10, array('d',binsx), 6, array('d',binsQ2))
Nout_Q2_x = ROOT.TH2F("Nout_Q2_x", "", 10, array('d',binsx), 6, array('d',binsQ2))
Nin_Q2_x  = ROOT.TH2F("Nin_Q2_x", "" , 10, array('d',binsx), 6, array('d',binsQ2))

##E vs Eta
Jet_eta_e = ROOT.TH2F("Jet_eta_e", "" , 100,-4.0,4.0, 100, 0, 150.0)

### Jet PT and Phi
minpt = 5
maxpt = 45
nbinspt = 10


##Diagonal
METMatrix = ROOT.TH2F("METMatrix", "Met Matrix", 60, minpt, maxpt, 60, minpt, maxpt)
JetMatrix = ROOT.TH2F("JetMatrix", "Jet Matrix", 100, minpt, 100.0, 100, minpt, 100.0)
ElectronMatrix = ROOT.TH2F("ElectronMatrix", "Electron Matrix", 100,10.0,40.0,100,10.0,40.0)

profile = {}
ResMatrix = {}
histo = {} 
distribution = {}


ResMatrix['x'] = ROOT.TH2F("ResMatrix_x", "" , 10, 0.0,1.0, 50,-1.0,1.0)
profile['x']   = ROOT.TProfile("profile_x", "", 10, 0.0,1.0, -1.0,1.0, 's')

##JET, phi, phi
ResMatrix['jetpt'] =    ROOT.TH2F("ResMatrix_jetpt",             "",  nbinspt, minpt, maxpt, 50, -1.0, 1.0)
profile['jetpt'] =      ROOT.TProfile("profile_jetpt",           "",  nbinspt, minpt, maxpt, -1.0,1.0,"s")
ResMatrix['jetphi']  =  ROOT.TH2F("ResMatrix_jetphi" ,           "",  20, 5, 100, 100, -0.3,0.3)  
profile['jetphi'] =     ROOT.TProfile("profile_jetphi",          "",  20, 5, 100, -0.3,0.3,"s")

for i in range(1,5):
    ResMatrix['jete_eta%i'%i] = ROOT.TH2F("ResMatrix_jete_eta%i"%i,        "",  20, 10, 200, 50, -1.0,1.0)
    profile['jete_eta%i'%i]   = ROOT.TProfile("profile_jete_eta%i"%i,        "",  20, 10, 200, -1.0, 1.0,"s")

ResMatrix['jete'] = ROOT.TH2F("ResMatrix_jete",        "",  20, 5, 100, 50, -1.0,1.0)
profile['jete']    = ROOT.TProfile("profile_jete", "", 20, 5, 100, -1.0,1.0,"s")
    
### MET PT and Phi
ResMatrix['met']       = ROOT.TH2F("ResMatrix_met",              "", nbinspt, minpt, maxpt, 50, -1.0, 1.0)
profile['met']         = ROOT.TProfile("profile_met",            "", nbinspt, minpt, maxpt, -1.0,1.0,"s")     
ResMatrix['metphi']    = ROOT.TH2F("ResMatrix_metphi",           "", nbinspt, minpt, maxpt, 100, -0.5, 0.5)
profile['metphi']      = ROOT.TProfile("profile_metphi",         "", nbinspt, minpt, maxpt, -0.5, 0.5, "s")

## MET performance no HCAL
ResMatrix['met_nobarrelHCAL']       = ROOT.TH2F("ResMatrix_met_nobarrelHCAL",              "", nbinspt, minpt, maxpt, 50, -1.0, 1.0)
profile['met_nobarrelHCAL']         = ROOT.TProfile("profile_met_nobarrelHCAL",            "", nbinspt, minpt, maxpt, -1.0,1.0,"s")
ResMatrix['metphi_nobarrelHCAL']    = ROOT.TH2F("ResMatrix_metphi_nobarrelHCAL",           "", nbinspt, minpt, maxpt, 100, -0.5, 0.5)
profile['metphi_nobarrelHCAL']      = ROOT.TProfile("profile_metphi_nobarrelHCAL",         "", nbinspt, minpt, maxpt, -0.5, 0.5, "s")


##Electron
ResMatrix['ept'] = ROOT.TH2F("ResMatrix_ept", "", 100, 10.0, 40.0, 100, -.20, .20)
profile['ept']       = ROOT.TProfile("profile_ept", "", 30, 10, 40, -.20,0.20,"s")


maxdphires = 0.5
if(isNC):
    maxdphires = 0.2
ResMatrix['dphi'] = ROOT.TH2F("ResMatrix_dph", "", 20, 5, 100.0, 100, -maxdphires, maxdphires)
profile['dphi']  = ROOT.TProfile("profile_dphi", "", 20, 5,100.0, -maxdphires,maxdphires,"s")
distribution['dphi_reco'] = ROOT.TH2F("distribution_dphi_reco", "",  6, 10.0, 40.0, 20, 2.8, ROOT.TMath.Pi())
distribution['dphi_gen'] = ROOT.TH2F("distribution_dphi_gen", "",  6, 10.0, 40.0, 20, 2.8, ROOT.TMath.Pi())     


ResMatrix['qtnormjet'] = ROOT.TH2F("ResMatrix_qtnormjet", "", 9, 10.0, 100.0, 100, -maxdphires, maxdphires)
profile['qtnormjet']   = ROOT.TProfile("profile_qtnormjet", "", 9, 10.0, 100.0, -maxdphires, maxdphires,"s")
distribution['qtnormjet_reco'] = ROOT.TH2F("distribution_qtnormjet_reco", "", 9, 10.0, 100.0, 20, 0,1.0)
distribution['qtnormjet_gen'] = ROOT.TH2F("distribution_qtnormjet_gen", "",  9, 10.0, 100.0, 20, 0,1.0)


ResMatrix['qt'] = ROOT.TH2F("ResMatrix_qt", "", 9, 10.0, 100.0, 100, -7, 7)
profile['qt']   = ROOT.TProfile("profile_qt", "", 9, 10.0, 100.0, -7, 7,"s")



histo['delta_reco'] = ROOT.TH1F("delta_reco","", 100,0.0,30.0)
histo['delta_reco_noel'] = ROOT.TH1F("delta_reco_noel","", 100,0.0,30.0)
histo['delta_reco_noBarrelHCAL'] = ROOT.TH1F("delta_reco_noBarrelHCAL","", 100, 0.0,30.0)
histo['delta_gen'] = ROOT.TH1F("delta_gen","", 100,0.0,30.0)

histo['Vratio']  = ROOT.TH1F('Vratio','', 100,0.0, 1.0)
histo['Vratio_truth']  = ROOT.TH1F('Vratio_truth','', 100,0.0, 1.0)
histo['Vratio_noBarrelHCAL']  = ROOT.TH1F('Vratio_BarrelHCAL','', 100,0.0, 1.0)

h_qt_reco = {}
h_qt_truth = {}

for i in range(1,5):
    h_qt_reco['bin%i'%i] = ROOT.TH1F("qt_reco_bin%i"%i, "qt reco #%i"%i, 20, 0,0.8)
    h_qt_truth['bin%i'%i] = ROOT.TH1F("qt_truth_bin%i"%i, "qt truth #%i"%i, 20, 0,0.8) 

    

ROOT.gStyle.SetOptStat("")
ROOT.gStyle.SetPalette(112)

c = ROOT.TCanvas("c", "c", 800,600)

JetMatrix.SetTitle("; generated E_{jet} [GeV]; reconstructed E_{jet} [GeV]")
METMatrix.SetTitle("; generated MET [GeV]; reconstructed MET [GeV]")



ResMatrix['x'].SetTitle("Response Matrix for x; x_{gen}  ; (x_{reco}-x_{gen})/x_{gen}")
ResMatrix['jetpt'].SetTitle("Jet response matrix; p_{T}^{gen} [GeV]; (p_{T}^{reco}-p_{T}^{gen})/p_{T}^{gen}")
ResMatrix['met'].SetTitle("Met response matrix; MET_{T}^{gen} [GeV]; (MET^{reco}-MET^{gen})/MET^{gen}")          
ResMatrix['metphi'].SetTitle(" MET_{T}^{gen} [GeV]; (MET^{reco}-MET^{gen})/MET^{gen}")

Q2_Matrix.SetTitle(" Q^{2} response matrix, JB method ; generated Q^{2} [GeV^{2}]; reconstructed Q^{2} [GeV^{2}]")
x_Matrix.SetTitle(" x response matrix, JB method; generated x; reconstructed x ")
y_Matrix.SetTitle(" y response matrix, JB method; generated y; reconstructed y")

# Loop over all events
for entry in range(0, numberOfEntries):
    if entry%10000==0:
        print 'event ' , entry
    #if entry>5000:
    #   break
    # Load selected branches with data from specified event
    treeReader.ReadEntry(entry)


    x = 0
    y = 0
    Q2 = 0
    W2 = 0
    
    if(isCC or isNC):  
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
    delta_track = 0            
    temp_p = ROOT.TVector3()
    for i in range(branchEFlowTrack.GetEntries()):
       track_mom = branchEFlowTrack.At(i).P4()
       delta_track += (track_mom.E() - track_mom.Pz())
       temp_p = temp_p + track_mom.Vect()
       #print track_mom.E() , ' ' , track_mom.Pz()

    if branchElectron.GetEntries()>0:
        e = branchElectron.At(0).P4()
        delta_track_noel = delta_track - (e.E() - e.Pz())   

    delta_photon = 0
    for i in range(branchEFlowPhoton.GetEntries()):
       pf_mom = branchEFlowPhoton.At(i).P4()        
       delta_photon += (pf_mom.E() - pf_mom.Pz())
       temp_p = temp_p + pf_mom.Vect()

    delta_neutral = 0
    delta_neutral_noBarrel = 0
    pt_noBarrel = ROOT.TVector3()
    
    for i in range(branchEFlowNeutralHadron.GetEntries()):
       pf_mom = branchEFlowNeutralHadron.At(i).P4()
       delta_neutral += (pf_mom.E() - pf_mom.Pz())
       temp_p = temp_p+ pf_mom.Vect()
       if abs(pf_mom.Eta())>1.0:
           delta_neutral_noBarrel += (pf_mom.E() - pf_mom.Pz())
           pt_noBarrel = pt_noBarrel + pf_mom.Vect()
    
    delta = delta_track + delta_photon + delta_neutral
    
    #delta_noel = delta_track_noel + delta_photon + delta_neutral
    #delta_noel_noBarrel = delta_track_noel + delta_photon + delta_neutral_noBarrel
    delta_noBarrel = delta_track + delta_photon + delta_neutral_noBarrel
    
    y_JB   = delta/(2.0*10.0)
    ptmiss = temp_p.Perp()
    Q2_JB  = (ptmiss*ptmiss)/(1-y_JB)
    s     = 4*10.0*275.0
    if(y_JB>0):
        x_JB  = Q2_JB/(s*y_JB)
    else:
        x_JB = -999

    pt_all = temp_p
    ## Compute VA and VAP
    VP = 0
    VAP = 0
    for i in range(branchEFlowTrack.GetEntries()):
        track_mom = branchEFlowTrack.At(i).P4()   
        dot =  track_mom.Vect().Dot(pt_all.Unit())
        if(dot>0):
            VP = VP + dot
        elif(dot<0):
            VAP = VAP - dot
    for i in range(branchEFlowPhoton.GetEntries()):
        pf_mom = branchEFlowPhoton.At(i).P4()
        dot =  pf_mom.Vect().Dot(pt_all.Unit())
        if(dot>0):
            VP = VP + dot
        elif(dot<0):
            VAP = VAP - dot


    VP_nobarrel = VP
    VAP_nobarrel = VAP
    for i in range(branchEFlowNeutralHadron.GetEntries()):
        pf_mom = branchEFlowNeutralHadron.At(i).P4()
        dot =  pf_mom.Vect().Dot(pt_all.Unit())
        if(dot>0):
            VP = VP + dot
        elif(dot<0):
            VAP = VAP - dot

        if( abs(pf_mom .Eta())<1.0):
             if(dot>0):
                 VP_nobarrel = VP_nobarrel + dot
             elif(dot<0):
                 VAP_nobarrel = VAP_nobarrel - dot



    #print 'VAP', VAP , ' ' , VAP_nobarrel
    #print 'VP' , VP ,  ' ' , VP_nobarrel
    if(VP>0):  
         histo['Vratio'].Fill(VAP/VP)
    if(VP_nobarrel>0):
        histo['Vratio_noBarrelHCAL'].Fill(VAP_nobarrel/VP_nobarrel)
    
    #delta at generator level
    ptgen_all = ROOT.TVector3()
    gendelta = 0
    for i in range(branchParticle.GetEntries()):
        particle = branchParticle.At(i)
        gen_mom = particle.P4()
        status = particle.Status                       
        if(status!=1): continue

        ptgen_all = gen_mom.Vect()
        gendelta = gendelta + (gen_mom.E() - gen_mom.Pz())
    ## Compute VAP and VP at the generator level
    VAP_gen = 0
    VP_gen = 0
    for i in range(branchParticle.GetEntries()):
        particle = branchParticle.At(i)
        gen_mom = particle.P4()
        status = particle.Status
        if(status!=1): continue
        dot =  gen_mom.Vect().Dot(ptgen_all.Unit())  
        if(dot>0):
            VP_gen = VP + dot
        elif(dot<0):
            VAP_gen = VAP - dot

    if(VP_gen>0):
        histo['Vratio_truth'].Fill(VAP_gen/VP_gen) 
    #print 'Gen delta ' , gendelta , ' reco delta ', delta 
    
    histo['delta_gen'].Fill(gendelta)
    histo['delta_reco'].Fill(delta)
    histo['delta_reco_noel'].Fill(0)#delta_noel)
    histo['delta_reco_noBarrelHCAL'].Fill(delta_noBarrel)

    
    ##Fill purity histograms y-x and Q2-x
    genbin  = Ngen_y_x.FindBin(y, x)
    recobin = Ngen_y_x.FindBin(y_JB,x_JB)
    Ngen_y_x.Fill(y,x)
    if(genbin!=recobin):
        Nout_y_x.Fill(y ,x) ### In a given generator bin, how many left out.
        Nin_y_x.Fill(y_JB, x_JB)
   

    genbin  = Ngen_Q2_x.FindBin(x,Q2)
    recobin = Ngen_Q2_x.FindBin(x_JB,Q2_JB)
    Ngen_Q2_x.Fill(x,Q2)
    if(genbin!=recobin):
        Nout_Q2_x.Fill(x,Q2) ### In a given generator bin, how many left out.
        Nin_Q2_x.Fill(x_JB, Q2_JB)

    x_Matrix.Fill(x, x_JB)
    if(x>0):
        ResMatrix['x'].Fill(x, (x_JB-x)/x)
        profile['x'].Fill(x, (x_JB-x)/x)
    
    Q2_Matrix.Fill(ROOT.TMath.Sqrt(Q2), ROOT.TMath.Sqrt(Q2_JB))
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
                    profile['ept'].Fill(gen_electron.PT, res)
                    ResMatrix['ept'].Fill(gen_electron.PT, res)   
                                                                                                                                                                                        
  ##Jet response matrix
    if branchJet.GetEntries() > 0:
        # Take first jet
        jet = branchJet.At(0)

        deltaR = 999
        ipar_best = 999
        for ipar in range(branchGenJet.GetEntries()):
             particle = branchGenJet.At(ipar) #branchParticle.At(ipar)
             genJetMomentum = particle.P4()
             if(genJetMomentum.Px() == 0 and genJetMomentum.Py() == 0): continue

             if(genJetMomentum.DeltaR(jet.P4()) < deltaR):
                 deltaR = genJetMomentum.DeltaR(jet.P4());
                 bestGenJetMomentum = genJetMomentum;
                 ipar_best = ipar
        #print ipar
        #print deltaR
        if (deltaR>0.3): continue
        #if(abs(genjet.Eta())>3.0): continue    
        genjet = bestGenJetMomentum #branchGenJet.At(0)
        if(abs(genjet.Eta())>3.0): continue 
        # Print jet transverse momentum
        #JetMatrix.Fill(genjet.Pt()*ROOT.TMath.CosH(genjet.Eta()), jet.PT*ROOT.TMath.CosH(jet.Eta))
        JetMatrix.Fill(genjet.E(), jet.P4().E())
        
        res = (jet.PT-genjet.Pt())/genjet.Pt()
        profile['jetpt'].Fill(genjet.Pt(), res)
        ResMatrix['jetpt'].Fill(genjet.Pt(),res)
        profile['jetphi'].Fill(genjet.E(),   jet.Phi-genjet.Phi())
        ResMatrix['jetphi'].Fill(genjet.E(), jet.Phi-genjet.Phi())


        
#       print genjet.PT*ROOT.TMath.CosH(genjet.Eta), ' ',  genjet.P4().E()
        Jet_eta_e.Fill(genjet.Eta(), genjet.Pt()*ROOT.TMath.CosH(genjet.Eta()))
            
        if(genjet.Pt()<5.0):
           continue
        if abs(genjet.Eta())<1.0:
            etabin = 1
        elif genjet.Eta()>1.0 and genjet.Eta()<2.0:
            etabin = 2
        elif genjet.Eta()>2.0 and genjet.Eta()<3.0:
            etabin = 3
        else:
            etabin=4
        ResMatrix['jete_eta%i'%etabin].Fill(genjet.E(), (jet.P4().E()-genjet.E())/genjet.E())
        profile['jete_eta%i'%etabin].Fill(genjet.E(), (jet.P4().E()-genjet.E())/genjet.E())  

        ResMatrix['jete'].Fill(genjet.E(), (jet.P4().E()-genjet.E())/genjet.E())
        profile['jete'].Fill(genjet.E(), (jet.P4().E()-genjet.E())/genjet.E())

    ##Met response matrix
    if branchMet.GetEntries() > 0:
        met = branchMet.At(0)
        genmet = branchGenMet.At(0)
        if(genmet.MET==0): continue
        #print 'MET', met.MET
        #print 'pt all ' , ptmiss_all.Perp()
        #print 'pt no barrel HCAL ', (ptmiss_all-ptmiss_noBarrel).Perp()
        #print pts_noBarrel.Perp()
        #print pt_all.Phi() , ' ' , met.Phi
        #print (-1.0*ptmiss_all).Phi() , ' ' , met.Phi
        #print 'GenMET', genmet.MET
        if branchGenMet.GetEntries()>0:
            METMatrix.Fill(genmet.MET, met.MET)
            res = (pt_all.Perp()-genmet.MET)/genmet.MET
            profile['met'].Fill(genmet.MET,res)
            ResMatrix['met'].Fill(genmet.MET, res)
            profile['metphi'].Fill(genmet.MET,met.Phi-genmet.Phi)
            ResMatrix['metphi'].Fill(genmet.MET, met.Phi-genmet.Phi)

            ##no hcal barrel
            res = ((pt_all-pt_noBarrel).Perp()-genmet.MET)/genmet.MET
            profile['met_nobarrelHCAL'].Fill(genmet.MET,res)
            ResMatrix['met_nobarrelHCAL'].Fill(genmet.MET, res)

            phi = (-1.0*(pt_all-pt_noBarrel)).Phi()
            profile['metphi_nobarrelHCAL'].Fill(genmet.MET,   phi-genmet.Phi)
            ResMatrix['metphi_nobarrelHCAL'].Fill(genmet.MET, phi-genmet.Phi)
            
            ##Fill purity matrices
            ## MET vs x
            genbin = Ngen_met_x.FindBin(genmet.MET, x)
            recobin = Ngen_met_x.FindBin(met.MET,x_JB)
            Ngen_met_x.Fill(genmet.MET,x)
            if(genbin!=recobin):
                Nout_met_x.Fill(genmet.MET,x) ### In a given generator bin, how many left out.
                Nin_met_x.Fill(met.MET, x_JB)
            
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
    ## dphi  
    if leptonOK and branchJet.GetEntries() > 0 and branchGenJet.GetEntries()>0:
        #jet = branchJet.At(0)
        #genjet = branchGenJet.At(0)
        genlepton_phi = genlepton.Phi()
        lepton_phi =    lepton.Phi()
      
        jet_phi =    jet.P4().Phi()
        genjet_phi=  genjet.Phi()

        deltaR = genjet.DeltaR(jet.P4());
        if(deltaR>0.3):
            print deltaR, '  delta R'
            continue
        dphi_gen = ROOT.TVector2.Phi_mpi_pi(genlepton_phi-genjet_phi)
        dphi_reco = ROOT.TVector2.Phi_mpi_pi(lepton_phi-jet_phi)

        qT_reco = ROOT.TVector2(jet.P4().Px() + lepton.Px(), jet.P4().Py() + lepton.Py() ).Mod()
        #qT_reco = qT_reco.Mod()/ROOT.TMath.Sqrt(Q2)
      
        qT_truth = ROOT.TVector2(genjet.Px() + genlepton.Px(), genjet.Py() + genlepton.Py() ).Mod()
        #qT_truth = qT_truth.Mod()/ROOT.TMath.Sqrt(Q2)   

        ResMatrix['dphi'].Fill(genjet.E(), dphi_reco-dphi_gen)
        profile['dphi'].Fill(genjet.E(), dphi_reco-dphi_gen)
        distribution['dphi_gen'].Fill(genjet.E(), dphi_gen)
        distribution['dphi_reco'].Fill(genjet.E(), dphi_reco)

        #print 'qt/jetpt' , ROOT.TMath.ACos(qT_truth/genlepton.Pt())
        #print 'qT/genpt' , ROOT.TMath.ASin(qT_truth/genjet.Pt())
        #print 'sin dphi', ROOT.TMath.Sin(dphi_gen + 3.14159)#qT_truth/genjet.Pt())
        #print 'tan dphi ',   ROOT.TMath.Tan(dphi_gen + 3.14159)
        #print 'cos dphi ',   ROOT.TMath.Cos(dphi_gen + 3.14159) 
        #print 'dphi ' , dphi_gen
        #print 'dphi ' , dphi_gen+3.14159
        #print ' ' 
        ResMatrix['qtnormjet'].Fill(genjet.E(), qT_reco/jet.PT -qT_truth/genjet.Pt())
        profile['qtnormjet'].Fill(genjet.E(), qT_reco/jet.PT -qT_truth/genjet.Pt()) 
        distribution['qtnormjet_reco'].Fill(genjet.E(), qT_reco/jet.PT)
        distribution['qtnormjet_gen'].Fill(genjet.E(), qT_truth/genjet.Pt())

        ResMatrix['qt'].Fill(genjet.E(), qT_reco-qT_truth)
        profile['qt'].Fill(genjet.E(),  qT_reco-qT_truth)
        
###Colors
for key in profile:
    profile[key].SetLineColor(2)
    profile[key].SetLineWidth(2)
    profile[key].SetMarkerColor(2)


##Spectrum
#for key in neutral_E.keys():
#    neutral_E[key].SetLineColor(2)

#for key in neutral_E.keys():                                                                                                                                                                                     #     photon_E[key].SetLineColor(4)

#for i in range(1,5):
#    neutral_E['etabin%i'%i].Draw()
#    photon_E['etabin%i'%i].Draw("same")
#    track_E['etabin%i'%i].Draw("same")
#    ROOT.gPad.SetLogy(1)
#    c.SaveAs("plots/Spectrum_etabin%i_%s.png"%(i,inputFile))
#    c.SaveAs("plots/Spectrum_etabin%i_%s.pdf"%(i,inputFile))

ROOT.gPad.SetLogy(0)


c.Clear()
METMatrix.Draw("colz")
c.SaveAs("plots/Diagonal_MET_%s.png"%inputFile)
c.SaveAs("plots/pdf/Diagonal_MET_%s.pdf"%inputFile)

c.Clear()
JetMatrix.Draw("colz")
c.SaveAs("plots/Diagonal_Jet_%s.png"%inputFile) 
c.SaveAs("plots/pdf/Diagonal_Jet_%s.pdf"%inputFile)


c.Clear()
y_Matrix.Draw("colz")
c.SaveAs("plots/y_response_%s.png"%inputFile)
c.SaveAs("plots/pdf/y_response_%s.pdf"%inputFile)

c.Clear()
x_Matrix.Draw("colz")
c.SaveAs("plots/x_response_%s.png"%inputFile)
c.SaveAs("plots/pdf/x_response_%s.pdf"%inputFile)

c.Clear()

Q2_Matrix.Draw("colz")
#c.SetLogy(1)
#c.SetLogx(1)
c.SaveAs("plots/Q2_response_%s.png"%inputFile)
c.SaveAs("plots/pdf/Q2_response_%s.pdf"%inputFile)
c.SetLogy(0)
c.SetLogx(0)

#input("Press enter to continue...")
c.Clear()
ResMatrix['jetpt'].SetTitle("; jet p_{T}^{gen}; (jet p_{T}^{reco} - jet p_{T}^{gen})/jet p_{T}^{gen}")   
ResMatrix['jetpt'].Draw("colz")
profile['jetpt'].Draw("same")
c.SaveAs("plots/jetpt_%s_profile.png"%inputFile)    
c.SaveAs("plots/pdf/jetpt_%s_profile.pdf"%inputFile)

##Jet E:
for i in range(1,4):
    c.Clear()
    ResMatrix['jete_eta%i'%i].SetTitle("; jet E^{gen}; (jet E^{reco} - jet E^{gen})/jet E^{gen}")
    ResMatrix['jete_eta%i'%i].Draw("colz")
    profile['jete_eta%i'%i].Draw("same")
    #c.SaveAs("plots/profile_jetE_eta%i_%s.png"%(i,inputFile))
    #c.SaveAs("plots/profile_jetE_eta%i_%s.pdf"%(i,inputFile))



#input("Press enter to continue...")
c.Clear()
ResMatrix['met'].SetTitle("; generated MET; (MET^{reco} - MET^{gen})/MET^{gen}")
ResMatrix['met'].Draw("colz")
profile['met'].Draw("same")
c.SaveAs("plots/MET_profile%s.png"%inputFile)
c.SaveAs("plots/pdf/MET_profile%s.pdf"%inputFile)

c.Clear()
ResMatrix['metphi'].SetTitle("; generated MET; MET #phi^{reco} - MET #phi^{gen} [rad]")
ResMatrix['metphi'].Draw("colz")
profile['metphi'].Draw("same")
c.SaveAs("plots/METphi_%s_profile.png"%inputFile)
c.SaveAs("plots/pdf/METphi_%s_profile.pdf"%inputFile) 


##no barrel hcal
c.Clear()
ResMatrix['met_nobarrelHCAL'].SetTitle("; generated MET; (MET^{reco} - MET^{gen})/MET^{gen}")
ResMatrix['met_nobarrelHCAL'].Draw("colz")
profile['met_nobarrelHCAL'].Draw("same")
#c.SaveAs("plots/profile_METnobarrelHCAL_%s.png"%inputFile)
#c.SaveAs("plots/profile_METnobarrelHCAL_%s.pdf"%inputFile)

c.Clear()
ResMatrix['metphi_nobarrelHCAL'].SetTitle("; generated MET; MET #phi^{reco} - MET #phi^{gen} [rad]")
ResMatrix['metphi_nobarrelHCAL'].Draw("colz")
profile['metphi_nobarrelHCAL'].Draw("same")
#c.SaveAs("plots/profile_METphi_nobarrelHCAL_%s.png"%inputFile)
#c.SaveAs("plots/profile_METphi_nobarrelHCAL_%s.pdf"%inputFile)




c.Clear()
ResMatrix['jetphi'].SetTitle("; jet E^{gen} [GeV]; jet #phi^{reco} - jet #phi^{gen} [rad]")  
ResMatrix['jetphi'].Draw("colz")
profile['jetphi'].Draw("same")
c.SaveAs("plots/jetphi_%s_profile.png"%inputFile) 
c.SaveAs("plots/pdf/jetphi_%s_profile.pdf"%inputFile)  

c.Clear()
if(isCC):
    ResMatrix['dphi'].SetTitle("; generated MET; #Delta #phi_{reco} - #Delta #phi_{gen}")
elif(isNC):
    ResMatrix['dphi'].SetTitle("; jet E^{gen} [GeV]; #Delta #phi_{reco} - #Delta #phi_{gen}")
    
ResMatrix['dphi'].Draw("colz")
profile['dphi'].Draw("same")
c.SaveAs("plots/dphi_%s_profile.png"%inputFile)
c.SaveAs("plots/pdf/dphi_%s_profile.pdf"%inputFile) 


c.Clear()
ResMatrix['qtnormjet'].Draw("colz")
profile['qtnormjet'].Draw("same")
c.SaveAs("plots/profile_qtnormjet_%s.png"%inputFile)
c.SaveAs("plots/pdf/profile_qtnormjet_%s.pdf"%inputFile)

c.Clear()
ResMatrix['qt'].Draw("colz")
profile['qt'].Draw("same")
c.SaveAs("plots/qt_%s_profile.png"%inputFile)
c.SaveAs("plots/pdf/qt_%s_profile.pdf"%inputFile)


c.Clear()
ResMatrix['x'].Draw("colz")
profile['x'].Draw("same")
c.SaveAs("plots/x_%s_profile.png"%inputFile)
if isCC:
    for i in range(1,7):
        print 'bin ' , i , ' ' , profile['metphi'].GetBinError(i)
        print profile['metphi'].GetNbinsX()

##projections
#    for i in range(1,7):
#        c.Clear()
#        h = ResMatrix['metphi'].ProjectionY("h",i,i)
#        h.Draw()
#        h.Fit("gaus")
#        f = h.GetFunction("gaus")
#        if(f):
#            f.SetLineColor(2)
#            f.Draw("same")
#            h.Draw("same")
#        c.SaveAs("plots/projection_metphi_%s_%i.png"%(inputFile,i))
#        c.SaveAs("plots/projection_metphi_%s_%i.pdf"%(inputFile,i))
#        h_nobarrel =ResMatrix['metphi_nobarrelHCAL'].ProjectionY("h_nobarrel",i,i)
#        h_nobarrel.SetLineColor(4)
#        h_nobarrel.Draw("same")
#        c.SaveAs("plots/projection_metphi_comparison%s_%i.png"%(inputFile,i))  

#qt
ytitle = {}
ytitle['qt'] = " #Delta q_{T} (GeV)"
ytitle['qtnormjet'] = " q_{T}/p_{T}^{jet} resolution "
ytitle['jetphi'  ]  = " #phi^{jet} resolution (radians)"
ytitle['jete']  = " Relative energy resolution "
ytitle['dphi']  = " #Delta #phi resolution (radians)"
ytitle['met']   = ' Relative MET resolution '
ytitle['metphi']  = 'MET #phi resolution (radians)'
ytitle['jetpt']  = 'Relative p_{T} resolution' 
ytitle['x'] = 'Relative x resolution'

xtitle = {}
xtitle['jete'] = 'generated jet energy [GeV]'
xtitle['dphi'] = 'generated jet energy [GeV]'
xtitle['jetphi'] = 'generated jet energy [GeV]'
xtitle['qt'] = 'generated jet energy [GeV]'
xtitle['qtnormjet'] = 'generated jet energy [GeV]'
xtitle['met']  = 'Generated MET'
xtitle['metphi'] = 'Generated MET'
xtitle['jetpt']  = 'generated jet p_{T} [GeV]'
xtitle['x' ] = 'generated x'


out_tree = ROOT.TFile("resolutions_%s.root"%inputFile,"RECREATE")

for variable in ['qt','qtnormjet','jetphi','jete','jetpt','dphi','met','metphi','x']:
    g_sigma = ROOT.TGraph()
    g_rms  = ROOT.TGraph()
    g_multi = ROOT.TMultiGraph()
    for i in range(1,ResMatrix[variable].GetNbinsX()+1):
        print 'jet resolution projections '
        c.Clear()
        x= ResMatrix[variable].GetXaxis().GetBinCenter(i)
        h =ResMatrix[variable].ProjectionY("h",i,i)
        h.Draw()
        h.Fit("gaus")
        f = h.GetFunction("gaus")
        
        if not (f):
            continue
        print 'Fitted Gaussian' , f.GetParameter(2)
        print 'Standard deviation '            , h.GetStdDev()    
        f.SetLineColor(2)
        f.Draw("same")
        
        g_sigma.SetPoint(i-1, x, f.GetParameter(2))
        g_rms.SetPoint(i-1, x, h.GetRMS())

        DrawText (0.6 ,0.80 , "#sigma = %2.2f"%(f.GetParameter(2)), 2)
        DrawText (0.6 ,0.74 , "Stdv =%2.2f"%(h.GetRMS()),1)
        minimo = ResMatrix[variable].GetXaxis().GetBinLowEdge(i)
        maximo = minimo + ResMatrix[variable].GetXaxis().GetBinWidth(i)
        DrawText (0.6 ,0.69 , "%2.2f-%2.2f GeV"%(minimo,maximo),1) 
        h.Draw("same")
        c.SaveAs("plots/projections/%s_%s_%i_projection.png"%(variable,inputFile,i))
    c.Clear()
    g_rms.SetMarkerColor(2)
    g_rms.SetLineColor(2) 
    g_multi.Add(g_sigma)
    g_multi.Add(g_rms)
    g_multi.Draw("APL")
    g_multi.SetTitle(" ; %s; %s" %(xtitle[variable], ytitle[variable]))
    g_multi.SetMinimum(0)

        #latex = ROOT.TLatex()
        #latex.SetNDC()
        #latex.SetTextSize(0.06)
    DrawText (0.5 ,0.80 , "Fitted Gaussian")
    DrawText (0.5 ,0.74 , "Standard deviation",2)  
    c.SaveAs("plots/resolution_%s_%s.png"%(variable,inputFile))
    c.SaveAs("plots/pdf/resolution_%s_%s.pdf"%(variable,inputFile))
    g_sigma.Write("resolution_sigma_%s"%(variable))
    g_rms.Write("resolution_stdv_%s"%(variable))

out_tree.Write()
out_tree.Close()


    
#for i in range(1,7):
#    c.Clear()
#    h =ResMatrix['met'].ProjectionY("h",i,i)
#    h.Draw()
#    h.Fit("gaus")
#    f = h.GetFunction("gaus")
#    if(f):
#        f.SetLineColor(2)
#        f.Draw("same")
#        h.Draw("same")
#    c.SaveAs("plots/projection_met_%s_%i.png"%(inputFile,i))

#    h_nobarrel =ResMatrix['met_nobarrelHCAL'].ProjectionY("h_nobarrel",i,i)
#    h_nobarrel.SetLineColor(4)
#    h_nobarrel.Draw("same")
#    c.SaveAs("plots/projection_met_comparison%s_%i.png"%(inputFile,i))

    
##Reco vs truth distributions
    
#for i in range(distribution['dphi_reco'].GetNbinsX()):
#    c.Clear()
#    hgen = distribution['dphi_gen'].ProjectionY("h",i,i).Clone('hgen')
#    hgen.SetTitle("; | #phi_{lepton} -#phi_{jet} | [rad] ; 1/#sigma d\sigma/d#Delta#phi")
#    hreco  = distribution['dphi_reco'].ProjectionY("h",i,i).Clone('hreco')
#    hgen.SetLineColor(4)    
#    hgen.DrawNormalized()
#    #hgen.SetTitle("; | #phi_{lepton} -#phi_{jet} | [rad] ; 1/#sigma d\sigma/d#Delta#phi")
#    hreco.SetLineColor(2)
#    hreco.DrawNormalized("same")
#    hreco.SetLineColor(2)
#    c.SaveAs("plots/projection_dphirecogen_%s_%i.png"%(inputFile,i))

#for i in range(distribution['qtnormjet_reco'].GetNbinsX()):
#    c.Clear()
#    hgen = distribution['qtnormjet_gen'].ProjectionY("h",i,i).Clone('hgen')
#    hreco  = distribution['qtnormjet_reco'].ProjectionY("h",i,i).Clone('hreco')
#    hgen.SetLineColor(4)
#    hgen.DrawNormalized()
#    hreco.DrawNormalized("same")
#    hreco.SetLineColor(2)
#    c.SaveAs("plots/projection_qtnormjetrecogen_%s_%i.png"%(inputFile,i))

    
#c.Clear()
#Ngen_met_x.Draw("colz")
#c.SaveAs("plots/Ngen_met_x_%s.png"%(inputFile))



c.Clear()

histo['delta_gen'].Draw()
histo['delta_gen'].SetTitle("; #delta [GeV] ; entries")   
histo['delta_gen'].SetLineColor(2)

histo['delta_reco'].Draw("same")
histo['delta_reco'].SetLineColor(1)
#histo['delta_reco_noel'].Draw("same")
#histo['delta_reco_noel'].SetLineColor(8)

histo['delta_reco_noBarrelHCAL'].Draw("same")
histo['delta_reco_noBarrelHCAL'].SetLineColor(4)
#c.SaveAs("plots/delta_%s.png"%(inputFile))
#c.SaveAs("plots/pdf/delta_%s.pdf"%(inputFile))

ROOT.gPad.SetLogy(1)
c.SaveAs("plots/delta_log%s.png"%(inputFile))
c.SaveAs("plots/pdf/delta_log%s.pdf"%(inputFile))
ROOT.gPad.SetLogy(1)



histo['delta_gen'].Scale(1.0/histo['delta_gen'].Integral())
histo['delta_reco'].Scale(1.0/histo['delta_reco'].Integral())
#histo['delta_reco_noel'].Scale(1.0/histo['delta_reco_noel'].Integral())
histo['delta_reco_noBarrelHCAL'].Scale(1.0/histo['delta_reco_noBarrelHCAL'].Integral())      

histo['delta_gen'].GetCumulative().Draw()
histo['delta_reco'].GetCumulative().Draw("same")
#histo['delta_reco_noel'].GetCumulative().Draw("same")
histo['delta_reco_noBarrelHCAL'].GetCumulative().Draw("same")

c.SaveAs("plots/delta_cumulative_log%s.png"%(inputFile))
c.SaveAs("plots/pdf/delta_cumulative_log%s.pdf"%(inputFile))
ROOT.gPad.SetLogy(0) 



##Purity

ROOT.gStyle.SetPalette(55)
#delta
purity_num = Ngen_met_x.Clone("purity_num")
purity_den = Ngen_met_x.Clone("purity_den")
print 'Ngen ', Ngen_met_x.GetBinContent(11)
print 'Nout ', Nout_met_x.GetBinContent(11)
print 'Nin  ', Nin_met_x.GetBinContent(11)
purity_num.Add(Nout_met_x,-1)
purity_den.Add(Nout_met_x,-1)
purity_den.Add(Nin_met_x)
purity = purity_num.Clone("purity")
purity.Divide(purity_den)
purity.SetMaximum(1.0)
purity.SetMinimum(0.0)  
print 'purity ', purity.GetBinContent(11)

c.Clear()
purity.Draw("colz")
purity.SetTitle("Purity; reco MET; reco x_{JB}")    
c.SaveAs("plots/purity_met_x%s.png"%(inputFile))
c.SaveAs("plots/pdf/purity_met_x%s.pdf"%(inputFile))




ROOT.gStyle.SetPalette(112)

c.Clear()
ROOT.gPad.SetLogy(0)
ROOT.gPad.SetLogx(0)
ResMatrix['jete'].Draw("colz")
ResMatrix['jete'].SetTitle("; jet E^{gen}; (jet E^{reco} - jet E^{gen})/jet E^{gen}")
profile['jete'].Draw("same")
c.SaveAs("plots/profile_jetE_%s.png"%(inputFile))
c.SaveAs("plots/pdf/profile_jetE_%s.pdf"%(inputFile))

c.Clear()
Jet_eta_e.Draw("colz")
Jet_eta_e.SetTitle("; jet #eta^{gen}; jet E^{gen}")
c.SaveAs("plots/profile_jetE_eta_%s.png"%(inputFile))
c.SaveAs("plots/pdf/profile_jetE_eta_%s.pdf"%(inputFile))





c.Clear()
histo['Vratio'].Draw()
histo['Vratio'].SetTitle("; V_{AP}/V_{P} ; Entries")
histo['Vratio'].SetLineColor(1)
histo['Vratio_truth'].SetLineColor(2)
histo['Vratio_truth'].Draw("same")
histo['Vratio_noBarrelHCAL'].SetLineColor(4)
histo['Vratio_noBarrelHCAL'].Draw("same")



ROOT.gPad.SetLogy(1)
DrawText (0.5 ,0.80 , "Reco")
DrawText (0.5 ,0.75 , "Reco, no HCAL",4)
DrawText (0.5 ,0.70 , "Truth",2)  
c.SaveAs("plots/Vratio_%s.png"%(inputFile))

c.SaveAs("plots/pdf/Vratio_%s.pdf"%(inputFile))
ROOT.gPad.SetLogy(0)






if(isCC):
    c.Clear()

##Purity y vs x 
    purity_num = Ngen_y_x.Clone("purity_num")
    purity_den = Ngen_y_x.Clone("purity_den")
    purity_num.Add(Nout_y_x,-1)
    purity_den.Add(Nout_y_x,-1)
    purity_den.Add(Nin_y_x)
    purity = purity_num.Clone("purity")
    purity.Divide(purity_den)

    purity.SetMaximum(1.0)
    purity.SetMinimum(0.0)
    purity.Draw("colz")
    purity.SetTitle("Purity; reco x_{JB}; reco y_{JB}")
    c.SaveAs("plots/purity_y_x_%s.png"%(inputFile))
    c.SaveAs("plots/pdf/purity_y_x_%s.pdf"%(inputFile))

    c.Clear()
##Purity Q2 vs x
    purity_num = Ngen_Q2_x.Clone("purity_num")
    purity_den = Ngen_Q2_x.Clone("purity_den")
    purity_num.Add(Nout_Q2_x,-1)
    purity_den.Add(Nout_Q2_x,-1)
    purity_den.Add(Nin_Q2_x)
    purity = purity_num.Clone("purity")

    purity.Divide(purity_den)
    purity.SetMaximum(1.0)
    purity.SetMinimum(0.0)
    ROOT.gPad.SetLogy(1)
    ROOT.gPad.SetLogx(1)
    purity.Draw("colz")
    purity.SetTitle("Purity; reco x_{JB}; reco Q^{2}_{JB} [GeV^{2}]")
    c.SaveAs("plots/purity_Q2_x_%s.png"%(inputFile))
    c.SaveAs("plots/pdf/purity_Q2_x_%s.pdf"%(inputFile))

