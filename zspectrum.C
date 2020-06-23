/*
This macro shows how to access the particle-level reference for reconstructed objects.
It is also shown how to loop over the jet constituents.

root -l examples/Example3.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif



//------------------------------------------------------------------------------

struct TestPlots
{
  TH1 *fElectronDeltaPT;

  TH1* h_zjet_truth;
  TH1* h_zjet_reco;
  TH1* h_zjet_res;

  TH1* h_jete_truth;
  TH1* h_jete_reco;
  TH1* h_jete_res;

  TH1* Ngen_z;
  TH1* Nin_z;
  TH1* Nout_z;
  TH1* purity_z;
};

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, TestPlots *plots)
{
  TLegend *legend;
  TPaveText *comment;

  plots->fElectronDeltaPT = result->AddHist1D(   "electron_delta_pt", "(p_{T}^{particle} - p_{T}^{electron})/p_{T}^{particle}",   "(p_{T}^{particle} - p_{T}^{electron})/p_{T}^{particle}", "number of electrons",
    100, -0.1, 0.1);
  
  plots->h_zjet_truth = result->AddHist1D("h_zjet_truth","", "","", 200,0,2.0);
  plots->h_zjet_reco = result->AddHist1D("h_zjet_reco","", "","", 200,0,2.0);
  plots->h_zjet_res = result->AddHist1D("h_zjet_res","", "","", 100,-1.0,1.0);

  plots->h_jete_truth = result->AddHist1D("h_jete_truth","", "","", 200,0,100);
  plots->h_jete_reco = result->AddHist1D("h_jete_reco","", "","", 200,0,100);
  plots->h_jete_res = result->AddHist1D("h_jete_res","", "","", 100,-1.0,1.0);



  plots->Ngen_z = result->AddHist1D("Ngen_z","", "","", 10,0.0,1.0);
  plots->Nin_z = result->AddHist1D("Nin_z","", "","", 10,0.0,1.0);
  plots->Nout_z = result->AddHist1D("Nout_z","", "","", 10,0.0,1.0);
  plots->purity_z =  result->AddHist1D("purity_z","", "","", 10,0.0,1.0);   

}

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");

  Long64_t allEntries = treeReader->GetEntries();
  
  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  Electron *electron;

  Track *track;
  Tower *tower;

  Jet *jet;
  Jet *genjet;
  Jet *matchedjet;
  TObject *object;

  TLorentzVector momentum;

  Float_t Eem, Ehad;
  Bool_t skip;

  Long64_t entry;

  Int_t i, j, pdgCode;

  // Loop over all events

  TF1 *f1 = new TF1("f1","1+0.05*sin(x)",-TMath::Pi(),TMath::Pi());
  int njets = 0;

  //GenParticle* pProton, pleptonIn, pleptonOut, pPhoton;

  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // Loop over all electrons in event
    int index_e = 0;
    double maxpt = 0;
    if(branchElectron->GetEntriesFast()==0) continue;
    //get max pt electron
    for(i = 0; i < branchElectron->GetEntriesFast(); ++i)
    {
      electron = (Electron*) branchElectron->At(i);
      if(electron->PT<5) continue;
      if(electron->PT>maxpt){
	index_e = i;
	maxpt = electron->PT;
      }
      
      particle = (GenParticle*) electron->Particle.GetObject();

      plots->fElectronDeltaPT->Fill((particle->PT - electron->PT)/particle->PT);
    }

    electron = (Electron*) branchElectron->At(index_e);
    auto genelectron = (GenParticle*) electron->Particle.GetObject();

    if(!genelectron){
      std::cout << "could not find gen electron " << std::endl;
      continue;
    }
    //std::cout << " reco electron pT " << electron->PT << std::endl;
    //std::cout << " gen electron pT " << genelectron->PT << std::endl;
    
    auto pProton      = ((GenParticle*)branchParticle->At(0))->P4();
    auto pleptonIn    = ((GenParticle*)branchParticle->At(3))->P4();
    auto pleptonOut   = genelectron->P4(); //((GenParticle*)branchParticle->At(5))->P4();
    auto pPhoton       = pleptonIn - pleptonOut;
    auto pPhoton_reco  = pleptonIn - electron->P4(); 
    //std::cout << " plepton out " << ((GenParticle*)branchParticle->At(5))->PT << std::endl;
      
    double Q2 = -pPhoton.M2();
    double W2 = (pProton + pPhoton).M2();
    double x = Q2 / (2. * pProton.Dot(pPhoton));
    double y = (pProton.Dot(pPhoton)) / (pProton.Dot(pleptonIn));

    double Q2_reco = -pPhoton_reco.M2();
    //std::cout << " Q2 " << Q2 << " Q2_reco " << Q2_reco << std::endl;
    
    //Vec4 p (boosted_proton.px(), boosted_proton.py(), boosted_proton.pz(), boosted_proton.e());
    //Vec4 q (boosted_gamma.px(), boosted_gamma.py(), boosted_gamma.pz(), boosted_gamma.e());        
    //std::cout <<"Q2" << Q2 << " x " << x << " y" << y << std::endl;
     
    //event selection
    if(Q2<100) continue;
    if(y<0.1) continue;
    if(y>0.85) continue;
    if(electron->PT<10) continue;
    // Loop over all jets in event
    for(int i = 0; i < branchJet->GetEntriesFast(); ++i)
    {
      jet = (Jet*) branchJet->At(i);
      //if(jet->PT<10.0) continue;
      
      float deltaR = 999;
      
      auto jetMomentum = jet->P4();

      deltaR = 999;
      int matched_index = -999;
      for(j = 0; j < branchGenJet->GetEntriesFast(); j++)
	{
	  genjet = (Jet*) branchGenJet->At(j);
	  auto genJetMomentum = genjet->P4();
 	  // this is simply to avoid warnings from initial state particle
	  // having infite rapidity ...
	  if(genJetMomentum.Px() == 0 && genJetMomentum.Py() == 0) continue;

	  // take the closest parton candidate
	  if(genJetMomentum.DeltaR(jetMomentum) < deltaR)
	    {
	      deltaR = genJetMomentum.DeltaR(jetMomentum);
	      auto matchedJetMomentum = genJetMomentum;
	      matched_index = j;
	    }
	}
          
      if(deltaR>0.3) continue;
      //cout << deltaR << endl;
      genjet = (Jet*) branchGenJet->At(matched_index);
      double z_jet_truth = pProton.Dot(genjet->P4())/pProton.Dot(pPhoton);

      double z_jet_reco = pProton.Dot(jet->P4())/pProton.Dot(pPhoton_reco);

      //std::cout << " z_jet_truth " << z_jet_truth <<  " z_jet_reco:" << z_jet_reco<< std::endl;
      plots->h_zjet_truth->Fill(z_jet_truth);
      plots->h_zjet_reco->Fill(z_jet_reco);
      plots->h_zjet_res->Fill((z_jet_reco-z_jet_truth)/z_jet_truth);

      //jet energy as well
      //std::cout << "energy " << genjet->P4().E() << std::endl;
      //std::cout << "energy " << jet->P4().E() << std::endl;      
      plots->h_jete_truth->Fill(genjet->P4().E());
      plots->h_jete_reco->Fill(jet->P4().E());
      plots->h_jete_res->Fill( (jet->P4().E() - genjet->P4().E())/genjet->P4().E());
      
      int genbin  = plots->Ngen_z->FindBin(z_jet_truth);
      int recobin = plots->Ngen_z->FindBin(z_jet_reco);
      plots->Ngen_z->Fill(z_jet_truth);
      if(genbin!=recobin){
          plots->Nout_z->Fill(z_jet_truth); //
          plots->Nin_z->Fill(z_jet_reco);
      }  

      // zjet = P.P_jet/P.q
      
      // genjet = bestGenJetMomentum #branchGenJet.At(0)

    }//loop over jets 
  }//loop over entries
  std::cout << plots->Ngen_z->GetBinContent(7) << std::endl;
  std::cout << plots->Nout_z->GetBinContent(7) << std::endl;
  std::cout << plots->Nin_z->GetBinContent(7) << std::endl;
  //compute purity
  //plots->purity_jetpt_z->SetMaximum(1.0);
  //plots->purity_jetpt_z->SetMinimum(0.0); 

    
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  result->Print("png");
}

//------------------------------------------------------------------------------

void zspectrum(const char *inputFile)
{
  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  TestPlots *plots = new TestPlots;

  BookHistograms(result, plots);

  AnalyseEvents(treeReader, plots);

  PrintHistograms(result, plots);

  result->Write("zspectrum_results.root");

  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}

//------------------------------------------------------------------------------
