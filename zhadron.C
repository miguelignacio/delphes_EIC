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
    
  TH2 *res_z_10_20;
  TH2 *res_z_20_40;
  TH2 *res_z_40_60;
  TH2 *res_z_60_80;


  TProfile* profile_z_10_20;
  TProfile* profile_z_20_40;
  TProfile* profile_z_40_60;
  TProfile* profile_z_60_80;
  
  TProfile* profile_p_10_20;
  TProfile* profile_p_20_40;
  TProfile* profile_p_40_60;
  TProfile* profile_p_60_80;
};

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, TestPlots *plots)
{
  TLegend *legend;
  TPaveText *comment;

  plots->res_z_10_20 = result->AddHist2D("res_z_10_20", "Collins Angle resolution, 10 <E_{jet}< 20 GeV", "generated z", "(z^{reco}-z^{gen})/z^{gen}", 10, 0.0, 1.0,  100,-1,1);
  plots->res_z_20_40 = result->AddHist2D("res_z_20_40", "Collins Angle resolution, 20 <E_{jet}< 40 GeV", "generated z", "(z^{reco}-z^{gen})/z^{gen}", 10,0.0,1.0, 100,-1,1);
  plots->res_z_40_60 = result->AddHist2D("res_z_40_60", "Collins Angle resolution, 40 <E_{jet}< 60 GeV", "generated z", "(z^{reco}-z^{gen})/z^{gen}", 10,0.0,1.0, 100,-1,1);
  plots->res_z_60_80 = result->AddHist2D("res_z_60_80", "Collins Angle resolution, 60 <E_{jet}< 80 GeV", "generated z", "(z^{reco}-z^{gen})/z^{gen}", 10,0.0,1.0, 100,-1,1);      
  plots->profile_z_10_20 = result->AddProfile("profile_z_10_20", "z resolution, 10 <E_{jet}< 20 GeV", "generated z", "(z^{reco}-z^{gen})/z^{gen}",
						   10,0.0,1.0, -1.0,1.0);

  plots->profile_z_20_40 = result->AddProfile("profile_z_20_40", "z resolution, 20 <E_{jet}< 40 GeV", "generated z", "(z^{reco}-z^{gen})/z^{gen}",
					      10,0.0,1.0, -1.0,1.0);
  plots->profile_z_40_60 = result->AddProfile("profile_z_40_60", "z resolution, 40 <E_{jet}< 60 GeV", "generated z", "(z^{reco}-z^{gen})/z^{gen}",
					      10,0.0,1.0, -1.0,1.0);

  plots->profile_z_60_80 = result->AddProfile("profile_z_60_80", "z resolution, 60 <E_{jet}< 80 GeV", "generated z", "(z^{reco}-z^{gen})/z^{gen}",                                                                                                               10,0.0,1.0, -1.0,1.0);
  
  plots->profile_p_10_20 = result->AddProfile("profile_p_10_20", "Momentum resolution 10 <E_{jet}< 20 GeV", "generated p", "(p^{reco}-p^{gen})/p^{gen}",
					     50,0.0,50.0, -1.0,1.0);
  plots->profile_p_20_40 = result->AddProfile("profile_p_20_40", "Momentum resolution 20 <E_{jet}< 40 GeV", "generated p", "(p^{reco}-p^{gen})/p^{gen}",
					      50,0.0,50.0, -1.0,1.0);
  plots->profile_p_40_60 = result->AddProfile("profile_p_40_60", "Momentum resolution 40 <E_{jet}< 60 GeV", "generated p", "(p^{reco}-p^{gen})/p^{gen}",
					      50,0.0,50.0, -1.0,1.0);
  plots->profile_p_60_80 = result->AddProfile("profile_p_60_80", "Momentum resolution 60 <E_{jet}< 80 GeV", "generated p", "(p^{reco}-p^{gen})/p^{gen}",
					      50,0.0,50.0, -1.0,1.0);
  //  TProfile *AddProfile(const char *name, const char *title,  const char *xlabel, const char *ylabel,    Int_t nxbins, Axis_t xmin, Axis_t xmax,    Int_t logx = 0, Int_t logy = 0);
  
  
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
  for(entry = 0; entry < allEntries; ++entry)
  // for(entry = 0; entry < 1000; ++entry)   
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // Loop over all electrons in event
    for(i = 0; i < branchElectron->GetEntriesFast(); ++i)
    {
      electron = (Electron*) branchElectron->At(i);
      particle = (GenParticle*) electron->Particle.GetObject();

    }

  
  
    // Loop over all jets in event
    for(int i = 0; i < branchJet->GetEntriesFast(); ++i)
    {
      jet = (Jet*) branchJet->At(i);
      if(jet->PT<5.0) continue;
      
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
      // genjet = bestGenJetMomentum #branchGenJet.At(0)
      if( genjet->PT<15 and genjet->PT>10){
      	njets = njets+1;
      }
      
      momentum.SetPxPyPzE(0.0, 0.0, 0.0, 0.0);

      // cout<<"Looping over jet constituents. Jet pt: "<<jet->PT<<", eta: "<<jet->Eta<<", phi: "<<jet->Phi<<endl;

      // Loop over all jet's constituents
      for(j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
      {
        object = jet->Constituents.At(j);
        // Check if the constituent is accessible
        if(object == 0) continue;
        if(object->IsA() == Track::Class())
        {
          track = (Track*) object;
	  if(track->PT<0.100) continue;  // if track has less than 200 MeV, then do not consider it
          //cout << "    Track pt: " << track->PT << ", eta: " << track->Eta << ", phi: " << track->Phi << endl;
          momentum += track->P4();
	  if(abs(track->PID)!=211) continue; //if not a charged pion, continue
	  //cout << " track charge " << track->Charge << " PID " << track->PID<<endl;
	  double pxh, pyh, pzh, cross;
         
	  pxh = track->P4().Px();
	  pyh = track->P4().Py();
	  pzh = track->P4().Pz();
          double ptrack = TMath::Sqrt(pxh*pxh + pyh*pyh + pzh*pzh);
	  
	  double pxj, pyj, pzj;
          pxj = jet->P4().Px();
	  pyj = jet->P4().Py();
	  pzj = jet->P4().Pz();
          double p = TMath::Sqrt(pxj*pxj + pyj*pyj + pzj*pzj);

	  TVector3 crossproduct = jet->P4().Vect().Unit().Cross(track->P4().Vect().Unit());
	  double sin = crossproduct.Mag()/(jet->P4().Vect().Mag());
	    
	  double z = jet->P4().Vect().Dot( track->P4().Vect() )/(jet->P4().P()*jet->P4().P());
          double r = TMath::Sqrt( pow(jet->P4().Phi() - track->P4().Phi(),2.0) + pow(jet->P4().Eta() - track->P4().Eta(),2.0));
	  TVector3 zaxis(0,0,1);
          TVector3 N = zaxis.Cross(jet->P4().Vect());
	  TVector3 S = N.Cross(jet->P4().Vect());
	  N = N.Unit();
	  S = S.Unit();
	  TVector3 jt  = track->P4().Vect().Dot(N)*N + track->P4().Vect().Dot(S)*S;
          double angle = track->P4().Vect().Dot(N)/track->P4().Vect().Dot(S);	  
          double phi_h = TMath::ATan(angle); 
	  
          //track truth
	  auto gentrack = (GenParticle*) track->Particle.GetObject();
          double p_truth = gentrack->P4().P();
	  double p_reco  = track->P4().P();
	  
          TVector3 gen_crossproduct = genjet->P4().Vect().Cross(gentrack->P4().Vect());
	  double genz = genjet->P4().Vect().Dot( gentrack->P4().Vect() )/(genjet->P4().P()*genjet->P4().P());
	  //if( !(genz>0.4 and genz<0.6)) continue;
	  double genr = TMath::Sqrt( pow(genjet->P4().Phi() - gentrack->P4().Phi(),2.0) + pow(genjet->P4().Eta() - gentrack->P4().Eta(),2.0));      
	  TVector3 genN = zaxis.Cross(genjet->P4().Vect());
	  TVector3 genS = genN.Cross(genjet->P4().Vect());
	  genN = genN.Unit();
	  genS = genS.Unit();
	  TVector3 genjt  = gentrack->P4().Vect().Dot(genN)*genN + gentrack->P4().Vect().Dot(genS)*genS;
	  double genphi_h = TMath::ATan(gentrack->P4().Vect().Dot(genN)/gentrack->P4().Vect().Dot(genS)); 

	  
     
	  double dphi = TVector2::Phi_mpi_pi(genphi_h-phi_h);
	  if(dphi> TMath::Pi()/2.0){
	    dphi = dphi- TMath::Pi();
	  }
	  else if(dphi<-TMath::Pi()/2.0){
	    dphi = dphi+ TMath::Pi();
	  }
	 
       
	  

	  if(genjet->P4().E()>10 and genjet->P4().E()<20){
	    plots->res_z_10_20->Fill(genz,(z-genz)/genz);
	    plots->profile_z_10_20->Fill(genz,(z-genz)/genz);
	    plots->profile_p_10_20->Fill(p_truth, (p_reco-p_truth)/p_truth);
	  }
	  else if(genjet->P4().E()>20 and genjet->P4().E()<40){
	    plots->res_z_20_40->Fill(genz, (z-genz)/genz);
	    plots->profile_z_20_40->Fill(genz,(z-genz)/genz);      
	    plots->profile_p_20_40->Fill(p_truth, (p_reco-p_truth)/p_truth);   
	  }
	  else if(genjet->P4().E()>40 and genjet->P4().E()<60){
	    plots->res_z_40_60->Fill(genz, (z-genz)/genz);
	    plots->profile_z_40_60->Fill(genz,(z-genz)/genz);      
	    plots->profile_p_40_60->Fill(p_truth, (p_reco-p_truth)/p_truth);   
	  }
	  else if(genjet->P4().E()>60 and genjet->P4().E()<80){
	    plots->res_z_60_80->Fill(genz, (z-genz)/genz);
	    plots->profile_z_60_80->Fill(genz,(z-genz)/genz);      
	    plots->profile_p_60_80->Fill(p_truth, (p_reco-p_truth)/p_truth);   
	  }  
	                

        }
        else if(object->IsA() == Tower::Class())
        {
          tower = (Tower*) object;
          //cout << "    Tower pt: " << tower->ET << ", eta: " << tower->Eta << ", phi: " << tower->Phi << endl;
          momentum += tower->P4();
        }
      }//end loop over constituents

      

    }//loop over jets
    
  }//loop over entries


}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  result->Print("png");
}

//------------------------------------------------------------------------------

void zhadron(const char *inputFile)
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

  result->Write("zh_results.root");

  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}

//------------------------------------------------------------------------------
