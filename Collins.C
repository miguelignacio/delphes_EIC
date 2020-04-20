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
  TH1 *fElectronDeltaEta;
  TH1 *fJetDeltaPT;
  TH1 *hjt;
  TH1 *hphi;
  TH1 *hz;
  TH1 *hr;
};

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, TestPlots *plots)
{
  TLegend *legend;
  TPaveText *comment;

  plots->fElectronDeltaPT = result->AddHist1D(
    "electron_delta_pt", "(p_{T}^{particle} - p_{T}^{electron})/p_{T}^{particle}",
    "(p_{T}^{particle} - p_{T}^{electron})/p_{T}^{particle}", "number of electrons",
    100, -0.1, 0.1);

  plots->fElectronDeltaEta = result->AddHist1D(
    "electron_delta_eta", "(#eta^{particle} - #eta^{electron})/#eta^{particle}",
    "(#eta^{particle} - #eta^{electron})/#eta^{particle}", "number of electrons",
    100, -0.1, 0.1);

  plots->fJetDeltaPT = result->AddHist1D(
    "jet_delta_pt", "(p_{T}^{jet} - p_{T}^{constituents})/p_{T}^{jet}",
    "(p_{T}^{jet} - p_{T}^{constituents})/p_{T}^{jet}", "number of jets",
    100, -1.0e-1, 1.0e-1);


  plots->hjt  = result->AddHist1D("jt distribution of hadrons in jet", "", "jt " , "entries",50,0,3.0);
  plots->hz   = result->AddHist1D("z distribution of hadrons in jet", ""," z " , "entries",50,0,1.0);
  plots->hphi = result->AddHist1D("phi distribution of hadrons in jet", "", "phi " , "entries",50,-1.0,1.0);
  plots->hr   = result->AddHist1D("r distribution of hadrons in het", "" , "z", "entries", 50, 0,1.0);
  
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

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  Electron *electron;

  Track *track;
  Tower *tower;

  Jet *jet;
  TObject *object;

  TLorentzVector momentum;

  Float_t Eem, Ehad;
  Bool_t skip;

  Long64_t entry;

  Int_t i, j, pdgCode;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // Loop over all electrons in event
    for(i = 0; i < branchElectron->GetEntriesFast(); ++i)
    {
      electron = (Electron*) branchElectron->At(i);
      particle = (GenParticle*) electron->Particle.GetObject();

      plots->fElectronDeltaPT->Fill((particle->PT - electron->PT)/particle->PT);
      plots->fElectronDeltaEta->Fill((particle->Eta - electron->Eta)/particle->Eta);
    }

    // Loop over all photons in event

    // Loop over all jets in event
    for(i = 0; i < branchJet->GetEntriesFast(); ++i)
    {
      jet = (Jet*) branchJet->At(i);
      if(jet->PT<10.0) continue;
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
          //cout << "    Track pt: " << track->PT << ", eta: " << track->Eta << ", phi: " << track->Phi << endl;
          momentum += track->P4();
	  if(abs(track->PID)==211) continue;
	  double pxh, pyh, pzh, cross;
         
	  pxh = track->P4().Px();
	  pyh = track->P4().Py();
	  pzh = track->P4().Pz();

	  double pxj, pyj, pzj;
          pxj = jet->P4().Px();
	  pyj = jet->P4().Py();
	  pzj = jet->P4().Pz();

	  TVector3 crossproduct = jet->P4().Vect().Cross(track->P4().Vect());
	  double z = jet->P4().Vect().Dot( track->P4().Vect() )/(jet->P4().P()*jet->P4().P());
          double r = TMath::Sqrt( pow(jet->P4().Phi() - track->P4().Phi(),2.0) + pow(jet->P4().Eta() - track->P4().Eta(),2.0));
	  TVector3 zaxis(0,0,1);
          TVector3 N = zaxis.Cross(jet->P4().Vect());
	  TVector3 S = N.Cross(jet->P4().Vect());
	  N = N.Unit();
	  S = S.Unit();
	  TVector3 jt  = track->P4().Vect().Dot(N)*N + track->P4().Vect().Dot(S)*S;
	  //cout << " alternative to jt " << jt_2.Mag() << endl;
	  double phi_h = jt.Unit().Dot(S);
	  plots->hphi->Fill(phi_h);
	  plots->hz->Fill(z);
	  plots->hjt->Fill(jt.Mag());
	  plots->hr->Fill(r);
	   
	  //cout << " angle " << phi_h << std::endl;
	  //cout << " r " << r << " jt " << jt << " z " << z << " " << jet->P4().P() <<  " " << jet->PT << endl;
        }
        else if(object->IsA() == Tower::Class())
        {
          tower = (Tower*) object;
          //cout << "    Tower pt: " << tower->ET << ", eta: " << tower->Eta << ", phi: " << tower->Phi << endl;
          momentum += tower->P4();
        }
      }
      plots->fJetDeltaPT->Fill((jet->PT - momentum.Pt())/jet->PT);
    }
  }
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  result->Print("png");
}

//------------------------------------------------------------------------------

void Collins(const char *inputFile)
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

  result->Write("results.root");

  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}

//------------------------------------------------------------------------------
