/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, and plot simple quantities such as the jet pt and the di-electron invariant
mass.

root -l examples/Example1.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

//------------------------------------------------------------------------------

void jetconstituent(const char *inputFile)
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  Track* track;
  TObject* object;
  Tower* tower;   
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis

  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");


  TClonesArray *branchEvent = treeReader->UseBranch("Event");

  // Book histograms
  TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 100, 0.0, 100.0);
  TH1 *histMass = new TH1F("mass", "M_{inv}(e_{1}, e_{2})", 100, 40.0, 140.0);
  TH2 *constituent_e_eta = new TH2F("constituent_e_eta","jet constituent energy vs eta", 100, 0, 30, 100,-4.,4.);
  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    //Float_t weight = event->Weight;

    // If event contains at least 1 jet
    if(branchJet->GetEntries() > 0)
    {
      // Take first jet
      Jet *jet = (Jet*) branchJet->At(0);

      // Plot jet transverse momentum
      histJetPT->Fill(jet->PT);

      // Print jet transverse momentum
      
      if(jet->PT<10 or jet->PT>15) continue;
      if(jet->Eta>2.0 or jet->Eta<1.5) continue;

      cout << "Jet pt: "<<jet->PT << endl;
      cout << "Jet eta: " << jet->Eta << endl;
      
      //loop over its constituents
      for(int j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
	{
	  
	  object = jet->Constituents.At(j);
	  // Check if the constituent is accessible
	  if(object == 0) continue;
	  //std::cout << 'as'<<std::endl;
	  if(object->IsA() == Track::Class())
	    {
	      track = (Track*) object;
	      if(track->PT<0.100) continue;  // if track has less than 200 MeV, then do not consider it
	      //cout << "    Track pt: " << track->PT << ", eta: " << track->Eta << ", phi: " << track->Phi << endl;
	      constituent_e_eta->Fill( track->P4().Vect().Mag(), track->Eta);
	    }

	  else if(object->IsA() == Tower::Class()){
	    tower = (Tower*) object;
	    cout << "Tower Energy " << tower->E << " Tower eta" << tower->Eta << endl;
	  }
	}
    }

   Jet *jet1, *jet2;

    // If event contains at least 2 electrons
    if(branchJet->GetEntries() > 1)
    {
      // Take first two electrons
      jet1 = (Jet *) branchJet->At(0);
      jet2 = (Jet *) branchJet->At(1);

      // Plot their invariant mass
      histMass->Fill(((jet1->P4()) + (jet2->P4())).M());
    }
  }

  // Show resulting histograms
  histJetPT->Draw();
  histMass->Draw();
  constituent_e_eta->Draw("colz");
}

