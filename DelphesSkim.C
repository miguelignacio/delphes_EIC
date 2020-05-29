#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TCut.h"
#include "TString.h"
#include <iostream>

void DelphesSkim(TString input_filename, TString output_filename, TString skimCuts = "")
{
  auto input_file = TFile::Open(input_filename);

   if (input_file == nullptr || input_file->IsZombie() || input_file->TestBit(TFile::kRecovered)) {
    std::cout << "The input ROOT file is corrupted. Exiting..." << std::endl;
    return;
  }

  auto input_tree = static_cast<TTree*>(input_file->Get("Delphes"));

  std::cout << "Input Delphes TTree contains " << input_tree->GetEntries() << " events..." << std::endl;

  auto output_file = TFile::Open(output_filename,"RECREATE");
  output_file->cd();
  auto output_tree = static_cast<TTree*>(input_tree->CopyTree(skimCuts));
  output_tree->Write();

  std::cout << "... skimmed output Delphes TTree contains " << output_tree->GetEntries() << " events." << std::endl;

  output_file->Close();
}
