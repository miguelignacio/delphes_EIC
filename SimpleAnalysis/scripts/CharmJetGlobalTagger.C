#include "TROOT.h"
#include "TChain.h"
#include "TFile.h"
#include "TEfficiency.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TCut.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TMath.h"
#include "TLine.h"
#include "TLatex.h"
#include "TRatioPlot.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"

#include <glob.h>
#include <iostream>
#include <iomanip>
#include <vector>

#include "PlotFunctions.h"


void CharmJetGlobalTagger(TString dir, TString input, TString filePattern = "*/out.root")
{
  // Global options
  gStyle->SetOptStat(0);

  // Create the TCanvas
  TCanvas *pad = new TCanvas("pad",
                             "",
                             800,
                             600);
  TLegend *legend    = nullptr;
  TH1F    *htemplate = nullptr;

  auto default_data = new TChain("tree");
  default_data->SetTitle(input.Data());
  auto files = fileVector(Form("%s/%s/%s", dir.Data(), input.Data(), filePattern.Data()));

  for (auto file : files)
  {
    default_data->Add(file.c_str());
  }

  // Create the signal and background trees

  auto signal_train = default_data->CopyTree("jet_flavor==4 && jet_n>0", "", TMath::Floor(default_data->GetEntries() / 1.0));
  std::cout << "Signal Tree (Training): " << signal_train->GetEntries() << std::endl;

  auto background_train = default_data->CopyTree("(jet_flavor<4||jet_flavor==21) && jet_n>0", "", TMath::Floor(default_data->GetEntries() / 1.0));
  std::cout << "Background Tree (Training): " << background_train->GetEntries() << std::endl;

  // auto u_background_train = default_data->CopyTree("(jet_flavor==2)", "",
  // TMath::Floor(default_data->GetEntries()/1.0));
  // std::cout << "u Background Tree (Training): " <<
  // u_background_train->GetEntries() << std::endl;

  // auto d_background_train = default_data->CopyTree("(jet_flavor==1)", "",
  // TMath::Floor(default_data->GetEntries()/1.0));
  // std::cout << "d Background Tree (Training): " <<
  // d_background_train->GetEntries() << std::endl;

  // auto g_background_train = default_data->CopyTree("(jet_flavor==21)", "",
  // TMath::Floor(default_data->GetEntries()/1.0));
  // std::cout << "g Background Tree (Training): " <<
  // g_background_train->GetEntries() << std::endl;

  // Create the TMVA tools
  TMVA::Tools::Instance();

  auto outputFile = TFile::Open("CharmJetGlobalTagger_Results.root", "RECREATE");

  TMVA::Factory factory("TMVAClassification",
                        outputFile,
                        "!V:ROC:!Correlations:!Silent:Color:DrawProgressBar:AnalysisType=Classification");

  //
  // Kaon Tagger
  //

  TMVA::DataLoader loader("dataset");

  // loader.AddVariable("jet_pt",  "Jet p_{T}",          "GeV", 'F',
  // 0.0,
  // 1000.0);
  // loader.AddVariable("jet_eta", "Jet Pseudorapidity", "",    'F',
  // -5.0, 5.0);

  loader.AddSpectator("jet_pt");
  loader.AddSpectator("jet_eta");
  loader.AddSpectator("jet_flavor");
  loader.AddSpectator("met_et");
  loader.AddVariable("jet_mlp_ip3dtagger");
  loader.AddVariable("jet_mlp_ktagger");
  loader.AddVariable("jet_mlp_eltagger");
  loader.AddVariable("jet_mlp_mutagger");
  loader.AddSignalTree(signal_train, 1.0);
  loader.AddBackgroundTree(background_train, 1.0);


  //  loader.AddVariable("jet_sip3dtag", "sIP3D Jet-Level Tag", "", 'B', -10,
  // 10);
  // loader.AddVariable("jet_charge");

  // loader.AddVariable("jet_ehadoveremratio");


  // loader.AddTree( signal_train, "strange_jets" );
  // loader.AddTree( u_background_train, "up jets" );
  // loader.AddTree( d_background_train, "down jets" );
  // loader.AddTree( g_background_train, "gluon jets" );

  loader.PrepareTrainingAndTestTree(TCut("jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && jet_flavor==4 && met_et > 10"),
				    TCut("jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && (jet_flavor<4||jet_flavor==21) && met_et > 10"),
				    "nTrain_Signal=50000:nTrain_Background=500000:nTest_Signal=50000:nTest_Background=500000:SplitMode=Random:NormMode=NumEvents:!V");
  
  // Declare the classification method(s)
  // factory.BookMethod(&loader,TMVA::Types::kBDT, "BDT",
  //
  //
  //  "!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20"
  // );
  factory.BookMethod(&loader,    TMVA::Types::kMLP, "CharmGlobalTagger",
                     "!H:!V:NeuronType=ReLU:VarTransform=Norm:NCycles=1000:HiddenLayers=N+8:TestRate=5:!UseRegulator");

  // Train
  factory.TrainAllMethods();

  // Test
  factory.TestAllMethods();
  factory.EvaluateAllMethods();

  // Plot a ROC Curve
  pad->cd();
  pad = factory.GetROCCurve(&loader);
  pad->Draw();

  pad->SaveAs("CharmJetGlobalTagger_ROC.pdf");
}
