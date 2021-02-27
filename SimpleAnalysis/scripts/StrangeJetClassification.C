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


void StrangeJetClassification(TString dir, TString input, TString filePattern = "*/out.root")
{
  // Global options
  gStyle->SetOptStat(0);

  // Create the TCanvas
  TCanvas* pad = new TCanvas("pad","",800,600);
  TLegend * legend = nullptr;
  TH1F* htemplate = nullptr;

  auto default_data = new TChain("tree");
  default_data->SetTitle(input.Data());
  auto files = fileVector(Form("%s/%s/%s", dir.Data(), input.Data(), filePattern.Data()));

  for (auto file : files) 
    {
      default_data->Add(file.c_str());
    }

  // Create the signal and background trees

  auto signal_train = default_data->CopyTree("jet_flavor==3", "", TMath::Floor(default_data->GetEntries()/1.0));
  std::cout << "Signal Tree (Training): " << signal_train->GetEntries() << std::endl;

  auto background_train = default_data->CopyTree("(jet_flavor<3||jet_flavor==21)", "", TMath::Floor(default_data->GetEntries()/1.0));
  std::cout << "Background Tree (Training): " << background_train->GetEntries() << std::endl;

  // auto u_background_train = default_data->CopyTree("(jet_flavor==2)", "", TMath::Floor(default_data->GetEntries()/1.0));
  // std::cout << "u Background Tree (Training): " << u_background_train->GetEntries() << std::endl;

  // auto d_background_train = default_data->CopyTree("(jet_flavor==1)", "", TMath::Floor(default_data->GetEntries()/1.0));
  // std::cout << "d Background Tree (Training): " << d_background_train->GetEntries() << std::endl;

  // auto g_background_train = default_data->CopyTree("(jet_flavor==21)", "", TMath::Floor(default_data->GetEntries()/1.0));
  // std::cout << "g Background Tree (Training): " << g_background_train->GetEntries() << std::endl;

  // Create the TMVA tools
  TMVA::Tools::Instance();

  auto outputFile = TFile::Open("StrangeJetClassification_Results.root", "RECREATE");

  TMVA::Factory factory("TMVAClassification", outputFile,
			"!V:ROC:!Correlations:!Silent:Color:DrawProgressBar:AnalysisType=Classification" );


  TMVA::DataLoader loader("dataset");

  loader.AddVariable("jet_pt");
  loader.AddVariable("jet_eta");
  loader.AddVariable("jet_nconstituents");
  loader.AddVariable("jet_Ks_leading_zhadron");
  loader.AddVariable("jet_K_leading_zhadron");
  loader.AddVariable("jet_charge");
  loader.AddVariable("jet_Ks_sumpt");
  loader.AddVariable("jet_K_sumpt");
  loader.AddVariable("jet_ehadoveremratio");
  loader.AddSpectator("jet_flavor");

  loader.AddSignalTree( signal_train, 1.0 );
  loader.AddBackgroundTree( background_train, 1.0 );
  // loader.AddTree( signal_train, "strange_jets" );
  // loader.AddTree( u_background_train, "up jets" );
  // loader.AddTree( d_background_train, "down jets" );
  // loader.AddTree( g_background_train, "gluon jets" );

  loader.PrepareTrainingAndTestTree(TCut("jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && jet_flavor==3"), 
				    TCut("jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && (jet_flavor<3||jet_flavor==21)"),
   				    "nTrain_Signal=50000:nTrain_Background=100000:nTest_Signal=50000:nTest_Background=100000:SplitMode=Random:NormMode=NumEvents:!V" );

  // Declare the classification method(s)
  factory.BookMethod(&loader,TMVA::Types::kBDT, "BDT",
		     "!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20" );

  // Train
  factory.TrainAllMethods();

  // Test
  factory.TestAllMethods();
  factory.EvaluateAllMethods();

  // Plot a ROC Curve
  pad->cd();
  pad = factory.GetROCCurve(&loader);
  pad->Draw();

  pad->SaveAs("StrangeJetClassification_ROC.pdf");
				    

}

