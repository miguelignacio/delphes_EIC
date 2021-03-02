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


void CharmJetClassification(TString dir, TString input, TString filePattern = "*/out.root")
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

  auto outputFile = TFile::Open("CharmJetClassification_Results.root", "RECREATE");

  TMVA::Factory factory("TMVAClassification",
                        outputFile,
                        "!V:ROC:!Correlations:!Silent:Color:DrawProgressBar:AnalysisType=Classification");

  //
  // Kaon Tagger
  //

  TMVA::DataLoader loader_ktagger("dataset_ktagger");
  TMVA::DataLoader loader_etagger("dataset_etagger");
  TMVA::DataLoader loader_mutagger("dataset_mutagger");
  TMVA::DataLoader loader_ip3dtagger("dataset_ip3dtagger");

  // loader_ktagger.AddVariable("jet_pt",  "Jet p_{T}",          "GeV", 'F',
  // 0.0,
  // 1000.0);
  // loader_ktagger.AddVariable("jet_eta", "Jet Pseudorapidity", "",    'F',
  // -5.0, 5.0);

  loader_ktagger.AddSpectator("jet_pt");
  loader_ktagger.AddSpectator("jet_eta");
  loader_ktagger.AddVariable("jet_k1_pt");
  loader_ktagger.AddVariable("jet_k1_sIP3D");
  loader_ktagger.AddVariable("jet_k2_pt");
  loader_ktagger.AddVariable("jet_k2_sIP3D");
  loader_ktagger.AddSpectator("jet_flavor");
  loader_ktagger.AddSpectator("met_et");
  loader_ktagger.AddSignalTree(signal_train, 1.0);
  loader_ktagger.AddBackgroundTree(background_train, 1.0);

  // loader_ktagger.AddVariable("jet_nconstituents");

  loader_etagger.AddSpectator("jet_pt");
  loader_etagger.AddSpectator("jet_eta");
  loader_etagger.AddVariable("jet_e1_pt");
  loader_etagger.AddVariable("jet_e1_sIP3D");
  loader_etagger.AddVariable("jet_e2_pt");
  loader_etagger.AddVariable("jet_e2_sIP3D");
  loader_etagger.AddSpectator("jet_flavor");
  loader_etagger.AddSpectator("met_et");
  loader_etagger.AddSignalTree(signal_train, 1.0);
  loader_etagger.AddBackgroundTree(background_train, 1.0);

  loader_mutagger.AddSpectator("jet_pt");
  loader_mutagger.AddSpectator("jet_eta");
  loader_mutagger.AddVariable("jet_mu1_pt");
  loader_mutagger.AddVariable("jet_mu1_sIP3D");
  loader_mutagger.AddVariable("jet_mu2_pt");
  loader_mutagger.AddVariable("jet_mu2_sIP3D");
  loader_mutagger.AddSpectator("jet_flavor");
  loader_mutagger.AddSpectator("met_et");
  loader_mutagger.AddSignalTree(signal_train, 1.0);
  loader_mutagger.AddBackgroundTree(background_train, 1.0);


  loader_ip3dtagger.AddSpectator("jet_pt");
  loader_ip3dtagger.AddSpectator("jet_eta");
  loader_ip3dtagger.AddVariable("jet_t1_pt");
  loader_ip3dtagger.AddVariable("jet_t1_sIP3D");
  loader_ip3dtagger.AddVariable("jet_t2_pt");
  loader_ip3dtagger.AddVariable("jet_t2_sIP3D");
  loader_ip3dtagger.AddVariable("jet_t3_pt");
  loader_ip3dtagger.AddVariable("jet_t3_sIP3D");
  loader_ip3dtagger.AddVariable("jet_t4_pt");
  loader_ip3dtagger.AddVariable("jet_t4_sIP3D");
  loader_ip3dtagger.AddSpectator("jet_flavor");
  loader_ip3dtagger.AddSpectator("met_et");
  loader_ip3dtagger.AddSignalTree(signal_train, 1.0);
  loader_ip3dtagger.AddBackgroundTree(background_train, 1.0);

  //  loader.AddVariable("jet_sip3dtag", "sIP3D Jet-Level Tag", "", 'B', -10,
  // 10);
  // loader.AddVariable("jet_charge");

  // loader.AddVariable("jet_ehadoveremratio");


  // loader.AddTree( signal_train, "strange_jets" );
  // loader.AddTree( u_background_train, "up jets" );
  // loader.AddTree( d_background_train, "down jets" );
  // loader.AddTree( g_background_train, "gluon jets" );

  loader_ktagger.PrepareTrainingAndTestTree(TCut("jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && jet_flavor==4 && met_et > 10"),
                                            TCut("jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && (jet_flavor<4||jet_flavor==21) && met_et > 10"),
                                            "nTrain_Signal=50000:nTrain_Background=500000:nTest_Signal=50000:nTest_Background=500000:SplitMode=Random:NormMode=NumEvents:!V");
  loader_etagger.PrepareTrainingAndTestTree(TCut("jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && jet_flavor==4 && met_et > 10"),
                                            TCut("jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && (jet_flavor<4||jet_flavor==21) && met_et > 10"),
                                            "nTrain_Signal=100000:nTrain_Background=1000000:nTest_Signal=100000:nTest_Background=1000000:SplitMode=Random:NormMode=NumEvents:!V");
  loader_mutagger.PrepareTrainingAndTestTree(TCut("jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && jet_flavor==4 && met_et > 10"),
                                             TCut("jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && (jet_flavor<4||jet_flavor==21) && met_et > 10"),
                                             "nTrain_Signal=100000:nTrain_Background=1000000:nTest_Signal=100000:nTest_Background=1000000:SplitMode=Random:NormMode=NumEvents:!V");
  loader_ip3dtagger.PrepareTrainingAndTestTree(TCut("jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && jet_flavor==4 && met_et > 10"),
                                               TCut("jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && (jet_flavor<4||jet_flavor==21) && met_et > 10"),
                                               "nTrain_Signal=10000:nTrain_Background=100000:nTest_Signal=10000:nTest_Background=100000:SplitMode=Random:NormMode=NumEvents:!V");

  // Declare the classification method(s)
  // factory.BookMethod(&loader,TMVA::Types::kBDT, "BDT",
  //
  //
  //  "!V:NTrees=1000:MinNodeSize=2.5%:MaxDepth=4:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20"
  // );
  factory.BookMethod(&loader_ktagger,    TMVA::Types::kMLP, "CharmKTagger",
                     "!H:!V:NeuronType=ReLU:VarTransform=Norm:NCycles=1000:HiddenLayers=N+12:TestRate=5:!UseRegulator");
  // factory.BookMethod(&loader_etagger,    TMVA::Types::kMLP, "CharmETagger",
  //                    "!H:!V:NeuronType=ReLU:VarTransform=Norm:NCycles=1000:HiddenLayers=N+12:TestRate=5:!UseRegulator");
  // factory.BookMethod(&loader_mutagger,   TMVA::Types::kMLP, "CharmMuTagger",
  //                    "!H:!V:NeuronType=ReLU:VarTransform=Norm:NCycles=1000:HiddenLayers=N+12:TestRate=5:!UseRegulator");
  // factory.BookMethod(&loader_ip3dtagger, TMVA::Types::kMLP, "CharmIP3DTagger",
  //                    "!H:!V:NeuronType=ReLU:VarTransform=Norm:NCycles=1000:HiddenLayers=N+16:TestRate=5:!UseRegulator");

  // factory.BookMethod(&loader, TMVA::Types::kMLP,
  // "CharmETagger","!H:!V:NeuronType=ReLU:VarTransform=Norm:NCycles=600:HiddenLayers=N+8:TestRate=5:!UseRegulator");
  // factory.BookMethod(&loader, TMVA::Types::kMLP,
  // "CharmMuTagger","!H:!V:NeuronType=ReLU:VarTransform=Norm:NCycles=600:HiddenLayers=N+8:TestRate=5:!UseRegulator");

  // Train
  factory.TrainAllMethods();

  // Test
  factory.TestAllMethods();
  factory.EvaluateAllMethods();

  // Plot a ROC Curve
  pad->cd();
  pad = factory.GetROCCurve(&loader_ktagger);
  // pad = factory.GetROCCurve(&loader_etagger);
  // pad = factory.GetROCCurve(&loader_mutagger);
  // pad = factory.GetROCCurve(&loader_ip3dtagger);
  pad->Draw();

  pad->SaveAs("CharmJetClassification_ROC.pdf");
}
