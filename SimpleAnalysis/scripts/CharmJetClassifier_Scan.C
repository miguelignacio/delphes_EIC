#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TTreeFormula.h"
#include <iostream>
#include "TMath.h"
#include "TH1.h"
#include "TCut.h"
#include <map>

void CharmJetClassifier_Scan(TString inputfile = "CharmJetClassification_Results.root",
			     TString foldername = "dataset",
			     TString taggername = "MLP",
			     Float_t bkg_eff = 0.004) {
  TFile *_file0    = TFile::Open(inputfile.Data());
  TTree *TestTree  = static_cast < TTree * > (_file0->Get(Form("%s/TestTree", foldername.Data())));
  TTree *TrainTree = static_cast < TTree * > (_file0->Get(Form("%s/TrainTree", foldername.Data())));

  Float_t nlight = TestTree->GetEntries("jet_flavor<4 || jet_flavor==22");
  Float_t ncharm = TestTree->GetEntries("jet_flavor==4");


  Float_t opt_cut  = 0.000;
  Float_t min_dist = 1000.0;

  for (Float_t mva_cut = 0.000; mva_cut < 1.000; mva_cut += 0.001) {
    Float_t nlight_pass = TestTree->GetEntries(Form("jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && met_et > 10.0 && (jet_flavor<4 || jet_flavor==22) && %s>%f", taggername.Data(), mva_cut));
    Float_t ncharm_pass = TestTree->GetEntries(Form("jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && met_et > 10.0 && (jet_flavor==4) && %s>%f", taggername.Data(), mva_cut));

    Float_t light_eff = nlight_pass / nlight;

    if (TMath::Abs(light_eff - bkg_eff) < min_dist) {
      min_dist = TMath::Abs(light_eff - bkg_eff);
      opt_cut  = mva_cut;
    }
  }

  // Print the final selection information
  Float_t nlight_pass = TestTree->GetEntries(Form("jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && met_et > 10.0 && (jet_flavor<4 || jet_flavor==22) && %s>%f", taggername.Data(), opt_cut));
  Float_t ncharm_pass = TestTree->GetEntries(Form("jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && met_et > 10.0 && (jet_flavor==4) && %s>%f", taggername.Data(), opt_cut));

  std::cout << Form("MLP cut %.3f", opt_cut) << " yields Eff_c = "
            << Form("%.3f", ncharm_pass / ncharm) << " and Eff_light = "
            << Form("%.3e", nlight_pass / nlight) << std::endl;


  /// Use the K-S Test to Compare the Training/Testing Shapes in Signal and
  // Background.
  std::cout << "Over-Training Study" << std::endl;
  std::cout << "=============================================" << std::endl;
  std::map < TString, TCut > samples;
  samples["Signal"]     = TCut("jet_flavor==4");
  samples["Background"] = TCut("(jet_flavor<4 || jet_flavor==22)");

  for (auto sample : samples) {
    TH1F *h_train = new TH1F("h_train",
                             "",
                             100,
                             0,
                             1);
    h_train->Sumw2();
    auto h_test = static_cast < TH1F * > (h_train->Clone("h_test"));

    TrainTree->Project(h_train->GetName(), Form("%s",taggername.Data()), sample.second);
    TestTree->Project(h_test->GetName(), Form("%s",taggername.Data()), sample.second);
    Float_t KS_Test = h_test->KolmogorovTest(h_train);
    std::cout << sample.first << ": K-S Test p-value = " << KS_Test << std::endl;

    if (h_test) delete h_test;

    if (h_train) delete h_train;
  }
}
