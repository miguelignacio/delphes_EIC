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


#include <glob.h>
#include <iostream>
#include <iomanip>
#include <vector>

#include "PlotFunctions.h"


void CharmJetMVAOptimization(TString filename = "CharmJetClassification_Results.root",
			     TString foldername = "dataset",
			     TString classifiername = "MLP")
{
  // Global options
  gStyle->SetOptStat(0);


  // Scan the output for the optimal performance point and quote the results
  float target_xsec = LookupCrossSection("CC_DIS");
  float eff_xsec_light = target_xsec * 0.983;
  float eff_xsec_charm = target_xsec * 0.020;
  float target_lumi = 100*u_fb;

  auto result_file = TFile::Open(filename.Data());
  auto result_tree = static_cast<TTree*>(result_file->Get(Form("%s/TestTree", foldername.Data())));

  float max_charm_yield_signifiance = 0.0;
  float optimal_mva_cut = 0.0;
	float optimal_charm_efficiency = 0.0;
  float optimal_light_efficiency = 0.0;
	float optimal_charm_yield = 0.0;
  float optimal_light_yield = 0.0;


  float n_light_all = static_cast<float>(result_tree->GetEntries("(jet_flavor==21 || jet_flavor < 4)"));
  float n_charm_all = static_cast<float>(result_tree->GetEntries("(jet_flavor == 4)"));

  float w_light = target_lumi * eff_xsec_light / n_light_all;
  float w_charm = target_lumi * eff_xsec_charm / n_charm_all;

  TH1F* h_light = new TH1F("h_light","",200,0,1.0);
  h_light->Sumw2();
  auto h_charm = static_cast<TH1F*>(h_light->Clone("h_charm"));

  result_tree->Project(h_light->GetName(), classifiername.Data(), "jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && (jet_flavor==21 || jet_flavor < 4) && met_et > 10");
  result_tree->Project(h_charm->GetName(), classifiername.Data(), "jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && (jet_flavor==4) && met_et > 10");


  Int_t nbins = h_light->GetNbinsX();

  for (Int_t i = 1; i <= h_light->GetNbinsX(); i++) {
    float min = h_light->GetBinLowEdge(i);
    std::cout << min << std::endl;

    float n_light = h_light->Integral(i, nbins);
    float n_charm = h_charm->Integral(i, nbins);

    float eff_light = n_light/n_light_all;
    float eff_charm = n_charm/n_charm_all;


    // Scale the number of light and charm yields to target_lumi
    n_light = n_light * w_light;
    n_charm = n_charm * w_charm;
    float n_total = n_light + n_charm;

    float err_light = TMath::Sqrt(n_light);
    float err_charm = TMath::Sqrt(n_charm);


    float err_yield = TMath::Sqrt(n_total + n_light);
    float yield_significance = n_charm / err_yield;

    if (yield_significance > max_charm_yield_signifiance) {
      max_charm_yield_signifiance = yield_significance;
      optimal_mva_cut = min;
			optimal_charm_efficiency = eff_charm;
      optimal_light_efficiency = eff_light;
			optimal_charm_yield = n_charm;
      optimal_light_yield = n_light;
    }

  }


  // for (float min = 0.0; min <= 1.0; min += 0.02) {
  //   std::cout << min << std::endl;

  //   float n_light = static_cast<float>(result_tree->GetEntries(Form("jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && (jet_flavor==21 || jet_flavor < 4) && met_et > 10 && MLP > %.5f", min)));
  //   float n_charm = static_cast<float>(result_tree->GetEntries(Form("jet_pt>5.0 && TMath::Abs(jet_eta) < 3.0 && (jet_flavor == 4) && met_et > 10 && MLP > %.5f", min)));


  //   float eff_light = n_light/n_light_all;
  //   float eff_charm = n_charm/n_charm_all;

  //   // Scale the number of light and charm yields to target_lumi
  //   n_light = n_light * w_light;
  //   n_charm = n_charm * w_charm;
  //   float n_total = n_light + n_charm;

  //   float err_light = TMath::Sqrt(n_light);
  //   float err_charm = TMath::Sqrt(n_charm);


  //   float err_yield = TMath::Sqrt(n_total + n_light);
  //   float yield_significance = n_charm / err_yield;

  //   if (yield_significance > max_charm_yield_signifiance) {
  //     max_charm_yield_signifiance = yield_significance;
  //     optimal_mva_cut = min;
  //     optimal_charm_efficiency = eff_charm;
  //     optimal_light_efficiency = eff_light;
  //   }


  // }

  std::cout << "Optimal MVA cut: " << optimal_mva_cut
	    << " with charm yield significance of " << max_charm_yield_signifiance
	    << " charm (light) efficiency " << optimal_charm_efficiency << " (" << optimal_light_efficiency << ")"
	    << " in 100/fb" << std::endl;
	std::cout << "Estimated charm (light) yield: " << optimal_charm_yield << "(" << optimal_light_yield << ")" << std::endl;

}
