#include "TROOT.h"
#include "TChain.h"
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


float target_xsec = LookupCrossSection("CC_DIS");

float target_lumi = 100 * u_fb;

float n_gen = -1;


void CharmJetClassification_Plots()
{
  // Global options
  gStyle->SetOptStat(0);

  // Create the TCanvas
  TCanvas *pad       = new TCanvas("pad", "", 800, 600);
  TLegend *legend    = nullptr;
  TH1F    *htemplate = nullptr;

  std::vector<TString> input_folders = { "dataset_ip3dtagger", "dataset_ktagger", "dataset_etagger", "dataset_mutagger" };
  std::vector<TString> input_vars    = { "jet_t1_pt", "jet_t1_sIP3D", "jet_t2_pt", "jet_t2_sIP3D", "jet_t3_pt", "jet_t3_sIP3D", "jet_t4_pt", "jet_t4_sIP3D", "jet_k1_pt", "jet_k1_sIP3D", "jet_k2_pt", "jet_k2_sIP3D", "jet_e1_pt", "jet_e1_sIP3D", "jet_e2_pt", "jet_e2_sIP3D", "jet_mu1_pt", "jet_mu1_sIP3D", "jet_mu2_pt", "jet_mu2_sIP3D" };

  for (auto input_folder : input_folders) {
    std::cout << input_folder.Data() << std::endl;

    for (auto input_var : input_vars) {
      std::cout << input_var.Data() << std::endl;

      TChain testing_data(Form("%s/TestTree", input_folder.Data()));
      TChain training_data(Form("%s/TrainTree", input_folder.Data()));

      testing_data.Add("CharmJetClassification_Results.root");
      training_data.Add("CharmJetClassification_Results.root");

      pad->Clear();

      // pad->Divide(1,2);
      pad->cd();

      TH1F *h_signal     = nullptr;
      TH1F *h_background = nullptr;

      if (input_var.Contains("pt") == kTRUE) {
        h_signal = new TH1F("h_signal", "signal distribution", 100, 0, 15);
        h_signal->Sumw2();
        pad->SetLogy(kFALSE);
      }

      if (input_var.Contains("sIP3D") == kTRUE) {
        h_signal = new TH1F("h_signal", "signal distribution", 200, -100, 100);
        h_signal->Sumw2();
        pad->SetLogy();
      }

      if (h_signal == nullptr) {
        continue;
      }

      h_background = static_cast<TH1F *>(h_signal->Clone("h_background"));
      training_data.Project(h_signal->GetName(),     input_var.Data(), "jet_flavor==4");
      training_data.Project(h_background->GetName(), input_var.Data(), "jet_flavor<4 || jet_flavor==22");

      if ((h_signal->GetEntries() == 0) || (h_background->GetEntries() == 0)) {
        if (h_signal) {
          delete h_signal;
        }

        if (h_background) {
          delete h_background;
        }
        continue;
      }

      auto h_signal_norm     = h_signal->DrawNormalized("HIST");
      auto h_background_norm = h_background->DrawNormalized("HIST SAME");

      h_signal_norm->SetLineWidth(2);
      h_signal_norm->SetLineColor(kBlue);

      h_background_norm->SetLineWidth(2);
      h_background_norm->SetLineColor(kRed);
      h_background_norm->SetLineStyle(kDashed);

      Float_t max_value = TMath::Max(h_signal_norm->GetMaximum(), h_background_norm->GetMaximum());

      Float_t min_value = 0.0;

      if (input_var.Contains("sIP3D") == kTRUE) {
        min_value = 1e-6;
      }

      h_signal_norm->SetMaximum(max_value * 1.25);
      h_background_norm->SetMaximum(max_value * 1.25);
      h_signal_norm->SetMinimum(min_value);
      h_background_norm->SetMinimum(min_value);

      // Label axes
      h_signal_norm->SetYTitle("Probability");

      if (input_var.Contains("sIP3D") == kTRUE) {
        h_signal_norm->SetXTitle("Track sIP_{3D}");
      }

      if (input_var.Contains("pt") == kTRUE) {
        h_signal_norm->SetXTitle("Track p_{T}");
      }

      pad->SetGrid(1, 1);

      // Legend
      TLegend *legend = smart_legend("upper right");
      legend->SetFillStyle(0);
      legend->SetBorderSize(0);
      legend->AddEntry(h_signal_norm,     "True Charm Jets", "lf");
      legend->AddEntry(h_background_norm, "True Light Jets", "lp");
      legend->Draw();


      pad->SaveAs(Form("%s.pdf", input_var.Data()));

      if (h_signal) {
        delete h_signal;
      }

      if (h_background) {
        delete h_background;
      }

      if (legend)
        delete legend;
    }
  }
}
