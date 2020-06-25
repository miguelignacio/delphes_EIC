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

void StrangeTagPlots(TString dir, TString input, TString xvar, TString filePattern = "*/out.root")
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


  // Setup plot configurations
  TString common_title = "NC-DIS, 10GeVx275GeV, Q^{2}>100 GeV^{2}";

  plot_config strangejets;
  strangejets.markerstyle = kFourSquaresX;
  strangejets.markercolor = kGreen+3;
  strangejets.linestyle   = kSolid;
  strangejets.linecolor   = kGreen+3;
  
  plot_config lightjets;
  lightjets.markerstyle = kFullCircle;
  lightjets.markercolor = kRed+1;
  lightjets.linestyle   = kSolid;
  lightjets.linecolor   = kRed+1;
  
  // Compare the z shapes for light and light jets
  strangejets.htemplate = new TH1F("ks_zhadron_template","",100,0,1.5);
  strangejets.xtitle = "K_{s}^{0} z_{hadron}";
  strangejets.ytitle = "K_{s}^{0} Candidates in Jets";

  lightjets.htemplate = strangejets.htemplate;
  lightjets.xtitle = strangejets.xtitle;
  lightjets.ytitle = strangejets.ytitle;

  auto strange_Ks_zhadron = GeneratePlot(strangejets, default_data, "jet_Ks_zhadron_strange", "jet_Ks_zhadron", TCut("jet_flavor==3"));
  auto light_Ks_zhadron   = GeneratePlot(lightjets, default_data, "jet_Ks_zhadron_light", "jet_Ks_zhadron", TCut("(jet_flavor<3 || jet_flavor==21)"));

  pad->cd();

  light_Ks_zhadron->DrawNormalized("HIST");
  strange_Ks_zhadron->DrawNormalized("E1 SAME");
  
  legend = smart_legend();
  legend->AddEntry(light_Ks_zhadron, "Light Jets", "lp");
  legend->AddEntry(strange_Ks_zhadron, "Strange Jets", "lp");
  legend->Draw();

  make_title(common_title);

  pad->SaveAs("Ks_zhadron.pdf");


  // Charged Kaons
  strangejets.htemplate = new TH1F("k_zhadron_template","",100,0,1.5);
  strangejets.xtitle = "K^{#pm} z_{hadron}";
  strangejets.ytitle = "K^{#pm} Candidates in Jets";

  lightjets.htemplate = strangejets.htemplate;
  lightjets.xtitle = strangejets.xtitle;
  lightjets.ytitle = strangejets.ytitle;

  auto strange_K_zhadron = GeneratePlot(strangejets, default_data, "jet_K_zhadron_strange", "jet_K_zhadron", TCut("jet_flavor==3"));
  auto light_K_zhadron   = GeneratePlot(lightjets, default_data, "jet_K_zhadron_light", "jet_K_zhadron", TCut("(jet_flavor<3 || jet_flavor==21)"));

  pad->cd();

  light_K_zhadron->DrawNormalized("HIST");
  strange_K_zhadron->DrawNormalized("E1 SAME");
  
  legend = smart_legend();
  legend->AddEntry(light_K_zhadron, "Light Jets", "lp");
  legend->AddEntry(strange_K_zhadron, "Strange Jets", "lp");
  legend->Draw();

  make_title(common_title);

  pad->SaveAs("K_zhadron.pdf");


  // Jet Charge
  strangejets.htemplate = new TH1F("jetcharge_template","",100,-1.5,1.5);
  strangejets.xtitle = "Jet Charge (#kappa = 0.5)";
  strangejets.ytitle = "Jet Candidates";

  lightjets.htemplate = strangejets.htemplate;
  lightjets.xtitle = strangejets.xtitle;
  lightjets.ytitle = strangejets.ytitle;

  auto ujets = lightjets;
  ujets.linestyle = kDashed;

  auto djets = lightjets;
  djets.linestyle = kSolid;

  auto gjets = lightjets;
  gjets.linestyle = kSolid;
  gjets.linecolor = kGray;

  


  auto strange_jetcharge = GeneratePlot(strangejets, default_data, "jet_charge_strange", "jet_charge", TCut("jet_flavor==3"));
  auto u_jetcharge   = GeneratePlot(ujets, default_data, "jet_charge_u", "jet_charge", TCut("(jet_flavor==2)"));
  auto d_jetcharge   = GeneratePlot(djets, default_data, "jet_charge_d", "jet_charge", TCut("(jet_flavor==1)"));
  auto g_jetcharge   = GeneratePlot(gjets, default_data, "jet_charge_g", "jet_charge", TCut("(jet_flavor==21)"));

  pad->cd();

  auto hist = g_jetcharge->DrawNormalized("HIST");
  u_jetcharge->DrawNormalized("HIST SAME");
  d_jetcharge->DrawNormalized("HIST SAME");
  strange_jetcharge->DrawNormalized("E1 SAME");
  
  hist->SetAxisRange(0.0, 0.035, "Y");

  legend = smart_legend();
  legend->AddEntry(g_jetcharge, "Gluon Jets", "lf");
  legend->AddEntry(u_jetcharge, "Up Jets", "lf");
  legend->AddEntry(d_jetcharge, "Down Jets", "lf");
  legend->AddEntry(strange_jetcharge, "Strange Jets", "lp");
  legend->Draw();

  make_title(common_title);


  pad->Update();
  pad->SaveAs("jet_charge.pdf");



}
