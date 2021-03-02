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
#include "TFile.h"

#include <glob.h>
#include <iostream>
#include <iomanip>
#include <vector>

#include "PlotFunctions.h"


float target_xsec = LookupCrossSection("CC_DIS");

float target_lumi = 100 * u_fb;

float n_gen = -1;


TEfficiency DifferentialTaggingEfficiency(TTree *data, plot_config draw_config, TString xvar, TString which = "all")
{
  std::cout << "DifferentialTaggingEfficiency: processing " << which.Data() << " ... " << std::endl;

  TH1F *tru_yield = static_cast<TH1F *>(draw_config.htemplate->Clone(Form("tru_yield_%d", histUID)));
  TH1F *tag_yield = static_cast<TH1F *>(draw_config.htemplate->Clone(Form("tag_yield_%d", histUID)));

  TCut *tru_selection = nullptr;
  TCut *tag_selection = nullptr;

  if (which == "light") {
    tru_selection = new TCut("jet_flavor < 4 || jet_flavor == 21");
    tag_selection = new TCut(*tru_selection && TCut(draw_config.cuts));
  } else if (which == "charm") {
    tru_selection = new TCut("jet_flavor == 4");
    tag_selection = new TCut(*tru_selection && TCut(draw_config.cuts));
  }


  data->Project(tru_yield->GetName(), xvar.Data(), *tru_selection);
  data->Project(tag_yield->GetName(), xvar.Data(), *tag_selection);

  std::cout << "  True Yield: " << tru_yield->Integral() << std::endl;
  std::cout << "  Tag  Yield: " << tag_yield->Integral() << std::endl;

  TEfficiency eff(*tag_yield, *tru_yield);

  if (tru_yield)
    delete tru_yield;

  if (tag_yield)
    delete tag_yield;

  histUID++;

  return eff;
}

TH1F* DifferentialTaggingYield(TTree *data, plot_config draw_config, TString xvar, TString which = "all")
{
  std::cout << "DifferentialTaggingYield: processing " << which.Data() << " ... " << std::endl;

  TH1F *tag_yield = static_cast<TH1F *>(draw_config.htemplate->Clone(Form("%s_%s_tag_yield_%d", which.Data(), xvar.Data(), histUID)));
  tag_yield->Sumw2();

  TCut *tag_selection = nullptr;

  if (which == "light") {
    tag_selection = new TCut(TCut("(jet_flavor < 4 || jet_flavor == 21)") && TCut(draw_config.cuts));
  } else if (which == "charm") {
    tag_selection = new TCut(TCut("jet_flavor == 4") && TCut(draw_config.cuts));
  }

  // Add other analysis-level cuts to the tag selection
  *tag_selection = (*tag_selection) && TCut("met_et > 10.0");


  data->Project(tag_yield->GetName(), xvar.Data(), *tag_selection);

  std::cout << "  Unnormalized Tag Yield: " << tag_yield->GetEntries() << std::endl;
  tag_yield->Scale(target_xsec * target_lumi / data->GetEntries());
  std::cout << "    Normalized Tag Yield: " << tag_yield->Integral() << std::endl;


  // auto tag_graph = new TGraphErrors(tag_yield);

  histUID++;
  return tag_yield;
}

void DrawTagEfficiencyPlot(TCanvas *pad, TTree *data, plot_config draw_config, TString xvar)
{
  auto light_eff = DifferentialTaggingEfficiency(data, draw_config, xvar, "light");
  auto charm_eff = DifferentialTaggingEfficiency(data, draw_config, xvar, "charm");


  // make a template histogram to fine-tune layout of the final plot
  TH1F *htemplate = new TH1F("htemplate", "", 1, draw_config.xlimits[0], draw_config.xlimits[1]);

  set_axis_range(htemplate, draw_config.ylimits[0], draw_config.ylimits[1], "Y");
  set_axis_range(htemplate, draw_config.xlimits[0], draw_config.xlimits[1], "X");

  // Configure plots
  configure_plot<TEfficiency>(&charm_eff, draw_config, "charm");
  configure_plot<TEfficiency>(&light_eff, draw_config, "light");

  // Draw Layout
  htemplate->Draw("HIST");
  charm_eff.Draw("P E1 SAME");
  light_eff.Draw("E1 P SAME");

  // Title
  TLatex plot_title = make_title();
  plot_title.Draw("SAME");

  // Configure the Pad
  pad->SetLogy(draw_config.logy);
  pad->SetLogx(draw_config.logx);
  pad->SetGrid(1, 1);


  // Configure the Legend
  TLegend *legend = smart_legend("lower right");
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(&charm_eff, "Charm Jets", "lp");
  legend->AddEntry(&light_eff, "Light Jets", "lp");
  legend->Draw();

  pad->Update();
  pad->SaveAs(Form("CharmTagPlot_tagging_efficiency_%s_%s.pdf", data->GetTitle(), xvar.Data()));
  pad->SaveAs(Form("CharmTagPlot_tagging_efficiency_%s_%s.C", data->GetTitle(), xvar.Data()));

  // Cleanup
  if (htemplate != nullptr)
    delete htemplate;

  if (legend != nullptr)
    delete legend;
}

void DrawTagYieldPlot(TCanvas *pad, TTree *data, plot_config draw_config, TString xvar)
{
  std::cout << "Drawing tag yield plot for variable " << xvar << std::endl;
  auto charm_yield = DifferentialTaggingYield(data, draw_config, xvar, "charm");
  auto light_yield = DifferentialTaggingYield(data, draw_config, xvar, "light");

  // Configure plots

  configure_plot<TH1F>(charm_yield, draw_config, "charm");
  configure_plot<TH1F>(light_yield, draw_config, "light");

  // generate template histogram
  TH1F *htemplate = new TH1F("htemplate", "", 1, draw_config.xlimits[0], draw_config.xlimits[1]);

  set_axis_range(htemplate, draw_config.ylimits[0], draw_config.ylimits[1], "Y");
  set_axis_range(htemplate, draw_config.xlimits[0], draw_config.xlimits[1], "X");
  htemplate->SetXTitle(draw_config.xtitle);
  htemplate->SetYTitle(draw_config.ytitle);

  pad->SetGrid(1, 1);
  pad->SetLogy(draw_config.logy);
  pad->SetLogx(draw_config.logx);

  htemplate->Draw("HIST");
  charm_yield->Draw("E1 P SAME");
  light_yield->Draw("E1 P SAME");

  TLatex plot_title = make_title();
  plot_title.Draw("SAME");


  // Legend
  TLegend *legend = smart_legend("upper right");
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->AddEntry(charm_yield, "Charm Jets", "lp");
  legend->AddEntry(light_yield, "Light Jets", "lp");
  legend->Draw();
  pad->Update();


  pad->SaveAs(Form("CharmTagPlot_tagging_yield_%s_%s.pdf", data->GetTitle(), xvar.Data()));
  pad->SaveAs(Form("CharmTagPlot_tagging_yield_%s_%s.C", data->GetTitle(), xvar.Data()));

  // Cleanup
  if (htemplate != nullptr)
    delete htemplate;

  if (legend != nullptr)
    delete legend;
}

float CutEfficiency(TTree *data, TCut cut, TString cutdesc = "")
{
  if (n_gen < 0) {
    n_gen = static_cast<float>(data->GetEntries());
  }

  float n_yield = data->GetEntries(cut.GetTitle());

  TString desc = cutdesc;

  if (desc == "") {
    desc = cut.GetTitle();
  }

  std::cout << setw(30) << desc << " efficiency: " << setprecision(4) << n_yield / n_gen * 100.0 << "%" << std::endl;

  return n_yield / n_gen;
}

void CharmTagPlots(TString dir, TString input, TString xvar, TString filePattern = "*/out.root")
{
  // Global options
  gStyle->SetOptStat(0);

  // Create the TCanvas
  TCanvas *pad       = new TCanvas("pad", "", 800, 600);
  TLegend *legend    = nullptr;
  TH1F    *htemplate = nullptr;

  auto default_data = new TChain("tree");
  default_data->SetTitle(input.Data());
  auto files = fileVector(Form("%s/%s/%s", dir.Data(), input.Data(), filePattern.Data()));

  for (auto file : files)
  {
    default_data->Add(file.c_str());
  }


  // Some basic cut efficiency information
  std::cout << "Inclusive Jet Efficiency Information: " << std::endl;
  TCut selection("jet_n>0");
  CutEfficiency(default_data, selection, "Jet Reconstructed");

  selection = selection && TCut("TMath::Abs(jet_eta) < 3.0");
  CutEfficiency(default_data, selection, "Jet Fiducial");

  selection = selection && TCut("met_et > 10.0");
  CutEfficiency(default_data, selection, "MET > 10 GeV");

  TCut main_preselection = selection;

  // Print charm jet MET cut efficiency
  std::cout << "Charm Jet Efficiency Information: " << std::endl;
  selection = TCut("jet_n > 0 && TMath::Abs(jet_eta) < 3.0 && jet_flavor == 4");
  CutEfficiency(default_data, selection, "Charm Jet Reconstructed");

  selection = selection && TCut("met_et > 10.0");
  CutEfficiency(default_data, selection, "MET > 10 GeV");


  TTree *default_data_selected = default_data->CopyTree(main_preselection.GetTitle());
  default_data_selected->SetTitle(input.Data());

  bool do_TagEffPlot    = kFALSE;
  bool do_TagYieldPlot  = kFALSE;  // kTRUE;
  bool do_100fbProjPlot = kTRUE; // kTRUE;
  bool do_Helicity      = kFALSE; // kTRUE;


  // plot configurations
  plot_config draw_config;

  if (xvar == "jet_pt") {
    float xbins[] = { 10, 12.5, 15, 17.5, 20, 25, 35, 60 };
    int   nbins   = sizeof(xbins) / sizeof(xbins[0]) - 1;
    draw_config.htemplate = new TH1F(xvar, "", nbins, xbins);
    draw_config.xlimits   = std::vector<float>();
    draw_config.xlimits.push_back(0);
    draw_config.xlimits.push_back(60);
    draw_config.ylimits = std::vector<float>();
    draw_config.ylimits.push_back(1e-3);
    draw_config.ylimits.push_back(1.0);
    draw_config.xtitle = "Reconstructed Jet p_{T} [GeV]";
    draw_config.ytitle = "#varepsilon_{tag}";
    draw_config.logy   = kTRUE;
    draw_config.logx   = kFALSE;
    draw_config.cuts   = "jet_sip3dtag==1";
  } else if (xvar == "bjorken_x") {
    float xbins[] = { 0.01, 0.043333, 0.076666, 0.1, 0.25, 0.5, 1.0 };
    int   nbins   = sizeof(xbins) / sizeof(xbins[0]) - 1;
    draw_config.htemplate = new TH1F(xvar, "", nbins, xbins);
    draw_config.xlimits   = std::vector<float>();
    draw_config.xlimits.push_back(0.01);
    draw_config.xlimits.push_back(1.0);
    draw_config.ylimits = std::vector<float>();
    draw_config.ylimits.push_back(1e-3);
    draw_config.ylimits.push_back(1.0);
    draw_config.xtitle = "Bjorken x";
    draw_config.ytitle = "#varepsilon_{tag}";
    draw_config.logy   = kTRUE;
    draw_config.logx   = kTRUE;
    draw_config.cuts   = "jet_sip3dtag==1";
  } else if (xvar == "jb_x") {
    float xbins[] = { 0.01, 0.043333, 0.076666, 0.1, 0.25, 0.5, 1.0 };
    int   nbins   = sizeof(xbins) / sizeof(xbins[0]) - 1;
    draw_config.htemplate = new TH1F(xvar, "", nbins, xbins);
    draw_config.xlimits   = std::vector<float>();
    draw_config.xlimits.push_back(0.01);
    draw_config.xlimits.push_back(1.0);
    draw_config.ylimits = std::vector<float>();
    draw_config.ylimits.push_back(1e-3);
    draw_config.ylimits.push_back(1.0);
    draw_config.xtitle = "Reconstructed x_{JB}";
    draw_config.ytitle = "#varepsilon_{tag}";
    draw_config.logy   = kTRUE;
    draw_config.logx   = kTRUE;
    draw_config.cuts   = "jet_sip3dtag==1";
  } else {
    do_TagEffPlot = kFALSE;
  }


  //
  // Efficiency Plot
  //

  if (do_TagEffPlot) {
    DrawTagEfficiencyPlot(pad, default_data_selected, draw_config, xvar);
  }

  //
  // Cleanup
  //

  pad->Clear();

  //
  // Yield plot
  //
  gStyle->SetErrorX(0.5);

  if (xvar == "jet_pt") {
    draw_config.ylimits = std::vector<float>();

    // draw_config.ylimits.push_back(0.1);
    // draw_config.ylimits.push_back(5000);
    draw_config.ylimits.push_back(0.0);
    draw_config.ylimits.push_back(2500);
    draw_config.logy = kFALSE;

    draw_config.xtitle = "Reconstructed Jet p_{T} [GeV]";
    draw_config.ytitle = "#varepsilon_{tag} #times #sigma_{CC-DIS} #times 100fb^{-1}";
    draw_config.cuts   = "jet_sip3dtag==1";
  } else if (xvar == "bjorken_x") {
    draw_config.ylimits = std::vector<float>();
    draw_config.ylimits.push_back(0.1);
    draw_config.ylimits.push_back(5000);
    draw_config.xtitle = "Bjorken x";
    draw_config.ytitle = "#varepsilon_{tag} #times #sigma_{CC-DIS} #times 100fb^{-1}";
    draw_config.cuts   = "jet_sip3dtag==1";
  } else if (xvar == "jb_x") {
    draw_config.ylimits = std::vector<float>();
    draw_config.ylimits.push_back(0.1);
    draw_config.ylimits.push_back(5000);
    draw_config.xtitle = "Reconstructed x_{JB}";
    draw_config.ytitle = "#varepsilon_{tag} #times #sigma_{CC-DIS} #times 100fb^{-1}";
    draw_config.cuts   = "jet_sip3dtag==1";
  } else if (xvar == "jet_mlp_ktagger") {
    draw_config.htemplate = new TH1F(xvar, "", 50, -0.00001, 1.00001);
    draw_config.xlimits   = std::vector<float>();
    draw_config.xlimits.push_back(0);
    draw_config.xlimits.push_back(1);
    draw_config.ylimits = std::vector<float>();
    draw_config.ylimits.push_back(0.1);
    draw_config.ylimits.push_back(1e7);
    draw_config.xtitle = "MLP Kaon-based Charm Jet Tagger Output";
    draw_config.ytitle = "#varepsilon_{tag} #times #sigma_{CC-DIS} #times 100fb^{-1}";
    draw_config.logy   = kTRUE;
    draw_config.logx   = kFALSE;
    draw_config.cuts   = "1==1";
  } else if (xvar == "jet_mlp_eltagger") {
    draw_config.htemplate = new TH1F(xvar, "", 50, -0.00001, 1.00001);
    draw_config.xlimits   = std::vector<float>();
    draw_config.xlimits.push_back(0);
    draw_config.xlimits.push_back(1);
    draw_config.ylimits = std::vector<float>();
    draw_config.ylimits.push_back(0.1);
    draw_config.ylimits.push_back(1e7);
    draw_config.xtitle = "MLP Electron-based Charm Jet Tagger Output";
    draw_config.ytitle = "#varepsilon_{tag} #times #sigma_{CC-DIS} #times 100fb^{-1}";
    draw_config.logy   = kTRUE;
    draw_config.logx   = kFALSE;
    draw_config.cuts   = "1==1";
  } else if (xvar == "jet_mlp_mutagger") {
    draw_config.htemplate = new TH1F(xvar, "", 50, -0.00001, 1.00001);
    draw_config.xlimits   = std::vector<float>();
    draw_config.xlimits.push_back(0);
    draw_config.xlimits.push_back(1);
    draw_config.ylimits = std::vector<float>();
    draw_config.ylimits.push_back(0.1);
    draw_config.ylimits.push_back(1e7);
    draw_config.xtitle = "MLP Muon-based Charm Jet Tagger Output";
    draw_config.ytitle = "#varepsilon_{tag} #times #sigma_{CC-DIS} #times 100fb^{-1}";
    draw_config.logy   = kTRUE;
    draw_config.logx   = kFALSE;
    draw_config.cuts   = "1==1";
  } else if (xvar == "jet_mlp_ip3dtagger") {
    draw_config.htemplate = new TH1F(xvar, "", 50, -0.00001, 1.00001);
    draw_config.xlimits   = std::vector<float>();
    draw_config.xlimits.push_back(0);
    draw_config.xlimits.push_back(1);
    draw_config.ylimits = std::vector<float>();
    draw_config.ylimits.push_back(0.1);
    draw_config.ylimits.push_back(1e7);
    draw_config.xtitle = "MLP sIP_{3D}-based Charm Jet Tagger Output";
    draw_config.ytitle = "#varepsilon_{tag} #times #sigma_{CC-DIS} #times 100fb^{-1}";
    draw_config.logy   = kTRUE;
    draw_config.logx   = kFALSE;
    draw_config.cuts   = "1==1";
  } else if (xvar == "jet_mlp_globaltagger") {
    draw_config.htemplate = new TH1F(xvar, "", 50, -0.00001, 1.00001);
    draw_config.xlimits   = std::vector<float>();
    draw_config.xlimits.push_back(0);
    draw_config.xlimits.push_back(1);
    draw_config.ylimits = std::vector<float>();
    draw_config.ylimits.push_back(0.1);
    draw_config.ylimits.push_back(1e7);
    draw_config.xtitle = "MLP Global Charm Jet Tagger Output";
    draw_config.ytitle = "#varepsilon_{tag} #times #sigma_{CC-DIS} #times 100fb^{-1}";
    draw_config.logy   = kTRUE;
    draw_config.logx   = kFALSE;
    draw_config.cuts   = "1==1";
  } else {
    do_TagYieldPlot = kFALSE;
  }


  if (do_TagYieldPlot) {
    DrawTagYieldPlot(pad, default_data_selected, draw_config, xvar);
  }

  //
  // Error Bands at 100/fb
  //

  if (xvar == "jet_pt") {
    draw_config.ylimits = std::vector<float>();
    draw_config.ylimits.push_back(0.1);
    draw_config.ylimits.push_back(5000);
    draw_config.xtitle = "Reconstructed Jet p_{T} [GeV]";
    draw_config.ytitle = "Relative Variation to Suppressed Strangeness";
  } else if (xvar == "bjorken_x") {
    draw_config.ytitle = "Relative Variation to Suppressed Strangeness";
  } else if (xvar == "jb_x") {
    draw_config.ytitle = "Relative Variation to Suppressed Strangeness";
  } else {
    do_100fbProjPlot = kFALSE;
  }

  bool do_BkgSubtraction = kTRUE;

  if (do_100fbProjPlot == kTRUE) {
    pad->Clear();
    gStyle->SetErrorX(0.5);

    auto charm_yield       = DifferentialTaggingYield(default_data_selected, draw_config, xvar, "charm");
    auto charm_yield_100fb = static_cast<TH1F *>(charm_yield->Clone("charm_yield_100fb"));

    for (Int_t i = 1; i <= charm_yield_100fb->GetNbinsX(); i++) {
      charm_yield_100fb->SetBinError(i, TMath::Sqrt(charm_yield_100fb->GetBinContent(i)));
    }

    TH1F *light_yield       = nullptr;
    TH1F *light_yield_100fb = nullptr;

    if (do_BkgSubtraction == kTRUE) {
      std::cout << "   For Yield Estimate, performing a naive background subtraction error estimation..." << std::endl;
      light_yield       = DifferentialTaggingYield(default_data_selected, draw_config, xvar, "light");
      light_yield_100fb = static_cast<TH1F *>(light_yield->Clone("light_yield_100fb"));

      for (Int_t i = 1; i <= light_yield_100fb->GetNbinsX(); i++) {
        light_yield_100fb->SetBinError(i, TMath::Sqrt(light_yield_100fb->GetBinContent(i)));
      }

      // inflate the charm yield errors to account for "background subtraction"
      for (Int_t i = 1; i <= light_yield_100fb->GetNbinsX(); i++) {
        float err_charm = charm_yield_100fb->GetBinError(i);
        float err_light = light_yield_100fb->GetBinError(i);
        float err_total = TMath::Sqrt(charm_yield_100fb->GetBinContent(i) + light_yield_100fb->GetBinContent(i));

        float new_charm_err = TMath::Sqrt(TMath::Power(err_total, 2.0) + TMath::Power(err_light, 2.0));


        charm_yield_100fb->SetBinError(i, new_charm_err);
      }
    }

    std::map<TString, TString> alt_samples;
    //alt_samples["CC_DIS_e10_p275_lha_21Rs2"] = "CT18ZNNLO with enhanced strangeness [R_{s}=0.863]";
    alt_samples["CC_DIS_e10_p275_lha_21Rs2"] = "Rs-High NNLO (enhanced strange)";
    alt_samples["CC_DIS_e10_p275_CT18ANNLO"] = "CT18A NNLO (intermediate strange)";
    //alt_samples["CC_DIS_e10_p275_lha_22Rs2"] = "CT18ZNNLO with intermediate strangeness (22Rs2)";
    //alt_samples["CC_DIS_e10_p275_lha_22Rs2ver3"] = "R_{s}^{INT}=0.594";
    //alt_samples["CC_DIS_e10_p275_MMHT2014nnlo68cl"] = "MMHT2014nnlo68cl";
    //alt_samples["CC_DIS_e10_p275_NNPDF31_nnlo_as_0118"] = "NNPDF31_nnlo_as_0118";

    std::map<TString, Int_t> alt_marker;
    alt_marker["CC_DIS_e10_p275_lha_21Rs2"] = kOpenDiamond;
    alt_marker["CC_DIS_e10_p275_CT18ANNLO"] = kOpenCrossX;
    alt_marker["CC_DIS_e10_p275_lha_22Rs2ver3" ] = kOpenCross;

    std::map<TString, Int_t> alt_color;
    alt_color["CC_DIS_e10_p275_lha_21Rs2"] = kBlue - 5;
    alt_color["CC_DIS_e10_p275_CT18ANNLO"] = kViolet - 5;
    alt_color["CC_DIS_e10_p275_lha_22Rs2ver3" ] = kSpring - 6;
    

    


    // Load alternative samples
    auto ccdis_20Rs2_data = new TChain("tree");
    files = fileVector(Form("%s/CC_DIS_e10_p275_lha_20Rs2/%s", dir.Data(), filePattern.Data()));
    for (auto file : files)
    {
      ccdis_20Rs2_data->Add(file.c_str());
    }

    std::cout << "Running pre-selection on alternative samples." << std::endl;
    std::cout << "... suppressed strangeness ..." << std::endl;
    TTree *ccdis_20Rs2_data_selected = ccdis_20Rs2_data->CopyTree(main_preselection.GetTitle());
    auto ccdis_20Rs2_yield = DifferentialTaggingYield(ccdis_20Rs2_data_selected, draw_config, xvar, "charm");
    auto ccdis_20Rs2_yield_100fb = static_cast<TH1F *>(ccdis_20Rs2_yield->Clone("20Rs2_yield_100fb"));

    for (Int_t i = 1; i <= charm_yield_100fb->GetNbinsX(); i++) {
      // charm_yield_100fb->SetBinError(i,
      // TMath::Sqrt(charm_yield_100fb->GetBinContent(i)) );
      ccdis_20Rs2_yield_100fb->SetBinError(i, TMath::Sqrt(ccdis_20Rs2_yield_100fb->GetBinContent(i)));

      if (do_BkgSubtraction == kTRUE) {
        float err_charm = ccdis_20Rs2_yield_100fb->GetBinError(i);
        float err_light = light_yield_100fb->GetBinError(i);
        float err_total = TMath::Sqrt(ccdis_20Rs2_yield_100fb->GetBinContent(i) + light_yield_100fb->GetBinContent(i));

        float new_charm_err = TMath::Sqrt(TMath::Power(err_total, 2.0) + TMath::Power(err_light, 2.0));
        ccdis_20Rs2_yield_100fb->SetBinError(i, new_charm_err);
      }

    }

    // Generate the error band histogram for the suppressed case
    TH1F *error_band_100fb = static_cast<TH1F *>(ccdis_20Rs2_yield_100fb->Clone("error_band_100fb"));

    for (Int_t i = 1; i <= error_band_100fb->GetNbinsX(); i++) {
      error_band_100fb->SetBinError(i, error_band_100fb->GetBinError(i) / error_band_100fb->GetBinContent(i));
      error_band_100fb->SetBinContent(i, 1.0);
    }


    // Now handle the suppressed-only error band
    error_band_100fb->SetLineColor(kGray);
    error_band_100fb->SetFillColor(kGray);
    error_band_100fb->SetMarkerSize(0.0001);
    error_band_100fb->SetMarkerColor(kGray);
    error_band_100fb->SetTitle("#splitline{Stat. Uncertainty}{[CT18NNLO, R_{s}=2s/(#bar{u} + #bar{d})= 0.325 (suppressed)]}");

    htemplate = new TH1F("htemplate", "", 1, draw_config.xlimits[0], draw_config.xlimits[1]);
    htemplate->SetXTitle(draw_config.xtitle);
    htemplate->SetYTitle(draw_config.ytitle);
    htemplate->GetYaxis()->SetTitleSize(0.04);
    set_axis_range(htemplate, 0.0,                    2.0,                    "Y");
    set_axis_range(htemplate, draw_config.xlimits[0], draw_config.xlimits[1], "X");

    htemplate->Draw("HIST");
    error_band_100fb->Draw("E2 SAME");


    // Legend
    legend = smart_legend("lower left", 0.75, 0.28);
    legend->SetTextSize(0.045);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    //legend->AddEntry(error_band_100fb, "Stat. Uncertainty [CT18NNLO, R_{s}=2s/(#bar{u} + #bar{d})= 0.325 (suppressed)]", "lf");
    legend->AddEntry(error_band_100fb, "#splitline{Stat. Uncertainty}{[Rs-Low NNLO, suppressed strange]}", "lf");


    // Other alternative samples

    Int_t plot_index = 0;
    for (auto& alt_sample : alt_samples) {
      auto alt_tree = new TChain("tree");
      files = fileVector(Form("%s/%s/%s", dir.Data(), alt_sample.first.Data(), filePattern.Data()));

      if (files.size() == 0) {
	std::cout << "Sample " << alt_sample.first.Data() << " appears to contain no files - check the name of the directory, etc." << std::endl;
	exit(0);
      }

      for (auto file : files)
	{
	  // Check file first
	  TFile* testfile = TFile::Open(file.c_str());
	  if (testfile != nullptr) {
	    alt_tree->Add(file.c_str());
	    delete testfile;
	  } else {
	    std::cout << "File " << file.c_str() << " is damaged!" << std::endl;
	  }
	}

      std::cout << "Generated events in sample " << alt_sample.first.Data() << ": " << alt_tree->GetEntries() << std::endl;

      TTree* alt_selected = alt_tree->CopyTree(main_preselection.GetTitle());
      auto yield = DifferentialTaggingYield(alt_selected, draw_config, xvar, "charm");

      delete alt_selected;

      auto yield_100fb = static_cast<TH1F *>(yield->Clone(alt_sample.second));
      

      for (Int_t i = 1; i <= charm_yield_100fb->GetNbinsX(); i++) {
	yield_100fb->SetBinError(i, TMath::Sqrt(yield_100fb->GetBinContent(i)));
      }
      
      // Draw the enhanced-suppressed range overlay
      std::cout << "Draw " << alt_sample.second << std::endl;
      TH1F *range_band_100fb = static_cast<TH1F *>(ccdis_20Rs2_yield_100fb->Clone(Form("range plot %s", alt_sample.second.Data())));
      range_band_100fb->Scale(-1.0);
      range_band_100fb->Add(yield_100fb);
      range_band_100fb->Divide(ccdis_20Rs2_yield_100fb);
      
      for (Int_t i = 1; i <= range_band_100fb->GetNbinsX(); i++) {
	range_band_100fb->SetBinError(i, 0.0);
	range_band_100fb->SetBinContent(i, range_band_100fb->GetBinContent(i) + 1.0);
      }
      TGraphErrors *range_plot_100fb = new TGraphErrors(range_band_100fb);
      range_plot_100fb->SetName(Form("range_plot_100fb_%d", plot_index));
      range_plot_100fb->SetMarkerColor(alt_color[alt_sample.first]);
      range_plot_100fb->SetLineColor(alt_color[alt_sample.first]);
      range_plot_100fb->SetMarkerSize(2.0);
      range_plot_100fb->SetMarkerStyle(alt_marker[alt_sample.first]);
      range_plot_100fb->SetTitle(alt_sample.second.Data());

      range_plot_100fb->Draw("E1 SAME P");
      legend->AddEntry(range_plot_100fb, alt_sample.second.Data(),                              "lp");

      plot_index++;
    }


    // auto ccdis_21Rs2_data = new TChain("tree");
    // files = fileVector(Form("%s/CC_DIS_e10_p275_lha_21Rs2/%s", dir.Data(), filePattern.Data()));

    // for (auto file : files)
    // {
    //   ccdis_21Rs2_data->Add(file.c_str());
    // }

    // auto ccdis_CT18ANNLO_data = new TChain("tree");
    // files = fileVector(Form("%s/CC_DIS_e10_p275_CT18ANNLO/%s", dir.Data(), filePattern.Data()));

    // for (auto file : files)
    // {
    //   ccdis_CT18ANNLO_data->Add(file.c_str());
    // }


    // std::cout << "... enhanced strangeness ..." << std::endl;
    // TTree *ccdis_21Rs2_data_selected = ccdis_21Rs2_data->CopyTree(main_preselection.GetTitle());
    // std::cout << "... CT18ANNLO ..." << std::endl;
    // TTree *ccdis_CT18ANNLO_data_selected = ccdis_CT18ANNLO_data->CopyTree(main_preselection.GetTitle());

    // auto ccdis_21Rs2_yield = DifferentialTaggingYield(ccdis_21Rs2_data_selected, draw_config, xvar, "charm");
    // auto ccdis_CT18ANNLO_yield = DifferentialTaggingYield(ccdis_CT18ANNLO_data_selected, draw_config, xvar, "charm");

    // auto ccdis_21Rs2_yield_100fb = static_cast<TH1F *>(ccdis_21Rs2_yield->Clone("ccdis_21Rs2_yield_100fb"));
    // auto ccdis_CT18ANNLO_yield_100fb = static_cast<TH1F *>(ccdis_CT18ANNLO_yield->Clone("ccdis_CT18ANNLO_yield_100fb"));




    // // Draw the CT18ANNLO-suppressed range overlay
    // std::cout << "Draw the Rs range band for CT18ANNLO vs. Suppressed Strangeness..." << std::endl;
    // TH1F *range_CT18ANNLO_band_100fb = static_cast<TH1F *>(ccdis_20Rs2_yield_100fb->Clone("range_CT18ANNLO_band_100fb"));
    // range_CT18ANNLO_band_100fb->Scale(-1.0);
    // range_CT18ANNLO_band_100fb->Add(ccdis_CT18ANNLO_yield_100fb);
    // range_CT18ANNLO_band_100fb->Divide(ccdis_20Rs2_yield_100fb);

    // for (Int_t i = 1; i <= range_CT18ANNLO_band_100fb->GetNbinsX(); i++) {
    //   range_CT18ANNLO_band_100fb->SetBinError(i, 0.0);
    //   range_CT18ANNLO_band_100fb->SetBinContent(i, range_CT18ANNLO_band_100fb->GetBinContent(i) + 1.0);
    // }
    // TGraphErrors *range_CT18ANNLO_plot_100fb = new TGraphErrors(range_CT18ANNLO_band_100fb);
    // range_CT18ANNLO_plot_100fb->SetMarkerColor(kGreen - 5);
    // range_CT18ANNLO_plot_100fb->SetLineColor(kGreen - 5);
    // range_CT18ANNLO_plot_100fb->SetMarkerSize(2.0);
    // range_CT18ANNLO_plot_100fb->SetMarkerStyle(kOpenTriangleUp);
    // range_CT18ANNLO_plot_100fb->SetTitle("CT18ANNLO");


    // range_plot_100fb->Draw("E1 SAME P");
    // range_CT18ANNLO_plot_100fb->Draw("E1 SAME P");

    TLine OneLine(draw_config.xlimits[0], 1.0, draw_config.xlimits[1], 1.0);
    OneLine.SetLineWidth(2);
    OneLine.SetLineColor(kBlack);
    OneLine.Draw("SAME");

    TLatex plot_title = make_title();
    plot_title.Draw("SAME");


    pad->SetLogy(kFALSE);

    if (xvar.Contains("_x"))
      pad->SetLogx(kTRUE);

    pad->SetGrid(1, 1);

    legend->Draw();

    pad->Update();

    pad->SaveAs(Form("CharmTagPlot_tagging_yield_100fb_%s_%s.pdf", input.Data(), xvar.Data()));
    pad->SaveAs(Form("CharmTagPlot_tagging_yield_100fb_%s_%s.C", input.Data(), xvar.Data()));


    if (htemplate != nullptr)
      delete htemplate;
  }


  //
  // Polarized strangeness sensitivity estimate plot
  //

  if (xvar == "jet_pt") {
    draw_config.ylimits = std::vector<float>();
    draw_config.ylimits.push_back(-100);
    draw_config.ylimits.push_back(100);
    draw_config.xtitle = "Reconstructed Jet p_{T} [GeV]";
    draw_config.ytitle = "Asymmetry Uncertainty [%]";
    draw_config.logy   = kFALSE;
    draw_config.logx   = kFALSE;
  } else {
    do_Helicity = kFALSE;
  }

  if (do_Helicity == kTRUE) {
    pad->Clear();
    gStyle->SetErrorX(0.5);

    auto charm_yield = DifferentialTaggingYield(default_data_selected, draw_config, xvar, "charm");

    // auto light_yield = DifferentialTaggingYield(default_data_selected,
    // draw_config, xvar, "light");
    Float_t polarization_e = 0.70;
    Float_t polarization_p = 0.70;
    Float_t polarization   = 0.70; // # single beam polarization


    // Generate the error band histogram for the suppressed case
    TH1F   *dA_band_100fb = static_cast<TH1F *>(charm_yield->Clone("dA_band_100fb"));
    Float_t n_gen         = static_cast<Float_t>(default_data_selected->GetEntries());
    std::cout << "N_gen: " << n_gen << std::endl;

    for (Int_t i = 1; i <= dA_band_100fb->GetNbinsX(); i++) {
      Float_t n_events = dA_band_100fb->GetBinContent(i);
      Float_t err_n    = TMath::Sqrt(n_events * (1.0 + polarization_e) / 2.0);
      Float_t dA       = 0.0;

      if (err_n > 0.0) {
        dA = 1.0 / (err_n * polarization_p) * 100.0;
      }
      dA_band_100fb->SetBinError(i, dA);
      dA_band_100fb->SetBinContent(i, 0.0);
    }

    htemplate = new TH1F("htemplate", "", 1, draw_config.xlimits[0], draw_config.xlimits[1]);
    htemplate->SetXTitle(draw_config.xtitle);
    htemplate->SetYTitle(draw_config.ytitle);
    set_axis_range(htemplate, draw_config.ylimits[0], draw_config.ylimits[1], "Y");
    set_axis_range(htemplate, draw_config.xlimits[0], draw_config.xlimits[1], "X");

    dA_band_100fb->SetFillColor(kBlue + 2);
    dA_band_100fb->SetLineColor(kBlue + 2);
    dA_band_100fb->SetMarkerSize(0.001);
    dA_band_100fb->SetMarkerColor(kBlue + 2);
    dA_band_100fb->SetFillColorAlpha(kBlue + 0.2, 0.5);

    htemplate->Draw("HIST");
    dA_band_100fb->Draw("E2 SAME");

    TLine ZeroLine(draw_config.xlimits[0], 0.0, draw_config.xlimits[1], 0.0);
    ZeroLine.SetLineWidth(2);
    ZeroLine.SetLineColor(kBlack);
    ZeroLine.Draw("SAME");

    TLatex plot_title = make_title();
    plot_title.Draw("SAME");


    pad->SetLogy(draw_config.logx);
    pad->SetLogx(draw_config.logy);

    pad->SetGrid(1, 1);

    // Legend
    TLegend *legend = smart_legend("upper right", 0.25, 0.10);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->AddEntry(dA_band_100fb, "p=70%", "f");
    legend->Draw();

    pad->Update();

    pad->SaveAs(Form("CharmTagPlot_dA_strange_helicity_100fb_%s_%s.pdf", input.Data(), xvar.Data()));
    pad->SaveAs(Form("CharmTagPlot_dA_strange_helicity_100fb_%s_%s.C", input.Data(), xvar.Data()));
  }


  //
  // Generator-level plots
  //

  pad->Clear();

  bool do_GenJetPlot = kTRUE;

  if (xvar == "GenJet.PT") {
    float xbins[] = { 0, 2.5, 5.0, 7.5, 10, 12.5, 15, 20, 25, 35, 60 };
    int   nbins   = sizeof(xbins) / sizeof(xbins[0]) - 1;
    draw_config.htemplate = new TH1F(xvar, "", nbins, xbins);

    draw_config.xlimits = std::vector<float>();
    draw_config.xlimits.push_back(0);
    draw_config.xlimits.push_back(60);

    draw_config.ylimits = std::vector<float>();
    draw_config.ylimits.push_back(1);
    draw_config.ylimits.push_back(2e5);

    draw_config.logy = kTRUE;
    draw_config.logx = kFALSE;

    draw_config.xtitle = "Generated Jet p_{T} [GeV]";
    draw_config.ytitle = "d#sigma/dp_{T} #times 100fb^{-1} [GeV^{-1}]";
  } else {
    do_GenJetPlot = kFALSE;
  }


  if (do_GenJetPlot == kTRUE) {
    pad = new TCanvas("pad", "", 900, 1200);

    // Load Delphes samples for this plot
    auto delphes_data = new TChain("Delphes");
    files = fileVector(Form("%s/../%s/%s", dir.Data(), input.Data(), filePattern.Data()));

    for (auto file : files)
    {
      delphes_data->Add(file.c_str());
    }

    TH1F *light_genjets = static_cast<TH1F *>(draw_config.htemplate->Clone(Form("GenJet_light_%s_%d", xvar.Data(), getHistUID())));
    light_genjets->Sumw2();

    TH1F *charm_genjets = static_cast<TH1F *>(draw_config.htemplate->Clone(Form("GenJet_charm_%s_%d", xvar.Data(), getHistUID())));
    charm_genjets->Sumw2();

    delphes_data->Project(light_genjets->GetName(), xvar.Data(), "");
    delphes_data->Project(charm_genjets->GetName(), xvar.Data(), "GenJet.Flavor == 4");

    // Transform these into differential cross-section distributions
    Float_t n_gen = static_cast<float>(delphes_data->GetEntries());

    for (Int_t ibin = 1; ibin <= light_genjets->GetNbinsX(); ibin++) {
      Float_t dydx        = light_genjets->GetBinContent(ibin) * target_lumi * target_xsec / n_gen / light_genjets->GetBinWidth(ibin);
      Float_t dydx_relerr = light_genjets->GetBinError(ibin) / light_genjets->GetBinContent(ibin);
      light_genjets->SetBinContent(ibin, dydx);
      light_genjets->SetBinError(ibin, dydx * dydx_relerr);

      dydx        = charm_genjets->GetBinContent(ibin) * target_lumi * target_xsec / n_gen / charm_genjets->GetBinWidth(ibin);
      dydx_relerr = charm_genjets->GetBinError(ibin) / charm_genjets->GetBinContent(ibin);
      charm_genjets->SetBinContent(ibin, dydx);
      charm_genjets->SetBinError(ibin, dydx * dydx_relerr);
    }

    configure_plot<TH1F>(charm_genjets, draw_config, "charm");
    configure_plot<TH1F>(light_genjets, draw_config, "light");

    // generate template histogram
    htemplate = new TH1F("htemplate", "", 1, draw_config.xlimits[0], draw_config.xlimits[1]);

    set_axis_range(htemplate, draw_config.ylimits[0], draw_config.ylimits[1], "Y");
    set_axis_range(htemplate, draw_config.xlimits[0], draw_config.xlimits[1], "X");
    htemplate->SetXTitle(draw_config.xtitle);
    htemplate->SetYTitle(draw_config.ytitle);

    set_axis_range(charm_genjets, draw_config.ylimits[0], draw_config.ylimits[1], "Y");
    set_axis_range(charm_genjets, draw_config.xlimits[0], draw_config.xlimits[1], "X");
    charm_genjets->SetXTitle(draw_config.xtitle);
    charm_genjets->SetYTitle(draw_config.ytitle);

    set_axis_range(light_genjets, draw_config.ylimits[0], draw_config.ylimits[1], "Y");
    set_axis_range(light_genjets, draw_config.xlimits[0], draw_config.xlimits[1], "X");
    light_genjets->SetXTitle(draw_config.xtitle);
    light_genjets->SetYTitle(draw_config.ytitle);


    pad->cd();
    pad->SetLogy(draw_config.logy);
    pad->SetLogx(draw_config.logx);
    pad->SetGrid(1, 1);
    pad->Update();

    auto rp = new TRatioPlot(charm_genjets, light_genjets);

    rp->SetH1DrawOpt("E1 P");
    rp->SetH2DrawOpt("E1 P");

    rp->Draw();
    rp->GetLowerRefYaxis()->SetTitle("Charm-to-All Ratio");
    rp->GetLowerRefYaxis()->SetLabelSize(0.022);
    rp->GetLowerRefYaxis()->SetTitleSize(0.027);
    rp->GetLowerRefYaxis()->SetTitleOffset(2);
    rp->SetLeftMargin(0.15);
    rp->SetRightMargin(0.05);
    rp->SetLowBottomMargin(0.55);

    rp->GetLowerRefYaxis()->SetRangeUser(0.0, 0.05);
    rp->GetLowerRefYaxis()->SetNdivisions(5, 1, 0, kFALSE);
    rp->GetLowYaxis()->SetRangeUser(-0.01, 0.051);
    rp->GetLowYaxis()->SetNdivisions(5, 1, 0, kTRUE);

    TLatex plot_title = make_title();
    plot_title.Draw("SAME");


    // Legend
    rp->GetUpperPad()->cd();
    TLegend *legend = smart_legend("lower left", 0.55, 0.15);
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->AddEntry(light_genjets, "All Jets [CT18NNLO]",   "lp");
    legend->AddEntry(charm_genjets, "Charm Jets [CT18NNLO]", "lp");
    legend->Draw();

    rp->GetUpperPad()->SetGrid(1, 1);

    pad->Update();
    pad->SaveAs(Form("CharmTagPlot_differential_xs_%s_%s.pdf", delphes_data->GetTitle(), xvar.Data()));
    pad->SaveAs(Form("CharmTagPlot_differential_xs_%s_%s.C", delphes_data->GetTitle(), xvar.Data()));

    // Cleanup
    if (htemplate != nullptr)
      delete htemplate;

    if (legend != nullptr)
      delete legend;
  }
}
