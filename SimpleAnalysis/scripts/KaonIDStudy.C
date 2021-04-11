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
#include "TString.h"

#include <glob.h>
#include <iostream>
#include <iomanip>
#include <vector>

#include "PlotFunctions.h"

void ConfigHistogram(TH1 *h, std::string which = "charm")
{
  if (which == "charm") {
    h->SetLineColor(kBlue);
    h->SetMarkerColor(kBlue);
    h->SetMarkerStyle(kOpenSquare);
    h->SetMarkerColor(kBlue);
  } else if (which == "light") {
    h->SetLineColor(kRed);
    h->SetMarkerColor(kRed);
    h->SetMarkerStyle(kFullCrossX);
    h->SetMarkerColor(kRed);
  }
}

void KaonIDStudy(TString dir,
                 TString input,
                 TString trackname   = "barrelDircTrack",
                 TString filePattern = "*/out.root")
{
  auto data = new TChain("tree");

  data->SetTitle(input.Data());
  auto files = fileVector(Form("%s/%s/%s", dir.Data(), input.Data(), filePattern.Data()));

  for (auto file : files)
  {
    data->Add(file.c_str());
  }


  // TString trackname = "CF4RICHTrack";
  // TString trackname = "mRICHTrack";
  // TString trackname = "barrelDircTrack";
  // TString trackname = "tofBarrelTrack";
  // TString trackname = "dRICHagTrack";

  Float_t etaMin = -8;
  Float_t etaMax = 8;

  TCut   *cut_fiducial = nullptr;
  Float_t xmin         = 0.0;
  Float_t xmax         = 55.0;
  Bool_t logplot = kFALSE;

  if (trackname == "mRICHTrack") {
    // mRICH
    etaMin = -3.5;
    etaMax = -1.0;
    xmax   = 12.0;
  } else if (trackname == "barrelDIRCTrack") {
    // barrel DIRC
    etaMin = -1.0;
    etaMax = 1.0;
    xmax   = 55.0;
  } else if (trackname == "CF4RICHTrack") {
    // CF4RICH
    etaMin = 1.0;
    etaMax = 4.0;
  } else if (trackname == "tofBarrelTrack") {
    // TOF Barrel Detector
    etaMin = -2.0;
    etaMax = 2.0;
  } else if (trackname == "dualRICHTrack") {
    // dualRICH, aerogel-based
    etaMin = 1.00;
    etaMax = 3.50;
    
    // etaMin       = 1.48;
    // etaMax       = 3.91;
    xmax = 50.0;
    logplot = kTRUE;
    
  }
  cut_fiducial = new TCut(Form("%0.1f < pid_track_eta && pid_track_eta < %0.1f && pid_track_pt>0.1", etaMin, etaMax));

  Float_t xbinsize = 0.5;
  Int_t   nbinsx   = TMath::Ceil((xmax - xmin) / xbinsize);

  TCut                        cut_truekaon("TMath::Abs(pid_track_true_pid) == 321");
  TCut                        cut_recokaon("TMath::Abs(pid_track_reco_pid) == 321");
  TCut                        cut_truepion("TMath::Abs(pid_track_true_pid) == 211");
  TCut                        cut_recopion("TMath::Abs(pid_track_reco_pid) == 211");

  plot_config draw_config;
  draw_config.htemplate = new TH1F("TrackPT",
                                   "",
                                   nbinsx,
                                   xmin,
                                   xmax);

  // draw_config.xlimits   = std::vector<float>();
  // draw_config.xlimits.push_back(0);
  // draw_config.xlimits.push_back(20);
  draw_config.ylimits = std::vector < float > ();
  draw_config.ylimits.push_back(0.0);
  draw_config.ylimits.push_back(1e4);
  draw_config.xtitle = "Track Momentum [GeV]";
  draw_config.ytitle = "Frequency";
  draw_config.logy   = kFALSE;
  draw_config.logx   = kFALSE;
  draw_config.cuts   = "";

  // Kaon -> Kaon ID
  auto true_kaon_pt = GeneratePlot(draw_config, data, "True Kaon PT", "pid_track_pt*TMath::CosH(pid_track_eta)", *cut_fiducial && cut_truekaon);
  auto reco_kaon_pt = GeneratePlot(draw_config, data, "Reconstructed Kaon PT", "pid_track_pt*TMath::CosH(pid_track_eta)", *cut_fiducial && cut_truekaon && cut_recokaon);


  TH1F *h_template = static_cast < TH1F * > (reco_kaon_pt->Clone("h_template"));
  h_template->SetAxisRange(0.0, 1.12, "Y");
  if (logplot)
    h_template->SetAxisRange(1.0e-5, 10.0, "Y");

  h_template->GetYaxis()->SetTitle("Efficiency");

  auto pid_efficiency = new TEfficiency(*reco_kaon_pt,
                                        *true_kaon_pt);
  pid_efficiency->SetMarkerColor(kBlue);
  pid_efficiency->SetFillColor(kBlue);
  pid_efficiency->SetLineColor(kBlue);
  pid_efficiency->SetMarkerStyle(kOpenSquare);

  // Pion -> Kaon ID
  auto true_pion_pt = GeneratePlot(draw_config, data, "True Pion PT", "pid_track_pt*TMath::CosH(pid_track_eta)", *cut_fiducial && cut_truepion);
  auto reco_pion_pt = GeneratePlot(draw_config, data, "Reconstructed Pion PT", "pid_track_pt*TMath::CosH(pid_track_eta)", *cut_fiducial && cut_truepion && cut_recokaon);


  auto misid_efficiency = new TEfficiency(*reco_pion_pt,
                                          *true_pion_pt);
  misid_efficiency->SetMarkerColor(kRed);
  misid_efficiency->SetFillColor(kRed);
  misid_efficiency->SetLineColor(kRed);
  misid_efficiency->SetMarkerStyle(kFullCrossX);

  // Draw!
  TCanvas c1("c1",
             "",
             800,
             600);

  c1.cd();
  c1.SetLogy(logplot);

  // true_kaon_pt->Draw("HIST");
  // reco_kaon_pt->Draw("E1 SAME");
  h_template->Draw("AXIS");
  pid_efficiency->Draw("E1 SAME");
  misid_efficiency->Draw("E1 SAME");

  TLatex etaRange((xmin + (xmax - xmin) * 0.67),
                  1.05,
                  Form("%0.1f < #eta < %0.1f", etaMin, etaMax));
  etaRange.Draw();

  // legend
  auto *legend = smart_legend("lower left");
  legend->AddEntry(pid_efficiency,   "K^{#pm} #rightarrow K^{#pm}",   "lp");
  legend->AddEntry(misid_efficiency, "#pi^{#pm} #rightarrow K^{#pm}", "lp");
  legend->Draw();

  c1.SaveAs(Form("%s_KaonIDStudy.pdf", trackname.Data()));


  return;

  // Plot the pT spectrum for true and reconstructed candidates from charm and
  // light jets

  // Kaon -> Kaon ID
  TCut               true_charm_mother("pid_track_jetmother == 4");
  TCut               true_light_mother("pid_track_jetmother > 0 && pid_track_jetmother < 4");

  if (cut_fiducial != nullptr) delete cut_fiducial;
  cut_fiducial = new TCut(Form("-4.0 < pid_track_eta && pid_track_eta < 4.0 && pid_track_pt>0.1 && pid_track_jet_pt > 10.0 && TMath::Abs(pid_track_jet_eta)<3.0"));


  auto charm_true_kaon_pt = GeneratePlot(draw_config, data, "True Kaon PT", "pid_track_pt*TMath::CosH(pid_track_eta)", *cut_fiducial && true_charm_mother && cut_truekaon);
  auto light_true_kaon_pt = GeneratePlot(draw_config, data, "True Kaon PT", "pid_track_pt*TMath::CosH(pid_track_eta)", *cut_fiducial && true_light_mother && cut_truekaon);

  c1.Clear();

  charm_true_kaon_pt->SetYTitle("Probability");

  ConfigHistogram(charm_true_kaon_pt, "charm");
  ConfigHistogram(light_true_kaon_pt, "light");

  charm_true_kaon_pt->Rebin(2);
  light_true_kaon_pt->Rebin(2);

  TH1 *charm_true_kaon = charm_true_kaon_pt->DrawNormalized("E1");
  TH1 *light_true_kaon = light_true_kaon_pt->DrawNormalized("E1 SAME");

  charm_true_kaon->SetAxisRange(1e-5, 1.0, "Y");

  c1.SetLogy();
  c1.SaveAs(Form("True_Kaon_Momentum.pdf"));


  auto charm_reco_kaon_pt = GeneratePlot(draw_config, data, "reco Kaon PT", "pid_track_pt*TMath::CosH(pid_track_eta)", *cut_fiducial && true_charm_mother && cut_truekaon && cut_recokaon);
  auto light_reco_kaon_pt = GeneratePlot(draw_config, data, "True Kaon PT", "pid_track_pt*TMath::CosH(pid_track_eta)", *cut_fiducial && true_light_mother && cut_truekaon && cut_recokaon);

  c1.Clear();

  charm_reco_kaon_pt->SetYTitle("Probability");

  ConfigHistogram(charm_reco_kaon_pt, "charm");
  ConfigHistogram(light_reco_kaon_pt, "light");

  charm_reco_kaon_pt->Rebin(2);
  light_reco_kaon_pt->Rebin(2);

  TH1 *charm_reco_kaon = charm_reco_kaon_pt->DrawNormalized("E1");
  TH1 *light_reco_kaon = light_reco_kaon_pt->DrawNormalized("E1 SAME");

  charm_reco_kaon->SetAxisRange(1e-5, 1.0, "Y");

  c1.SetLogy();
  c1.SaveAs(Form("Reco_Kaon_Momentum.pdf"));

  //
  // Zoom in on the low-momentum region, where most of the true kaons in CC DIS
  // live.
  //

  plot_config kaon_low_p_config;
  kaon_low_p_config.htemplate = new TH1F("TrackMomentum_Low",
                                         "",
                                         100,
                                         0.0,
                                         10.0);

  // kaon_low_p_config.xlimits = std::vector < float > ();
  // kaon_low_p_config.xlimits.push_back(0);
  // kaon_low_p_config.xlimits.push_back(10);
  kaon_low_p_config.ylimits = std::vector < float > ();
  kaon_low_p_config.ylimits.push_back(0.0);
  kaon_low_p_config.ylimits.push_back(1.0e6);
  kaon_low_p_config.xtitle = "Track Momentum [GeV]";
  kaon_low_p_config.ytitle = "Frequency";
  kaon_low_p_config.logy   = kFALSE;
  kaon_low_p_config.logx   = kFALSE;
  kaon_low_p_config.cuts   = "";

  // auto reco_kaon_pt = GeneratePlot(draw_config, data, "Reconstructed Kaon
  // PT", "pid_track_pt*TMath::CosH(pid_track_eta)", *cut_fiducial &&
  // cut_truekaon && cut_recokaon);

  auto charm_true_kaon_p_low = GeneratePlot(kaon_low_p_config, data, "Charm True Kaon P", "pid_track_pt*TMath::CosH(pid_track_eta)", *cut_fiducial && true_charm_mother && cut_truekaon);
  auto light_true_kaon_p_low = GeneratePlot(kaon_low_p_config, data, "Light True Kaon P", "pid_track_pt*TMath::CosH(pid_track_eta)", *cut_fiducial && true_light_mother && cut_truekaon);

  std::cout << charm_true_kaon_p_low->GetNbinsX() << std::endl;
  std::cout << light_true_kaon_p_low->GetNbinsX() << std::endl;

  c1.Clear();
  c1.SetLogy(kaon_low_p_config.logy);
  c1.SetLogx(kaon_low_p_config.logx);

  charm_true_kaon_p_low->SetYTitle("Probability");

  ConfigHistogram(charm_true_kaon_p_low, "charm");
  ConfigHistogram(light_true_kaon_p_low, "light");

  charm_true_kaon = charm_true_kaon_p_low->DrawNormalized("E1");
  light_true_kaon = light_true_kaon_p_low->DrawNormalized("E1 SAME");

  charm_true_kaon->SetAxisRange(0.0, 0.04, "Y");

  c1.SaveAs(Form("True_Kaon_Momentum_Low.pdf"));

  // Zoomed, reconstructed momentum
  auto charm_reco_kaon_p_low = GeneratePlot(kaon_low_p_config, data, "Charm True Kaon P", "pid_track_pt*TMath::CosH(pid_track_eta)", *cut_fiducial && true_charm_mother && cut_truekaon && cut_recokaon);
  auto light_reco_kaon_p_low = GeneratePlot(kaon_low_p_config, data, "Light True Kaon P", "pid_track_pt*TMath::CosH(pid_track_eta)", *cut_fiducial && true_light_mother && cut_truekaon && cut_recokaon);

  std::cout << charm_reco_kaon_p_low->GetNbinsX() << std::endl;
  std::cout << light_reco_kaon_p_low->GetNbinsX() << std::endl;

  c1.Clear();
  c1.SetLogy(kaon_low_p_config.logy);
  c1.SetLogx(kaon_low_p_config.logx);

  charm_reco_kaon_p_low->SetYTitle("Probability");

  ConfigHistogram(charm_reco_kaon_p_low, "charm");
  ConfigHistogram(light_reco_kaon_p_low, "light");

  charm_reco_kaon = charm_reco_kaon_p_low->DrawNormalized("E1");
  light_reco_kaon = light_reco_kaon_p_low->DrawNormalized("E1 SAME");

  charm_reco_kaon->SetAxisRange(0.0, 0.04, "Y");

  c1.SaveAs(Form("Reco_Kaon_Momentum_Low.pdf"));
}
