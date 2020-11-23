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

void              KaonIDStudy(TString dir, TString input, TString trackname = "barrelDircTrack", TString filePattern = "*/out.root")
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

  if (trackname == "mRICHTrack") {
    // mRICH
    etaMin       = -4.0;
    etaMax       = -1.0;
    xmax         = 12.0;
  } else if (trackname == "barrelDircTrack") {
    // barrel DIRC
    etaMin       = -1.0;
    etaMax       = 1.0;
    xmax         = 55.0;
  } else if (trackname == "CF4RICHTrack") {
    // CF4RICH
    etaMin       = 1.0;
    etaMax       = 4.0;
  } else if (trackname == "tofBarrelTrack") {
    // TOF Barrel Detector
    etaMin       = -2.0;
    etaMax       = 2.0;
  } else if (trackname == "dRICHTrack") {
    // dualRICH, aerogel-based
    etaMin       = 1.48;
    etaMax       = 3.91;
    xmax         = 75.0;
  }
  cut_fiducial = new TCut(Form("%0.1f < pid_track_eta && pid_track_eta < %0.1f && pid_track_pt>0.1", etaMin, etaMax));

  Float_t xbinsize = 0.5;
  Int_t   nbinsx   = TMath::Ceil((xmax - xmin) / xbinsize);

  TCut                        cut_truekaon("((pid_track_pid & 0xffff0000) >> 16) == 321");
  TCut                        cut_recokaon("(pid_track_pid & 0xffff) == 321");
  TCut                        cut_truepion("((pid_track_pid & 0xffff0000) >> 16) == 211");
  TCut                        cut_recopion("(pid_track_pid & 0xffff) == 211");

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
}
