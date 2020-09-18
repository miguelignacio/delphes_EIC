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

//
// Some universal constants for code to use
//

// Units, Cross-Sections
// standard baseline units: fb, GeV, mm, seconds

float u_fb = 1.0;
float u_mb = 1.0e12 * u_fb; // convert to fb

float u_GeV = 1.0;
float u_MeV = 1.0e3 * u_GeV; // convert to GeV

float u_mm = 1.0;
float u_cm = 1e2 * u_mm;


struct plot_config {
  TString xtitle = "";
  TString ytitle = "";
  TH1F* htemplate = nullptr;
  std::vector<float> xlimits;
  std::vector<float> ylimits;
  bool logy = kFALSE;
  bool logx = kFALSE;
  TString cuts = "";

  int markercolor = kBlack;
  int markerstyle = kDot;
  int linecolor   = kBlack;
  int linestyle   = kSolid;
  int linewidth   = 3;
  int fillcolor   = 0;
  int fillstyle   = 0;

};




static Int_t histUID = 0;

inline Int_t getHistUID()
{
  Int_t thisID = histUID;
  histUID++;
  return thisID;
}


inline float smart_legend_x(float x_trial, float x_width)
{
  float excess = 0.92 - (x_trial + x_width);
  float x_new = x_trial;

  if (excess < 0.0) {
    std::cout << "smart_legend_x(): x width overshoots range - shifting by " << excess << std::endl;
    x_new = x_trial + excess;
  }

  return x_new;

}

inline float smart_legend_y(float y_trial, float y_height)
{
  float excess = 0.92 - (y_trial + y_height);
  float y_new = y_trial;

  if (excess < 0.0) {
    std::cout << "smart_legend_y(): y height overshoots range - shifting by " << excess << std::endl;
    y_new = y_trial + excess;
  }

  return y_new;

}

inline TLegend* smart_legend(std::string where = "upper right", float legend_width = 0.30, float legend_height = 0.10)
{
  TLegend* legend = nullptr;

  float legend_x = 0.65;
  float legend_y = 0.85;

  if (where == "upper right") {
    legend_x = smart_legend_x(0.65, legend_width);
    legend_y = smart_legend_y(0.82, legend_height);
  } else if (where == "center right") {
    legend_x = smart_legend_x(0.65, legend_width);
    legend_y = smart_legend_y(0.55, legend_height);
  } else if (where == "lower right") {
    legend_x = smart_legend_x(0.65, legend_width);
    legend_y = smart_legend_y(0.20, legend_height);
  } else if (where == "upper left") {
    legend_x = smart_legend_x(0.15, legend_width);
    legend_y = smart_legend_y(0.82, legend_height);
  } else if (where == "center left") {
    legend_x = smart_legend_x(0.15, legend_width);
    legend_y = smart_legend_y(0.55, legend_height);
  } else if (where == "lower left") {
    legend_x = smart_legend_x(0.15, legend_width);
    legend_y = smart_legend_y(0.20, legend_height);
  } else {
    std::cout << "You specified a placement of " << where << " that is unknown..." << std::endl;

  }

  legend = new TLegend(legend_x, legend_y, legend_x + legend_width, legend_y + legend_height);
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);

  return legend;

}

inline void set_axis_range(TH1F* hist, float min, float max, std::string which = "X") {
  hist->SetAxisRange(min, max, which.c_str());
  std::cout << "set_axis_range: setting " << which << " axis range to: [" << min << ", " << max << "]" << std::endl;
  if (which == "X") {
    hist->GetXaxis()->SetLimits(min, max);
    hist->GetXaxis()->SetRangeUser(min, max);
  } else if (which == "Y") {
    hist->GetYaxis()->SetLimits(min, max);
    hist->GetYaxis()->SetRangeUser(min, max);
  }

}

template <class T> void configure_plot(T* object, plot_config options, std::string which = "" )
{

  if (which == "") {
    return;
  } else if (which == "charm") {
    object->SetLineColor(kBlue+1);
    object->SetMarkerColor(kBlue+1);
    object->SetMarkerStyle(kOpenDiamond);
  } else if (which == "light") {
    object->SetLineColor(kRed+1);
    object->SetMarkerColor(kRed+1);
    object->SetMarkerStyle(kFullCircle);
  } else if (which == "errorband") {

  }

  object->SetLineWidth(2);
  object->SetMarkerSize(2);

    // object->SetXTitle( options.xtitle );
    // object->SetYTitle( options.ytitle );
}


inline TLatex make_title(TString text = "CC-DIS, 10GeVx275GeV, Q^{2}>100 GeV^{2}") {

  TLatex plot_title;
  plot_title.SetTextSize(0.035);
  plot_title.SetTextAlign(22); // center-center
  plot_title.DrawLatexNDC(0.55, 0.97, text.Data());

  return plot_title;
}


inline std::vector<std::string> fileVector(const std::string& pattern){
  glob_t glob_result;
  glob(pattern.c_str(),GLOB_TILDE,NULL,&glob_result);
  std::vector<std::string> files;
  for(unsigned int i=0;i<glob_result.gl_pathc;++i){
    files.push_back(std::string(glob_result.gl_pathv[i]));
  }
  globfree(&glob_result);
  return files;
}


float LookupCrossSection(TString samplename)
{
  float xsection = 0.0;
  if (samplename.Contains("CC") && samplename.Contains("DIS")) {
    // Q^2 > 100 GeV^2
    xsection = 1.47637e-08*u_mb;
  } else if (samplename.Contains("NC") && samplename.Contains("DIS")) {
    // Q^2 > 50 GeV^2
    xsection = 7.444e-06*u_mb;
  }
  return xsection;
}



TH1F* GeneratePlot(plot_config draw_config,
		   TTree* data,
		   TString name,
		   TString x,
		   TCut selection)
{
  TH1F* plot = static_cast<TH1F*>(draw_config.htemplate->Clone(Form("%s_%d", name.Data(), getHistUID())));

  data->Project(plot->GetName(), x.Data(), selection);

  plot->SetLineColor(draw_config.linecolor);
  plot->SetLineStyle(draw_config.linestyle);

  plot->SetFillStyle(draw_config.fillstyle);
  plot->SetFillColor(draw_config.fillcolor);

  plot->SetMarkerColor(draw_config.markercolor);
  plot->SetMarkerStyle(draw_config.markerstyle);

  plot->SetXTitle( draw_config.xtitle );
  plot->SetYTitle( draw_config.ytitle );

  return plot;
}
