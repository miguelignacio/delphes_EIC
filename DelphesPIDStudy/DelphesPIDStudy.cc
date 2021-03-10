#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1.h>
#include <TString.h>
#include <TObjString.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "TInterpreter.h"
#include "TEfficiency.h"

#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <getopt.h>
#include <glob.h>
#include <vector>
#include <map>
#include <any>

#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"

static std::string input_dir = "";
static std::string output_file = "";
static int nevents = -1;


// HELPER METHODS

void PrintHelp()
{
  std::cout <<
    "--input_dir <i>:       Directory containing all the ROOT files you want to process\n"
    "--output_file <o>:     Output ROOT file to store results\n"
    "--nevents <n>:         The total number of events to process, starting from the zeroth event in the input.\n"
    "--help:                Show this helpful message!\n";
  exit(1);
}

std::vector<std::string> fileVector(const std::string& pattern){
  glob_t glob_result;
  glob(pattern.c_str(),GLOB_TILDE,NULL,&glob_result);
  std::vector<std::string> files;
  for(unsigned int i=0;i<glob_result.gl_pathc;++i){
    files.push_back(std::string(glob_result.gl_pathv[i]));
  }
  globfree(&glob_result);
  return files;
}

// MAIN FUNCTION


int main(int argc, char *argv[])
{
  std::cout <<
    "===================== DELPHESPIDSTUDY =====================" << std::endl;


  // Handle complex TTree data storage types by defining them for ROOT
  gInterpreter->GenerateDictionary("std::vector<std::vector<float>>","vector");

	    

  if (argc <= 1) {
    PrintHelp();
  }

  const char* const short_opts = "i:o:n:h";
  const option long_opts[] = {
    {"input_dir", required_argument, nullptr, 'i'},
    {"output_file", required_argument, nullptr, 'o'},
    {"nevents", optional_argument, nullptr, 'n'},
    {"help", no_argument, nullptr, 'h'},
    {nullptr, no_argument, nullptr, 0}
  };

  while (true)
    {
      const auto opt = getopt_long(argc, argv, short_opts, long_opts, nullptr);
      
      if (-1 == opt)
	break;
      
      switch (opt)
        {
        case 'i':
	  input_dir = optarg;
	  std::cout << "Input Directory: " << input_dir << std::endl;
	  break;

        case 'o':
	  output_file = optarg;
	  std::cout << "Output File: " << output_file << std::endl;
	  break;

        case 'n':
	  nevents = std::stoi(optarg);
	  std::cout << "Number of events to process: " << nevents << std::endl;
	  break;
	  

        case 'h': // -h or --help
        case '?': // Unrecognized option
	  PrintHelp();
	  break;
        default:
	  PrintHelp();
	  break;
        }
    }
  

  auto data = new TChain("Delphes");

  auto files = fileVector(input_dir);

  for (auto file : files) 
    {
      data->Add(file.c_str());
    }

  ExRootTreeReader *treeReader = new ExRootTreeReader(data);
  int n_entries = data->GetEntries();

  std::cout 
    << "The provided data set contains the following number of events: " << std::endl
    << n_entries
    << std::endl;

  // Load object pointers
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchPIDSystemsTrack = treeReader->UseBranch("PIDSystemsTrack");


  if (nevents < 0) {
    std::cout
      << "Processing all events in the sample..." << std::endl;
  } else {
    std::cout
      << "Processing "<< nevents << " events in the sample..." << std::endl;
  }


  // top-level map is truth ID, lower-level maps are PID designation pointing to TH1F
  std::map<int, std::map<int, TH1F*>*> PIDmap;

  std::vector<int> species;
  species.push_back(321);
  species.push_back(211);
  species.push_back(2212);
  species.push_back(11);

  for(int i=0; i < n_entries; ++i) {
    // event number printout
    if(i%1000==0) {
      std::cout << "Processing Event " << i << std::endl;
    }

    if (nevents >= 0 && i >= nevents)
      break;

    // read the data for i-th event
    // data->GetEntry(i);
    // Load selected branches with data from specified event
    treeReader->ReadEntry(i);


    // Put tracks in vectors for easier use
    std::vector<Track *> true_tracks;

    for (int itrk = 0; itrk < branchEFlowTrack->GetEntries(); itrk++) {
      auto track = (Track *)branchEFlowTrack->At(itrk);

      if (track->PT < 0.1 || TMath::Abs(track->Eta) > 3.5)
	continue;

      true_tracks.push_back(track);
    }

    std::vector<Track *> pid_tracks;

    for (int itrk = 0; itrk < branchPIDSystemsTrack->GetEntries(); itrk++) {
      auto track = (Track *)branchPIDSystemsTrack->At(itrk);

      if (track->PT < 0.1 || TMath::Abs(track->Eta) > 3.5)
	continue;

      pid_tracks.push_back(track);
    }


    // in each event, loop over the EFlowTracks and find all particles of a certain type.
    // then loop over the PID tracks and match the tracks. Use the ID in the PID system
    // to assign efficiency/misid

    for (auto specie : species) {
      if (PIDmap.find(specie) == PIDmap.end()) {
	PIDmap[specie] = new std::map<int, TH1F*>();
      }
      std::map<int, TH1F*>* specie_pid = PIDmap[specie];

      std::vector<Track*> true_specie_tracks;
      for (auto track : true_tracks) {
	if (TMath::Abs(track->PID) == specie)
	  true_specie_tracks.push_back(track);
      }
      
      for (auto track : true_specie_tracks) {
	for (auto pid_track : pid_tracks) {
	  //if (pid_track->P4().DeltaR( track->P4() ) < 1.0e-5) {
	  if (pid_track->Particle == track->Particle) {
	    // particle-level match!
	    int pid_id = TMath::Abs(pid_track->PID);

	    if (specie_pid->find(pid_id) == specie_pid->end()) {
	      (*specie_pid)[pid_id] = new TH1F(Form("h1_%d_%d",specie, pid_id), "", 100, 0, 50.0);
	      (*specie_pid)[pid_id]->Sumw2();

	    }
	    (*specie_pid)[pid_id]->Fill(pid_track->P4().Vect().Mag());
	  }
	}
      }
     

    }
    

  }

  TFile ResultsFile(output_file.c_str(),"RECREATE");
  ResultsFile.cd();


  for (auto specie : species) {
    std::map<int, TH1F*>* specie_pid = PIDmap[specie];

    // Now generate efficiency plots
    TH1F* h_all_species = nullptr;
    for ( auto[pid_id, histo] : *specie_pid) {
      std::cout << pid_id << std::endl;
      if (h_all_species == nullptr)
	h_all_species = static_cast<TH1F*>(histo->Clone("h_all_species"));
      else
	h_all_species->Add(histo);
    }
  
    for (auto[pid_id, histo] : *specie_pid) {
      TEfficiency eff_plot(*histo, *h_all_species);
      eff_plot.SetName(Form("eff_%d_%d", specie, pid_id));
      eff_plot.SetTitle(Form("Efficiency (%d #rightarrow %d) vs. Momentum;Track Momentum [GeV];Efficiency",specie, pid_id));
      eff_plot.Write();
    }
  }    
  
  ResultsFile.Close();


  std::cout <<
    "========================== FINIS =========================" << std::endl;

  exit(EXIT_SUCCESS);  
}
