#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TString.h>
#include <TObjString.h>
#include "TInterpreter.h"

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

#include "ModuleHandler.h"
#include "TreeHandler.h"

static std::string input_dir = "";
static std::string output_file = "";
static std::string module_sequence = "";
static int nevents = -1;

ModuleHandler *ModuleHandler::instance = 0;
TreeHandler *TreeHandler::instance = 0;


// HELPER METHODS

void PrintHelp()
{
  std::cout <<
    "--input_dir <i>:       Directory containing all the ROOT files you want to process\n"
    "--output_file <o>:     Output ROOT file to store results\n"
    "--module_sequence <s>: A string comma-separated list of modules to load; order is preserved in execution.\n"
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
    "===================== SIMPLEANALYSIS =====================" << std::endl;


  // Handle complex TTree data storage types by defining them for ROOT
  gInterpreter->GenerateDictionary("std::vector<std::vector<float>>","vector");

	    

  if (argc <= 1) {
    PrintHelp();
  }

  const char* const short_opts = "i:o:h";
  const option long_opts[] = {
    {"input_dir", required_argument, nullptr, 'i'},
    {"output_file", required_argument, nullptr, 'o'},
    {"module_sequence", required_argument, nullptr, 's'},
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

        case 's':
	  module_sequence = optarg;
	  std::cout << "Module sequence: " << module_sequence << std::endl;
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
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
  TClonesArray *branchGenParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchRawTrack = treeReader->UseBranch("Track");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchMET = treeReader->UseBranch("MissingET");


  // Setup the module handler
  ModuleHandler *module_handler = module_handler->getInstance(treeReader);
  auto module_list = TString(module_sequence).Tokenize(",");

  // Setup the output storage
  TreeHandler *tree_handler = tree_handler->getInstance(output_file.c_str(), "tree");
  tree_handler->initialize();

  for (int i = 0; i < module_list->GetEntries(); i++) {
    auto name = static_cast<TObjString*>(module_list->At(i))->GetString().Data();
    std::cout << "   Appending module " << name << std::endl;
    module_handler->addModule(name);
  }

  for (auto module : module_handler->getModules()) {
    module->initialize();
  }

  if (nevents < 0) {
    std::cout
      << "Processing all events in the sample..." << std::endl;
  } else {
    std::cout
      << "Processing "<< nevents << " events in the sample..." << std::endl;
  }

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

    std::map<std::string, std::any> DataStore;


    for (auto module : module_handler->getModules()) {
      module->setJets(branchJet);
      module->setGenJets(branchGenJet);
      module->setEFlowTracks(branchEFlowTrack);
      module->setTracks(branchRawTrack);
      module->setGenParticles(branchGenParticle);
      module->setPhotons(branchPhoton);
      module->setElectrons(branchElectron);
      module->setNeutralHadrons(branchNeutralHadron);
      module->setMET(branchMET);

      bool result = module->execute(&DataStore);
      if (result == false) 
	break;
    }

    tree_handler->execute();
    

    // if (DataStore.find("CharmJets") != DataStore.end()) {
    //   std::vector<Jet*> charm_jets = std::any_cast<std::vector<Jet*>>(DataStore["CharmJets"]);
    // }

  }

  for (auto module : module_handler->getModules()) {
    module->finalize();
  }
  tree_handler->finalize();


  std::cout <<
    "========================== FINIS =========================" << std::endl;

  exit(EXIT_SUCCESS);  
}
