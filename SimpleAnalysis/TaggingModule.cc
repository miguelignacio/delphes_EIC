#include "TaggingModule.h"

#include "TClonesArray.h"

#include "AnalysisFunctions.cc"

#include <iostream>
#include <iomanip>  
#include <fstream>

#include "TreeHandler.h"

TaggingModule::TaggingModule(ExRootTreeReader* data)
  : Module(data)
{

}

TaggingModule::~TaggingModule()
{

}

void TaggingModule::initialize()
{

}

void TaggingModule::finalize()
{

}
bool TaggingModule::execute(std::map<std::string, std::any>* DataStore)
{
  auto data = getData();

  // Construct the general vector of jets
  std::vector<Jet*> all_jets;
  for (int ijet = 0; ijet < getJets()->GetEntries(); ijet++) 
    {
      Jet *jet = (Jet*) getJets()->At(ijet);
      all_jets.push_back(jet);
    }

  // select jets well in the fiducial region
  std::vector<Jet*> fiducial_jets = SelectorFcn<Jet>(all_jets, [](Jet* j){ return (TMath::Abs(j->Eta) < 3.0 && j->PT > 5.0); });


  // Retrieve the general tracks list
  auto tracks = getEFlowTracks();

  // Produce a list of tagged jets
  std::vector<Jet*> charmtagged_jets;
  for (auto jet : fiducial_jets) {
    if (Tagged_sIP3D(jet, *tracks, 3.75, 1.00, 2.0) == true) {
      charmtagged_jets.push_back(jet);
    }
  }
  (*DataStore)["CharmJets"] = charmtagged_jets;


  


  return true;
}

