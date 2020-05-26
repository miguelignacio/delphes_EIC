#include "CharmJetModule.h"

#include "TClonesArray.h"

#include "classes/DelphesClasses.h"

#include "AnalysisFunctions.cc"

#include <iostream>

CharmJetModule::CharmJetModule(ExRootTreeReader* data)
  : Module(data)
{

}

CharmJetModule::~CharmJetModule()
{

}



bool CharmJetModule::execute(std::map<std::string, std::any>* DataStore)
{
  auto data = getData();

  // If event contains at least 1 jet
  std::vector<Jet*> all_jets;
  for (int ijet = 0; ijet < getJets()->GetEntries(); ijet++) 
    {
      // Take first jet
      Jet *jet = (Jet*) getJets()->At(ijet);
      all_jets.push_back(jet);
    }

  std::vector<Jet*> fiducial_jets = SelectorFcn<Jet>(all_jets, [](Jet* j){ return (TMath::Abs(j->Eta) < 3.0 && j->PT > 5.0); });
  std::vector<Jet*> charmJets = SelectorFcn<Jet>(fiducial_jets, [](Jet* j){ return (j->Flavor == 4); });


  (*DataStore)["CharmJets"] = charmJets;

  // for (int igjet = 0; igjet < getGenJets()->GetEntries(); igjet++) 
  // 	{
  // 	  Jet *tru_jet = (Jet*) getGenJets()->At(igjet);
  // 	  // if (jet->P4().DeltaR(tru_jet->P4()) < 0.5 &) {
	    
  // 	  // }
  // 	}

  // auto jet_constituents = jet->Constituents;

  // for (int iconst = 0; iconst < jet_constituents.GetEntries(); iconst++) {
  // 	auto constituent = jet_constituents.At(iconst);

  // 	if(constituent == 0) continue;

  // 	if (constituent->IsA() == Track::Class()) {
  // 	  auto track = static_cast<Track*>(constituent);
  // 	}
  // }
      


  return true;
}

