#ifndef TAGGINGSTUDYMODULE_HH
#define TAGGINGSTUDYMODULE_HH

#include "classes/DelphesClasses.h"

#include "Module.h"

class TaggingStudyModule : public Module {

 public:

  TaggingStudyModule(ExRootTreeReader* data);

  ~TaggingStudyModule();

  void initialize() override;
  bool execute(std::map<std::string, std::any>* DataStore) override;
  void finalize() override;


 private:
  // errd0 and errz0 are in millimeters
  // bool Tagged(Jet* jet, float minSignif=2.2, float minPT = 1.0, int minTracks = 3, float errd0 = -1, float errz0 = -1);

  // Branch variables for storage to disk
  float _jet_pt;
  float _jet_eta;
  float _jet_flavor;
  float _jet_tagged;
  float _jet_btag;


  // study map linking variations to numbers and kinds of jet

  std::map<std::string, std::map<std::string, int>> study_variations;
  std::vector<int> minTrackVar{2, 3, 4};
  std::vector<float> minSignifVar{1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.00, 4.25, 4.50, 4.75, 5.00, 6.00, 7.00, 8.00, 9.00,10.00,12.00,14.00, 16.00};
  std::vector<float> minTrkPTVar{0.10, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5};

};

#endif
