#pragma once

#include <iostream>
#include <string>
// ROOT
#include "TChain.h"
#include "TBranch.h"
// GENIE
#include "Framework/EventGen/EventRecord.h"
#include "Framework/GHEP/GHepParticle.h"
#include "Framework/GHEP/GHepUtils.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/Ntuple/NtpMCEventRecord.h"
// sbnanaobj
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/SRGlobal.h"
#include "sbnanaobj/StandardRecord/Flat/FlatRecord.h"
// nusystematics
#include "nusystematics/utility/response_helper.hh"
// sbnnusyst
#include "sbnnusyst/utility/Utilities.h"

namespace sbnnusyst{

class WeightUpdater{

public:

  WeightUpdater(
    std::string caftreename,
    std::string srname,
    std::string globaltreename,
    std::string srglobalname,
    std::string genietreename,
    std::string genierecname
  );
  ~WeightUpdater();

  nusyst::response_helper* fRH;
  void SetResponseHelper(std::string fclname);

  std::string fCAFTreeName;
  std::string fSRName;
  std::string fGlobalTreeName;
  std::string fSRGlobalName;
  std::string fGENIETreeName;
  std::string fGENIERecName;
  size_t NProcessedCAFEvents;
  size_t NMaxCAFEventsToProcess;
  void SetNMaxCAFEventsToProcess(size_t nmax);
  size_t GlobalGENIEEventCounter;
  void ProcessFile(std::string inputfile);
  size_t NProcessedFiles;

  TFile *fOutputFile;
  TTree *fOutputCAFTree;
  TTree *fOutputGENIETree;
  TTree *fOutputGlobalTree;
  caf::FlatStandardRecord* fOutputFlatSR;
  genie::NtpMCEventRecord *fOutputGENIENtp;

  void SetOutputFileName(std::string FileName);
  void CreateMetadataTree(); // TODO
  void CreateGlobalTree(caf::SRGlobal* input_srglobal);

  void SetOutputPOTHistName(std::string name);
  void SetOutputLivetimeHistName(std::string name);
  std::string fPOTHistName;
  std::string fLivetimeHistName;
  TH1D *fOutputPOT;
  TH1D *fOutputLivetime;
  bool AddPOTHist(TH1D *h_input);
  bool AddLivetimeHist(TH1D *h_input);

  void Save();

  bool DoDebug;

};

} // END namespace sbnnusyst
