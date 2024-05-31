// std
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>
// ROOT
#include "TObjString.h"
#include "TChain.h"
#include "TFile.h"
// systematicstools
#include "systematicstools/interface/ISystProviderTool.hh"
#include "systematicstools/interface/SystMetaData.hh"
#include "systematicstools/interface/types.hh"
#include "systematicstools/utility/ParameterAndProviderConfigurationUtility.hh"
#include "systematicstools/utility/md5.hh"
#include "systematicstools/utility/printers.hh"
#include "systematicstools/utility/string_parsers.hh"
// sbnanaobj
#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/SRGlobal.h"
#include "sbnanaobj/StandardRecord/Flat/FlatRecord.h"

namespace cliopts {
  std::string input_filename = "";
  size_t NMax = std::numeric_limits<size_t>::max();
  size_t NSkip = 0;
} // namespace cliopts

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n" << std::endl;
  std::cout << "\t-?|--help        : Show this message.\n"
               "\t-i <ghep.root>   : GENIE TChain descriptor to read events\n"
               "\t                   from. (n.b. quote wildcards).\n"
               "\t-N <NMax>        : Maximum number of events to process.\n"
               "\t-s <NSkip>       : Number of events to skip.\n"
            << std::endl;
}

void HandleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-i") {
      cliopts::input_filename = argv[++opt];
    } else if (std::string(argv[opt]) == "-N") {
      cliopts::NMax = systtools::str2T<size_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-s") {
      cliopts::NSkip = systtools::str2T<size_t>(argv[++opt]);
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

int main(int argc, char const *argv[]) {

  HandleOpts(argc, argv);
  if (!cliopts::input_filename.size()) {
    std::cout << "[ERROR]: Expected to be passed a -i option." << std::endl;
    SayUsage(argv);
    return 1;
  }

  std::string fCAFTreeName = "recTree";
  std::string fSRName = "rec";
  std::string fGlobalTreeName = "globalTree";
  std::string fSRGlobalName = "global";

  // Input file
  TFile *f_input = new TFile(cliopts::input_filename.c_str());

  // Global tree
  TTree *t_input_global = (TTree *)f_input->Get(fGlobalTreeName.c_str());
  caf::SRGlobal* srglobal = nullptr;
  t_input_global->SetBranchAddress(fSRGlobalName.c_str(), &srglobal);
  t_input_global->GetEntry(0);

  std::cout << srglobal->wgts.size() << " parameter sets:" << std::endl;
  for(unsigned int i = 0; i < srglobal->wgts.size(); ++i){
    const caf::SRWeightPSet& pset = srglobal->wgts[i];
    std::cout << "  " << i << ": " << pset.name << ", type " << pset.type << ", " << pset.nuniv << " universes, adjusted parameters:" << std::endl;

    for(const caf::SRWeightMapEntry& entry: pset.map){
      std::cout << "    " << entry.param.name << std::endl;
      for(const auto& val: entry.vals){
        std::cout << val << " ";
      }
      std::cout << std::endl;
    }
  }

  // CAF tree
  TTree *fInputCAFTree = (TTree *)f_input->Get(fCAFTreeName.c_str());
  size_t ThisNCAFEvents = fInputCAFTree->GetEntries();
  size_t NToRead = std::min(ThisNCAFEvents, cliopts::NMax);
  printf("@@ Number of CAF events = %ld\n", ThisNCAFEvents);
  printf("@@ Number of CAF events to read = %ld\n", NToRead);
 
  const caf::CAFType caftype = caf::GetCAFType(fInputCAFTree);

  // - SRProxy to access record
  caf::SRSpillProxy* srproxy = new caf::SRSpillProxy(fInputCAFTree, fSRName.c_str());

  // Loop over CAFTree
  for (size_t cafev_it = 0; cafev_it < NToRead; ++cafev_it) {

    fInputCAFTree->GetEntry(cafev_it);

    printf("- %ld-th event\n", cafev_it);
    std::cout << "  * sourceName = " << srproxy->hdr.sourceName.GetValue() << std::endl;
    std::cout << "  * sourceNameHash = " << srproxy->hdr.sourceNameHash << std::endl;
    //printf("  * sourceName = %s\n",sourceName.c_str());
    //printf("  * sourceNameHash = %zd\n",sourceNameHash);
    printf("  * Neutrino loop *\n");
    int i_nu = 0;
    for(auto& nu: srproxy->mc.nu){
      printf("  - %d-th neutrino\n", i_nu);
      int i_w=0;
      for(auto& wgt: nu.wgt){
        printf("    - %d-th PSet (%s)\n", i_w, srglobal->wgts[i_w].name.c_str());
        for(int i_univ=0; i_univ<wgt.univ.size(); i_univ++){
          std::cout << "      - univ: " << i_univ << ", weight = " << wgt.univ[i_univ] << std::endl;
        }
        i_w++;
      }
      i_nu++;
    }

/*
    printf("  * Slice loop *\n");
    int i_slc=0;
    for(auto& slc: srproxy->slc){
      printf("  - %d-th slice\n", i_slc);
      if(slc.truth.index>=0){
        printf("    - Matched to neutrino\n");

        int i_w=0;
        for(auto& wgt: slc.truth.wgt){
          printf("      - %d-th PSet (%s)\n", i_w, srglobal->wgts[i_w].name.c_str());
          for(int i_univ=0; i_univ<wgt.univ.size(); i_univ++){
            std::cout << "        - univ: " << i_univ << ", weight = " << wgt.univ[i_univ] << std::endl;
          }
          i_w++;
        }

      }
      else{
        //printf("    - Not matched to neutrino\n");
      }
      i_slc++;
    }
*/


  }

}
