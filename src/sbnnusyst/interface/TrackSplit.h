#pragma once

#include <iostream>
#include <string>
// ROOT
#include "TH2.h"
#include "TSpline.h"
#include "TProfile.h"
// sbnanaobj
#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"
#include "sbnanaobj/StandardRecord/StandardRecord.h"
// sbnnusyst
#include "sbnnusyst/utility/Utilities.h"

namespace sbnnusyst{

class TrackSplit{

public:

  TrackSplit();
  ~TrackSplit();

  void Init(std::string tracksplit_filename, std::string chi2_filename);

  bool DoDebug;

  struct TrkChi2Results { ///< determined particle ID
    Float_t chi2_kaon, chi2_muon, chi2_pion, chi2_proton, pida;
    Int_t pid_ndof;
  };
  struct TrkMomentumResults { ///< range-based momentum results
    Float_t p_muon, p_pion, p_proton;
  };
  struct PreservedInitialTrkResults { ///< stuff from the unchanged initial track that I want to copy to the secondary track
    Float_t dir_end_x, dir_end_y, dir_end_z, end_x, end_y, end_z;
    unsigned producer;
  };

  // split probabilities
  TH2D *splitProb[2][2]; // In ICARUS, idx1: 0 = East and 1 = West, idx2: 0 = cathode, 1 = z=0 point
  void ReadTrackSplitMap(std::string filepath);

  void ApplyTrackSplit(caf::SRSlice& slc);
  caf::SRCaloPoint FillCaloPointFrom( const caf::SRCaloPoint& inCaloPt );
  caf::SRPFP FillPtrPFP( const caf::SRPFP& inPfp );

  // chi2
  TSpline3 KEvsR_spline3;

  TProfile *dedx_range_pro;   ///< proton template
  TProfile *dedx_range_ka;    ///< kaon template
  TProfile *dedx_range_pi;    ///< pion template
  TProfile *dedx_range_mu;    ///< muon template
  TrkChi2Results CalculateChi2(const caf::SRTrackCalo& calo);
  TrkMomentumResults CalculateMomenta(const float length);

};

} // END namespace sbnnusyst
