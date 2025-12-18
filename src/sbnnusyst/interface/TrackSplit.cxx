#include "WeightUpdater.h"
#include "TROOT.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TGraph.h"

namespace sbnnusyst{

TrackSplit::TrackSplit(){

  DoDebug = false;

}

TrackSplit::~TrackSplit(){

}

void TrackSplit::Init(std::string tracksplit_filename, std::string chi2_filename){

  // Track split 

  std::cout << "[TrackSplit::Init] Track split probability from " << tracksplit_filename.c_str() << std::endl;

  TFile* file_tracksplitdata = TFile::Open(tracksplit_filename.c_str());

  splitProb[0][0] = (TH2D*)file_tracksplitdata->Get("z_cathodeangle_east");
  splitProb[0][1] = (TH2D*)file_tracksplitdata->Get("x_zgapangle_east");
  splitProb[1][0] = (TH2D*)file_tracksplitdata->Get("z_cathodeangle_west");
  splitProb[1][1] = (TH2D*)file_tracksplitdata->Get("x_zgapangle_west");

  // Get the dEdx templates as in Calo Syst

  std::cout << "[TrackSplit::Init] Chi2-pid template from " << chi2_filename.c_str() << std::endl;

  TFile *file_Chi2Template = TFile::Open(chi2_filename.c_str());
  dedx_range_pro = (TProfile*)file_Chi2Template->Get("dedx_range_pro");
  dedx_range_ka  = (TProfile*)file_Chi2Template->Get("dedx_range_ka");
  dedx_range_pi  = (TProfile*)file_Chi2Template->Get("dedx_range_pi");
  dedx_range_mu  = (TProfile*)file_Chi2Template->Get("dedx_range_mu");

  // Set up the spline used for momentum calculation, as in
  // https://github.com/LArSoft/larreco/blob/LARSOFT_SUITE_v09_72_00/larreco/RecoAlg/TrackMomentumCalculator.cxx
  std::array<float, 29> Range_grampercm{
    {9.833E-1/1.396, 1.786E0/1.396, 3.321E0/1.396, 6.598E0/1.396, 1.058E1/1.396, 3.084E1/1.396, 4.250E1/1.396, 6.732E1/1.396, 1.063E2/1.396, 1.725E2/1.396,
     2.385E2/1.396,  4.934E2/1.396, 6.163E2/1.396, 8.552E2/1.396, 1.202E3/1.396, 1.758E3/1.396, 2.297E3/1.396, 4.359E3/1.396, 5.354E3/1.396, 7.298E3/1.396,
     1.013E4/1.396,  1.469E4/1.396, 1.910E4/1.396, 3.558E4/1.396, 4.326E4/1.396, 5.768E4/1.396, 7.734E4/1.396, 1.060E5/1.396, 1.307E5/1.396}};

  constexpr std::array<float, 29> KE_MeV{
    {10,    14,    20,    30,    40,     80,     100,    140,    200,   300,
    400,   800,   1000,  1400,  2000,   3000,   4000,   8000,   10000, 14000,
    20000, 30000, 40000, 80000, 100000, 140000, 200000, 300000, 400000}};
  TGraph const KEvsR{29, Range_grampercm.data(), KE_MeV.data()};
  KEvsR_spline3 = TSpline3("KEvsRS", &KEvsR);

  std::cout << "[TrackSplit::Init] Done" << std::endl;

}

void TrackSplit::ApplyTrackSplit(caf::SRSlice& slc){

  std::map< unsigned int, std::vector<std::vector<caf::SRCaloPoint>> > scdryInitIdxToPfpCaloPts;
  std::map< unsigned int, PreservedInitialTrkResults > scdryInitIdxToPreservedInfo;
  unsigned int expectedScdryTrks = 0;

  if(DoDebug){
    std::cout << "[TrackSplit::ApplyTrackSplit] Called, slc.reco.pfp.size() = " << slc.reco.pfp.size() << std::endl;
  }

  int pfpIdx = -1;
  // In this first pfp loop, we will "update" the crossing track if it splits
  // We'll take care of the secondary track later
  for( auto& pfp : slc.reco.pfp ) {

    pfpIdx+=1;

    bool crossesCath = false;
    bool crossesMidZ = false;
    if ( std::isnan(pfp.trk.len) || std::isnan(pfp.trk.start.x) ||
         std::isnan(pfp.trk.start.z) || std::isnan(pfp.trk.end.x) || std::isnan(pfp.trk.end.z) ) {
      continue;
    }
    if ( ((fabs(pfp.trk.start.x)-210.2)*(fabs(pfp.trk.end.x)-210.2)) < 0. )
      crossesCath = true;
    if ( (pfp.trk.start.z * pfp.trk.end.z) < 0. )
      crossesMidZ = true;

    if ( crossesCath && crossesMidZ )
      crossesMidZ = false;

    if ( !crossesCath && !crossesMidZ ) {
      continue;
    }

    // Should we split the track?
    double pToSplit = 0.;
    TVector3 trkAvgDir( pfp.trk.end.x - pfp.trk.start.x, pfp.trk.end.y - pfp.trk.start.y, pfp.trk.end.z - pfp.trk.start.z );
    double nSlopes = 0.;
    if ( crossesCath ) {
      bool isWest = pfp.trk.start.x > 0;
      nSlopes = (isWest ? (210.2-pfp.trk.start.x)/trkAvgDir.X() : (pfp.trk.start.x+210.2)/trkAvgDir.X() );
      double zAtCathode = pfp.trk.start.z + nSlopes*trkAvgDir.Z();
      TVector3 vNormPlane( 1., 0., 0. );
      double ThCathode = trkAvgDir.Angle(vNormPlane);
      ThCathode = std::min(ThCathode, TMath::Pi()-ThCathode);
      double cosThCathode = TMath::Cos( ThCathode );

      const int xBin = splitProb[1][0]->GetXaxis()->FindBin(zAtCathode);
      const int yBin = splitProb[1][0]->GetYaxis()->FindBin(cosThCathode);
      if (xBin == 0 || xBin == splitProb[1][0]->GetXaxis()->GetNbins() + 1 ||
          yBin == 0 || yBin == splitProb[1][0]->GetYaxis()->GetNbins() + 1) {
        if(DoDebug){
          std::cout << "Note that z = " << zAtCathode << " and cosTh = " << cosThCathode << " for cathode." << std::endl;
          std::cout << "BINS: x=" << xBin << ", y=" << yBin << std::endl;
          std::cout << "   start: (" << pfp.trk.start.x << ", " << pfp.trk.start.y << ", " << pfp.trk.start.z << ")" << std::endl;
          std::cout << "   end: (" << pfp.trk.end.x << ", " << pfp.trk.end.y << ", " << pfp.trk.end.z << ")" << std::endl;
        }
        pToSplit = 0;
      }
      else{
/*
        // Using histogram
        pToSplit = (isWest ? splitProb[1][0]->GetBinContent(xBin,yBin) :
                       splitProb[0][0]->GetBinContent(xBin,yBin) );
*/
        // Using fixed number from BNB
        pToSplit = (isWest ? 0.02 : 0.06);
      }
    } // END if crossesCath
    else {
      bool isWest = pfp.trk.start.x > 0;
      bool isNorth = pfp.trk.start.z > 0;
      nSlopes = (!isNorth ? (0.-pfp.trk.start.z)/trkAvgDir.Z() : (pfp.trk.start.z)/trkAvgDir.Z() );
      double xAtMiddleZ = pfp.trk.start.x + nSlopes*trkAvgDir.X();
      TVector3 vNormPlane( 0., 0., 1. );
      double ThMiddleZ = trkAvgDir.Angle(vNormPlane);
      ThMiddleZ = std::min(ThMiddleZ, TMath::Pi()-ThMiddleZ);
      double cosThMiddleZ = TMath::Cos( ThMiddleZ );

      const int xBin = (isWest ? splitProb[1][1]->GetXaxis()->FindBin(xAtMiddleZ) :
                                 splitProb[0][1]->GetXaxis()->FindBin(xAtMiddleZ));
      const int yBin = splitProb[1][1]->GetYaxis()->FindBin(cosThMiddleZ);
      if (xBin == 0 || xBin == splitProb[1][1]->GetXaxis()->GetNbins() + 1 ||
          yBin == 0 || yBin == splitProb[1][1]->GetYaxis()->GetNbins() + 1) {
        if(DoDebug){
          std::cout << "Note that x = " << xAtMiddleZ << " and cosTh = " << cosThMiddleZ << " for z=0 gap." << std::endl;
          std::cout << "BINS: x=" << xBin << ", y=" << yBin << std::endl;
          std::cout << "   start: (" << pfp.trk.start.x << ", " << pfp.trk.start.y << ", " << pfp.trk.start.z << ")" << std::endl;
          std::cout << "   end: (" << pfp.trk.end.x << ", " << pfp.trk.end.y << ", " << pfp.trk.end.z << ")" << std::endl;
        }
        pToSplit = 0;
      }
      else{
/*
        // Using histogram
        pToSplit = (isWest ? splitProb[1][1]->GetBinContent(xBin,yBin) :
                                splitProb[0][1]->GetBinContent(xBin,yBin) );
*/
        // Using fixed number from BNB
        pToSplit = 0.19;

      }
    } // END if NOT crossesCath

    if ( DoDebug ){
      std::cout << "This track crosses a boundary and its p for splitting is: " << pToSplit << std::endl;
    }


    // Set up the random engine
    // Set a seed using track variables
    std::uint32_t seed = int(pfp.trk.len) +
                         10 * abs(pfp.trk.start.x) +
                         100 * abs(pfp.trk.start.z) +
                         1000 * abs(pfp.trk.end.x) +
                         10000 * abs(pfp.trk.end.z);
    TRandom3 tRand(seed);
    double sigma = 1.0;

    // JSKIMDEBUG TODO
    //sigma = 0.4;

    if ( tRand.Rndm() > pToSplit*sigma ) {
      if( DoDebug ){
        std::cout << "... but we did not split it." << std::endl;
      }
      continue;
    }

    // DEBUGGING!
    if( DoDebug ){
      std::cout << "Where? ==> " << (crossesCath ? "Cathode" : "MiddleZ") << std::endl;
      std::cout << "Trk Start (x,y,z) = (" << pfp.trk.start.x << "," << pfp.trk.start.y << "," << pfp.trk.start.z << ")" << std::endl;
      std::cout << "Trk End (x,y,z) = (" << pfp.trk.end.x << "," << pfp.trk.end.y << "," << pfp.trk.end.z << ")" << std::endl;
      std::cout << "Trk original length = " << pfp.trk.len << std::endl;
      std::cout << "... and we've decided to split the track." << std::endl;
    }

    // Find the calo point at which to break the track on each plane, and save the remaining calo points
    caf::SRCaloPoint splitPoint[3];
    caf::SRCaloPoint firstPoint[3];

    std::vector< std::vector<caf::SRCaloPoint> > savedCaloPoints;
    std::vector< std::vector<caf::SRCaloPoint> > savedScdryCaloPoints;

    PreservedInitialTrkResults PreservedTrkResults;

    for ( unsigned int idxPlane=0; idxPlane < 3; ++idxPlane ){
      savedCaloPoints.push_back({});
      savedScdryCaloPoints.push_back({});
      double minDist = std::numeric_limits<double>::max();
      double largestRR = 0.;
      unsigned int nRRThatAreNaN = 0;
      unsigned int nPoints = 0;

      // For each plane, find 
      // 1) split point
      // 2) First point (i.e., larged RR)
      for ( auto const& caloPt : pfp.trk.calo[idxPlane].points ) {
        // TODO: REALLY THESE SHOULD GRAB CLOSEST POINT *ON* THE SIDE OF THE SPLIT. THIS IS JUST CLOSEST POINT OVERALL...
        if ( crossesCath ) {
/*
          if ( DoDebug ){
            std::cout << caloPt.x << "(" << fabs(caloPt.x)-210.2 << ")";
          }
*/
          if ( fabs(fabs(caloPt.x)-210.2) < minDist ){
            minDist = fabs(fabs(caloPt.x)-210.2);
            splitPoint[idxPlane] = FillCaloPointFrom( caloPt );
/*
            if ( DoDebug ){
              std::cout << "[" << splitPoint[idxPlane].x << " " << caloPt.x << "] ";
            }
*/
          }
        }
        else {
/*
          if ( DoDebug ){
            std::cout << caloPt.z << "(" << fabs(caloPt.z) << ")";
          }
*/
          if ( fabs(caloPt.z) < minDist ){
            minDist = fabs(caloPt.z);
            splitPoint[idxPlane] = FillCaloPointFrom( caloPt );
/*
            if ( DoDebug ){
              std::cout << "[" << splitPoint[idxPlane].z << " " << caloPt.z << "] ";
            }
*/
          }
        }

        if ( !std::isnan(caloPt.rr) && caloPt.rr > largestRR ) {
          largestRR = caloPt.rr;
          firstPoint[idxPlane] = FillCaloPointFrom( caloPt );
        }
        else if ( std::isnan(caloPt.rr) ) {
          nRRThatAreNaN+=1;
        }
      } // END Loop calao points

      if ( DoDebug ){
        std::cout << std::endl;
      }

      // Now for the first track, shift its RRs to the split point
      // For the second track, RRs remain the same, and copy to savedScdryCaloPoints
      for ( auto const& caloPt : pfp.trk.calo[idxPlane].points ) {
        caf::SRCaloPoint savedCaloPoint = FillCaloPointFrom(caloPt);

        if ( caloPt.rr >= splitPoint[idxPlane].rr ) {
          // this is before the split
          savedCaloPoint.rr = caloPt.rr - splitPoint[idxPlane].rr;
          savedCaloPoints[idxPlane].push_back(savedCaloPoint);
          if ( savedCaloPoint.dedx < 1000. ) nPoints+=1;
        }
        else {
          // beyond split point, save to savedScdryCaloPoints
          savedScdryCaloPoints[idxPlane].push_back( savedCaloPoint );
        }
      } // END Loop calao points

      // Update current track to the first track
      // Update the calo info as needed
      pfp.trk.calo[idxPlane].nhit = nPoints;
      pfp.trk.calo[idxPlane].ke = caf::kSignalingNaN;
      pfp.trk.calo[idxPlane].charge = caf::kSignalingNaN;
      pfp.trk.calo[idxPlane].points = savedCaloPoints[idxPlane];

      // Re-calc chi2 for the first track
      TrkChi2Results chi2output = CalculateChi2(pfp.trk.calo[idxPlane]);
      pfp.trk.chi2pid[idxPlane].chi2_proton = chi2output.chi2_proton;
      pfp.trk.chi2pid[idxPlane].chi2_kaon = chi2output.chi2_kaon;
      pfp.trk.chi2pid[idxPlane].chi2_pion = chi2output.chi2_pion;
      pfp.trk.chi2pid[idxPlane].chi2_muon = chi2output.chi2_muon;
      pfp.trk.chi2pid[idxPlane].pida = chi2output.pida;
      pfp.trk.chi2pid[idxPlane].pid_ndof = chi2output.pid_ndof;

      // DEBUGS
      if( DoDebug ){
        std::cout << "... PLANE DEBUG INFO: " << idxPlane << std::endl;
        //std::cout << "    initial points = " << initialPoints << std::endl;
        std::cout << "    --> Post systematic shift" << std::endl;
        std::cout << "    nRRThatAreNaN = " << nRRThatAreNaN << std::endl;
        std::cout << "    nhit = " << nPoints << std::endl;
        std::cout << "    first point (x,y,z) = (" << firstPoint[idxPlane].x << "," << firstPoint[idxPlane].y << "," << firstPoint[idxPlane].z << ")" << std::endl;
        std::cout << "    first point rr = " << firstPoint[idxPlane].rr << std::endl;
        std::cout << "    split point (x,y,z) = (" << splitPoint[idxPlane].x << "," << splitPoint[idxPlane].y << "," << splitPoint[idxPlane].z << ")" << std::endl;
        std::cout << "    split point rr = " << splitPoint[idxPlane].rr << std::endl;
      }
    } // loop planes

    // Save the Preserved Track Info
    PreservedTrkResults.producer = pfp.trk.producer;
    PreservedTrkResults.end_x = pfp.trk.end.x;
    PreservedTrkResults.end_y = pfp.trk.end.y;
    PreservedTrkResults.end_z = pfp.trk.end.z;
    PreservedTrkResults.dir_end_x = pfp.trk.dir_end.x;
    PreservedTrkResults.dir_end_y = pfp.trk.dir_end.y;
    PreservedTrkResults.dir_end_z = pfp.trk.dir_end.z;

    // Update the track info
    //pfp.trk.len = ((firstPoint[2].rr-splitPoint[2].rr) + (firstPoint[1].rr-splitPoint[1].rr))/2.;
    pfp.trk.len = std::max(firstPoint[0].rr, std::max(firstPoint[1].rr, firstPoint[2].rr)) - std::min(splitPoint[0].rr, std::min(splitPoint[1].rr, splitPoint[2].rr));
    pfp.trk.npts = (pfp.trk.calo[2].nhit + pfp.trk.calo[1].nhit + pfp.trk.calo[0].nhit);
    pfp.trk.bestplane = ( pfp.trk.calo[2].nhit >= pfp.trk.calo[1].nhit ? (pfp.trk.calo[2].nhit >= pfp.trk.calo[0].nhit ? caf::Plane_t(2) : caf::Plane_t(0) ) :
                                                                         (pfp.trk.calo[0].nhit >= pfp.trk.calo[1].nhit ? caf::Plane_t(0) : caf::Plane_t(1) ));

    if ( DoDebug ){
      std::cout << "... and its length post-split is now: " << pfp.trk.len << std::endl;
    }

    // updated end point: 24 May 2024
    double endRR = std::numeric_limits<double>::max();
    double endPtX = -9999.;
    double endPtY = -9999.;
    double endPtZ = -9999.;
    for ( unsigned int idxPt = 0; idxPt < 3; ++idxPt ) {
      if ( std::isnan(splitPoint[idxPt].rr) || std::isinf(splitPoint[idxPt].rr) ) continue;
      if ( splitPoint[idxPt].rr < endRR ) {
        endRR = splitPoint[idxPt].rr;
        endPtX = splitPoint[idxPt].x;
        endPtY = splitPoint[idxPt].y;
        endPtZ = splitPoint[idxPt].z;
      }
    }
    caf::SRVector3D trkEnd( endPtX, endPtY, endPtZ );
    pfp.trk.end = trkEnd;

    caf::SRVector3D trkEndDir( caf::kSignalingNaN, caf::kSignalingNaN, caf::kSignalingNaN );
    pfp.trk.dir_end = trkEndDir;

    TrkMomentumResults p_output = CalculateMomenta((float)pfp.trk.len);
    pfp.trk.rangeP.p_muon = p_output.p_muon;
    pfp.trk.rangeP.p_proton = p_output.p_proton;
    pfp.trk.rangeP.p_pion = p_output.p_pion;

    // NaN some of the remaining variables we don't want the user to use with this shift...
    // TODO: maybe think of rerunning MCS with calo point info?
    caf::SRTrkMCS newMCS;
    pfp.trk.mcsP = newMCS;

    caf::SRTrackScatterClosestApproach newScatter;
    pfp.trk.scatterClosestApproach = newScatter;

    caf::SRTrackStoppingChi2Fit newStopChi2;
    pfp.trk.stoppingChi2Fit = newStopChi2;

    caf::SRTrackDazzle newDazzle;
    pfp.trk.dazzle = newDazzle;

    // Update shower variables... for this one, we ONLY update the shower length to be 90% the track length...
    // -- We may want to update others later...
    pfp.shw.len = 0.9*pfp.trk.len;

    // Save the secondary track info
    unsigned int pfpIdxUI = (unsigned int)pfpIdx;
    scdryInitIdxToPfpCaloPts[ pfpIdxUI ] = savedScdryCaloPoints;
    scdryInitIdxToPreservedInfo[ pfpIdxUI ] = PreservedTrkResults;
    expectedScdryTrks+=1;

  } // END Loop pfp

  // Now take care of the secondary track
  if ( expectedScdryTrks > 0 ){
    if ( DoDebug ){
      std::cout << "... Dealing with " << expectedScdryTrks << " secondary tracks." << std::endl;
    }

    for( auto const& [initIdx, caloPtsVec] : scdryInitIdxToPfpCaloPts ) {

      caf::SRPFP newPfp;

      // PFP Elements:
      newPfp.id = -999; // Some default...
      newPfp.parent = slc.reco.pfp[initIdx].id;
      newPfp.parent_is_primary = false;
      newPfp.trackScore = slc.reco.pfp[initIdx].trackScore;
      // --> also the score vars...
      newPfp.pfochar.chgendfrac = slc.reco.pfp[initIdx].pfochar.chgendfrac;
      newPfp.pfochar.chgfracspread = slc.reco.pfp[initIdx].pfochar.chgfracspread;
      newPfp.pfochar.linfitdiff = slc.reco.pfp[initIdx].pfochar.linfitdiff;
      newPfp.pfochar.linfitlen = slc.reco.pfp[initIdx].pfochar.linfitlen;
      newPfp.pfochar.linfitgaplen = slc.reco.pfp[initIdx].pfochar.linfitgaplen;
      newPfp.pfochar.linfitrms = slc.reco.pfp[initIdx].pfochar.linfitrms;
      newPfp.pfochar.openanglediff = slc.reco.pfp[initIdx].pfochar.openanglediff;
      newPfp.pfochar.pca2ratio = slc.reco.pfp[initIdx].pfochar.pca2ratio;
      newPfp.pfochar.pca3ratio = slc.reco.pfp[initIdx].pfochar.pca3ratio;
      newPfp.pfochar.vtxdist = slc.reco.pfp[initIdx].pfochar.vtxdist;
      /////////////////////////////
      newPfp.slcID = slc.reco.pfp[initIdx].slcID;
      newPfp.t0 = (std::isnan(slc.reco.pfp[initIdx].t0) ? caf::kSignalingNaN : (float)slc.reco.pfp[initIdx].t0);

      // Track Elements:
      // Things unchanged from initial track...
      newPfp.trk.producer = scdryInitIdxToPreservedInfo[ initIdx ].producer;

      caf::SRVector3D trkEnd( scdryInitIdxToPreservedInfo[ initIdx ].end_x,
                              scdryInitIdxToPreservedInfo[ initIdx ].end_y,
                              scdryInitIdxToPreservedInfo[ initIdx ].end_z );
      newPfp.trk.end = trkEnd;

      caf::SRVector3D trkEndDir( scdryInitIdxToPreservedInfo[ initIdx ].dir_end_x,
                                 scdryInitIdxToPreservedInfo[ initIdx ].dir_end_y,
                                 scdryInitIdxToPreservedInfo[ initIdx ].dir_end_z );
      newPfp.trk.dir_end = trkEndDir;

      // Things changed from initial track...
      std::vector< std::map<double, caf::SRCaloPoint> > caloPointRRMap;

      for ( unsigned int idxPlane=0; idxPlane < 3; ++idxPlane ){
        caloPointRRMap.push_back( std::map<double, caf::SRCaloPoint>() );
        unsigned int nRRThatAreNaN = 0;
        unsigned int nPoints = 0;

        for ( auto const& caloPt : caloPtsVec[idxPlane] ) {
          if ( !std::isnan(caloPt.rr) ) {
            caloPointRRMap[idxPlane][ caloPt.rr ] = caloPt;
          }
          else if ( std::isnan(caloPt.rr) ) {
            nRRThatAreNaN+=1;
          }
        } // Get the ordering of points to be able to find "first" 5

        for ( auto const& caloPt : caloPtsVec[idxPlane] ) {
          if ( caloPt.dedx < 1000. ) nPoints+=1;
        } // calo points to save

        // Update the calo info as needed
        newPfp.trk.calo[idxPlane].nhit = nPoints;
        newPfp.trk.calo[idxPlane].ke = caf::kSignalingNaN;
        newPfp.trk.calo[idxPlane].charge = caf::kSignalingNaN;
        newPfp.trk.calo[idxPlane].points = caloPtsVec[idxPlane];

        // The Chi2PID
        TrkChi2Results scdryChi2output = CalculateChi2(newPfp.trk.calo[idxPlane]);
        newPfp.trk.chi2pid[idxPlane].chi2_proton = scdryChi2output.chi2_proton;
        newPfp.trk.chi2pid[idxPlane].chi2_kaon = scdryChi2output.chi2_kaon;
        newPfp.trk.chi2pid[idxPlane].chi2_pion = scdryChi2output.chi2_pion;
        newPfp.trk.chi2pid[idxPlane].chi2_muon = scdryChi2output.chi2_muon;
        newPfp.trk.chi2pid[idxPlane].pida = scdryChi2output.pida;
        newPfp.trk.chi2pid[idxPlane].pid_ndof = scdryChi2output.pid_ndof;
      } // loop planes

/*
      if ( DoDebug ) {
        std::cout << "Secondary track RRs: ";
        for ( auto const& [rrs, pt] : caloPointRRMap[2] ) {
          std::cout << rrs << " ";
        }
        std::cout << std::endl;
      }
*/
      // Attempt to get the start direction
      std::vector< std::vector<double> > startPoints;
      std::vector< std::vector<double> > dirs;
      std::vector<double> scdryLengths;
      for ( unsigned int idxPlane=1; idxPlane < 3; ++idxPlane ) {
        unsigned int nPoints = 0;
        std::map<double, caf::SRCaloPoint>::iterator it = caloPointRRMap[idxPlane].end();
        if ( caloPointRRMap[idxPlane].size() == 0 ) continue; // if no points then we skip this plane
        --it; // get away from .end()
        scdryLengths.push_back( it->second.rr );
        std::vector<double> startPoint = { it->second.x, it->second.y, it->second.z };
        startPoints.push_back(startPoint);
        while ( it!=caloPointRRMap[idxPlane].begin() ) {
          double xyz_start[3] = { it->second.x, it->second.y, it->second.z };
          //std::cout << "xyz: " << xyz_start[0] << " " << xyz_start[1] << " " << xyz_start[2] << std::endl;
          --it;
          double xyz_end[3] = { it->second.x, it->second.y, it->second.z };
          nPoints+=1;
          std::vector<double> dir = { xyz_end[0]-xyz_start[0], xyz_end[1]-xyz_start[1], xyz_end[2]-xyz_start[2] };
          dirs.push_back(dir);
          if ( nPoints==5 ) break;
        }
      }

      std::vector< double > ave_dir = {0., 0., 0.};
      for ( auto const& dir: dirs ) {
        TVector3 dir_tv3( dir[0], dir[1], dir[2] );
        dir_tv3 = dir_tv3.Unit();
        ave_dir[0] += dir_tv3.X();
        ave_dir[1] += dir_tv3.Y();
        ave_dir[2] += dir_tv3.Z();
      }
      if ( dirs.size() != 0 ) {
        ave_dir[0] /= double(dirs.size());
        ave_dir[1] /= double(dirs.size());
        ave_dir[2] /= double(dirs.size());
      }

      TVector3 ave_dir_tv3( ave_dir[0], ave_dir[1], ave_dir[2] );
      if ( dirs.size() != 0 ) ave_dir_tv3 = ave_dir_tv3.Unit();

      caf::SRVector3D trkStartDir( ave_dir_tv3.X(),
                                   ave_dir_tv3.Y(),
                                   ave_dir_tv3.Z() );
      newPfp.trk.dir = trkStartDir;

      // Length
      double scdryLength = 0.;
      // Commenting out to use this updated length definition (24 march 2024)
      //for ( auto const& scdryLenVal : scdryLengths ) scdryLength+=scdryLenVal;
      //if( scdryLengths.size() != 0 ) scdryLength /= double(scdryLengths.size());
      for ( auto const& scdryLenVal : scdryLengths ) {
        if ( scdryLenVal > scdryLength ) scdryLength=scdryLenVal;
      }
      newPfp.trk.len = scdryLength;

      // Start point (updated 24 may 2024)
      double startRR = 0.;
      double startPtX = caf::kSignalingNaN;
      double startPtY = caf::kSignalingNaN;
      double startPtZ = caf::kSignalingNaN;
      for ( unsigned int idxPt = 0 ; idxPt < startPoints.size(); ++idxPt ) {
        if ( scdryLengths[idxPt] > startRR ) {
          startRR = scdryLengths[idxPt];
          startPtX = startPoints[idxPt][0];
          startPtY = startPoints[idxPt][1];
          startPtZ = startPoints[idxPt][2];
        }
      }
      caf::SRVector3D trkStartLoc( startPtX,
                                   startPtY,
                                   startPtZ );
      newPfp.trk.start = trkStartLoc;

      // N points & best plane
      newPfp.trk.npts = (newPfp.trk.calo[2].nhit + newPfp.trk.calo[1].nhit + newPfp.trk.calo[0].nhit);
      newPfp.trk.bestplane = ( newPfp.trk.calo[2].nhit >= newPfp.trk.calo[1].nhit ? (newPfp.trk.calo[2].nhit >= newPfp.trk.calo[0].nhit ? caf::Plane_t(2) : caf::Plane_t(0) ) :
                                                                                    (newPfp.trk.calo[0].nhit >= newPfp.trk.calo[1].nhit ? caf::Plane_t(0) : caf::Plane_t(1) ));

      TrkMomentumResults p_output = CalculateMomenta(newPfp.trk.len);
      newPfp.trk.rangeP.p_muon = p_output.p_muon;
      newPfp.trk.rangeP.p_proton = p_output.p_proton;
      newPfp.trk.rangeP.p_pion = p_output.p_pion;

      // Shower Elements:
      // for our shower cut we check shw.start.x, shw.len, shw.conversion_gap
      // Update length as above
      newPfp.shw.len = 0.9*newPfp.trk.len;
      // Make start X, Y, Z be the same as the track start
      newPfp.shw.start.x = !std::isnan(newPfp.trk.start.x) ? newPfp.trk.start.x : caf::kSignalingNaN;
      newPfp.shw.start.y = !std::isnan(newPfp.trk.start.y) ? newPfp.trk.start.y : caf::kSignalingNaN;
      newPfp.shw.start.z = !std::isnan(newPfp.trk.start.z) ? newPfp.trk.start.z : caf::kSignalingNaN;
      // Make conversion gap be the distance between the slice vertex and this point...
      if ( std::isnan(newPfp.shw.start.x) || std::isnan(newPfp.shw.start.y) || std::isnan(newPfp.shw.start.z) ) {
        newPfp.shw.conversion_gap = caf::kSignalingNaN;
      }
      else {
        newPfp.shw.conversion_gap = std::hypot( newPfp.shw.start.x - slc.vertex.x,
                                                newPfp.shw.start.y - slc.vertex.y,
                                                newPfp.shw.start.z - slc.vertex.z );
      }

      // TODO add the pfp list

      if( DoDebug ){
        std::cout << "... Secondary track length is " << newPfp.trk.len << std::endl;
      }

      slc.reco.pfp.push_back( newPfp );

    } // END loop over scdryInitIdxToPfpCaloPts

  } // END if expectedScdryTrks>0

  if( DoDebug ){
    std::cout << "[TrackSplit::ApplyTrackSplit] Done; now slc.reco.pfp.size() = " << slc.reco.pfp.size() << std::endl;
  }

}

caf::SRCaloPoint TrackSplit::FillCaloPointFrom( const caf::SRCaloPoint& inCaloPt ) {

  caf::SRCaloPoint outCaloPt;
  outCaloPt.rr = inCaloPt.rr;
  outCaloPt.dqdx = inCaloPt.dqdx;
  outCaloPt.dedx = inCaloPt.dedx;
  outCaloPt.pitch = inCaloPt.pitch;
  outCaloPt.t = inCaloPt.t;
  outCaloPt.x = inCaloPt.x;
  outCaloPt.y = inCaloPt.y;
  outCaloPt.z = inCaloPt.z;
  outCaloPt.integral = inCaloPt.integral;
  outCaloPt.sumadc = inCaloPt.sumadc;
  outCaloPt.width = inCaloPt.width;
  outCaloPt.mult = inCaloPt.mult;
  outCaloPt.wire = inCaloPt.wire;
  outCaloPt.tpc = inCaloPt.tpc;
  outCaloPt.start = inCaloPt.start;
  outCaloPt.end = inCaloPt.end;
  outCaloPt.channel = inCaloPt.channel;

  caf::SRTrueCaloPoint outTruth;
  outTruth.h_nelec = inCaloPt.truth.h_nelec;
  outTruth.h_e = inCaloPt.truth.h_nelec;
  outTruth.p_nelec = inCaloPt.truth.h_nelec;
  outTruth.h_e = inCaloPt.truth.h_nelec;
  outTruth.x = inCaloPt.truth.h_nelec;
  outTruth.y = inCaloPt.truth.h_nelec;
  outTruth.z = inCaloPt.truth.h_nelec;
  outTruth.rr = inCaloPt.truth.h_nelec;
  outTruth.pitch = inCaloPt.truth.h_nelec;

  outCaloPt.truth = outTruth;

  return outCaloPt;
}

caf::SRPFP TrackSplit::FillPtrPFP( const caf::SRPFP& inPfp ){

  caf::SRPFP ret;

  ret.id = inPfp.id;
  ret.ndaughters = inPfp.ndaughters;
  // Set ret.daughters = inPfp.daughters;
  for ( auto const& daughter : inPfp.daughters ) ret.daughters.push_back( int(daughter) );

  ret.parent = inPfp.parent;
  ret.parent_is_primary = inPfp.parent_is_primary;

  ret.trackScore = inPfp.trackScore;
  // Set ret.pfochar = inPfp.pfochar;
  ret.pfochar.chgendfrac = inPfp.pfochar.chgendfrac;
  ret.pfochar.chgfracspread = inPfp.pfochar.chgfracspread;
  ret.pfochar.linfitdiff = inPfp.pfochar.linfitdiff;
  ret.pfochar.linfitlen = inPfp.pfochar.linfitlen;
  ret.pfochar.linfitgaplen = inPfp.pfochar.linfitgaplen;
  ret.pfochar.linfitrms = inPfp.pfochar.linfitrms;
  ret.pfochar.openanglediff = inPfp.pfochar.openanglediff;
  ret.pfochar.pca2ratio = inPfp.pfochar.pca2ratio;
  ret.pfochar.pca3ratio = inPfp.pfochar.pca3ratio;
  ret.pfochar.vtxdist = inPfp.pfochar.vtxdist;

  ret.slcID = inPfp.slcID;
  ret.t0 = inPfp.t0;

  // Set ret.trk = inPfp.trk;
  ret.trk.producer = inPfp.trk.producer;
  ret.trk.npts = inPfp.trk.npts;
  ret.trk.len = inPfp.trk.len;
  ret.trk.costh = inPfp.trk.costh;
  ret.trk.phi = inPfp.trk.phi;
  ret.trk.dir.x = inPfp.trk.dir.x;
  ret.trk.dir.y = inPfp.trk.dir.y;
  ret.trk.dir.z = inPfp.trk.dir.z;
  ret.trk.dir_end.x = inPfp.trk.dir_end.x;
  ret.trk.dir_end.y = inPfp.trk.dir_end.y;
  ret.trk.dir_end.z = inPfp.trk.dir_end.z;
  ret.trk.start.x = inPfp.trk.start.x;
  ret.trk.start.y = inPfp.trk.start.y;
  ret.trk.start.z = inPfp.trk.start.z;
  ret.trk.end.x = inPfp.trk.end.x;
  ret.trk.end.y = inPfp.trk.end.y;
  ret.trk.end.z = inPfp.trk.end.z;
  ret.trk.bestplane = inPfp.trk.bestplane;
  // Trk Chi2PID
  ret.trk.chi2pid[0].pdg = inPfp.trk.chi2pid[0].pdg;
  ret.trk.chi2pid[0].pid_ndof = inPfp.trk.chi2pid[0].pid_ndof;
  ret.trk.chi2pid[0].chi2_muon = inPfp.trk.chi2pid[0].chi2_muon;
  ret.trk.chi2pid[0].chi2_pion = inPfp.trk.chi2pid[0].chi2_pion;
  ret.trk.chi2pid[0].chi2_kaon = inPfp.trk.chi2pid[0].chi2_kaon;
  ret.trk.chi2pid[0].chi2_proton = inPfp.trk.chi2pid[0].chi2_proton;
  ret.trk.chi2pid[0].pida = inPfp.trk.chi2pid[0].pida;
  ret.trk.chi2pid[1].pdg = inPfp.trk.chi2pid[1].pdg;
  ret.trk.chi2pid[1].pid_ndof = inPfp.trk.chi2pid[1].pid_ndof;
  ret.trk.chi2pid[1].chi2_muon = inPfp.trk.chi2pid[1].chi2_muon;
  ret.trk.chi2pid[1].chi2_pion = inPfp.trk.chi2pid[1].chi2_pion;
  ret.trk.chi2pid[1].chi2_kaon = inPfp.trk.chi2pid[1].chi2_kaon;
  ret.trk.chi2pid[1].chi2_proton = inPfp.trk.chi2pid[1].chi2_proton;
  ret.trk.chi2pid[1].pida = inPfp.trk.chi2pid[1].pida;
  ret.trk.chi2pid[2].pdg = inPfp.trk.chi2pid[2].pdg;
  ret.trk.chi2pid[2].pid_ndof = inPfp.trk.chi2pid[2].pid_ndof;
  ret.trk.chi2pid[2].chi2_muon = inPfp.trk.chi2pid[2].chi2_muon;
  ret.trk.chi2pid[2].chi2_pion = inPfp.trk.chi2pid[2].chi2_pion;
  ret.trk.chi2pid[2].chi2_kaon = inPfp.trk.chi2pid[2].chi2_kaon;
  ret.trk.chi2pid[2].chi2_proton = inPfp.trk.chi2pid[2].chi2_proton;
  ret.trk.chi2pid[2].pida = inPfp.trk.chi2pid[2].pida;
  // Trk Calo
  ret.trk.calo[0].nhit = inPfp.trk.calo[0].nhit;
  ret.trk.calo[0].ke = inPfp.trk.calo[0].ke;
  ret.trk.calo[0].charge = inPfp.trk.calo[0].charge;
  for ( auto const& caloPt : inPfp.trk.calo[0].points ) ret.trk.calo[0].points.push_back( FillCaloPointFrom(caloPt) );
  ret.trk.calo[1].nhit = inPfp.trk.calo[1].nhit;
  ret.trk.calo[1].ke = inPfp.trk.calo[1].ke;
  ret.trk.calo[1].charge = inPfp.trk.calo[1].charge;
  for ( auto const& caloPt : inPfp.trk.calo[1].points ) ret.trk.calo[1].points.push_back( FillCaloPointFrom(caloPt) );
  ret.trk.calo[2].nhit = inPfp.trk.calo[2].nhit;
  ret.trk.calo[2].ke = inPfp.trk.calo[2].ke;
  ret.trk.calo[2].charge = inPfp.trk.calo[2].charge;
  for ( auto const& caloPt : inPfp.trk.calo[2].points ) ret.trk.calo[2].points.push_back( FillCaloPointFrom(caloPt) );
  // MCS P
  ret.trk.mcsP.fwdP_muon = inPfp.trk.mcsP.fwdP_muon;
  ret.trk.mcsP.fwdP_pion = inPfp.trk.mcsP.fwdP_pion;
  ret.trk.mcsP.fwdP_kaon = inPfp.trk.mcsP.fwdP_kaon;
  ret.trk.mcsP.fwdP_proton = inPfp.trk.mcsP.fwdP_proton;
  ret.trk.mcsP.fwdP_err_muon = inPfp.trk.mcsP.fwdP_err_muon;
  ret.trk.mcsP.fwdP_err_pion = inPfp.trk.mcsP.fwdP_err_pion;
  ret.trk.mcsP.fwdP_err_kaon = inPfp.trk.mcsP.fwdP_err_kaon;
  ret.trk.mcsP.fwdP_err_proton = inPfp.trk.mcsP.fwdP_err_proton;
  ret.trk.mcsP.bwdP_muon = inPfp.trk.mcsP.bwdP_muon;
  ret.trk.mcsP.bwdP_pion = inPfp.trk.mcsP.bwdP_pion;
  ret.trk.mcsP.bwdP_kaon = inPfp.trk.mcsP.bwdP_kaon;
  ret.trk.mcsP.bwdP_proton = inPfp.trk.mcsP.bwdP_proton;
  ret.trk.mcsP.bwdP_err_muon = inPfp.trk.mcsP.bwdP_err_muon;
  ret.trk.mcsP.bwdP_err_pion = inPfp.trk.mcsP.bwdP_err_pion;
  ret.trk.mcsP.bwdP_err_kaon = inPfp.trk.mcsP.bwdP_err_kaon;
  ret.trk.mcsP.bwdP_err_proton = inPfp.trk.mcsP.bwdP_err_proton;
  ret.trk.mcsP.is_bwd_muon = inPfp.trk.mcsP.is_bwd_muon;
  ret.trk.mcsP.is_bwd_pion = inPfp.trk.mcsP.is_bwd_pion;
  ret.trk.mcsP.is_bwd_kaon = inPfp.trk.mcsP.is_bwd_kaon;
  ret.trk.mcsP.is_bwd_proton = inPfp.trk.mcsP.is_bwd_proton;
  // Range P
  ret.trk.rangeP.p_muon = inPfp.trk.rangeP.p_muon;
  ret.trk.rangeP.p_pion = inPfp.trk.rangeP.p_pion;
  ret.trk.rangeP.p_proton = inPfp.trk.rangeP.p_proton;
  // Trk Truth
  ret.trk.truth.visEintrk = inPfp.trk.truth.visEintrk;
  ret.trk.truth.eff = inPfp.trk.truth.eff;
  ret.trk.truth.eff_cryo = inPfp.trk.truth.eff_cryo;
  ret.trk.truth.pur = inPfp.trk.truth.pur;
  ret.trk.truth.nmatches = inPfp.trk.truth.nmatches;
  // --> Leave the "matches" empty...
  ret.trk.truth.bestmatch.G4ID = inPfp.trk.truth.bestmatch.G4ID;
  ret.trk.truth.bestmatch.energy = inPfp.trk.truth.bestmatch.energy;
  ret.trk.truth.bestmatch.hit_completeness = inPfp.trk.truth.bestmatch.hit_completeness;
  ret.trk.truth.bestmatch.hit_purity = inPfp.trk.truth.bestmatch.hit_purity;
  ret.trk.truth.bestmatch.energy_completeness = inPfp.trk.truth.bestmatch.energy_completeness;
  ret.trk.truth.bestmatch.energy_purity = inPfp.trk.truth.bestmatch.energy_purity;
  ret.trk.truth.p.genE = inPfp.trk.truth.p.genE;
  // --> Leave the plane info NaN'ed/zero'ed out...
  ret.trk.truth.p.startE = inPfp.trk.truth.p.startE;
  ret.trk.truth.p.endE = inPfp.trk.truth.p.endE;
  ret.trk.truth.p.genT = inPfp.trk.truth.p.genT;
  ret.trk.truth.p.startT = inPfp.trk.truth.p.startT;
  ret.trk.truth.p.endT = inPfp.trk.truth.p.endT;
  ret.trk.truth.p.length = inPfp.trk.truth.p.length;
  ret.trk.truth.p.genp.x = inPfp.trk.truth.p.genp.x;
  ret.trk.truth.p.genp.y = inPfp.trk.truth.p.genp.y;
  ret.trk.truth.p.genp.z = inPfp.trk.truth.p.genp.z;
  ret.trk.truth.p.startp.x = inPfp.trk.truth.p.startp.x;
  ret.trk.truth.p.startp.y = inPfp.trk.truth.p.startp.y;
  ret.trk.truth.p.startp.z = inPfp.trk.truth.p.startp.z;
  ret.trk.truth.p.endp.x = inPfp.trk.truth.p.endp.x;
  ret.trk.truth.p.endp.y = inPfp.trk.truth.p.endp.y;
  ret.trk.truth.p.endp.z = inPfp.trk.truth.p.endp.z;
  ret.trk.truth.p.gen.x = inPfp.trk.truth.p.gen.x;
  ret.trk.truth.p.gen.y = inPfp.trk.truth.p.gen.y;
  ret.trk.truth.p.gen.z = inPfp.trk.truth.p.gen.z;
  ret.trk.truth.p.start.x = inPfp.trk.truth.p.start.x;
  ret.trk.truth.p.start.y = inPfp.trk.truth.p.start.y;
  ret.trk.truth.p.start.z = inPfp.trk.truth.p.start.z;
  ret.trk.truth.p.end.x = inPfp.trk.truth.p.end.x;
  ret.trk.truth.p.end.y = inPfp.trk.truth.p.end.y;
  ret.trk.truth.p.end.z = inPfp.trk.truth.p.end.z;
  ret.trk.truth.p.wallin = inPfp.trk.truth.p.wallin;
  ret.trk.truth.p.wallout = inPfp.trk.truth.p.wallout;
  ret.trk.truth.p.cont_tpc = inPfp.trk.truth.p.cont_tpc;
  ret.trk.truth.p.crosses_tpc = inPfp.trk.truth.p.crosses_tpc;
  ret.trk.truth.p.contained = inPfp.trk.truth.p.contained;
  ret.trk.truth.p.pdg = inPfp.trk.truth.p.pdg;
  ret.trk.truth.p.G4ID = inPfp.trk.truth.p.G4ID;
  ret.trk.truth.p.interaction_id = inPfp.trk.truth.p.interaction_id;
  ret.trk.truth.p.cryostat = inPfp.trk.truth.p.cryostat;
  for ( auto const& daughter : inPfp.trk.truth.p.daughters ) ret.trk.truth.p.daughters.push_back( (unsigned)daughter );
  ret.trk.truth.p.parent = inPfp.trk.truth.p.parent;
  ret.trk.truth.p.generator = inPfp.trk.truth.p.generator;
  ret.trk.truth.p.start_process = inPfp.trk.truth.p.start_process;
  ret.trk.truth.p.end_process = inPfp.trk.truth.p.end_process;
  ret.trk.truth.p.gstatus = inPfp.trk.truth.p.gstatus;
  // --> For now leave other stuff NaN'ed/zero'ed out...
  // Set ret.shw = inPfp.shw;
  ret.shw.bestplane = inPfp.shw.bestplane;
  ret.shw.bestplane_dEdx = inPfp.shw.bestplane_dEdx;
  ret.shw.bestplane_energy = inPfp.shw.bestplane_energy;
  ret.shw.conversion_gap = inPfp.shw.conversion_gap;
  ret.shw.density = inPfp.shw.density;
  ret.shw.len = inPfp.shw.len;
  ret.shw.open_angle = inPfp.shw.open_angle;
  ret.shw.dir.x = inPfp.shw.dir.x;
  ret.shw.dir.y = inPfp.shw.dir.y;
  ret.shw.dir.z = inPfp.shw.dir.z;
  ret.shw.start.x = inPfp.shw.start.x;
  ret.shw.start.y = inPfp.shw.start.y;
  ret.shw.start.z = inPfp.shw.start.z;
  ret.shw.end.x = inPfp.shw.end.x;
  ret.shw.end.y = inPfp.shw.end.y;
  ret.shw.end.z = inPfp.shw.end.z;
  ret.shw.cosmicDist = inPfp.shw.cosmicDist;
  ret.shw.producer = inPfp.shw.producer;
  // Shw Truth
  ret.shw.truth.visEintrk = inPfp.shw.truth.visEintrk;
  ret.shw.truth.eff = inPfp.shw.truth.eff;
  ret.shw.truth.eff_cryo = inPfp.shw.truth.eff_cryo;
  ret.shw.truth.pur = inPfp.shw.truth.pur;
  ret.shw.truth.nmatches = inPfp.shw.truth.nmatches;
  // --> Leave the "matches" empty...
  ret.shw.truth.bestmatch.G4ID = inPfp.shw.truth.bestmatch.G4ID;
  ret.shw.truth.bestmatch.energy = inPfp.shw.truth.bestmatch.energy;
  ret.shw.truth.bestmatch.hit_completeness = inPfp.shw.truth.bestmatch.hit_completeness;
  ret.shw.truth.bestmatch.hit_purity = inPfp.shw.truth.bestmatch.hit_purity;
  ret.shw.truth.bestmatch.energy_completeness = inPfp.shw.truth.bestmatch.energy_completeness;
  ret.shw.truth.bestmatch.energy_purity = inPfp.shw.truth.bestmatch.energy_purity;
  ret.shw.truth.p.genE = inPfp.shw.truth.p.genE;
  // --> Leave the plane info NaN'ed/zero'ed out...
  ret.shw.truth.p.startE = inPfp.shw.truth.p.startE;
  ret.shw.truth.p.endE = inPfp.shw.truth.p.endE;
  ret.shw.truth.p.genT = inPfp.shw.truth.p.genT;
  ret.shw.truth.p.startT = inPfp.shw.truth.p.startT;
  ret.shw.truth.p.endT = inPfp.shw.truth.p.endT;
  ret.shw.truth.p.length = inPfp.shw.truth.p.length;
  ret.shw.truth.p.genp.x = inPfp.shw.truth.p.genp.x;
  ret.shw.truth.p.genp.y = inPfp.shw.truth.p.genp.y;
  ret.shw.truth.p.genp.z = inPfp.shw.truth.p.genp.z;
  ret.shw.truth.p.startp.x = inPfp.shw.truth.p.startp.x;
  ret.shw.truth.p.startp.y = inPfp.shw.truth.p.startp.y;
  ret.shw.truth.p.startp.z = inPfp.shw.truth.p.startp.z;
  ret.shw.truth.p.endp.x = inPfp.shw.truth.p.endp.x;
  ret.shw.truth.p.endp.y = inPfp.shw.truth.p.endp.y;
  ret.shw.truth.p.endp.z = inPfp.shw.truth.p.endp.z;
  ret.shw.truth.p.gen.x = inPfp.shw.truth.p.gen.x;
  ret.shw.truth.p.gen.y = inPfp.shw.truth.p.gen.y;
  ret.shw.truth.p.gen.z = inPfp.shw.truth.p.gen.z;
  ret.shw.truth.p.start.x = inPfp.shw.truth.p.start.x;
  ret.shw.truth.p.start.y = inPfp.shw.truth.p.start.y;
  ret.shw.truth.p.start.z = inPfp.shw.truth.p.start.z;
  ret.shw.truth.p.end.x = inPfp.shw.truth.p.end.x;
  ret.shw.truth.p.end.y = inPfp.shw.truth.p.end.y;
  ret.shw.truth.p.end.z = inPfp.shw.truth.p.end.z;
  ret.shw.truth.p.wallin = inPfp.shw.truth.p.wallin;
  ret.shw.truth.p.wallout = inPfp.shw.truth.p.wallout;
  ret.shw.truth.p.cont_tpc = inPfp.shw.truth.p.cont_tpc;
  ret.shw.truth.p.crosses_tpc = inPfp.shw.truth.p.crosses_tpc;
  ret.shw.truth.p.contained = inPfp.shw.truth.p.contained;
  ret.shw.truth.p.pdg = inPfp.shw.truth.p.pdg;
  ret.shw.truth.p.G4ID = inPfp.shw.truth.p.G4ID;
  ret.shw.truth.p.interaction_id = inPfp.shw.truth.p.interaction_id;
  ret.shw.truth.p.cryostat = inPfp.shw.truth.p.cryostat;
  for ( auto const& daughter : inPfp.shw.truth.p.daughters ) ret.shw.truth.p.daughters.push_back( (unsigned)daughter );
  ret.shw.truth.p.parent = inPfp.shw.truth.p.parent;
  ret.shw.truth.p.generator = inPfp.shw.truth.p.generator;
  ret.shw.truth.p.start_process = inPfp.shw.truth.p.start_process;
  ret.shw.truth.p.end_process = inPfp.shw.truth.p.end_process;
  ret.shw.truth.p.gstatus = inPfp.shw.truth.p.gstatus;
  // --> Leave plane, razzle, and selVars info NaN'ed/zero'ed out...

  return ret;

}

TrackSplit::TrkChi2Results TrackSplit::CalculateChi2(const caf::SRTrackCalo& calo){

  // copied&modified from https://github.com/LArSoft/larana/blob/develop/larana/ParticleIdentification/Chi2PIDAlg.cxx#L60

  int npt = 0;
  double chi2pro = 0;
  double chi2ka = 0;
  double chi2pi = 0;
  double chi2mu = 0;
  double PIDA = 0; //by Bruce Baller
  std::vector<double> vpida;

  int used_trkres = 0;
  for(unsigned i = 0; i < calo.points.size(); ++i) { //hits
    const auto& pt = calo.points[i];
    double hit_dedx = pt.dedx;
    double hit_rr = pt.rr;

    //ignore the first and the last point
    if (i == 0 || i == calo.points.size() - 1) continue;
    if (hit_rr < 30) {
      PIDA += hit_dedx * pow(hit_rr, 0.42);
      vpida.push_back(hit_dedx * pow(hit_rr, 0.42));
      used_trkres++;
    }
    if (hit_dedx > 1000) continue; //protect against large pulse height
    int bin = dedx_range_pro->FindBin(hit_rr);
    if (bin >= 1 && bin <= dedx_range_pro->GetNbinsX()) {
      double bincpro = dedx_range_pro->GetBinContent(bin);
      if (bincpro < 1e-6) { //for 0 bin content, using neighboring bins
        bincpro =
          (dedx_range_pro->GetBinContent(bin - 1) + dedx_range_pro->GetBinContent(bin + 1)) / 2;
      }
      double bincka = dedx_range_ka->GetBinContent(bin);
      if (bincka < 1e-6) {
        bincka =
          (dedx_range_ka->GetBinContent(bin - 1) + dedx_range_ka->GetBinContent(bin + 1)) / 2;
      }
      double bincpi = dedx_range_pi->GetBinContent(bin);
      if (bincpi < 1e-6) {
        bincpi =
          (dedx_range_pi->GetBinContent(bin - 1) + dedx_range_pi->GetBinContent(bin + 1)) / 2;
      }
      double bincmu = dedx_range_mu->GetBinContent(bin);
      if (bincmu < 1e-6) {
        bincmu =
          (dedx_range_mu->GetBinContent(bin - 1) + dedx_range_mu->GetBinContent(bin + 1)) / 2;
      }
      double binepro = dedx_range_pro->GetBinError(bin);
      if (binepro < 1e-6) {
        binepro =
          (dedx_range_pro->GetBinError(bin - 1) + dedx_range_pro->GetBinError(bin + 1)) / 2;
      }
      double bineka = dedx_range_ka->GetBinError(bin);
      if (bineka < 1e-6) {
        bineka = (dedx_range_ka->GetBinError(bin - 1) + dedx_range_ka->GetBinError(bin + 1)) / 2;
      }
      double binepi = dedx_range_pi->GetBinError(bin);
      if (binepi < 1e-6) {
        binepi = (dedx_range_pi->GetBinError(bin - 1) + dedx_range_pi->GetBinError(bin + 1)) / 2;
      }
      double binemu = dedx_range_mu->GetBinError(bin);
      if (binemu < 1e-6) {
        binemu = (dedx_range_mu->GetBinError(bin - 1) + dedx_range_mu->GetBinError(bin + 1)) / 2;
      }
      //double errke = 0.05*hit_dedx;   //5% KE resolution
      double errdedx = 0.04231 + 0.0001783 * hit_dedx * hit_dedx; //resolution on dE/dx
      errdedx *= hit_dedx;
      chi2pro += pow((hit_dedx - bincpro) / std::sqrt(pow(binepro, 2) + pow(errdedx, 2)), 2);
      chi2ka += pow((hit_dedx - bincka) / std::sqrt(pow(bineka, 2) + pow(errdedx, 2)), 2);
      chi2pi += pow((hit_dedx - bincpi) / std::sqrt(pow(binepi, 2) + pow(errdedx, 2)), 2);
      chi2mu += pow((hit_dedx - bincmu) / std::sqrt(pow(binemu, 2) + pow(errdedx, 2)), 2);
      //std::cout<<i<<" "<<hit_dedx<<" "<<hit_rr<<" "<<bincpro<<std::endl;
      ++npt;
    }
  }

  // Making output
  TrkChi2Results output;
  if (npt) {
    chi2pro /= npt;
    chi2ka /= npt;
    chi2pi /= npt;
    chi2mu /= npt;
  }
  static bool fUseMedian = true; // I think default fcl in ICARUS is true
  if (used_trkres > 0) {
    if (fUseMedian) {
      PIDA = TMath::Median(vpida.size(), &vpida[0]);
    }
    else { // use mean
      PIDA /= used_trkres;
    }
  }

  output.chi2_kaon = chi2ka;
  output.chi2_muon = chi2mu;
  output.chi2_pion = chi2pi;
  output.chi2_proton = chi2pro;
  output.pida = PIDA;
  output.pid_ndof = npt;

  return output;

}

TrackSplit::TrkMomentumResults TrackSplit::CalculateMomenta(const float length){

  // largely copying and modifying from
  // https://github.com/LArSoft/larreco/blob/develop/larreco/RecoAlg/TrackMomentumCalculator.cxx
  // from LARSOFT_SUITE_v09_72_00. The track momentum calculator is used by:
  // https://github.com/SBNSoftware/sbncode/blob/develop/sbncode/LArRecoProducer/RangePAllPID_module.cc
  // and we copy its implementation of pion momentum (the larsoft module covers muon and proton...)

  // Re-implementing pieces to store the values, but the calculation bits should be the same.

  if (length < 0 || std::isnan(length)) {
    //std::cout << "Invalid track range " << length << " return -1" << std::endl;
    TrkMomentumResults output;
    output.p_proton = -1.;
    output.p_pion = -1.;
    output.p_muon = -1.;

    return output;
  }

  constexpr double Muon_M = 105.7, Proton_M = 938.272;
  constexpr double Pion_M = 139.57; // https://pdg.lbl.gov/2018/listings/rpp2018-list-pi-plus-minus.pdf

  TrkMomentumResults output;

  // Muon
  double KE_Muon = KEvsR_spline3.Eval(length);
  double p_muon = -999.;
  if ( KE_Muon >= 0 ) p_muon = std::sqrt((KE_Muon * KE_Muon) + (2 * Muon_M * KE_Muon));
  output.p_muon = p_muon / 1000;

  // Pion, as implemented in sbncode RangeP code
  double KE_Pion_FromMuCalc = KEvsR_spline3.Eval(length*Muon_M/Pion_M);
  double p_pion = -999.;
  if ( KE_Pion_FromMuCalc >= 0 )
    p_pion = std::sqrt((KE_Pion_FromMuCalc * KE_Pion_FromMuCalc) + (2 * Muon_M * KE_Pion_FromMuCalc)) * (Pion_M/Muon_M);
  output.p_pion = p_pion / 1000;

  // Proton
  double KE_Proton = -999.;
  if (length > 0 && length <= 80)
    KE_Proton = 29.9317 * std::pow(length, 0.586304);
  else if (length > 80 && length <= 3.022E3)
    KE_Proton = 149.904 + (3.34146 * length) + (-0.00318856 * length * length) +
                (4.34587E-6 * length * length * length) +
                (-3.18146E-9 * length * length * length * length) +
                (1.17854E-12 * length * length * length * length * length) +
                (-1.71763E-16 * length * length * length * length * length * length);
  double p_proton = -999.;
  if ( KE_Proton >= 0 ) p_proton = std::sqrt((KE_Proton * KE_Proton) + (2 * Proton_M * KE_Proton));
  output.p_proton = p_proton / 1000;

  return output;

}

} // END namespace sbnnusyst
