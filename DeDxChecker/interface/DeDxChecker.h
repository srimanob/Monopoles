#ifndef Monopoles_DeDxChecker_H
#define Monopoles_DeDxChecker_H

#include "DataFormats/TrackReco/interface/Track.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
//#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"

#include "DataFormats/Common/interface/Handle.h"

//#include "RecoTracker/DeDx/interface/UnbinnedLikelihoodFit.h"

//#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

//#include <TLinearFitter.h>
#include <TFile.h>
#include <TTree.h>

#include <vector>
#include <string>
#include <set>

using namespace std;
using namespace edm;

//typedef std::vector<Trajectory> TrajectoryCollection;


class DeDxChecker: public edm::EDAnalyzer{
  public:
    DeDxChecker(const edm::ParameterSet&);
    virtual ~DeDxChecker();

    virtual void beginJob();
    virtual void endJob();
    void Init(const edm::EventSetup&) ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
  private:
    void GetDeDxStrips(const reco::Track &Track);

    std::string _Source;
    std::string _Output;

    Handle<reco::TrackCollection> _hTracks;
    //Handle<TrajectoryCollection> _hTrajectories;
    //Handle<TrajTrackAssociationCollection> _hTrajTrackAssociations;
    Handle<ValueMap<reco::DeDxData> > _hDeDx;
    //edm::ESHandle<GlobalTrackingGeometry> _TrackingGeom;

    float _Pt, _Eta, _D0, _Z0, _Phi0, _DeDx;
    int _NTracks, _SatStrips, _TotStrips;
    vector<int> _HitSatStrips, _HitTotStrips;

    TFile *_OutputFile;
    TTree *_Tree;
};
#endif
