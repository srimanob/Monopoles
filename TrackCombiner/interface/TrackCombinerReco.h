#ifndef Monopoles_TrackCombinerReco_H
#define Monopoles_TrackCombinerReco_H

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
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>

#include <vector>
#include <string>
#include <set>

using namespace std;
using namespace edm;

//typedef std::vector<Trajectory> TrajectoryCollection;


class TrackCombinerReco: public edm::EDAnalyzer{
  public:
    TrackCombinerReco(const edm::ParameterSet&);
    virtual ~TrackCombinerReco();

    virtual void beginJob();
    virtual void endJob();
    void Init(const edm::EventSetup&) ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
  private:
    int AddPoints(const reco::Track &Track);
    void RemovePoints(int n);
    void AddMoreTracks(vector<int> &Group);
    void FitXY(vector<int> &Group);
    void FitRZ(bool Debug=false);
    void FitDeDx();
    void AverageIso(vector<int> &Group);
    void Save(vector<int> &Group);
    void Clear();

    std::string _Source;
    std::string _Output;
    float _PhiCut, _Chi2Cut, _PtCut, _DeDxCut, _DefaultError, _ErrorFudge;
    float _MeVperADCPixel, _MeVperADCStrip;

    Handle<reco::TrackCollection> _hTracks;
    //Handle<TrajectoryCollection> _hTrajectories;
    //Handle<TrajTrackAssociationCollection> _hTrajTrackAssociations;
    Handle<ValueMap<reco::DeDxData> > _hDeDx;
    edm::ESHandle<GlobalTrackingGeometry> _TrackingGeom;

    map<uint, float> _NormMap;

    set<int> _Used;

    vector<GlobalPoint> _Points;
    vector<GlobalError> _Errors;
    vector<float>       _Charges;

    //TLinearFitter *_RZFitter;
    TF1 *_RZFunc;
    TF1 *_XYFunc;
    //UnbinnedLikelihoodFit _DeDxFitter;
    //TF1 *_DeDxFunc;

    float _XYPar[3], _XYErr[3], _RZPar[3], _RZErr[3];
    float _Chi2XY, _Chi2RZ;
    int _NdofXY, _NdofRZ;
    float _DeDx;
    float _Iso;

    TFile *_OutputFile;
    TTree *_Tree;
    vector<float> _vXYPar0, _vXYPar1, _vXYPar2, _vXYErr0, _vXYErr1, _vXYErr2, _vRZPar0, _vRZPar1, _vRZPar2, _vRZErr0, _vRZErr1, _vRZErr2, _vChi2XY, _vChi2RZ, _vNdofXY, _vNdofRZ, _vDeDx, _vIso;
    vector<string> _vGroup;

    bool _TrackHitOutput;
    TTree *_TrackHitTree;
    vector<int> _vTHTrack;
    vector<float> _vTHX, _vTHY, _vTHZ, _vTHErrX, _vTHErrY, _vTHErrZ;

};
#endif
