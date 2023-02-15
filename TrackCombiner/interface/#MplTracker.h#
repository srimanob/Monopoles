#ifndef Monopoles_MplTracker_H
#define Monopoles_MplTracker_H

//#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"

#include "DataFormats/Common/interface/Handle.h"

//#include "RecoTracker/DeDx/interface/UnbinnedLikelihoodFit.h"

//#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

#include "Monopoles/MonoAlgorithms/interface/MonoTrack.h"
#include "Monopoles/MonoAlgorithms/interface/MonoTrackMatcher.h"
#include "Monopoles/MonoAlgorithms/interface/MonoEcalObs0.h"

//#include <TLinearFitter.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFile.h>
#include <TTree.h>

#include <vector>
#include <string>
#include <set>

using namespace std;
//using namespace edm;

typedef std::vector<Trajectory> TrajectoryCollection;

enum EcalClustID {
  fEBClean=0,
  fEBUnclean,
  fEBCombined,
  fEEClean,
  fEEUnclean,
  fEECombined
};

namespace reco {
  class CaloCluster;
}

class MplTracker{
 public:
  MplTracker(edm::ParameterSet const& iPS, edm::ConsumesCollector && iC);
  ~MplTracker();
  
  void beginJob(TTree *Tree);
  void endJob();
  void Init(const edm::EventSetup&) ;
  void analyze(const edm::Event&, const edm::EventSetup&);
  
  void getTracks(std::vector<Mono::MonoTrack> &) const;
  void doMatch(unsigned,const Mono::MonoEcalCluster *,const Mono::EBmap &);
  void doMatch(unsigned,const reco::CaloCluster **,const EcalClustID);
  
  inline const std::vector<int> & getSubHits() const { return _vSubHits; }   
  inline const std::vector<int> & getSatSubHits() const { return _vSatSubHits; } 
  inline const std::vector<float> & getTIso() const { return _vIso; }
  inline const std::vector<float> & getNdof() const { return _vNdof; }
  
 private:
  int  AddPoints(const reco::Track &Track);
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
  
  edm::EDGetTokenT< reco::TrackCollection > m_TrackTag; 
  edm::EDGetTokenT<std::vector<Trajectory> > m_TrajectoryTag;
  edm::EDGetTokenT<TrajTrackAssociationCollection> m_assoMapTag;
  
  edm::Handle<reco::TrackCollection> _hTracks;
  edm::Handle<TrajectoryCollection> _hTrajectories;
  edm::Handle<TrajTrackAssociationCollection> _hTrajTrackAssociations;
  edm::Handle<edm::ValueMap<reco::DeDxData> > _hDeDx;
  edm::ESHandle<GlobalTrackingGeometry> _TrackingGeom;
  
  map<uint, float> _NormMap;
  
  set<int> _Used;
  
  vector<GlobalPoint> _Points;
  vector<GlobalError> _Errors;
  vector<float>       _Charges, _HighHits, _SumHits;
  
  //TLinearFitter *_RZFitter;
  TF1 *_RZFunc;
  TF1 *_XYFunc;
  //UnbinnedLikelihoodFit _DeDxFitter;
  //TF1 *_DeDxFunc;
  
  float _XYPar[3], _XYErr[3], _RZPar[3], _RZErr[3];
  float _Chi2XY, _Chi2RZ;
  int _Ndof;
  float _Iso;
  
  //TFile *_OutputFile;
  TTree *_Tree;
  vector<float> _vXYPar0, _vXYPar1, _vXYPar2, _vXYErr0, _vXYErr1, _vXYErr2, _vRZPar0, _vRZPar1, _vRZPar2, _vRZErr0, _vRZErr1, _vRZErr2, _vChi2XY, _vChi2RZ, _vNdof, _vIso;
  vector<int> _vHits, _vSatHits, _vSubHits, _vSatSubHits;
  vector<string> _vGroup;
  vector<int> _clustMatchEB; // match to ecal cluster
  vector<int> _clustMatchEBClean; // match to ecal cluster
  vector<int> _clustMatchEBUnclean; // match to ecal cluster
  vector<int> _clustMatchEE; // match to ecal cluster
  vector<int> _clustMatchEEClean; // match to ecal cluster
  vector<int> _clustMatchEEUnclean; // match to ecal cluster
  vector<double> _clustDistEB; // distance to ecal cluster
  vector<double> _clustDistEBClean; // distance to ecal cluster
  vector<double> _clustDistEBUnclean; // distance to ecal cluster
  vector<double> _clustDistEE; // distance to ecal cluster
  vector<double> _clustDistEEClean; // distance to ecal cluster
  vector<double> _clustDistEEUnclean; // distance to ecal cluster
  
  bool _TrackHitOutput;
  //TTree *_TrackHitTree;
  vector<int> _vTHTrack, _vTHStrips, _vTHSatStrips;
  vector<float> _vTHX, _vTHY, _vTHZ, _vTHErrX, _vTHErrY, _vTHErrZ;
  
  bool _FillSelf;
};
#endif
