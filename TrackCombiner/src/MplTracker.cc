#include "Monopoles/TrackCombiner/interface/MplTracker.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

#include "DataFormats/CaloRecHit/interface/CaloCluster.h"

#include "RecoTracker/DeDx/interface/DeDxTools.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"

#include "TFitResult.h"
#include "TVirtualFitter.h"

#include <sstream>

//#define DEBUG

using namespace std; using namespace edm;

/// Constructor
MplTracker::MplTracker(edm::ParameterSet const& parameterSet, edm::ConsumesCollector && iC)
{
  if(CLHEP::electron_charge==0) std::cout << "asdf" << std::endl;
  
  //_Source = parameterSet.getParameter<std::string>("TrackSource");
  m_TrackTag = iC.consumes< reco::TrackCollection >( parameterSet.getParameter< edm::InputTag >( "TrackTag" ) );
  m_TrajectoryTag = iC.consumes< std::vector< Trajectory > >( parameterSet.getParameter<edm::InputTag>("TrackTag") );
  m_assoMapTag = iC.consumes< TrajTrackAssociationCollection >( parameterSet.getParameter<edm::InputTag>("TrackTag") );

  //_Output = parameterSet.getParameter<std::string>("Output");
  _PhiCut = parameterSet.getUntrackedParameter<double>("TrackPhiCut", 0.5);
  _Chi2Cut = parameterSet.getUntrackedParameter<double>("TrackChi2Cut", 5.0);
  _PtCut = parameterSet.getUntrackedParameter<double>("TrackPtCut", 3.0);
  _DeDxCut = parameterSet.getUntrackedParameter<double>("TrackDeDxCut", 5.0);
  _DefaultError = pow(parameterSet.getUntrackedParameter<double>("TrackDefaultError", 0.05), 2);
  _ErrorFudge = pow(parameterSet.getUntrackedParameter<double>("TrackErrorFudge", 0.02), 2);

  _MeVperADCPixel = parameterSet.getUntrackedParameter<double>("TrackMeVperADCPixel", 3.61e-6);
  _MeVperADCStrip = parameterSet.getUntrackedParameter<double>("TrackMeVperADCStrip", 3.61e-6*265);

  _TrackHitOutput = parameterSet.getUntrackedParameter<bool>("TrackHitOutput", false);

  _RZFunc = new TF1("RZFunc", "[0] + [1]*x + [2]*x^2", 0, 200);
  // track starts along x axis, fit to semicircle:
  //_XYFunc = new TF1("XYFunc", "[0] + (sqrt(1 - ((x-[1])*[2])^2)-1)/[2]", 0, 200);
  _XYFunc = new TF1("XYFunc", "[0] + sqrt([2]^2 - (x-[1])^2)*TMath::Sign(1,[2]) - [2]", 0, 200);
  //_DeDxFunc = new TF1("DeDxFunc", "TMath::Landau(x,[0],[1],1)", 0, 255);
  //_DeDxFitter.setFunction(_DeDxFunc);

  _FillSelf = false;
}
 

/// Destructor
MplTracker::~MplTracker(){
}

void MplTracker::beginJob(TTree *Tree=NULL){
  //_OutputFile = new TFile(_Output.c_str(), "recreate");
  if(Tree) _Tree = Tree;
  else{
    _Tree = new TTree("MplTrackSets", "MplTrackSets");
    _FillSelf = true;
  }

  _Tree->Branch("Track_Group", &_vGroup);

  _Tree->Branch("Track_XYPar0", &_vXYPar0);
  _Tree->Branch("Track_XYPar1", &_vXYPar1);
  _Tree->Branch("Track_XYPar2", &_vXYPar2);
  _Tree->Branch("Track_XYErr0", &_vXYErr0);
  _Tree->Branch("Track_XYErr1", &_vXYErr1);
  _Tree->Branch("Track_XYErr2", &_vXYErr2);
  _Tree->Branch("Track_RZPar0", &_vRZPar0);
  _Tree->Branch("Track_RZPar1", &_vRZPar1);
  _Tree->Branch("Track_RZPar2", &_vRZPar2);
  _Tree->Branch("Track_RZErr0", &_vRZErr0);
  _Tree->Branch("Track_RZErr1", &_vRZErr1);
  _Tree->Branch("Track_RZErr2", &_vRZErr2);

  _Tree->Branch("Track_Chi2XY", &_vChi2XY);
  _Tree->Branch("Track_Chi2RZ", &_vChi2RZ);
  _Tree->Branch("Track_Ndof", &_vNdof);

  _Tree->Branch("Track_Hits", &_vHits);
  _Tree->Branch("Track_SatHits", &_vSatHits);
  _Tree->Branch("Track_SubHits", &_vSubHits);
  _Tree->Branch("Track_SatSubHits", &_vSatSubHits);
  _Tree->Branch("Track_Iso", &_vIso);

  _Tree->Branch("Track_clustMatchEB",&_clustMatchEB);
  _Tree->Branch("Track_clustDistEB",&_clustDistEB);
  _Tree->Branch("Track_clustMatchEE",&_clustMatchEE);
  _Tree->Branch("Track_clustDistEE",&_clustDistEE);

  _Tree->Branch("Track_clustMatchEBClean",&_clustMatchEBClean);
  _Tree->Branch("Track_clustDistEBClean",&_clustDistEBClean);
  _Tree->Branch("Track_clustMatchEBUnclean",&_clustMatchEBUnclean);
  _Tree->Branch("Track_clustDistEBUnclean",&_clustDistEBUnclean);

  _Tree->Branch("Track_clustMatchEEClean",&_clustMatchEEClean);
  _Tree->Branch("Track_clustDistEEClean",&_clustDistEEClean);
  _Tree->Branch("Track_clustMatchEEUnclean",&_clustMatchEEUnclean);
  _Tree->Branch("Track_clustDistEEUnclean",&_clustDistEEUnclean);

  if(_TrackHitOutput){
    _Tree->Branch("TrackHit_Track", &_vTHTrack);
    _Tree->Branch("TrackHit_X", &_vTHX);
    _Tree->Branch("TrackHit_Y", &_vTHY);
    _Tree->Branch("TrackHit_Z", &_vTHZ);
    _Tree->Branch("TrackHit_ErrX", &_vTHErrX);
    _Tree->Branch("TrackHit_ErrY", &_vTHErrY);
    _Tree->Branch("TrackHit_ErrZ", &_vTHErrZ);
    _Tree->Branch("TrackHit_Strips", &_vTHStrips);
    _Tree->Branch("TrackHit_SatStrips", &_vTHSatStrips);
  }
}

void MplTracker::endJob(){
  //_OutputFile->cd();
  _Tree->Write();
  //if(_TrackHitOutput) _TrackHitTree->Write();
  //_OutputFile->Close();
}

void MplTracker::Init(const edm::EventSetup& iSetup) {
  
#ifdef DEBUG
  std::cout<<"Start: MplTracker::Init"<<std::endl;
#endif

  iSetup.get<GlobalTrackingGeometryRecord>().get(_TrackingGeom);

  // Fill Normalization map
  edm::ESHandle<TrackerGeometry> tkGeom;
  iSetup.get<TrackerDigiGeometryRecord>().get( tkGeom );

  const vector<const GeomDet*> Det = tkGeom->dets();
  for(uint i=0; i<Det.size(); i++){
    DetId Detid = Det[i]->geographicalId();

    const StripGeomDetUnit* StripDetUnit = dynamic_cast<const StripGeomDetUnit*> (Det[i]);
    const PixelGeomDetUnit* PixelDetUnit = dynamic_cast<const PixelGeomDetUnit*> (Det[i]);

    if(StripDetUnit){
      _NormMap[Detid.rawId()] = _MeVperADCStrip / StripDetUnit->surface().bounds().thickness();
    }else if(PixelDetUnit){
      _NormMap[Detid.rawId()] = _MeVperADCPixel / PixelDetUnit->surface().bounds().thickness();
    }
  }

#ifdef DEBUG
  std::cout<<"End: MplTracker::Init"<<std::endl;
#endif
}

void MplTracker::analyze(const Event& event, const EventSetup& setup){

#ifdef DEBUG
  std::cout<<"Start: MplTracker::analyze"<<std::endl;
#endif

  if(_NormMap.size() == 0) Init(setup);

  //event.getByLabel(_Source, _hTracks);
  event.getByToken(m_TrackTag, _hTracks);
  event.getByToken(m_TrajectoryTag, _hTrajectories);
  event.getByToken(m_assoMapTag, _hTrajTrackAssociations);
  //event.getByLabel("dedxHarmonic2", _hDeDx);

#ifdef DEBUG
  std::cout<<"   + Pass:getByToken for new track collection"<<endl;
#endif

  _Used.clear();
  Clear();

  // Loop over the tracks
  for(uint i = 0; i!=_hTracks->size(); i++){
    //for(std::vector<Trajectory>::const_iterator Traj = _hTrajectories->begin(); Traj!=_hTrajectories->end(); Traj++){
    //edm::Ref<std::vector<Trajectory> > TrajRef(_hTrajectories, i);
    //edm::RefToBase<reco::Track> TrackRef ( (*_hTrajTrackAssociations.product())[TrajRef] );
    edm::Ref<std::vector<reco::Track> > TrackRef(_hTracks, i);

    if(TrackRef->pt() < _PtCut){
#ifdef DEBUG
      cout << "Track " << i << " Pt " << TrackRef->pt() << " too low." << endl;
#endif
      continue;
    }

    // Output Track hits:
    for (trackingRecHit_iterator iHit=TrackRef->recHitsBegin(); iHit!=TrackRef->recHitsEnd(); iHit++){
      _vTHTrack.push_back(i);

      const TrackingRecHit *Hit = *iHit;

      if(!Hit->isValid()) continue;

      LocalPoint LPos = Hit->localPosition();
      LocalError LErr = Hit->localPositionError();

      const GeomDet *Detector = _TrackingGeom->idToDet(Hit->geographicalId());

      GlobalPoint GPos = Detector->toGlobal(LPos);

      _vTHX.push_back(GPos.x());
      _vTHY.push_back(GPos.y());
      _vTHZ.push_back(GPos.z());

      if(LErr.valid()){
        GlobalError GErr = ErrorFrameTransformer::transform( LErr, Detector->surface() );
        //GlobalError GErr = GlobalError(LErr.xx(), LErr.xy(), LErr.yy(), 0, 0, 0);
        //GErr = GlobalError(Detector->GErr.matrix());
        _vTHErrX.push_back(sqrt(GErr.cxx()));
        _vTHErrY.push_back(sqrt(GErr.cyy()));
        _vTHErrZ.push_back(sqrt(GErr.czz()));
      }else{
        _vTHErrX.push_back(0);
        _vTHErrY.push_back(0);
        _vTHErrZ.push_back(0);
      }

      int Strips=0, SatStrips=0;
      if(const SiStripMatchedRecHit2D* matchedHit=dynamic_cast<const SiStripMatchedRecHit2D*>(Hit)){
        const vector<uint8_t>& Ampls = matchedHit->stereoCluster().amplitudes();
        Strips += Ampls.size();
        for(uint i=0; i<Ampls.size(); i++){
	  if(Ampls[i] >= 254) SatStrips++;
        }
      }else if(const ProjectedSiStripRecHit2D* projectedHit=dynamic_cast<const ProjectedSiStripRecHit2D*>(Hit)) {
        auto OrigHit=projectedHit->originalHit();
        const vector<uint8_t>& Ampls = OrigHit.stripCluster().amplitudes();
        Strips += Ampls.size();
        for(uint i=0; i<Ampls.size(); i++){
	  if(Ampls[i] >= 254) SatStrips++;
        }
      }else if(const SiStripRecHit2D* singleHit=dynamic_cast<const SiStripRecHit2D*>(Hit)){
        const vector<uint8_t>& Ampls = singleHit->stripCluster().amplitudes();
        Strips += Ampls.size();
        for(uint i=0; i<Ampls.size(); i++){
	  if(Ampls[i] >= 254) SatStrips++;
        }
      }else if(const SiStripRecHit1D* single1DHit=dynamic_cast<const SiStripRecHit1D*>(Hit)){
        const vector<uint8_t>& Ampls = single1DHit->stripCluster().amplitudes();
        Strips += Ampls.size();
        for(uint i=0; i<Ampls.size(); i++){
	  if(Ampls[i] >= 254) SatStrips++;
        }
// don't use pixels for now (since the standard algorithms don't use them)
      //}else if(const SiPixelRecHit* pixelHit=dynamic_cast<const SiPixelRecHit*>(Hit)){ 
//      Charge = pixelHit->cluster()->charge();
      }
      //Added dedx cut:
      if(SatStrips>=18 && SatStrips>=Strips-5) continue;

      _vTHStrips.push_back(Strips);
      _vTHSatStrips.push_back(SatStrips);
    }

    // skip this if it's already used:
    if (_Used.count(i) > 0) continue;

/*    float DeDx = (*_hDeDx.product())[TrackRef].dEdx();

    if( DeDx > 0 && DeDx < _DeDxCut){
      //cout << "Track " << i << " DeDx too low." << endl;
      continue;
    }*/

    vector<int> Group;
    Group.push_back(i);

    _Points.clear();
    _Errors.clear();
    _Charges.clear();
    _SumHits.clear();
    _HighHits.clear();

    AddPoints(*TrackRef);

    AddMoreTracks(Group);

    FitXY(Group);
    FitRZ();
    //FitDeDx();
    AverageIso(Group);

    //cout << "Just Fitted.  Sizes: " << _Points.size() << " " << _Errors.size() << " " << _Charges.size() << endl;


    Save(Group);

    for (uint j=0; j<Group.size(); j++)
      _Used.insert(Group[j]);
  }

  if(_FillSelf) _Tree->Fill();
  //if(_TrackHitOutput) _TrackHitTree->Fill();

#ifdef DEBUG
  std::cout<<"End: MplTracker::analyze"<<std::endl;
#endif
}

void MplTracker::AddMoreTracks(vector<int> &Group){
  edm::Ref<std::vector<reco::Track> > InitTrack(_hTracks, Group[0]);
  //edm::RefToBase<reco::Track> InitTrack ( (*_hTrajTrackAssociations.product())[InitTrajRef] );

  //cout << "Checking: " << Group[0] << " " << InitTrack->phi() << endl;

  for(uint i = Group.back()+1; i!=_hTracks->size(); i++){
    if(_Used.count(i) > 0) continue;

    edm::Ref<std::vector<reco::Track> > ThisTrack(_hTracks, i);
    //edm::RefToBase<reco::Track> ThisTrack ( (*_hTrajTrackAssociations.product())[ThisTrajRef] );

    if(ThisTrack->pt() < _PtCut){
      //cout << " Pt too low: " << ThisTrack->pt() << endl;
      continue;
    }

/*    float DeDx = (*_hDeDx.product())[ThisTrack].dEdx();

    if( DeDx > 0 && DeDx < _DeDxCut){
      //cout << "Track " << i << " DeDx too low." << endl;
      continue;
    }*/


    //cout << " " << i << " " << ThisTrack->phi() << endl;

    if(fabs(ThisTrack->phi() - InitTrack->phi()) > _PhiCut) continue;

    //cout << " Passed Phi Cut" << endl;

    int NPoints = AddPoints(*ThisTrack);
    if(NPoints == 0) continue;

    FitXY(Group);
    //cout << " Chi2XY: " << _Chi2XY / _NdofXY << endl;
    if(_Ndof > 0 && _Chi2XY / _Ndof > _Chi2Cut){
      RemovePoints(NPoints);
      continue;
    }

    FitRZ();
    //cout << " Chi2RZ: " << _Chi2RZ / _NdofRZ << endl;
    if(_Ndof > 0 && _Chi2RZ / _Ndof > _Chi2Cut){
      RemovePoints(NPoints);
      continue;
    }

    //cout << " Good!" << endl;
    Group.push_back(i);
  }
}


int MplTracker::AddPoints(const reco::Track &Track){
  int NPoints = 0;

  // Find the right Trajectory
  const Trajectory *Traj=NULL;
  for(uint i=0; i < _hTrajectories->size(); i++){
    edm::Ref<std::vector<Trajectory> > TrajRef(_hTrajectories, i);
    edm::RefToBase<reco::Track> TrackRef ( (*_hTrajTrackAssociations.product())[TrajRef] );
    if(TrackRef->pt() == Track.pt() && TrackRef->eta() == Track.eta() && TrackRef->phi() == Track.phi()){
      Traj = TrajRef.get();
      break;
    }
  }
  if (Traj == NULL) cout << "Incoming!!!" << endl;

  for (trackingRecHit_iterator iHit=Track.recHitsBegin(); iHit!=Track.recHitsEnd(); iHit++){

    const TrackingRecHit *Hit = *iHit;

    if(!Hit->isValid()) continue;

    LocalPoint LPos = Hit->localPosition();
    LocalError LErr = Hit->localPositionError();

    const GeomDet *Detector = _TrackingGeom->idToDet(Hit->geographicalId());

    GlobalPoint GPos = Detector->toGlobal(LPos);

    _Points.push_back(GPos);

    double r2 = GPos.perp2();

    if(LErr.valid()){
      GlobalError GErr = ErrorFrameTransformer::transform( LErr, Detector->surface() );
      //GlobalError GErr = GlobalError(LErr.xx(), LErr.xy(), LErr.yy(), 0, 0, 0);
      //GErr = GlobalError(Detector->GErr.matrix());

      // if errors are huge, treat them as invalid
      if(GErr.cxx() > 100 || GErr.cyy() > 100 || GErr.czz() > 100)
	_Errors.push_back(GlobalError(_DefaultError*r2, 0, _DefaultError*r2, 0, 0, _DefaultError*r2) );
      else
	_Errors.push_back(GErr + GlobalError(_ErrorFudge*r2, 0, _ErrorFudge*r2, 0, 0, _ErrorFudge*r2) ); //error bars are too small.
    }else{
      _Errors.push_back(GlobalError(_DefaultError*r2, 0, _DefaultError*r2, 0, 0, _DefaultError*r2) );
    }

    // add the dedx information: different accessor for every type of hit
    TrajectoryStateOnSurface State(Traj->geometricalInnermostState().globalParameters(), Detector->surface());
    LocalVector Direction = State.localDirection();
    double cosine = Direction.z()/Direction.mag();
    float Charge = 0;
    int HighHits = 0, SumHits = 0;

    if(const SiStripMatchedRecHit2D* matchedHit=dynamic_cast<const SiStripMatchedRecHit2D*>(Hit)){
      const vector<uint8_t>& Ampls = matchedHit->stereoCluster().amplitudes();
      SumHits += Ampls.size();
      for(uint i=0; i<Ampls.size(); i++){
	Charge += Ampls[i];
	if(Ampls[i] >= 254) HighHits++;
      }
    }else if(const ProjectedSiStripRecHit2D* projectedHit=dynamic_cast<const ProjectedSiStripRecHit2D*>(Hit)) {
      auto OrigHit=projectedHit->originalHit();
      const vector<uint8_t>& Ampls = OrigHit.stripCluster().amplitudes();
      SumHits += Ampls.size();
      for(uint i=0; i<Ampls.size(); i++){
	Charge += Ampls[i];
	if(Ampls[i] >= 254) HighHits++;
      }
    }else if(const SiStripRecHit2D* singleHit=dynamic_cast<const SiStripRecHit2D*>(Hit)){
      const vector<uint8_t>& Ampls = singleHit->stripCluster().amplitudes();
      SumHits += Ampls.size();
      for(uint i=0; i<Ampls.size(); i++){
	Charge += Ampls[i];
	if(Ampls[i] >= 254) HighHits++;
      }
    }else if(const SiStripRecHit1D* single1DHit=dynamic_cast<const SiStripRecHit1D*>(Hit)){
      const vector<uint8_t>& Ampls = single1DHit->stripCluster().amplitudes();
      SumHits += Ampls.size();
      for(uint i=0; i<Ampls.size(); i++){
	Charge += Ampls[i];
	if(Ampls[i] >= 254) HighHits++;
      }
// don't use pixels for now (since the standard algorithms don't use them)
    //}else if(const SiPixelRecHit* pixelHit=dynamic_cast<const SiPixelRecHit*>(Hit)){ 
//      Charge = pixelHit->cluster()->charge();
      //Charge = -1;
    }

    //Added dedx cut:
    if(HighHits>=18 && HighHits>=SumHits-5) continue;


    //cout << "Adding Charge: " << Charge << " * " << abs(cosine) << " * " << _NormMap[Hit->geographicalId().rawId()] << endl;

    float NormCharge = _NormMap[Hit->geographicalId().rawId()] * Charge * abs(cosine);

    if(NormCharge > 0 && NormCharge < _DeDxCut){
      _Points.pop_back();
      _Errors.pop_back();
      continue;
    }

    _Charges.push_back(NormCharge);
    _SumHits.push_back(SumHits);
    _HighHits.push_back(HighHits);

    NPoints++;
  }

  //cout << "Just Added " << NPoints << ".  Sizes: " << _Points.size() << " " << _Errors.size() << " " << _Charges.size() << endl;

  return NPoints;
}

void MplTracker::RemovePoints(int n){
  for(int i=0; i<n; i++){
    _Points.pop_back();
    _Errors.pop_back();
    _Charges.pop_back();
    _SumHits.pop_back();
    _HighHits.pop_back();
  }

  //cout << "Just Removed " << n << ".  Sizes: " << _Points.size() << " " << _Errors.size() << " " << _Charges.size() << endl;
}

void MplTracker::Clear(){
  _vGroup.clear();

  _vXYPar0.clear();
  _vXYPar1.clear();
  _vXYPar2.clear();
  _vXYErr0.clear();
  _vXYErr1.clear();
  _vXYErr2.clear();
  _vRZPar0.clear();
  _vRZPar1.clear();
  _vRZPar2.clear();
  _vRZErr0.clear();
  _vRZErr1.clear();
  _vRZErr2.clear();

  _vChi2XY.clear();
  _vChi2RZ.clear();
  _vNdof.clear();

  _vHits.clear();
  _vSatHits.clear();
  _vSubHits.clear();
  _vSatSubHits.clear();
  _vIso.clear();

  _vTHTrack.clear();
  _vTHX.clear();
  _vTHY.clear();
  _vTHZ.clear();
  _vTHErrX.clear();
  _vTHErrY.clear();
  _vTHErrZ.clear();

  _clustMatchEB.clear();
  _clustDistEB.clear();
  _clustMatchEBClean.clear();
  _clustDistEBClean.clear();
  _clustMatchEBUnclean.clear();
  _clustDistEBUnclean.clear();

  _clustMatchEE.clear();
  _clustDistEE.clear();
  _clustMatchEEClean.clear();
  _clustDistEEClean.clear();
  _clustMatchEEUnclean.clear();
  _clustDistEEUnclean.clear();
}

void MplTracker::Save(vector<int> &Group){
  /*cout << "Outputting a track..." << endl << "Constituents: ";
  for(uint i=0; i<Group.size(); i++){
    cout << Group[i] << " ";
  }
  cout << endl;
  cout << "XYPars: " << _XYPar[0] << " +/- " << _XYErr[0] << " " << _XYPar[1] << " +/- " << _XYErr[1] << " " << _XYPar[2] << " +/- " << _XYErr[2] << endl;
  cout << "XYChi2: " << _Chi2XY << " / " << _NdofXY << endl;
  cout << "RZPars: " << _RZPar[0] << " +/- " << _RZErr[0] << " " << _RZPar[1] << " +/- " << _RZErr[1] << " " << _RZPar[2] << " +/- " << _RZErr[2] << endl;
  cout << "RZChi2: " << _Chi2RZ << " / " << _NdofRZ << endl;*/

  ostringstream csv;
  csv << Group[0];
  for(uint i=1; i<Group.size(); i++) csv << "," << Group[i];
  _vGroup.push_back(csv.str());

  _vXYPar0.push_back(_XYPar[0]);
  _vXYPar1.push_back(_XYPar[1]);
  _vXYPar2.push_back(_XYPar[2]);
  _vXYErr0.push_back(_XYErr[0]);
  _vXYErr1.push_back(_XYErr[1]);
  _vXYErr2.push_back(_XYErr[2]);
  _vRZPar0.push_back(_RZPar[0]);
  _vRZPar1.push_back(_RZPar[1]);
  _vRZPar2.push_back(_RZPar[2]);
  _vRZErr0.push_back(_RZErr[0]);
  _vRZErr1.push_back(_RZErr[1]);
  _vRZErr2.push_back(_RZErr[2]);

  _vChi2XY.push_back(_Chi2XY);
  _vChi2RZ.push_back(_Chi2RZ);
  _vNdof.push_back(_Ndof);

  _vIso.push_back(_Iso);

  int Hits = 0, SatHits=0, SubHits = 0, SatSubHits=0;

  for(uint i=0; i<_SumHits.size(); i++){
    Hits++;
    if(2*_HighHits[i]>=_SumHits[i]) SatHits++;

    SubHits += _SumHits[i];
    SatSubHits += _HighHits[i];
  }

  _vHits.push_back(Hits);
  _vSatHits.push_back(SatHits);
  _vSubHits.push_back(SubHits);
  _vSatSubHits.push_back(SatSubHits);

  //for(uint i=0; i<_Charges.size(); i++)
    //cout << "Charge: " << _Charges[i] << endl;
}

void MplTracker::FitXY(vector<int> &Group){
  int NumPoints = _Points.size();
  TGraphErrors XYGraph(NumPoints);

  float AvePt = 0;
  for(uint i=0; i<Group.size(); i++){
    edm::Ref<std::vector<reco::Track> > ThisTrack(_hTracks, Group[i]);
    AvePt += ThisTrack->pt() * ThisTrack->charge();
  }
  AvePt /= Group.size();

  //cout << endl;

  float XMax = 0;

  // rotate so the initial path is along the x axis
  float RStart = _Points[0].perp2();
  float REnd = _Points[NumPoints-1].perp2();
  float Phi0;
  if(RStart > REnd) Phi0 = _Points[0].phi();
  else Phi0 = _Points[NumPoints-1].phi();

  for(int i=0; i<NumPoints; i++){
    //cout << "XYFitter: " << Points[i].x() << " " << Points[i].y() << " "
    // << sqrt(Errors[i].cxx()) << " " << sqrt(Errors[i].cyy()) << endl;
    float NewX = _Points[i].x()*cos(Phi0)+_Points[i].y()*sin(Phi0);
    float NewY = -_Points[i].x()*sin(Phi0)+_Points[i].y()*cos(Phi0);
    XYGraph.SetPoint(i, NewX, NewY);
    float ErrorX = sqrt(_Errors[i].cxx()*cos(Phi0)*cos(Phi0) + _Errors[i].cyy()*sin(Phi0)*sin(Phi0) + 2*_Errors[i].cyx()*sin(Phi0)*cos(Phi0));
    float ErrorY = sqrt(_Errors[i].cxx()*sin(Phi0)*sin(Phi0) + _Errors[i].cyy()*cos(Phi0)*cos(Phi0) - 2*_Errors[i].cyx()*sin(Phi0)*cos(Phi0));
    XYGraph.SetPointError(i, ErrorX, ErrorY);
    //cout << "Trans: " << NewX << " " << NewY << " " << ErrorX << " " << ErrorY << endl;
    if(NewX + ErrorX > XMax) XMax = NewX + ErrorX;
  }

  _XYFunc->SetParameters(1, 1, AvePt/0.0114);
  //_XYFunc->SetParLimits(2, -XMax, XMax);
  TFitResultPtr Result = XYGraph.Fit(_XYFunc, "Q S B");

  _XYPar[0] = sqrt(pow(Result->Parameter(0)-Result->Parameter(2),2)+pow(Result->Parameter(1),2)) - fabs(Result->Parameter(2)); // unsigned d0
  _XYErr[0] = Result->ParError(0); // not correct right now
  _XYPar[1] = Phi0 - atan(Result->Parameter(1)/(Result->Parameter(2)-Result->Parameter(0))); // phi0
  _XYErr[1] = Result->ParError(1); // not correct right now
  _XYPar[2] = Result->Parameter(2); // radius
  _XYErr[2] = Result->ParError(2);
  _Chi2XY = Result->Chi2();
  _Ndof = Result->Ndf();

#ifdef DEBUG
  cout << "DEBUG: " << _XYPar[0] << " " << _XYPar[1] << " " << _XYPar[2] << endl;
#endif
}

void MplTracker::FitRZ(bool Debug){ /* _RZFitter->ClearPoints();

  int NumPoints = Points.size();

  cout << endl;

  for(int i=0; i<NumPoints; i++){
    double x = sqrt(Points[i].x()*Points[i].x() + Points[i].y()*Points[i].y());
    cout << "RZFitter: " << Points[i].x() << " " << Points[i].y() << " " << Points[i].z() << " "
			 << sqrt(Errors[i].cxx()) << " " << sqrt(Errors[i].cyy()) << " " << sqrt(Errors[i].czz()) << endl;
    _RZFitter->AddPoint(&x, Points[i].z(), sqrt(Errors[i].czz()) );
  }

  _RZFitter->Eval();

  for(int i=0; i<3; i++){
    _RZPar[i] = _RZFitter->GetParameter(i);
    _RZErr[i] = _RZFitter->GetParError(i);
  }
  _Chi2RZ = _RZFitter->GetChisquare();
  _NdofRZ = NumPoints - 3;
  if(_NdofRZ <= 0) _NdofRZ = 1; */
  int NumPoints = _Points.size();
  TGraphErrors RZGraph(NumPoints);

  if(Debug) cout << endl;

  for(int i=0; i<NumPoints; i++){
    if(Debug){
      cout << _Points[i].x() << " " << _Points[i].y() << " " << _Points[i].z() << " "
	   << sqrt(_Errors[i].cxx()) << " " << sqrt(_Errors[i].cyy()) << " " << sqrt(_Errors[i].czz()) << endl;
    }
    RZGraph.SetPoint(i, sqrt(_Points[i].x()*_Points[i].x() + _Points[i].y()*_Points[i].y()), _Points[i].z());
    RZGraph.SetPointError(i, sqrt(_Errors[i].rerr(_Points[i])), sqrt(_Errors[i].czz()));
  }

  _RZFunc->SetParameters(0,0,0); //1,-1,0.01);
  TFitResultPtr Result = RZGraph.Fit(_RZFunc, "Q S");

  for(int i=0; i<3; i++){
    _RZPar[i] = Result->Parameter(i);
    _RZErr[i] = Result->ParError(i);
  }
  _Chi2RZ = Result->Chi2();
  //_NdofRZ = Result->Ndf();

}

void MplTracker::FitDeDx(){
  //use the UnbinnedFit method for now.
/*
  int NumPoints = _Charges.size();
  if(NumPoints > 100) NumPoints=100;

  cout << endl;

  // ignore negative charges (they're from ignored hits)
  double x[100];
  int iOld=0, iNew=0;
  for(iOld=0; iOld<NumPoints; iOld++){
    if(_Charges[iOld] < 0) continue;
    x[iNew] = _Charges[iOld];
    cout << "Charges: " << x[iNew] << endl;
    iNew++;
  }
  _DeDxFitter.setData(iNew, x);
  _DeDxFunc->SetParameters(3, 0.3);
  _DeDxFitter.fit();

  _DeDx = _DeDxFunc->GetParameter(0);*/

  // use the median DeDx:

  /*vector<float> SortedCharges;
  for(uint i=0; i<_Charges.size(); i++)
    if(_Charges[i]>0) SortedCharges.push_back(_Charges[i]);

  if(SortedCharges.size() == 0){
    _DeDx = 0;
    return;
  }

  sort(SortedCharges.begin(), SortedCharges.end());

  _DeDx = SortedCharges[SortedCharges.size()/2];
  */
  /*int SumHits1 = 0, SumHits2=0, HighHits1 = 0, HighHits2=0;

  for(uint i=0; i<_SumHits.size(); i++){
    SumHits1 += _SumHits[i];
    HighHits1 += _HighHits[i];
    SumHits2++;
    if(_HighHits[i]>0) HighHits2++;
  }
  _HighDeDx1 = (float)HighHits1/SumHits1;
  _HighDeDx2 = (float)HighHits2/SumHits2;
*/
}

void MplTracker::AverageIso(vector<int> &Group){
  edm::Ref<std::vector<reco::Track> > InitTrack(_hTracks, Group[0]);

  float InitEta = InitTrack->eta();
  float InitPhi = InitTrack->phi();
  float IsoPt = 0;

  for(uint i = 0; i!=_hTracks->size(); i++){
    if(count(Group.begin(), Group.end(), i) > 0) continue;

    edm::Ref<std::vector<reco::Track> > ThisTrack(_hTracks, i);

    //float dR2 = pow(ThisTrack->eta() - InitEta, 2) + pow(ThisTrack->phi() - InitPhi, 2);
    float dR2 = reco::deltaR(ThisTrack->eta(), ThisTrack->phi(), InitEta, InitPhi);

    if(dR2 > .16) continue; // dR cut of 0.4

    IsoPt += ThisTrack->pt();
  }

  _Iso = IsoPt;
}


void MplTracker::getTracks(std::vector<Mono::MonoTrack> &tracks) const
{
  const unsigned nTracks = _vXYPar0.size();
  tracks.resize(nTracks);

  for ( unsigned t=0; t != nTracks; t++ ) {
    tracks[t] = Mono::MonoTrack(_vXYPar0[t],_vXYPar1[t],_vXYPar2[t],_vRZPar0[t],_vRZPar1[t],_vRZPar2[t]);
  }

}


void MplTracker::doMatch(unsigned nClusters, const Mono::MonoEcalCluster *clusters,const Mono::EBmap &ecalMap)
{

  Mono::MonoTrackMatcher matcher(50.);

  std::vector<Mono::MonoTrack> tracks;
  this->getTracks(tracks);

  const unsigned nTracks = tracks.size();
  matcher.match(nClusters,clusters,ecalMap,nTracks,&tracks[0],_clustMatchEB,_clustDistEB);
 

}



void MplTracker::doMatch(unsigned nClusters, const reco::CaloCluster **clusters,const EcalClustID id=fEBCombined)
{

  Mono::MonoTrackMatcher matcher(50.);

  std::vector<Mono::MonoTrack> tracks;
  this->getTracks(tracks);

  const unsigned nTracks = tracks.size();
  if ( id == fEBCombined )
    matcher.match(nClusters,clusters,nTracks,&tracks[0],_clustMatchEB,_clustDistEB,true);
  else if ( id == fEBClean )
    matcher.match(nClusters,clusters,nTracks,&tracks[0],_clustMatchEBClean,_clustDistEBClean,true);
  else if ( id == fEBUnclean )
    matcher.match(nClusters,clusters,nTracks,&tracks[0],_clustMatchEBUnclean,_clustDistEBUnclean,true);
  else if ( id == fEECombined )
    matcher.match(nClusters,clusters,nTracks,&tracks[0],_clustMatchEE,_clustDistEE,false);
  else if ( id == fEEClean )
    matcher.match(nClusters,clusters,nTracks,&tracks[0],_clustMatchEEClean,_clustDistEEClean,false);
  else if ( id == fEEUnclean )
    matcher.match(nClusters,clusters,nTracks,&tracks[0],_clustMatchEEUnclean,_clustDistEEUnclean,false);
 

}
