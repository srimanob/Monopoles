#include "Monopoles/TrackCombiner/interface/TrackCombinerReco.h"

//#include "TrackingTools/TrackRefitter/interface/TrackTransformer.h"
//#include "TrackingTools/TrackRefitter/interface/TrackTransformerForGlobalCosmicMuons.h" 
//#include "TrackingTools/TrackRefitter/interface/TrackTransformerForCosmicMuons.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Common/interface/Handle.h"

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

#include "RecoTracker/DeDx/interface/DeDxTools.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"

#include <TFitResult.h>
#include <TVirtualFitter.h>

#include <sstream>

using namespace std; using namespace edm;

/// Constructor
TrackCombinerReco::TrackCombinerReco(const ParameterSet& parameterSet){
  _Source = parameterSet.getParameter<std::string>("Source");
  _Output = parameterSet.getParameter<std::string>("Output");
  _PhiCut = parameterSet.getUntrackedParameter<double>("PhiCut", 0.5);
  _Chi2Cut = parameterSet.getUntrackedParameter<double>("Chi2Cut", 5.0);
  _PtCut = parameterSet.getUntrackedParameter<double>("PtCut", 3.0);
  _DeDxCut = parameterSet.getUntrackedParameter<double>("DeDxCut", 5.0);
  _DefaultError = pow(parameterSet.getUntrackedParameter<double>("DefaultError", 0.05), 2);
  _ErrorFudge = pow(parameterSet.getUntrackedParameter<double>("ErrorFudge", 0.02), 2);

  _MeVperADCPixel = parameterSet.getUntrackedParameter<double>("MeVperADCPixel", 3.61e-6);
  _MeVperADCStrip = parameterSet.getUntrackedParameter<double>("MeVperADCStrip", 3.61e-6*265);

  _TrackHitOutput = parameterSet.getUntrackedParameter<bool>("TrackHitOutput", false);

  _RZFunc = new TF1("RZFunc", "[0] + [1]*x + [2]*x^2", 0, 200);
  // track starts along x axis, fit to semicircle:
  //_XYFunc = new TF1("XYFunc", "[0] + (sqrt(1 - ((x-[1])*[2])^2)-1)/[2]", 0, 200);
  _XYFunc = new TF1("XYFunc", "[0] + sqrt([2]^2 - (x-[1])^2)*TMath::Sign(1,[2]) - [2]", 0, 200);
  //_DeDxFunc = new TF1("DeDxFunc", "TMath::Landau(x,[0],[1],1)", 0, 255);
  //_DeDxFitter.setFunction(_DeDxFunc);
}
 

/// Destructor
TrackCombinerReco::~TrackCombinerReco(){
}

void TrackCombinerReco::beginJob(){
  _OutputFile = new TFile(_Output.c_str(), "recreate");
  _Tree = new TTree("MplTrackSets", "MplTrackSets");
  _Tree->Branch("Group", &_vGroup);

  _Tree->Branch("XYPar0", &_vXYPar0);
  _Tree->Branch("XYPar1", &_vXYPar1);
  _Tree->Branch("XYPar2", &_vXYPar2);
  _Tree->Branch("XYErr0", &_vXYErr0);
  _Tree->Branch("XYErr1", &_vXYErr1);
  _Tree->Branch("XYErr2", &_vXYErr2);
  _Tree->Branch("RZPar0", &_vRZPar0);
  _Tree->Branch("RZPar1", &_vRZPar1);
  _Tree->Branch("RZPar2", &_vRZPar2);
  _Tree->Branch("RZErr0", &_vRZErr0);
  _Tree->Branch("RZErr1", &_vRZErr1);
  _Tree->Branch("RZErr2", &_vRZErr2);

  _Tree->Branch("Chi2XY", &_vChi2XY);
  _Tree->Branch("Chi2RZ", &_vChi2RZ);
  _Tree->Branch("NdofXY", &_vNdofXY);
  _Tree->Branch("NdofRZ", &_vNdofRZ);

  _Tree->Branch("DeDx", &_vDeDx);
  _Tree->Branch("Iso", &_vIso);

  if(_TrackHitOutput){
    _TrackHitTree = new TTree("TrackHits", "TrackHits");
    _TrackHitTree->Branch("Track", &_vTHTrack);
    _TrackHitTree->Branch("X", &_vTHX);
    _TrackHitTree->Branch("Y", &_vTHY);
    _TrackHitTree->Branch("Z", &_vTHZ);
    _TrackHitTree->Branch("ErrX", &_vTHErrX);
    _TrackHitTree->Branch("ErrY", &_vTHErrY);
    _TrackHitTree->Branch("ErrZ", &_vTHErrZ);
  }
}

void TrackCombinerReco::endJob(){
  _OutputFile->cd();
  _Tree->Write();
  if(_TrackHitOutput) _TrackHitTree->Write();
  _OutputFile->Close();
}

void TrackCombinerReco::Init(const edm::EventSetup& iSetup) {
  iSetup.get<GlobalTrackingGeometryRecord>().get(_TrackingGeom);

  // Fill Normalization map
  edm::ESHandle<TrackerGeometry> tkGeom;
  iSetup.get<TrackerDigiGeometryRecord>().get( tkGeom );

  const std::vector<const GeomDet*> Det = tkGeom->dets();
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
}

void TrackCombinerReco::analyze(const Event& event, const EventSetup& setup){
  if(_NormMap.size() == 0) Init(setup);

  event.getByLabel(_Source, _hTracks);
  //event.getByLabel(_Source, _hTrajectories);
  //event.getByLabel(_Source, _hTrajTrackAssociations);

  event.getByLabel("dedxHarmonic2", _hDeDx);

  _Used.clear();
  Clear();

  // Loop over the trajectories
  for(uint i = 0; i!=_hTracks->size(); i++){
  //for(std::vector<Trajectory>::const_iterator Traj = _hTrajectories->begin(); Traj!=_hTrajectories->end(); Traj++){

    //edm::Ref<std::vector<Trajectory> > TrajRef(_hTrajectories, i);
    //edm::RefToBase<reco::Track> TrackRef ( (*_hTrajTrackAssociations.product())[TrajRef] );
    edm::Ref<std::vector<reco::Track> > TrackRef(_hTracks, i);

    if(TrackRef->pt() < _PtCut){
      //cout << "Track " << i << " Pt " << TrackRef->pt() << " too low." << endl;
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
      _vTHZ.push_back(GPos.y());

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

    AddPoints(*TrackRef);

    AddMoreTracks(Group);

    FitXY(Group);
    FitRZ();
    FitDeDx();
    AverageIso(Group);

    //cout << "Just Fitted.  Sizes: " << _Points.size() << " " << _Errors.size() << " " << _Charges.size() << endl;


    Save(Group);

    for (uint j=0; j<Group.size(); j++)
      _Used.insert(Group[j]);
  }

  _Tree->Fill();
  if(_TrackHitOutput) _TrackHitTree->Fill();
}

void TrackCombinerReco::AddMoreTracks(vector<int> &Group){
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
    if(_NdofXY > 0 && _Chi2XY / _NdofXY > _Chi2Cut){
      RemovePoints(NPoints);
      continue;
    }

    FitRZ();
    //cout << " Chi2RZ: " << _Chi2RZ / _NdofRZ << endl;
    if(_NdofRZ > 0 && _Chi2RZ / _NdofRZ > _Chi2Cut){
      RemovePoints(NPoints);
      continue;
    }

    //cout << " Good!" << endl;
    Group.push_back(i);
  }
}


int TrackCombinerReco::AddPoints(const reco::Track &Track){
  int NPoints = 0;

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
    //TrajectoryStateOnSurface State=Meas->updatedState();
    ///LocalVector Direction = State.localDirection();
    double cosine = 1.0; //Direction.z()/Direction.mag();
    float Charge = 0;

    if(const SiStripMatchedRecHit2D* matchedHit=dynamic_cast<const SiStripMatchedRecHit2D*>(Hit)){
      const vector<uint8_t>& Ampls = matchedHit->stereoCluster().amplitudes();
      for(uint i=0; i<Ampls.size(); i++) Charge += Ampls[i];
    }else if(const ProjectedSiStripRecHit2D* projectedHit=dynamic_cast<const ProjectedSiStripRecHit2D*>(Hit)) {
      auto OrigHit=projectedHit->originalHit();
      const vector<uint8_t>& Ampls = OrigHit.stripCluster().amplitudes();
      for(uint i=0; i<Ampls.size(); i++) Charge += Ampls[i];
    }else if(const SiStripRecHit2D* singleHit=dynamic_cast<const SiStripRecHit2D*>(Hit)){
      const vector<uint8_t>& Ampls = singleHit->stripCluster().amplitudes();
      for(uint i=0; i<Ampls.size(); i++) Charge += Ampls[i];
    }else if(const SiStripRecHit1D* single1DHit=dynamic_cast<const SiStripRecHit1D*>(Hit)){
      const vector<uint8_t>& Ampls = single1DHit->stripCluster().amplitudes();
      for(uint i=0; i<Ampls.size(); i++) Charge += Ampls[i];
// don't use pixels for now (since the standard algorithms don't use them)
    //}else if(const SiPixelRecHit* pixelHit=dynamic_cast<const SiPixelRecHit*>(Hit)){ 
//      Charge = pixelHit->cluster()->charge();
      //Charge = -1;
    }

    //cout << "Adding Charge: " << Charge << " * " << abs(cosine) << " * " << _NormMap[Hit->geographicalId().rawId()] << endl;

    float NormCharge = _NormMap[Hit->geographicalId().rawId()] * Charge * abs(cosine);

    if(NormCharge > 0 && NormCharge < _DeDxCut){
      _Points.pop_back();
      _Errors.pop_back();
      continue;
    }

    _Charges.push_back(NormCharge);

    NPoints++;
  }

  //cout << "Just Added " << NPoints << ".  Sizes: " << _Points.size() << " " << _Errors.size() << " " << _Charges.size() << endl;

  return NPoints;
}

void TrackCombinerReco::RemovePoints(int n){
  for(int i=0; i<n; i++){
    _Points.pop_back();
    _Errors.pop_back();
    _Charges.pop_back();
  }

  //cout << "Just Removed " << n << ".  Sizes: " << _Points.size() << " " << _Errors.size() << " " << _Charges.size() << endl;
}

void TrackCombinerReco::Clear(){
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
  _vNdofXY.clear();
  _vNdofRZ.clear();

  _vDeDx.clear();
  _vIso.clear();

  _vTHTrack.clear();
  _vTHX.clear();
  _vTHY.clear();
  _vTHZ.clear();
  _vTHErrX.clear();
  _vTHErrY.clear();
  _vTHErrZ.clear();
}

void TrackCombinerReco::Save(vector<int> &Group){
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
  _vNdofXY.push_back(_NdofXY);
  _vNdofRZ.push_back(_NdofRZ);

  _vDeDx.push_back(_DeDx);
  _vIso.push_back(_Iso);
}

void TrackCombinerReco::FitXY(vector<int> &Group){
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

  _XYPar[0] = Result->Parameter(0); // d0 (to first order)
  _XYErr[0] = Result->ParError(0); // not correct right now
  _XYPar[1] = Phi0 - atan(Result->Parameter(1)/(Result->Parameter(2)-Result->Parameter(0))); // phi0
  _XYErr[1] = Result->ParError(1); // not correct right now
  _XYPar[2] = Result->Parameter(2); // radius
  _XYErr[2] = Result->ParError(2);
  _Chi2XY = Result->Chi2();
  _NdofXY = Result->Ndf();

  //cout << _XYPar[0] << " " << _XYPar[1] << " " << _XYPar[2] << endl;
}

void TrackCombinerReco::FitRZ(bool Debug){ /* _RZFitter->ClearPoints();

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
  _NdofRZ = Result->Ndf();

}

void TrackCombinerReco::FitDeDx(){
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

  vector<float> SortedCharges;
  for(uint i=0; i<_Charges.size(); i++)
    if(_Charges[i]>0) SortedCharges.push_back(_Charges[i]);

  if(SortedCharges.size() == 0){
    _DeDx = 0;
    return;
  }

  sort(SortedCharges.begin(), SortedCharges.end());

  _DeDx = SortedCharges[SortedCharges.size()/2];
}

void TrackCombinerReco::AverageIso(vector<int> &Group){
  edm::Ref<std::vector<reco::Track> > InitTrack(_hTracks, Group[0]);

  float InitEta = InitTrack->eta();
  float InitPhi = InitTrack->phi();
  float IsoPt = 0;

  for(uint i = 0; i!=_hTracks->size(); i++){
    if(count(Group.begin(), Group.end(), i) > 0) continue;

    edm::Ref<std::vector<reco::Track> > ThisTrack(_hTracks, i);

    float dR2 = pow(ThisTrack->eta() - InitEta, 2) + pow(ThisTrack->phi() - InitPhi, 2);

    if(dR2 > .16) continue; // dR cut of 0.4

    IsoPt += ThisTrack->pt();
  }

  _Iso = IsoPt;
}

