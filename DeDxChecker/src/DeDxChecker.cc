#include "Monopoles/DeDxChecker/interface/DeDxChecker.h"

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

#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"

//#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
//#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
//#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
//#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
//#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"

#include "RecoTracker/DeDx/interface/DeDxTools.h"

//#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"

#include <TFitResult.h>
#include <TVirtualFitter.h>

#include <sstream>

using namespace std; using namespace edm;

/// Constructor
DeDxChecker::DeDxChecker(const ParameterSet& parameterSet){
  _Source = parameterSet.getParameter<std::string>("Source");
  _Output = parameterSet.getParameter<std::string>("Output");
}

/// Destructor
DeDxChecker::~DeDxChecker(){
}

void DeDxChecker::beginJob(){
  _OutputFile = new TFile(_Output.c_str(), "recreate");
  _Tree = new TTree("tracks", "tracks");

  _Tree->Branch("NTracks", &_NTracks);
  _Tree->Branch("Pt", &_Pt);
  _Tree->Branch("Eta", &_Eta);
  _Tree->Branch("Phi0", &_Phi0);
  _Tree->Branch("D0", &_D0);
  _Tree->Branch("Z0", &_Z0);
  _Tree->Branch("DeDx", &_DeDx);
  _Tree->Branch("SatStrips", &_SatStrips);
  _Tree->Branch("TotStrips", &_TotStrips);

  _Tree->Branch("HitSatStrips", &_HitSatStrips);
  _Tree->Branch("HitTotStrips", &_HitTotStrips);
}

void DeDxChecker::endJob(){
  _OutputFile->cd();
  _Tree->Write();
  _OutputFile->Close();
}

/*
void DeDxChecker::Init(const edm::EventSetup& iSetup) {
  iSetup.get<GlobalTrackingGeometryRecord>().get(_TrackingGeom);

  // Fill Normalization map
  edm::ESHandle<TrackerGeometry> tkGeom;
  iSetup.get<TrackerDigiGeometryRecord>().get( tkGeom );

  vector<GeomDet*> Det = tkGeom->dets();
  for(uint i=0; i<Det.size(); i++){
    DetId Detid = Det[i]->geographicalId();

    StripGeomDetUnit* StripDetUnit = dynamic_cast<StripGeomDetUnit*> (Det[i]);
    PixelGeomDetUnit* PixelDetUnit = dynamic_cast<PixelGeomDetUnit*> (Det[i]);

    if(StripDetUnit){
      _NormMap[Detid.rawId()] = _MeVperADCStrip / StripDetUnit->surface().bounds().thickness();
    }else if(PixelDetUnit){
      _NormMap[Detid.rawId()] = _MeVperADCPixel / PixelDetUnit->surface().bounds().thickness();
    }
  }
}*/

void DeDxChecker::analyze(const Event& event, const EventSetup& setup){
  //if(_NormMap.size() == 0) Init(setup);

  event.getByLabel("generalTracks", _hTracks);
  //event.getByLabel(_Source, _hTrajectories);
  //event.getByLabel(_Source, _hTrajTrackAssociations);

  event.getByLabel("dedxHarmonic2", _hDeDx);

  // Loop over the tracks
  _NTracks = _hTracks->size();
  for(int i = 0; i!=_NTracks; i++){
    edm::Ref<std::vector<reco::Track> > TrackRef(_hTracks, i);

    _Pt = TrackRef->pt();
    _Eta = TrackRef->eta();
    _Phi0 = TrackRef->phi();
    _D0 = TrackRef->d0();
    _Z0 = TrackRef->dsz();

    _DeDx = (*_hDeDx.product())[TrackRef].dEdx();

    if(_DeDx < 4) continue;

    GetDeDxStrips(*TrackRef);

    _Tree->Fill();
  }
}

void DeDxChecker::GetDeDxStrips(const reco::Track &Track){
  _SatStrips = 0;
  _TotStrips = 0;

  _HitSatStrips.clear();
  _HitTotStrips.clear();

  for (trackingRecHit_iterator iHit=Track.recHitsBegin(); iHit!=Track.recHitsEnd(); iHit++){

    const TrackingRecHit *Hit = *iHit;

    if(!Hit->isValid()) continue;

    // add the dedx information: different accessor for every type of hit

    if(const SiStripMatchedRecHit2D* matchedHit=dynamic_cast<const SiStripMatchedRecHit2D*>(Hit)){
      const vector<uint8_t>& Ampls = matchedHit->stereoCluster().amplitudes();
      _TotStrips += Ampls.size();
      _HitTotStrips.push_back(Ampls.size());

      int Saturated=0;
      for(uint i=0; i<Ampls.size(); i++){
	if(Ampls[i] >= 254) Saturated++;
      }
      _SatStrips += Saturated;
      _HitSatStrips.push_back(Saturated);
    }else if(const ProjectedSiStripRecHit2D* projectedHit=dynamic_cast<const ProjectedSiStripRecHit2D*>(Hit)) {
      const vector<uint8_t>& Ampls = projectedHit->originalHit().stripCluster().amplitudes();
      _TotStrips += Ampls.size();
      _HitTotStrips.push_back(Ampls.size());

      int Saturated=0;
      for(uint i=0; i<Ampls.size(); i++){
	if(Ampls[i] >= 254) Saturated++;
      }
      _SatStrips += Saturated;
      _HitSatStrips.push_back(Saturated);
    }else if(const SiStripRecHit2D* singleHit=dynamic_cast<const SiStripRecHit2D*>(Hit)){
      const vector<uint8_t>& Ampls = singleHit->stripCluster().amplitudes();
      _TotStrips += Ampls.size();
      _HitTotStrips.push_back(Ampls.size());

      int Saturated=0;
      for(uint i=0; i<Ampls.size(); i++){
	if(Ampls[i] >= 254) Saturated++;
      }
      _SatStrips += Saturated;
      _HitSatStrips.push_back(Saturated);
    }else if(const SiStripRecHit1D* single1DHit=dynamic_cast<const SiStripRecHit1D*>(Hit)){
      const vector<uint8_t>& Ampls = single1DHit->stripCluster().amplitudes();
      _TotStrips += Ampls.size();
      _HitTotStrips.push_back(Ampls.size());

      int Saturated=0;
      for(uint i=0; i<Ampls.size(); i++){
	if(Ampls[i] >= 254) Saturated++;
      }
      _SatStrips += Saturated;
      _HitSatStrips.push_back(Saturated);
// don't use pixels for now (since the standard algorithms don't use them)
    //}else if(const SiPixelRecHit* pixelHit=dynamic_cast<const SiPixelRecHit*>(Hit)){ 
//      Charge = pixelHit->cluster()->charge();
    }

    //cout << "Adding Charge: " << Charge << " * " << abs(cosine) << " * " << _NormMap[Hit->geographicalId().rawId()] << endl;
  }
}

