

#include "Monopoles/MonoAlgorithms/interface/EcalMapper.h"

#include <iostream>
#include <cmath>


// FWCore includes
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// DataFormats
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/Common/interface/SortedCollection.h"


// Ecal includes
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"


// root includes
#include "TH2D.h"
#include "TLatex.h"



namespace Mono {


// constructor
EcalMapper::EcalMapper(const edm::EventSetup &es )
{

  // get calorimetry geometry
  edm::ESHandle<CaloGeometry> calo;
  es.get<CaloGeometryRecord>().get(calo);
  m_caloGeo = (const CaloGeometry*)calo.product();


}



// destructor
EcalMapper::~EcalMapper()
{ }



// process over the event
void EcalMapper::fillMap(const edm::Event &ev, TFileDirectory *dir  )
{ 


  TH2D * eMap;
  TH2D * tMap;

  char eMapName[50];
  sprintf(eMapName,"eMap_%d_%d_%llu",ev.id().run(),ev.id().luminosityBlock(),ev.id().event());

  char tMapName[50];
  sprintf(tMapName,"tMap_%d_%d_%llu",ev.id().run(),ev.id().luminosityBlock(),ev.id().event());


  // get Ecal EB geometry
  const CaloSubdetectorGeometry *geom = m_caloGeo->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);
  EBDetId did;

  const unsigned N_ETA = 2*did.MAX_IETA;
  const double min_eta = -did.MAX_IETA*did.crystalUnitToEta;
  const double max_eta = did.MAX_IETA*did.crystalUnitToEta;
  const unsigned N_PHI = did.MAX_IPHI;
  const double min_phi = -M_PI;
  const double max_phi = M_PI;


  // if TFileDirectory is passed put histograms in it
  // otherwise, use TFileService to dump histograms
  if ( dir ) {

    eMap = dir->make<TH2D>(eMapName,eMapName,N_ETA,min_eta,max_eta,N_PHI,min_phi,max_phi);    
    tMap = dir->make<TH2D>(tMapName,tMapName,N_ETA,min_eta,max_eta,N_PHI,min_phi,max_phi);

  } else {

    edm::Service<TFileService> m_fs;

    eMap = m_fs->make<TH2D>(eMapName,eMapName,N_ETA,min_eta,max_eta,N_PHI,min_phi,max_phi); 
    tMap = m_fs->make<TH2D>(tMapName,tMapName,N_ETA,min_eta,max_eta,N_PHI,min_phi,max_phi);
  }

  eMap->GetXaxis()->SetTitle("#eta");
  eMap->GetYaxis()->SetTitle("#phi");
  
  tMap->GetXaxis()->SetTitle("#eta");
  tMap->GetYaxis()->SetTitle("#phi");

  // EcalRecHits
  edm::InputTag TagEcalEB_RecHits("ecalRecHit","EcalRecHitsEB");
  edm::Handle<EBRecHitCollection > EcalRecHits;
  ev.getByLabel(TagEcalEB_RecHits,EcalRecHits);




  // loop over EcalRecHits
  EBRecHitCollection::const_iterator itHit = EcalRecHits->begin();
  for ( ; itHit != EcalRecHits->end(); itHit++ ) {

    EBDetId detId( (*itHit).id() );
    //PHAT fix
    //const CaloCellGeometry *cell = geom->getGeometry( detId );
    //if ( cell )  {
    //const unsigned etaBin = eMap->GetXaxis()->FindBin( cell->getPosition().eta() );
    //const unsigned phiBin = eMap->GetYaxis()->FindBin( cell->getPosition().phi() );
    //eMap->SetBinContent(etaBin,phiBin,(*itHit).energy() );
    //tMap->SetBinContent(etaBin,phiBin,(*itHit).time() );
    //}
    if (geom->getGeometry(detId)){
      const unsigned etaBin = eMap->GetXaxis()->FindBin(geom->getGeometry(detId)->getPosition().eta());
      const unsigned phiBin = eMap->GetYaxis()->FindBin(geom->getGeometry(detId)->getPosition().phi());
      eMap->SetBinContent(etaBin,phiBin,(*itHit).energy());
      tMap->SetBinContent(etaBin,phiBin,(*itHit).time());
    }
  }





}


// setMarker
void EcalMapper::setMarker(const double eta, const double phi, const char * text)
{

  marker t(eta,phi,text);
  m_markers.push_back( t );

}



} // end Mono namespace
