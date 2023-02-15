#ifndef Monopoles_MonoAlgorithms_NPVHelper_h
#define Monopoles_MonoAlgorithms_NPVHelper_h

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

class NPVHelper{ 
 public:
  NPVHelper(edm::ParameterSet const& iPS, edm::ConsumesCollector && iC);
  ~NPVHelper();

  inline unsigned int getNPV(const edm::Event & iEvent, const edm::EventSetup &iSetup)
  {
    // PrimaryVertex analysis
    edm::EDGetTokenT< reco::VertexCollection > m_PVTag;
    m_PVTag = consumes< reco::VertexCollection>(edm::InputTag("offlinePrimaryVertices",""));
    //edm::Handle<reco::VertexCollection> handlePV;
    //iEvent.getByToken(InputTag("offlinePrimaryVertices"),handlePV);
    //edm::InputTag m_PVTag = edm::InputTag("offlinePrimaryVertices","");
    edm::Handle<reco::VertexCollection> handlePV;
    iEvent.getByToken(m_PVTag,handlePV);


    int totalNPV = 0.;
    
    reco::VertexCollection::const_iterator pv = handlePV->begin();
    for ( ; pv != handlePV->end(); pv++ ) {
      if ( !pv->isFake() && pv->ndof() > 4.0 ) {
	++totalNPV;
      }
    }
    
    return totalNPV;
  }
  
} // end of NPVHelper class

#endif
