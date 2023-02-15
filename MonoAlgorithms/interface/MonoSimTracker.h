#ifndef MONOALGORITHMS_MONOSIMTRACKER_H
#define MONOALGORITHMS_MONOSIMTRACKER_H

////////////////////////////////////////////////////////
// C S Cowden				23 February 2012
// Helper class to trace monopole simulation hits through
// tracker subdetectors
////////////////////////////////////////////////////////


// std includes
#include <vector>
#include <iostream>


// FWCore includes
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"

// Data Formats
#include "SimDataFormats/CaloHit/interface/PCaloHitContainer.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

// Geometries
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"


// Monopole includes
#include "Monopoles/MonoAlgorithms/interface/MonoDefs.h"




// forward class declarations
namespace edm { class Event; class EventSetup; }


namespace Mono {


// enumerator to describe sub detectors
enum subdet {
  EcalEB=1,
  EcalEE,
  PixelEBLowTof
};



static const char * collNames[] = { ""  	// simTracks
    , "EcalHitsEB" 				// EcalEB
    , "TrackerHitsPixelBarrelLowTof" 		//PixelEBLowTof
};



struct SumStruct {
  double sum;
  unsigned N;
  SumStruct(): sum(0.),N(0) {}
  inline void operator+=(double add) {
    N++;
    sum += add;
  }

};



template <class T,class S>
class MonoSimTracker {

  public:
  MonoSimTracker(const edm::Event &, const edm::EventSetup &, const subdet);
  ~MonoSimTracker();


  // access methods
  inline const unsigned size(const MonoEnum m) const { return m_hits[m].size(); }

  // geometry access
  const double eta(const MonoEnum m,const unsigned i) const;
  const double phi(const MonoEnum m,const unsigned i) const;
  const double x(const MonoEnum m,const unsigned i) const;
  const double y(const MonoEnum m,const unsigned i) const;
  const double z(const MonoEnum m,const unsigned i) const;

  // geometry by detId
  const double eta(const unsigned i) const;
  const double phi(const unsigned i) const;


  // hit access (could be skipped and access the hit directory through hit()
  const double energy(const MonoEnum m,const unsigned i) const;
  const double time(const MonoEnum m,const unsigned i) const;

  const T * hit(const MonoEnum m,const unsigned i) const;
  inline const S* geo() const { return m_geo; }

  inline const std::map<unsigned,SumStruct> & idSumMap(const MonoEnum m) const { return m_idMapSum[m]; }
 


  private:

  void fillHits(const edm::Event &, const edm::EventSetup &);

  // static const members
  /*static const char * collNames[] = { ""  	// simTracks
    , "EcalHitsEB" 				// EcalEB
    , "TrackerHitsPixelBarrelLowTof" 		//PixelEBLowTof
  };*/

  // subdetector
  const subdet m_subdet;

  // geometry
  const S * m_geo;
  
  // vector of hits
  std::vector<T> m_hits[2];

  // detector id map to summed energy
  std::map<unsigned,SumStruct> m_idMapSum[2];


};  // end MonoSimTracker class declaration




// templated class icc file include
#include "Monopoles/MonoAlgorithms/interface/MonoSimTracker.icc"


} // end Mono namespace





#endif


