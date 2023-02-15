#ifndef MONOTRUTHSNOOP_H
#define MONOTRUTHSNOOP_H

///////////////////////////////////////////////////
// C S Cowden				16 February 2012
// Helper class to find stable monopoles in HepMC event record.
///////////////////////////////////////////////////

//std includes
#include <vector>

// DataFormat includes
#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

// Monopole includes
#include "Monopoles/MonoAlgorithms/interface/MonoDefs.h"

//
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//forward class declarations
namespace edm { class Event; class EventSetup; class ConsumesCollector;}
namespace HepMC { class GenParticle; }

namespace Mono {

class MonoTruthSnoop {

 public: 
  MonoTruthSnoop( const edm::Event &, const edm::EventSetup &);
  ~MonoTruthSnoop();
  
  // access methods
  // return the monopole corresponding to the MonoEnum
  inline const HepMC::GenParticle * const  mono(MonoEnum pole) const { return m_mono[pole];}
  
  // return the hard (unstable status == 3) monopole from hard scatter process.
  // the daughters and stable monopoles result from this guy's end vertex.
  inline const HepMC::GenParticle * const monoGen(MonoEnum pole) const { return m_monoGen[pole]; }
  
  // return a vector of the duaghters of the monopole corresponding to MonoEnum
  inline const std::vector<const HepMC::GenParticle *> & monoDaughter(MonoEnum pole) const { return m_monoDaughter[pole];}
  
  // static consts
  static const int MONO_PID = MONOID;

  private:

  edm::InputTag hepmcCollectionToken_;
  
  // clear member data
  void clear();

  // search function
  bool snoop( const edm::Event &, const edm::EventSetup &);

  // member data
  const HepMC::GenParticle *m_mono[2];
  const HepMC::GenParticle *m_monoGen[2];
  std::vector<const HepMC::GenParticle *> m_monoDaughter[2];


}; // end MonoTruthSnoop class declaration



} // end Mono namespace


#endif
