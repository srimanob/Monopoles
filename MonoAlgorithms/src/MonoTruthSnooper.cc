#include "Monopoles/MonoAlgorithms/interface/MonoTruthSnooper.h"

// std includes
#include <cassert>

// FW includes
#include "FWCore/Framework/interface/Event.h"

// dataformats includes
#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

namespace Mono {

  MonoTruthSnoop::MonoTruthSnoop( const edm::Event &ev, const edm::EventSetup &es)
  { 
    hepmcCollectionToken_ = edm::InputTag("generator");
    m_mono[monopole] = NULL;
    m_mono[anti_monopole] = NULL;
    //iC.consumes<edm::HepMCProduct>(hepmcCollectionToken_);
    snoop(ev,es); 
  }

  MonoTruthSnoop::~MonoTruthSnoop()
  { 
  }

  void MonoTruthSnoop::clear() 
  {
    m_mono[monopole] = NULL;
    m_mono[anti_monopole] = NULL;
    m_monoDaughter[monopole].clear();
    m_monoDaughter[anti_monopole].clear();
  }

  bool MonoTruthSnoop::snoop( const edm::Event &ev, const edm::EventSetup &es) 
  { 
    using namespace edm;
    
    // Get the MC event
    //std::cout<<"Start MonoTruthSnoop"<<std::endl;
    //iC.consumes<edm::HepMCProduct>(hepmcCollectionToken_);    

    //iC.consumes<edm::HepMCProduct>(hepmcCollectionToken_);
    //edm::EDGetTokenT<edm::HepMCProduct> hepmcCollectionToken_;
    //hepmcCollectionToken_ = consumes < edm::HepMCProduct >( edm::InputTag("generator") );
    edm::InputTag hepmcCollectionToken_ = edm::InputTag("generatorSmeared","");
    
    edm::Handle<edm::HepMCProduct> mcproduct;
    ev.getByLabel(hepmcCollectionToken_, mcproduct);
    
    const HepMC::GenEvent* mc = mcproduct->GetEvent();
    assert(mc);
    //std::cout<<"Pass edm::Handle"<<std::endl;
    
    // print the event record
    // mc->print();  // just do this for debug purposes
    
    // Cycle over MC particles
    const HepMC::GenEvent::particle_const_iterator end = mc->particles_end();
    for (HepMC::GenEvent::particle_const_iterator p = mc->particles_begin();
	 p != end; ++p)
      {
	const HepMC::GenParticle* particle = *p;
	const reco::Candidate::LorentzVector p4( particle->momentum());
	
	// I only care about monopoles
	if ( abs(particle->pdg_id()) != MONO_PID ) continue;
	
	// if monopole is part of the hard process
	if ( particle->status() == 3 ) {
	  
	  // loop over decay products
	  const HepMC::GenVertex * vert = particle->end_vertex();
   	HepMC::GenVertex::particles_out_const_iterator daughter = vert->particles_out_const_begin();
        for ( ; daughter != vert->particles_out_const_end(); daughter++ ) {
	  const reco::Candidate::LorentzVector d4( (*daughter)->momentum() );
	  if ( abs((*daughter)->pdg_id()) != MONO_PID && particle->pdg_id() > 0 )  {
	    m_monoDaughter[monopole].push_back( (*daughter) );
	  }
	  else if ( abs((*daughter)->pdg_id()) != MONO_PID && particle->pdg_id() < 0 ) {
	    m_monoDaughter[anti_monopole].push_back( (*daughter) );
	  }
	  
	}
	
	if ( particle->pdg_id() > 0 ) {
	  m_monoGen[monopole] = particle;
	} else if ( particle->pdg_id() < 0 ) {
	  m_monoGen[anti_monopole] = particle;
	}
	
	
	// normally this should be a status of 1
      } else if ( particle->status() == 1 ) {
	
	
	// write out the two stable monopole's kinematic information
	if ( particle->pdg_id() > 0 ) {
	  m_mono[monopole] = particle;
	} else if ( particle->pdg_id() < 0 ) {
	  m_mono[anti_monopole] = particle;
	}
      }  
      }
  return true; 
}





}
