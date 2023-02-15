// -*- C++ -*-
//
// Package:    MonoHack
// Class:      MonoHack
// 
/**\class MonoHack MonoHack.cc Monopoles/MonoHack/src/MonoHack.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Christopher Cowden
//         Created:  Tue Feb  7 14:25:42 CST 2012
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


//data formats
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

// Monopole includes
#include "Monopoles/MonoAlgorithms/interface/MonoDefs.h"



const int g_mono_id = MONOID;


//
// class declaration
//

class MonoHack : public edm::EDAnalyzer {
   public:
      explicit MonoHack(const edm::ParameterSet&);
      ~MonoHack();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

      // ----------member data ---------------------------
  const edm::InputTag m_mc_label;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MonoHack::MonoHack(const edm::ParameterSet& iConfig)
  :m_mc_label(iConfig.getParameter<edm::InputTag>("MC_Label"))
{
   //now do what ever initialization is needed

}


MonoHack::~MonoHack()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
MonoHack::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

  edm::Handle<edm::HepMCProduct> mcproduct;
  iEvent.getByLabel(m_mc_label, mcproduct);
  HepMC::GenEvent *mc = const_cast<HepMC::GenEvent *>(mcproduct->GetEvent());
  assert(mc);


  HepMC::GenEvent::particle_iterator p = mc->particles_begin();
  for ( ; p != mc->particles_end(); ++p ) {

    if ( abs( (*p)->pdg_id()) == 17 ) {
      //std::cout << "Monopole: " << (*p)->pdg_id() << " " ;
      (*p)->set_pdg_id( g_mono_id*(*p)->pdg_id()/abs((*p)->pdg_id()) );
      //std::cout << (*p)->pdg_id() << std::endl;
    }

  } 
  
  

   
}


// ------------ method called once each job just before starting event loop  ------------
void 
MonoHack::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MonoHack::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
MonoHack::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
MonoHack::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MonoHack::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MonoHack::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MonoHack::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MonoHack);
