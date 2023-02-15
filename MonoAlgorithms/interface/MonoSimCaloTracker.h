#ifndef MONOALGORITHMS_MONOSIMCALOTRACKER_H
#define MONOALGORITHMS_MONOSIMCALOTRACKER_H

//////////////////////////////////////////////////////////
// C S Cowden				23 February 2012
// Helper class to trace monopole simulation hits throught
// calorimeters
///////////////////////////////////////////////////////////


// std includes
#include <vector>

// Monopole includes
#include "Monopoles/MonoAlgorithms/interface/MonoDefs.h"

// forward class declarations
namespace edm { class Event; class EventSetup; }


namespace Mono {

class MonoSimCaloTracker {

  public:
  MonoSimCaloTracker( const edm::EventSetup &);
  ~MonoSimCaloTracker();

  // access methods

  private:



}


} // end Mono namespace

#endif
