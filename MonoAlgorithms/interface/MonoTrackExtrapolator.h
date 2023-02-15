#ifndef Monopoles_MonoAlgorithms_MonoTrackExtrapolator_h
#define Monopoles_MonoAlgorithms_MonoTrackExtrapolator_h

//////////////////////////////////
// C S Cowden 18 March 2013
// Extrapolate monopole tracks to 
// Ecal.
//////////////////////////////////

#include <cmath>

#include "Monopoles/MonoAlgorithms/interface/MonoDefs.h"

namespace Mono {

class MonoTrackExtrapolator {

public:
  inline MonoTrackExtrapolator() { };
  inline virtual ~MonoTrackExtrapolator() { };

  // extrapolate to specified radius value
  // function arguments:
  // fit par0, par1, par2, r
  inline double zVr(double p0, double p1, double p2, double r)
    { return p0+p1*r+p2*r*r; }

  // extrapolate to specified radius value, returns phi
  // function arguments:
  // fit par0, par1, par2, r
  inline double phiVr(double p0, double p1, double p2, double r)
    { return p1-asin( (r*r-p0*p2)/(2*r*(p2-p0)) ); }

  // extrapolate to a specified Z value
  // function arguments:
  // fit par0, par1, par2, z
  inline double rVz(double p0, double p1, double p2, double z) { return -1; }

  // find eta from z and r
  inline double eta(double z, double r) 
    { return asinh( z/r ); }


};

} // end mono namespace

#endif
