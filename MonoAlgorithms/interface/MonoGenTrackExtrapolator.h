#ifndef MonoAlgorithms_MonoGenTrackExtrapolator_h
#define MonoAlgorithms_MonoGenTrackExtrapolator_h
///////////////////////////////////////////////
//  MonoGenTrackExtrapolator
//  A helper class to extrapolate the generator level monopole
//  This works in MKS units (so distances are in meters).
///////////////////////////////////////////////

#include <cmath>

#include "Monopoles/MonoAlgorithms/interface/MonoDefs.h"

namespace HepMC {
  class GenParticle;
  class FourVector;
}


namespace Mono {

class MonoGenTrackExtrapolator {

public:
  inline MonoGenTrackExtrapolator():
    s_charge(2*M_PI*s_eps0*s_hbar*s_c*s_c/s_e)  // conversion from fundamental unit to SI using Ampere's convention in the Dirac quantization condition
    ,s_mass(1e9*s_e/s_c/s_c)  // conversion from GeV to kg
    ,s_mom(1e9*s_e/s_c)
    ,m_x(0.),m_y(0.),m_z(0.)
    ,m_Bz(4.) { }

  virtual inline ~MonoGenTrackExtrapolator () { }

  // set the gen particle whose track we wish to extrapolate
  inline void setMonopole(const HepMC::GenParticle &part) {
    const int pdg_id = part.pdg_id();
    assert( abs(pdg_id) == MONOID );

    if ( pdg_id > 0 ) m_charge = s_charge;
    else m_charge = -1.*s_charge;

    const HepMC::FourVector & mom = part.momentum();

    m_m = mom.m()*s_mass;
    
    m_px = mom.px()*s_mom;
    m_py = mom.py()*s_mom;
    m_pz = mom.pz()*s_mom;
    m_pT = mom.perp()*s_mom;

    m_eta = mom.eta();
    m_phi = mom.phi();

    m_a = 0.5*m_charge*m_Bz*m_m/m_pT/m_pT;
    m_b = m_pz/m_pT;
    m_c = m_z/m_m;

  }

  // set the initial conditions
  inline void setVert(const double x, const double y, const double z) {
    m_x = x;
    m_y = y;
    m_z = z;
  }

  // accessor methods
  // return numerical values of interest as a funciton of input values

  // phi taken as a constant from generator
  inline double phi() const { return m_phi; }

  // position in Z as a function of r
  inline double zVr(double r) const { return m_a*r*r+m_b*r+m_c; }

  // position in r as a function of Z
  // returns -1 if no solution exists
  inline double rVz(double z) const;

  // position in eta as a function of z
  inline double etaVz(double z) const { 
    return eta(z, rVz(z));
  }

  // position in eta as a function of r
  inline double etaVr(double r) const { 
    return eta( zVr(r), r);
  }

  // time as a function of z 
  inline double tVz(double z) const { return tVr(rVz(z)); }

  // time as a function of r
  inline double tVr(double r) const { return m_m*r/m_pT; }


  // access class data
  inline double charge() const { return m_charge; }
  inline double mass() const { return m_m; }

  inline double a() const { return m_a; }
  inline double b() const { return m_b; }
  inline double c() const { return m_c; }


  // helper functions
  // calculate eta from given z and r coorinates
  // it does not use any member data, so it can be used in general
  inline double eta(double z, double r) const {
    const double tanTheta = r/fabs(z);
    const double tanTheta2 = tanTheta/(1.+sqrt(1.+tanTheta*tanTheta));

    if ( z >= 0. ) return -log(tanTheta2);
    else return log(tanTheta2);

  }


  

private:

  static constexpr double s_eps0 = 8.85418782e-12;
  static constexpr double s_hbar= 1.054571726e-34;
  static constexpr double s_c = 299792458.;
  static constexpr double s_e = 1.60217646e-19;
  const double s_charge;
  const double s_mass;
  const double s_mom;

  // monopole data
  double m_charge;
  double m_m;

  // initial momentum
  double m_px;
  double m_py;
  double m_pz;
  double m_pT;

  // initial eta and phi
  double m_phi;
  double m_eta;

  // initial coordinates
  double m_x;
  double m_y;
  double m_z;

  // magnetic field data
  double m_Bz;  

  // extrapolation parameters
  double m_a;
  double m_b;
  double m_c;


};


inline double MonoGenTrackExtrapolator::rVz(double z) const 
{
  // check we aren't dividing by zero
  assert( fabs(m_a) > 1e-10 );

  // check if soluation exists
  if ( m_a*z < 0 && m_b*m_b < -4.*m_a*z ) {
    return -1.;
  } 

  // return the radius
  return fabs((-m_b+sqrt(m_b*m_b+4*m_a*z))/(2*m_a)); 
  //return 1.;

}

} // end Mono namespace

#endif
