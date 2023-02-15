#ifndef MONOALGORITHMS_ENERGYFLOWFUNCTOR_H
#define MONOALGORITHMS_ENERGYFLOWFUNCTOR_H

#include <cmath>

namespace Mono {

class EnergyFlowFunctor {

public:
  inline EnergyFlowFunctor() { }
  inline EnergyFlowFunctor(const EnergyFlowFunctor &func)
    :m_N(func.m_N)
    ,m_sigmaEta(func.m_sigmaEta)
    ,m_sigmaPhi(func.m_sigmaPhi)
    { }

  inline virtual ~EnergyFlowFunctor() { }


  // set the function's parameters
  inline virtual void setParameters ( unsigned N, const double *pars ) {
    assert( N == 3 );
    m_N = pars[0];
    m_sigmaEta = pars[1];
    m_sigmaPhi = pars[2];
  }


  // compute the value of the function
  inline virtual double operator()(const double i, const double j) const {
    return m_N*exp(-pow(i/m_sigmaEta,2)/2.-pow(j/m_sigmaPhi,2)/2.);
  }



private:
  double m_N;
  double m_sigmaEta;
  double m_sigmaPhi;


};

} // end Mono namespace

#endif
