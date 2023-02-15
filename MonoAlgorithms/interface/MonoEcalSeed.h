#ifndef  MONOALGORITHMS_MONOECALSEED_H
#define  MONOALGORITHMS_MONOECALSEED_H
/////////////////////////////////////////////////////////////
// C S Cowden                                       20 June 2012
//Seed class for monopole Ecal observables
//////////////////////////////////////////////////////////////


namespace Mono {


class MonoEcalSeed {

public:

  inline MonoEcalSeed():m_seedLength(3U),m_iEta(0U),m_iPhi(0U),m_energy(0.) { } 

  inline MonoEcalSeed(const unsigned L, const unsigned iEta, const unsigned iPhi,const double E)
    :m_seedLength(L)
    ,m_iEta(iEta)
    ,m_iPhi(iPhi)
    ,m_energy(E)
  { }

  inline virtual ~MonoEcalSeed() { }

  //accessor methods
  inline const unsigned seedLength() const { return m_seedLength; }
  inline const unsigned ieta() const { return m_iEta; }
  inline const unsigned iphi() const { return m_iPhi; }
  inline const double   energy() const { return m_energy; }

private:
  // length of the seed in eta (cells long)
  unsigned m_seedLength;
  // eta bin
  unsigned m_iEta;
  // phi bin
  unsigned m_iPhi;
  // energy
  double m_energy;

};


} // end Mono namespace


#endif
