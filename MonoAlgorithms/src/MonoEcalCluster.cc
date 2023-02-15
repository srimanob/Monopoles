
#include "Monopoles/MonoAlgorithms/interface/MonoEcalCluster.h"

#include <cmath>
#include <cstdlib>

#include "Monopoles/MonoAlgorithms/interface/MonoEcalObs0.h"


namespace Mono {


const double MonoEcalCluster::energy(int etaDiff,int phiDiff, const EBmap &map) const
{

  const unsigned nEta = map.nEta();
  const unsigned nPhi = map.nPhi();

  assert( abs(etaDiff) < (int)nEta );

  unsigned newEta = 0U;
  unsigned newPhi = 0U;

  if ( etaDiff < 0 ) { 
    assert( abs(etaDiff) <= (int)m_iEta );
  } 
  newEta = m_iEta + etaDiff;

  assert( newEta < nEta );

  while ( phiDiff >= (int)nPhi ) phiDiff -= nPhi;
  while ( phiDiff < 0 ) phiDiff += nPhi;

  if ( phiDiff > 0 ) {
    newPhi = phiDiff + m_iPhi;
    if ( newPhi >= nPhi ) newPhi -= nPhi;
  } else {
    if ( abs(phiDiff) > (int)m_iPhi ) newPhi = nPhi - (m_iPhi + phiDiff);
    else newPhi = m_iPhi + phiDiff;
  } 

#ifdef DEBUG
  std::cout << "looking for energy at: " << newEta << " " << newPhi << " " << newPhi*nEta+newEta << std::endl;
  std::cout.flush();
#endif 

  assert( newPhi < nPhi );
  assert( newEta < nEta );

  const unsigned loc = newPhi*nEta+newEta;
  assert( loc < map.nCells() );
  return map[loc];

}

const double MonoEcalCluster::time(int etaDiff,int phiDiff, const EBmap &map) const
{

  const unsigned nEta = map.nEta();
  const unsigned nPhi = map.nPhi();

  assert( abs(etaDiff) < (int)nEta );

  unsigned newEta = 0U;
  unsigned newPhi = 0U;

  if ( etaDiff < 0 ) { 
    assert( abs(etaDiff) <= (int)m_iEta );
  } 
  newEta = m_iEta + etaDiff;

  assert( newEta < nEta );

  while ( phiDiff >= (int)nPhi ) phiDiff -= nPhi;
  while ( phiDiff < 0 ) phiDiff += nPhi;

  if ( phiDiff > 0 ) {
    newPhi = phiDiff + m_iPhi;
    if ( newPhi >= nPhi ) newPhi -= nPhi;
  } else {
    if ( abs(phiDiff) > (int)m_iPhi ) newPhi = nPhi - (m_iPhi + phiDiff);
    else newPhi = m_iPhi + phiDiff;
  } 

#ifdef DEBUG
  std::cout << "looking for energy at: " << newEta << " " << newPhi << " " << newPhi*nEta+newEta << std::endl;
  std::cout.flush();
#endif 

  assert( newPhi < nPhi );
  assert( newEta < nEta );

  const unsigned loc = newPhi*nEta+newEta;
  assert( loc < map.nCells() );
  return map.time(loc);

}

const EcalRecHit * MonoEcalCluster::getRecHit(int etaDiff, int phiDiff, const EBmap &map) const
{

  const unsigned nEta = map.nEta();
  const unsigned nPhi = map.nPhi();

  assert( abs(etaDiff) < (int)nEta );

  unsigned newEta = 0U;
  unsigned newPhi = 0U;

  if ( etaDiff < 0 ) { 
    assert( abs(etaDiff) <= (int)m_iEta );
  } 
  newEta = m_iEta + etaDiff;

  assert( newEta < nEta );

  while ( phiDiff >= (int)nPhi ) phiDiff -= nPhi;
  while ( phiDiff < 0 ) phiDiff += nPhi;

  if ( phiDiff > 0 ) {
    newPhi = phiDiff + m_iPhi;
    if ( newPhi >= nPhi ) newPhi -= nPhi;
  } else {
    if ( abs(phiDiff) > (int)m_iPhi ) newPhi = nPhi - (m_iPhi + phiDiff);
    else newPhi = m_iPhi + phiDiff;
  } 

#ifdef DEBUG
  std::cout << "looking for energy at: " << newEta << " " << newPhi << " " << newPhi*nEta+newEta << std::endl;
  std::cout.flush();
#endif 

  assert( newPhi < nPhi );
  assert( newEta < nEta );

  const unsigned loc = newPhi*nEta+newEta;
  assert( loc < map.nCells() );
  return map.getRecHit(loc);

}

}
