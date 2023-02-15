
#include "Monopoles/MonoAlgorithms/interface/MonoEcalObs0.h"

#include <iostream>

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "Monopoles/MonoAlgorithms/interface/MonoEcalCalibReader.h"

#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "CLHEP/Matrix/Matrix.h"



#ifdef __cplusplus
extern "C" {
#endif

void dgetrf_(int *M,int *N, double *A, int *lda, int *IPIV, int *INFO);
void dgetri_(int *N,double *A,int *lda,int *IPIV,double *WORK, int*lwork, int *INFO);

#ifdef __cplusplus
} // extern "C" ending
#endif



namespace Mono {


int nanChecker(const unsigned N, const double *data) {
  if(CLHEP::electron_charge==0) std::cout << "useless command" << std::endl;
  for ( unsigned i=0; i != N; i++ ) {
    if ( data[i] != data[i] ) return 1;
    else if ( isinf(data[i]) ) return -1;
  }
  return 0;
}


// --------------------------- EBmap member functions ---------------------------------

void EBmap::fillMap(const edm::Event &ev)
{
  // clear the maps
  clear();

  // fill Ecal Map with Rec Hit energy
  edm::InputTag tagEcalEB_RecHits("ecalRecHit","EcalRecHitsEB");
  edm::Handle<EBRecHitCollection> ecalRecHits;
  ev.getByLabel(tagEcalEB_RecHits,ecalRecHits);

  double maxE = 0.;

  const unsigned nHits = ecalRecHits->size();
  for ( unsigned i=0; i != nHits; i++ ) {
    EBDetId detId( (*ecalRecHits)[i].id() );
    const unsigned loc = findBin(detId);
    const double energy = (*ecalRecHits)[i].energy();
    if ( energy > maxE ) maxE = energy;
    m_ecalMap[loc] = energy;
    m_ecalTMap[loc] = (*ecalRecHits)[i].time();
    m_ecalRecHitMap[loc] = &(*ecalRecHits)[i];
  }  

}

void EBmap::clear()
{
  for ( unsigned i=0; i != m_nCells; i++ ) {
    m_ecalMap[i] = 0.;
    m_ecalTMap[i] = 0.;
  }
}

void EBmap::constructGeo(const edm::EventSetup & es)
{
  // get a handle on the Ecal Barrel subdetector geometry
  edm::ESHandle<CaloGeometry> calo;
  es.get<CaloGeometryRecord>().get(calo);
  const CaloGeometry * caloGeo = (const CaloGeometry*)calo.product();
  m_geom = caloGeo->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);
}

const unsigned EBmap::findBin( const EBDetId &did ) const
{
  
  //const CaloCellGeometry *cell = m_geom->getGeometry( did );
  //const double eta = cell->getPosition().eta();
  //const double phi = cell->getPosition().phi();
  const double eta = m_geom->getGeometry(did)->getPosition().eta();
  const double phi = m_geom->getGeometry(did)->getPosition().phi();
  return findBinEtaPhi(eta,phi);

}

const unsigned EBmap::findBinEtaPhi(const double eta, const double phi ) const
{

  if ( eta > m_maxEta || eta < m_minEta || phi > m_maxPhi || phi < m_minPhi )
    throw cms::Exception("eta,phi out of range") << "The eta and/or phi given to EBmap::findBinEtaPhi is out of range";

  assert(phi <= M_PI && phi >= -M_PI );
 
  unsigned phiBin = (phi - m_minPhi)/m_phiWidth;
  unsigned etaBin = (eta-m_minEta)/m_etaWidth;

  const unsigned bin = phiBin*m_nEta+etaBin;
  assert( bin < m_nCells );

  return bin;
 
} 


// --------------------------- StripSeedFinder member functions ------------------------

void StripSeedFinder::constructGeo(const edm::EventSetup &es)
{

  // get a handle on the Ecal Barrel subdetector geometry
  edm::ESHandle<CaloGeometry> calo;
  es.get<CaloGeometryRecord>().get(calo);
  const CaloGeometry * caloGeo = (const CaloGeometry*)calo.product();
  m_geom = caloGeo->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);


}


bool StripSeedFinder::find(const edm::Event &ev, const EBmap &ecalMap) 
{

  const unsigned nCells = ecalMap.nCells();
  assert( m_seeds.size() == nCells );

  // clear old seeds if any
  clear();

  const unsigned nEta = ecalMap.nEta();
  const unsigned nPhi = ecalMap.nPhi();


  // start looking for seeds
  const unsigned etaSearch = nEta - m_seedLength + 1U;
  for ( unsigned i=0; i != nPhi; i++ ) {
    for ( unsigned j=0; j != etaSearch; j++ ) {
      const unsigned loc = i*nEta+j;
      double energySum = ecalMap[loc];
      for ( unsigned s=1; s != m_seedLength; s++ ) 
	energySum += ecalMap[loc+s];

      if ( energySum > m_threshold ) addSeed(j,i,energySum);
    }
  } 

  // search in phi direction (cross-check of signature)
  /*const unsigned phiSearch = nPhi - m_seedLength + 1U;
  for ( unsigned i=0; i != phiSearch; i++ ) {
    for ( unsigned j=0; j != nEta; j++ ) {
      const unsigned loc = i*nEta+j;
      double energySum = ecalMap[loc];
      for ( unsigned s=1; s != m_seedLength; s++ )
	energySum += ecalMap[(i+s)*nEta+j];

      if ( energySum > m_threshold ) addSeed(j,i,energySum);
    }
  }*/

  if ( m_nSeeds > 1 ) mergeSeeds(ecalMap);
  centerSeeds(ecalMap);
  if ( m_nSeeds > 1 ) mergeSeeds(ecalMap);

  setLength(ecalMap);

  //for ( unsigned s=0; s != m_nSeeds; s++ )
  //  assert( m_seeds[s].seedLength() == m_clustLength );

  // check for seeds under threshold
  for ( unsigned s=0; s != m_nSeeds; s++ )
    assert( m_seeds[s].energy() >= 0. );

  return 1;
}



void StripSeedFinder::mergeSeeds(const EBmap &map)
{

  const unsigned nEta = map.nEta();

  // stage of indeces into m_seeds
  unsigned m_stageSize = 0U;

  // cycle over seeds
  // if no overlap of seed in stage, add to stage
  // else merge into seed in stage 
  m_stage[m_stageSize++] = 0U;
  for ( unsigned s = 1U; s != m_nSeeds; s++ ) {
    const MonoEcalSeed & seed = m_seeds[s];
    bool didMerge = false;

    for ( unsigned t = 0; t != m_stageSize; t++ ) {
      const unsigned seedIndex = m_stage[t];

      if ( overlapInEta(s,seedIndex) ) {
       	const unsigned sEta = seed.ieta();
     	const unsigned sPhi = seed.iphi();
        const MonoEcalSeed & teed = m_seeds[seedIndex];
    	const unsigned tEta = teed.ieta();
	const unsigned tPhi = teed.iphi();
      	const unsigned sLength = seed.seedLength();
	const unsigned tLength = seed.seedLength();

	if ( sEta > tEta ) {
	  const unsigned newLength = sLength + tLength - (tEta+tLength-sEta);
	  double newEnergy = teed.energy();
	  for ( unsigned i=tLength; i != newLength; i++ )
	    newEnergy += map[tPhi*nEta+tEta+i];
	  m_seeds[seedIndex] = MonoEcalSeed(newLength,tEta,tPhi,newEnergy);
	} else if ( tEta > sEta ) {
	  const unsigned newLength = tLength + sLength - (sEta+sLength-tEta);	
	  double newEnergy = seed.energy();
	  for ( unsigned i=sLength; i != newLength; i++ )
	    newEnergy += map[sPhi*nEta+sEta+i];
	  m_seeds[seedIndex] = MonoEcalSeed(newLength,sEta,sPhi,newEnergy);
	}

      	didMerge = true;	
      }
    }

    if ( !didMerge ) m_stage[m_stageSize++] = s;

  }

  // copy stage into m_seeds
  for ( unsigned s=0; s != m_stageSize; s++ ) 
    m_seeds[s] = m_seeds[m_stage[s]];
  m_nSeeds = m_stageSize;

  m_stageSize = 0U;

}


void StripSeedFinder::centerSeeds(const EBmap &map)
{

  const unsigned nEta = map.nEta();

  // cycle over seeds
  for ( unsigned s=0; s != m_nSeeds; s++ ) {
    const MonoEcalSeed & seed = m_seeds[s];
    const unsigned length = seed.seedLength();
    const unsigned loc = seed.iphi()*nEta+seed.ieta();
    unsigned hotSpot = 0U;
    double hotPeak = 0.;
    for ( unsigned i=0U; i != length; i++ ) {
      const double E = map[loc+i];
      if ( E > hotPeak ) {
	hotPeak = E;
	hotSpot = i;
      }
    }


    // is the hotspot centered?
    // if not...buffer the short side 
    const int middle = length/2;
    const int dist = hotSpot - middle;
    unsigned newLength = length;
    unsigned newEta = seed.ieta();
    if ( length % 2 == 0 ) {
      if ( dist > 1 ) {
	 newLength += dist;
      } else if ( dist < -1 ) {
	newLength -= dist;
	(unsigned)abs(dist) < newEta ? newEta -= (unsigned)abs(dist) : newEta = 0U;
      }
      else continue; // no need to do anything to the seed since the hot spot is centered
    } else if ( dist != 0 ) { 
      newLength += abs(dist);
      if ( dist < 0 ) (unsigned)abs(dist) < newEta ? newEta -= (unsigned)abs(dist) : newEta = 0U;
    } else continue; // no need to do anything to the seed since the hot spot is centered

    if ( newEta+newLength-1U >= nEta ) newLength = nEta-newEta; 

    assert( newEta < nEta );
    assert( newLength-1U + newEta < nEta );

    // find the new energy
    const unsigned newLoc = seed.iphi()*nEta+newEta;
    double newE = 0.;
    for ( unsigned i=0; i != newLength; i++ )
      newE += map[newLoc+i];
  
    // update seed list   
    m_seeds[s] = MonoEcalSeed(newLength,newEta,seed.iphi(),newE);
    
  }

}


void StripSeedFinder::setLength(const EBmap &map)
{

  // cycle over seeds find those that need to be truncated or extended
  for ( unsigned s=0; s != m_nSeeds; s++ ) {
    const MonoEcalSeed & seed = m_seeds[s];
    const unsigned length = seed.seedLength();

    if ( length > m_clustLength ) truncSeed(map,s);
    else if ( length < m_clustLength ) extendSeed(map,s);
    
  }
   

}

void StripSeedFinder::truncSeed(const EBmap &map, const unsigned s)
{
  assert(s < m_nSeeds);

  const unsigned nEta = map.nEta();
  
  const MonoEcalSeed & seed = m_seeds[s]; 
  const unsigned sLength = seed.seedLength();
  if ( sLength < m_clustLength ) return;

  // truncate the seed to the desired length and maximise the energy
  // in that length
  const unsigned iPhi = seed.iphi();
  std::vector<std::pair<unsigned,double> > subSeeds;
  unsigned curLoc = seed.ieta();
  const unsigned sEnd = curLoc + sLength;
  unsigned curEnd = curLoc+m_clustLength;
  while ( curEnd != sEnd ) {
    double energy = 0.;
    unsigned loc = iPhi*nEta+curLoc; 
    for ( unsigned i=0; i != m_clustLength; i++ ) {
      energy += map[loc+i]; 
    }
    subSeeds.push_back( std::pair<unsigned,double>(curLoc,energy) );
    curLoc++;
    curEnd = curLoc+m_clustLength; 
  }

  // find the largest energy in subSeeds and keep as the seed
  const unsigned nSubSeeds = subSeeds.size();
  unsigned maxEl = 0U;
  double maxVal = 0.;
  for ( unsigned i=0; i != nSubSeeds; i++ ) {
    const unsigned el = subSeeds[i].first;
    const double val = subSeeds[i].second;
    if ( val > maxVal ) {
      maxVal = val;
      maxEl = el;
    }
  }

  assert( maxEl+m_clustLength < nEta );
  assert( maxEl < 200 );
  assert( maxVal > 0.);

  // put replace seed s by new seed
  m_seeds[s] = MonoEcalSeed(m_clustLength,maxEl,iPhi,maxVal); 
  

}

void StripSeedFinder::extendSeed(const EBmap &map, const unsigned s)
{
  assert( s < m_nSeeds );

  const unsigned nEta = map.nEta();

  const MonoEcalSeed & seed = m_seeds[s];
  const unsigned sLength = seed.seedLength();
  if ( sLength > m_clustLength ) return;

  // extend the seed to desired length by symmetrically adding cells to either side
  const unsigned iPhi = seed.iphi();
  const unsigned diff = m_clustLength - sLength;
  const unsigned leftAdd = diff/2U;
  const unsigned rightAdd = diff-leftAdd;

  const unsigned curEta = seed.ieta();


  assert( leftAdd+rightAdd == diff );

  unsigned newEta = 9999;
  // make sure we don't hang off the end of the barrel
  if ( leftAdd <= curEta && rightAdd+sLength-1U+curEta < nEta ) {
    newEta = curEta-leftAdd; 
  } else if ( leftAdd > curEta ) {
    newEta=0U;
  } else if ( rightAdd+sLength-1U+curEta >= nEta ) {
    newEta = nEta-m_clustLength;
  } else {
    // don't know how we got here
    std::cerr << "StripSeedFinder::extendSeed unkown seed location" << std::endl;
  }
  
  
  double energy=0;
  const unsigned loc = iPhi*nEta+newEta;
  for ( unsigned i=0; i != m_clustLength; i++ ) 
    energy += map[loc+i];


  assert( newEta+m_clustLength-1U<nEta);
  assert( newEta < nEta ) ;
  assert( energy > 0. );

  m_seeds[s] = MonoEcalSeed(m_clustLength,newEta,iPhi,energy);
    

}


// --------------------------- ClusterBuilder member functions ---------------------

void ClusterBuilder::buildClusters(const unsigned nSeeds, const MonoEcalSeed *seeds, const EBmap &map)
{

  assert( seeds );
  m_nClusters = 0U;
  if ( nSeeds > m_clusters.size() ) m_clusters.resize(nSeeds);

  const unsigned nEta = map.nEta();
  const unsigned nPhi = map.nPhi();

  // add N eta strips to each side of the seed
  const unsigned N=2;

  for ( unsigned s=0; s != nSeeds; s++ ) {
    const MonoEcalSeed & seed = seeds[s];
    const unsigned sEta = seed.ieta();
    const unsigned sPhi = seed.iphi();
    const double sEnergy = seed.energy();
    const unsigned length = seed.seedLength();

    double energy = sEnergy;
    // add below in phi
    for ( unsigned i=N; i != 0; i-- ) {
      int newPhiS = (int)sPhi-i;
      unsigned newPhi = 0U;
      if ( newPhiS < 0 ) newPhi = nPhi+newPhiS;
      else newPhi = newPhiS;
      assert( newPhi < nPhi );
      for ( unsigned i=0; i != length; i++ ) {
	energy += map[newPhi*nEta+sEta+i];
      } 
    }
    // add above in phi
    for ( unsigned i=1; i != N+1; i++ ) {
      unsigned newPhi = sPhi+i;
      if ( newPhi >= nPhi ) newPhi = newPhi-nPhi;
      assert ( newPhi < nPhi );
      for ( unsigned i=0; i != length; i++ ) {
	energy += map[newPhi*nEta+sEta+i];
      }
    }

    m_clusters[m_nClusters++] = MonoEcalCluster(length,2*N+1U,sEta,sPhi,energy,seed); 
  }
  
}

// --------------------------- GenMonoClusterTagger member functions ----------------

void GenMonoClusterTagger::initialize( const edm::Event &ev, const edm::EventSetup &es)
{
  
  // get the generator product
  MonoTruthSnoop snooper(ev,es);
  const HepMC::GenParticle *mono = snooper.mono(Mono::monopole);
  const HepMC::GenParticle *anti = snooper.mono(Mono::anti_monopole);

  if ( mono ) {
    m_extrap.setMonopole(*mono);
    m_monoPID.push_back( mono->pdg_id() );
    m_monoPt.push_back( mono->momentum().perp() );
    m_monoPhi.push_back( m_extrap.phi() );
    
    if ( m_tagEB ) {
      m_monoEta.push_back( m_extrap.etaVr(s_EcalR) );
      m_monoTime.push_back( m_extrap.tVr(s_EcalR) );
    } else {
      double z = s_EEz;
      double tmp = m_extrap.rVz(z);
      if ( tmp == -1 ) {
      	z = -z;
	tmp = m_extrap.rVz(tmp);
      }
      m_monoEta.push_back( m_extrap.eta(z,tmp) );
      m_monoTime.push_back( m_extrap.tVz(z) );
    }
  } 

  if ( anti ) {
    m_extrap.setMonopole(*anti);
    m_monoPID.push_back( anti->pdg_id() );
    m_monoPt.push_back( anti->momentum().perp() );
    m_monoPhi.push_back( m_extrap.phi() );

    if ( m_tagEB ) {
      m_monoEta.push_back( m_extrap.etaVr(s_EcalR) );
      m_monoTime.push_back( m_extrap.tVr(s_EcalR) );
    } else {
      double z = s_EEz;
      double tmp = m_extrap.rVz(z);
      if ( tmp == -1 ) {
      	z = -z;
	tmp = m_extrap.rVz(tmp);
      }
      m_monoEta.push_back( m_extrap.eta(z,tmp) );
      m_monoTime.push_back( m_extrap.tVz(z) );
    }
  }


}

void GenMonoClusterTagger::tag(const unsigned nClusters, const MonoEcalCluster *clusters, const EBmap &map)
{
  assert( clusters );
  const unsigned nMonopoles = m_monoPID.size();

  if ( nClusters > m_tagged.size() ) {
    m_tagged.resize(nClusters);
    m_dR.resize(nClusters);
    m_matchTime.resize(nClusters);
    m_matchPt.resize(nClusters);
    m_monoMatch.resize(nClusters);
  }

  for ( unsigned c=0; c != nClusters; c++ ) {
    double minDR = DBL_MAX;
    double minTime = DBL_MAX;
    double minPt = DBL_MAX;
    int minPID = 0;
    const MonoEcalCluster & cluster = clusters[c];
    const unsigned iEta = cluster.ieta();
    const unsigned iPhi = cluster.iphi();
    const double eta = map.eta(iEta);
    const double phi = map.phi(iPhi);
    for ( unsigned m=0; m != nMonopoles; m++ ) {
      const double dR = reco::deltaR(m_monoEta[m],m_monoPhi[m],eta,phi);
      if ( dR < minDR ) {
	minDR = dR;
  	minTime = m_monoTime[m];
	minPt = m_monoPt[m];
	minPID = m_monoPID[m];
      }
    }
    m_dR[c] = minDR;
    m_matchTime[c] = minTime;
    m_matchPt[c] = minPt;
    m_monoMatch[c] = minPID;
    const bool tagged = minDR < m_dRcut;
    m_tagged[c] = tagged;
    m_nClusters++;
  } 

}

// --------------------------- MonoEcalObs0 member functions ------------------------

double MonoEcalObs0::calculate(const edm::EventSetup &es, const edm::Event &ev
  , std::vector<double> *betas, std::vector<double> *betaTs)
{

  assert( betas );
  assert( betaTs );

  m_ecalMap.constructGeo(es);
  m_seedFinder.constructGeo(es);

  // fill ecal map
  m_ecalMap.fillMap(ev);

  // run seed finder
  m_seedFinder.find(ev,m_ecalMap);

  // run the cluster builders
  m_clusterBuilder.buildClusters(m_seedFinder.nSeeds(),m_seedFinder.seeds(),m_ecalMap);

  // cycle over found clusters in event
  /*const unsigned nClusters = m_clusterBuilder.nClusters();
  betas->resize( nClusters );
  const MonoEcalCluster * clusters = m_clusterBuilder.clusters();
  for ( unsigned c=0; c != nClusters; c++ ) {
    const unsigned length = clusters[c].clusterLength();
    const unsigned width = clusters[c].clusterWidth();
    const ClustCategorizer cat(length,width);
    const unsigned side = length*width;
    const unsigned size = side*side;
   
    std::vector<double> & hMat = m_hMatMap[cat];
    std::vector<double> & eMat = m_Eavg[cat];
    // assertion here on size, do better error checking and handling in the future
    // what to do if encounter un-categorized cluster size
    //assert( hMat.size() );
    if ( !hMat.size() ) {
      std::cerr << "MonoEcalObs0 encountered an unknown cluster size: " << length << " " << width
      << " in event: " << ev.id().event() << std::endl;
      continue;
    }

    if ( m_wsSize < side ) {
      m_wsSize = side;
      m_workspace.resize(m_wsSize);
    }

    // clear workspace
    for ( unsigned i=0; i != m_wsSize; i++ ) m_workspace[i] = 0.;

  
    const double eTot = clusters[c].clusterEnergy();
    // fill workspace
    for ( unsigned i=0; i != width; i++ ) {
      int ki = (int)i-(int)width/2;
      for ( unsigned j=0; j != length; j++ ) {
	unsigned num = i*length+j;
	//m_workspace[num] = clusters[c].energy(j,ki,m_ecalMap)/eTot-m_functor(j-(int)length/2,ki-(int)width/2);
	m_workspace[num] = clusters[c].energy(j,ki,m_ecalMap)/eTot-eMat[i*length+j];
      }
    }

    // calculate beta for this cluster
    double & beta = (*betas)[c];
    beta = 0.;
    for ( unsigned i=0; i != side; i++ ) {
      const double eli = m_workspace[i];
      for ( unsigned j=0; j != side; j++ ) {
	const double term = eli*hMat[i*side+j]*m_workspace[j];
	beta += term;
      }
    }

    
  }*/


  return 0.;
}



// ------------------------- MonoEcalOb0Calibrator ---------------------------------

void MonoEcalObs0Calibrator::calculateMijn(const edm::EventSetup &es, const edm::Event &ev)
{

  m_ecalMap.constructGeo(es);
  m_seedFinder.constructGeo(es);

  // fill ecal map
  m_ecalMap.fillMap(ev);

  // run seed finder
  m_seedFinder.find(ev,m_ecalMap);

  // run the cluster builders
  m_clusterBuilder.buildClusters(m_seedFinder.nSeeds(),m_seedFinder.seeds(),m_ecalMap);


  //const unsigned nEta = m_ecalMap.nEta();

  // cycle over clusters and fill appropriate map
  const unsigned nClusters = m_clusterBuilder.nClusters();
  const MonoEcalCluster * clusters = m_clusterBuilder.clusters();
  for ( unsigned i=0; i != nClusters; i++ ) {
    const unsigned length = clusters[i].clusterLength();
    const unsigned width = clusters[i].clusterWidth();

    const ClustCategorizer cat(length,width);

    // check if catergory exists in map
    std::map<ClustCategorizer,std::vector<std::vector<double> > >::iterator iter = m_Mijn.find(cat);
    if ( iter == m_Mijn.end() ) {
      std::vector<std::vector<double> > & vec = m_Mijn[cat];
      vec.resize( width*width*length*length );
    }

    std::vector<std::vector<double> > & vec = m_Mijn[cat];



//    const unsigned ieta = clusters[i].ieta();
//    const unsigned iphi = clusters[i].iphi();
    const unsigned size = width*length;
    const double eTot = clusters[i].clusterEnergy();

    if ( size <= m_wsSize ) {
      m_wsSize = size;
      m_workspace.resize(m_wsSize);
    }

    // clear workspace
    for ( unsigned j=0; j != m_wsSize; j++ ) m_workspace[j]=0.;

    // fill workspace
    for ( unsigned k=0; k != width; k++ ) {
      int ki = (int)k-(int)width/2;
      for ( unsigned j=0; j != length; j++ ) {
	unsigned num = k*length+j;
	m_workspace[num] = clusters[i].energy(j,ki,m_ecalMap)/eTot-m_functor(j-(int)length/2,ki-(int)width/2);
      }
    }
  

    // calculate M_ij^n
    assert( size*size == vec.size() );
    for ( unsigned j=0; j != size; j++ ) {
      const double E1 = m_workspace[j];
      for ( unsigned k=0; k != size; k++ ) {
	vec[j*size+k].push_back(E1*m_workspace[k]);
      }
    }
	 
  }

}

void MonoEcalObs0Calibrator::fillClust(const edm::EventSetup &es, const edm::Event &ev)
{

  m_ecalMap.constructGeo(es);
  m_seedFinder.constructGeo(es);

  // fill ecal map
  m_ecalMap.fillMap(ev);

  // run seed finder
  m_seedFinder.find(ev,m_ecalMap);

  // run the cluster builders
  m_clusterBuilder.buildClusters(m_seedFinder.nSeeds(),m_seedFinder.seeds(),m_ecalMap);


  //const unsigned nEta = m_ecalMap.nEta();

  // cycle over clusters and fill appropriate map
  const unsigned nClusters = m_clusterBuilder.nClusters();
  const MonoEcalCluster * clusters = m_clusterBuilder.clusters();
  for ( unsigned i=0; i != nClusters; i++ ) {
    const unsigned length = clusters[i].clusterLength();
    const unsigned width = clusters[i].clusterWidth();

    const ClustCategorizer cat(length,width);

    // check if catergory exists in map
    std::map<ClustCategorizer,std::vector<std::vector<double> > >::iterator iter = m_Eclusts.find(cat);
    if ( iter == m_Eclusts.end() ) {
      std::vector<std::vector<double> > & vec = m_Eclusts[cat];
      vec.resize( width*length );
    }

    std::map<ClustCategorizer,std::vector<std::vector<double> > >::iterator tter = m_Tclusts.find(cat);
    if ( tter == m_Tclusts.end() ) {
      std::vector<std::vector<double> > & tvec = m_Tclusts[cat];
      tvec.resize( width*length );
    }

    std::vector<std::vector<double> > & vec = m_Eclusts[cat];
    std::vector<std::vector<double> > & tvec = m_Tclusts[cat];



    //const unsigned ieta = clusters[i].ieta();
    //const unsigned iphi = clusters[i].iphi();
    //const unsigned size = width*length;
    const double eTot = clusters[i].clusterEnergy();


    // fill workspace
    for ( unsigned k=0; k != width; k++ ) {
      int ki = (int)k-(int)width/2;
      for ( unsigned j=0; j != length; j++ ) {
	unsigned num = k*length+j;
	vec[num].push_back( clusters[i].energy(j,ki,m_ecalMap)/eTot );
	tvec[num].push_back( clusters[i].time(j,ki,m_ecalMap) );
      }
    }
  
  }
}

void MonoEcalObs0Calibrator::calculateHij()
{

  computeMij();

  assert( m_Mij.size() );
  assert( m_MTij.size() );

  // let's invert M
  MIJType::iterator mijiter = m_Mij.begin();
  MIJType::iterator mijEnd = m_Mij.end();

  //unsigned count=0;
  for( ; mijiter != mijEnd; mijiter++ ) {
    std::vector<double> & hVec = m_hij[mijiter->first];
    std::vector<double> & mVec = mijiter->second;
    const unsigned size = mVec.size();
    hVec.resize(size);

    const unsigned side = std::sqrt(size);

    //double minEl = DBL_MAX;
    //double maxEl = -(DBL_MAX-1.);

    // decompose mVec with lapack
    int n = side;
    int IPIV = side+1;
    int * lda = new int[IPIV];
    int INFO;
    dgetrf_(&n,&n,&mVec[0],&n,lda,&INFO); 

    bool infoTest = INFO==0;
    std::cout << "Decomposition result: " << INFO << std::endl;
    //assert(!INFO);

    double * work = new double[size];
    int lwork = size;
    dgetri_(&n,&mVec[0],&n,lda,work,&lwork,&INFO);

    infoTest = infoTest && INFO==0;
    std::cout << "Inversion result: " << INFO << std::endl;
    //assert(!INFO);

    delete [] lda;
    delete [] work;

    if ( !infoTest ) {
	hVec.clear();
	continue;
    }
    for ( unsigned i=0; i != size; i++ )
      hVec[i] = mVec[i];

    // check for nan elements
    const int nc = nanChecker(hVec.size(),&hVec[0]);
    if ( nc != 0 ) {
      throw cms::Exception("Found nan/inf Hij element") << " nc = " << nc;
    }

    /*CLHEP::HepMatrix mat(side,side);
    CLHEP::HepMatrix check(side,side);
    for ( unsigned i=0; i != side; i++ ) {
      for ( unsigned j=0; j != side; j++ ) {
	double el = mVec[i*side+j];
    	mat[i][j] = el;
    	check[i][j] = el;
	if ( el > maxEl ) maxEl = el;
	if ( el < minEl ) minEl = el;
	assert( el == mVec[j*side+i] );
      }
    }

    const double det = mat.determinant();
    assert( det );

    int success;
    mat.invert(success);
    assert(success);
    if ( success ) {
      for ( unsigned i=0; i != side; i++ ) {
   	for ( unsigned j=0; j != side; j++ ) {
	  hVec[i*side+j] = mat[i][j];
    	}
      }
    }

    CLHEP::HepMatrix res = check * mat;
    for ( unsigned i=0; i != side; i++ ) {
      std::cout << "i " << i << " " << res[i][i] << " " << check[i][i] << " " << mat[i][i] << std::endl;
    }
  
    count++; */

  }



  // let's invert M
  MIJType::iterator tijiter = m_MTij.begin();
  MIJType::iterator tijEnd = m_MTij.end();

  for( ; tijiter != tijEnd; tijiter++ ) {
    std::vector<double> & hVec = m_hTij[tijiter->first];
    std::vector<double> & mVec = tijiter->second;
    const unsigned size = mVec.size();
    hVec.resize(size);

    const unsigned side = std::sqrt(size);

    //double minEl = DBL_MAX;
    //double maxEl = -(DBL_MAX-1.);

    // decompose mVec with lapack
    int n = side;
    int IPIV = side+1;
    int * lda = new int[IPIV];
    int INFO;
    dgetrf_(&n,&n,&mVec[0],&n,lda,&INFO); 

    bool infoTest = INFO==0;
    std::cout << "Decomposition result: " << INFO << std::endl;
    //assert(!INFO);

    double * work = new double[size];
    int lwork = size;
    dgetri_(&n,&mVec[0],&n,lda,work,&lwork,&INFO);

    infoTest = infoTest && INFO==0;
    std::cout << "Inversion result: " << INFO << std::endl;
    //assert(!INFO);

    delete [] lda;
    delete [] work;

    if ( !infoTest ) {
	hVec.clear();
	continue;
    }
    for ( unsigned i=0; i != size; i++ )
      hVec[i] = mVec[i];


  }
  

}


void MonoEcalObs0Calibrator::computeMij()
{

  assert( m_Eclusts.size() );
  assert( m_Tclusts.size() );

  // find average E in each cluster bin
  MIJNType::iterator eijn = m_Eclusts.begin();
  MIJNType::iterator eijnEnd = m_Eclusts.end();
  for ( ; eijn != eijnEnd; eijn++ ) {
    std::vector<std::vector<double> > & vec = eijn->second;
    const unsigned size = vec.size();
    assert(size);
    const unsigned N = vec[0].size();
    assert(N);
    std::vector<double> & avgs = m_Eavg[eijn->first];
    avgs.clear();
    avgs.resize(size); 

    // for loop over cluster elements
    for ( unsigned i=0; i != size; i++ ) {
      long double mean = 0.L;
      // loop over cluster entries
      for ( unsigned n=0; n != N; n++ ) {
	mean += vec[i][n];	
      }
      mean /= N;
      avgs[i] = mean;
    }
  }
  
  // compute Mij
  eijn = m_Eclusts.begin();
  for ( ; eijn != eijnEnd; eijn++ ) {
    std::vector<std::vector<double> > & vec = eijn->second;
    const unsigned size = vec.size();
    assert(size);
    const unsigned N = vec[0].size();
    assert(N);
    std::vector<double> & avgs = m_Eavg[eijn->first];

    std::vector<double> & mijVec = m_Mij[eijn->first];
    mijVec.clear();
    mijVec.resize(size*size);

    // loop over cluster elements
    for ( unsigned i=0; i != size; i++ ) {
      for ( unsigned j=0; j != size; j++ ) {
	long double El = 0.L;
	for ( unsigned n=0; n != N; n++ ) {
	  El += (vec[i][n]-avgs[i])*(vec[j][n]-avgs[j]);
	}		
	El /= N;
	mijVec[i*size+j] = El;
      }
    }

    // check for nan values
    const int nc = nanChecker(mijVec.size(),&mijVec[0]);
    if ( nc != 0 ) {
      throw cms::Exception("Fond nan/inf in Mij") << " nc = " << nc;
    }
    
  }



  // find average T in each cluster bin
  MIJNType::iterator tijn = m_Tclusts.begin();
  MIJNType::iterator tijnEnd = m_Tclusts.end();
  for ( ; tijn != tijnEnd; tijn++ ) {
    std::vector<std::vector<double> > & vec = tijn->second;
    const unsigned size = vec.size();
    assert(size);
    const unsigned N = vec[0].size();
    assert(N);
    std::vector<double> & avgs = m_Tavg[tijn->first];
    avgs.clear();
    avgs.resize(size); 

    // for loop over cluster elements
    for ( unsigned i=0; i != size; i++ ) {
      long double mean = 0.L;
      // loop over cluster entries
      for ( unsigned n=0; n != N; n++ ) {
	mean += vec[i][n];	
      }
      mean /= N;
      avgs[i] = mean;
    }
  }
  
  // compute Mij
  tijn = m_Tclusts.begin();
  for ( ; tijn != tijnEnd; tijn++ ) {
    std::vector<std::vector<double> > & vec = tijn->second;
    const unsigned size = vec.size();
    assert(size);
    const unsigned N = vec[0].size();
    assert(N);
    std::vector<double> & avgs = m_Tavg[tijn->first];

    std::vector<double> & mijVec = m_MTij[tijn->first];
    mijVec.clear();
    mijVec.resize(size*size);

    // loop over cluster elements
    for ( unsigned i=0; i != size; i++ ) {
      for ( unsigned j=0; j != size; j++ ) {
	long double El = 0.L;
	for ( unsigned n=0; n != N; n++ ) {
	  El += (vec[i][n]-avgs[i])*(vec[j][n]-avgs[j]);
	}		
	El /= N;
	mijVec[i*size+j] = El;
      }
    }
    
  }




  /*MIJNType::iterator mijn = m_Mijn.begin();
  MIJNType::iterator mijnend = m_Mijn.end();
  for ( ; mijn != mijnend; mijn++ ) {
    std::vector<std::vector<double> > & data = mijn->second;
    const unsigned size = data.size();
    assert(size);
    const unsigned N = data[0].size();
    assert(N);
    std::vector<double> & vec = m_Mij[mijn->first];
    vec.clear();
    vec.resize(size);

    for ( unsigned i=0; i != size; i++ ) {
      long double mij = 0.;
      for ( unsigned n=0; n != N; n++ ) 
	mij += data[i][n];

      mij /= N;
      vec[i] = (double)mij;
    }
    
  }*/

}


} //end Mono namespace
