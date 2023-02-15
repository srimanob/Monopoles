#ifndef MONOALGORITHMS_MONOECALOBS0_H
#define MONOALGORITHMS_MONOECALOBS0_H
/////////////////////////////////////////////////////////
// C S Cowden                       20 June 2012
// First monopole physics observable for CMS ECAL
/////////////////////////////////////////////////////////

#include <iostream>
#include <vector>
#include <string>
#include <cassert>
#include <cmath>
#include <cfloat>
#include <cstring>

#include "Monopoles/MonoAlgorithms/interface/MonoDefs.h"
#include "Monopoles/MonoAlgorithms/interface/MonoEcalSeed.h"
#include "Monopoles/MonoAlgorithms/interface/MonoEcalCluster.h"
#include "Monopoles/MonoAlgorithms/interface/ClustCategorizer.h"
#include "Monopoles/MonoAlgorithms/interface/EnergyFlowFunctor.h"
#include "Monopoles/MonoAlgorithms/interface/MonoEcalCalibReader.h"
#include "Monopoles/MonoAlgorithms/interface/MonoTruthSnooper.h"
#include "Monopoles/MonoAlgorithms/interface/MonoGenTrackExtrapolator.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"

// forward declarations
namespace edm {
  class Event;
  class EventSetup;
}

namespace Mono {


// ------------------------------------------------------------------
// Some helper functions

// check a matrix for nan elements
// returns 1 if nan is first encountered
// returns -1 if inf is first encountered
// it returns at the first encounter of either nan or inf
// otherwise it returns 0
int nanChecker(unsigned N, const double *data);


// ---------------------------------------------------------------------
// Ecal barrel map class
class EBmap {

public:
  inline EBmap() { constructMap(); }

  virtual inline ~EBmap() { }

  // construct the map
  inline void constructMap() 
    {
      const EBDetId did;
      m_nEta = 2*did.MAX_IETA;
      m_nPhi = did.MAX_IPHI;
      m_nCells = m_nPhi*m_nEta;

      m_minEta = -did.MAX_IETA*did.crystalUnitToEta;
      m_maxEta = did.MAX_IETA*did.crystalUnitToEta;
      m_minPhi = -M_PI;
      m_maxPhi = M_PI;

      m_etaWidth = did.crystalUnitToEta;
      m_phiWidth = 2.*M_PI/m_nPhi;

      m_ecalMap.resize(m_nCells);
      m_ecalTMap.resize(m_nCells);
      m_ecalRecHitMap.resize(m_nCells);
    }

  // fill map
  void fillMap(const edm::Event &ev);

  // clear map
  void clear();

  // construct Ecal geometry stuff
  void constructGeo(const edm::EventSetup &);


  // accessor methods
  const inline unsigned nCells() const { return m_nCells; }
  const inline unsigned nEta() const { return m_nEta; }
  const inline unsigned nPhi() const { return m_nPhi; }

  // return the energy of the given bin
  const double inline operator[](const unsigned bin) const
   { 
     assert(m_ecalMap.size());
     assert(bin < m_nCells );
     return m_ecalMap[bin];
   }

  // return the time in the given bin
  const double inline time(const unsigned bin ) const
    {
      assert(m_ecalTMap.size());
      assert(bin < m_nCells);
      return m_ecalTMap[bin];
    }

  // return the RecHit in the given bin
  inline const EcalRecHit * getRecHit(const unsigned bin ) const
    {
      assert(m_ecalRecHitMap.size());
      assert(bin < m_nCells);
      return m_ecalRecHitMap[bin];
    }


  // return the ecalMap bin of an EBdetId
  const unsigned findBin(const EBDetId & ) const;
  const unsigned findBinEtaPhi(double,double ) const;


  // return eta and phi given eta/phi bin (not combined bin)
  const inline  double eta(unsigned i) const
    {
      assert(i < m_nEta );
      return i*m_etaWidth+m_minEta;
    }
  const inline  double phi(unsigned i) const
    {
      assert( i < m_nPhi );
      return i*m_phiWidth+m_minPhi;
    }


private:
  
  // ecalMap
  unsigned m_nCells;
  unsigned m_nEta;
  unsigned m_nPhi;
  double m_minEta;
  double m_maxEta;
  double m_minPhi;
  double m_maxPhi;
  double m_etaWidth;
  double m_phiWidth;

  // map of energy
  std::vector<double> m_ecalMap;
  // map of time
  std::vector<double> m_ecalTMap;
  // map of rechits
  std::vector<const EcalRecHit*> m_ecalRecHitMap;


  // calorimetry geometry
  const CaloSubdetectorGeometry *m_geom;

};


// ----------------------------------------------------------------------
// Seed finder class based on strips in eta
class StripSeedFinder {

public:
  inline StripSeedFinder()
    :m_seedLength(3U)
    ,m_clustLength(5U)
    ,m_threshold(20.)
    ,m_nSeeds(0U) 
    { }

  inline StripSeedFinder(const unsigned seedLength,const unsigned clustLength,const double threshold, const unsigned cells)
    :m_seedLength(seedLength)
    ,m_clustLength(clustLength)
    ,m_threshold(threshold)
    ,m_nSeeds(0U)
    ,m_maxSeeds(cells)
    { }


  inline virtual ~StripSeedFinder()
   { }

  // initialize the memory for seed array
  inline void initialize() { m_seeds.resize(m_maxSeeds); m_stage.resize(m_maxSeeds); }

  // construct the seed array and obtain geometry related stuff
  void constructGeo(const edm::EventSetup &);

  // run the finder
  // return false of an error occurred
  // This function expects the following arguments:
  // edm::Event, EBmap reference
  bool find(const edm::Event &, const EBmap &);

  // clear seeds
  inline void clear() { m_nSeeds = 0U; }


  // accessor methods
  inline const unsigned nSeeds() const { return m_nSeeds; }
  inline const MonoEcalSeed * seeds() const { return &m_seeds[0]; }
  inline const unsigned seedLength() const { return m_seedLength; }

  inline const bool adjacentInPhi(const unsigned i,const unsigned j) const
    {
      assert ( i < m_nSeeds );
      assert ( j < m_nSeeds );

      const MonoEcalSeed & iSeed = m_seeds[i];
      const MonoEcalSeed & jSeed = m_seeds[j];

      const unsigned iEta = iSeed.ieta();
      const unsigned jEta = jSeed.ieta();
      const unsigned iLength = iSeed.seedLength();
      const unsigned jLength = jSeed.seedLength();

      if ( iEta - jEta > jLength ) return false;
      else if ( jEta - iEta > iLength ) return false;
  
      //(unused) const unsigned iPhi = m_seeds[i].iphi(); 
      //(unused) const unsigned jPhi = m_seeds[j].iphi();
      //Phat modify; to check
      return true; //abs(iPhi-jPhi) == 1;
    }

  const bool adjacentInEta(const unsigned i,const unsigned j) const
    {
      assert ( i < m_nSeeds );
      assert ( j < m_nSeeds );

      const MonoEcalSeed & iSeed = m_seeds[i];
      const MonoEcalSeed & jSeed = m_seeds[j];

      if ( iSeed.iphi() != jSeed.iphi() ) return false;

      const unsigned iEta = iSeed.ieta();
      const unsigned jEta = jSeed.ieta();
      const unsigned iLength = iSeed.seedLength();
      const unsigned jLength = jSeed.seedLength();

      bool res = iEta+iLength+1U == jEta;
      res = res && (jEta+jLength+1U == iEta);
      return res;
    }

  const bool overlapInEta(const unsigned i,const unsigned j) const
    {
      assert ( i < m_nSeeds );
      assert ( j < m_nSeeds );

      const MonoEcalSeed & iSeed = m_seeds[i];
      const MonoEcalSeed & jSeed = m_seeds[j];

      if ( iSeed.iphi() != jSeed.iphi() ) return false;

      const unsigned iEta = iSeed.ieta();
      const unsigned jEta = jSeed.ieta();
      const unsigned iLength = iSeed.seedLength();
      const unsigned jLength = jSeed.seedLength();

      bool res = iEta+iLength >= jEta;
      res = res && ( jEta+jLength >= iEta );
      return res;
    }

private:

  // add seed to seed list
  inline void addSeed(const unsigned ieta, const unsigned iphi, const double E )
   {
      m_seeds[m_nSeeds++] = MonoEcalSeed(m_seedLength,ieta,iphi,E);
   }

  // merge seeds
  void mergeSeeds(const EBmap &);

  // center hot spot of seeds in strip
  void centerSeeds(const EBmap &);

  // set the length to the desired cluster length
  void setLength(const EBmap &);

  // truncate seed if too long, pass the seed index to truncate
  void truncSeed(const EBmap &,unsigned);

  // extend seed if too short, pass the seed index to extend
  void extendSeed(const EBmap &,unsigned);

  // length of the seeds to find (number of cells in eta)
  unsigned m_seedLength;

  // length desired for cluster
  unsigned m_clustLength;

  // energy threshold of seed
  double m_threshold;
 
  // number of seeds found 
  unsigned m_nSeeds;

  // seed array
  unsigned m_maxSeeds;
  std::vector<MonoEcalSeed> m_seeds;
  unsigned m_stageSize;
  std::vector<unsigned> m_stage;

  // calorimetry geometry
  const CaloSubdetectorGeometry *m_geom;

};

// ----------------------------------------------------------
// Simply extends seeds in eta strips
class SimplePathFinder {

};


// ---------------------------------------------------------
// Build clusters from seeds
class ClusterBuilder {

public:
  inline ClusterBuilder(): m_nClusters(0U) { m_clusters.resize(50U); }

  inline virtual ~ClusterBuilder () { }

  // build clusters around seeds
  void buildClusters(unsigned,const MonoEcalSeed *, const EBmap &);

  // accessor methods
  inline const unsigned nClusters() const { return m_nClusters; }
  inline const MonoEcalCluster * clusters() const { return &m_clusters[0]; }

private:

  // number of clusters
  unsigned m_nClusters;

  std::vector<MonoEcalCluster> m_clusters;

};


// ----------------------------------------------------------
// GenMonopole Cluster tagger
class GenMonoClusterTagger {

public:
  inline GenMonoClusterTagger(double dRcut,bool tagEB = true) 
    :m_dRcut(dRcut) 
    ,m_tagEB(tagEB)
    ,m_nClusters(0U)
    { 
      m_tagged.resize(50U);
      m_dR.resize(50U);
      m_matchTime.resize(50U);
      m_matchPt.resize(50U);
      m_monoMatch.resize(50U);
    }


  inline virtual ~GenMonoClusterTagger() { }

  // initialize 
  // setup geometry, find monopoles in generator 
  void initialize(const edm::Event &ev, const edm::EventSetup &es);

  // tag
  void tag(unsigned nClusters, const MonoEcalCluster *clusters, const EBmap &map);

  // generic tag routine
  template<class T>
  inline void tag(unsigned nObs, const T *obs);

  // add gen monopole
  inline void addMonopole(const HepMC::GenParticle part) {
    m_monoPID.push_back(part.pdg_id());
    m_monoPt.push_back(part.momentum().perp());
    m_extrap.setMonopole(part);

    if ( m_tagEB ) {
      m_monoEta.push_back(m_extrap.etaVr(s_EcalR));
      m_monoPhi.push_back(m_extrap.phi());
      m_monoTime.push_back(m_extrap.tVr(s_EcalR));
    } else {
      double z = s_EEz;
      double tmp = m_extrap.rVz(z);
      if ( tmp == -1 ) {
	z = -s_EEz;
	tmp = m_extrap.rVz(z);
      }
      m_monoEta.push_back(m_extrap.eta(z,tmp));
      m_monoTime.push_back(m_extrap.tVz(z));
    }
  } 

  // clear member data
  inline void clear() {
    clearTags();
    clearMonopoles();
  }

  // just clear the tag and cluster information
  inline void clearTags() {
    const unsigned vecSize = m_dR.size();
    m_nClusters = 0U;
    for ( unsigned i=0; i != vecSize; i++ ) {
      m_dR[i] = 0.;
      m_matchTime[i] = 0.;
      m_matchPt[i] = 0.;
      m_tagged[i] = 0;
      m_monoMatch[i] = 0;
    }
  }

  // clear the generator particles
  inline void clearMonopoles() {
    m_monoEta.clear();
    m_monoPhi.clear();
    m_monoTime.clear();
    m_monoPt.clear();
    m_monoPID.clear();
  }

  // accessor methods
  inline unsigned nClusters() const { return m_nClusters; }

  // closest match (not one to one matching)
  inline const double * matchDR() const { return &m_dR[0]; }

  // time estimate of monopole arrival
  inline const double * matchTime() const { return &m_matchTime[0]; }

  // pt (at gen level) of monopole
  inline const double * matchPt() const { return &m_matchPt[0]; }

  //  PID of closest match monopole (a check for monopole or anti-monopole)
  inline const int * matchPID() const { return &m_monoMatch[0]; } 
  
  // tag results 
  inline const int * tagResult() const { return &m_tagged[0]; }


private:

  MonoGenTrackExtrapolator m_extrap;

  std::vector<double> m_dR;
  std::vector<double> m_matchTime; // time estimate of monopole arrival
  std::vector<double> m_matchPt;
  std::vector<int>    m_monoMatch;  // PID to distinguish monopole from anti-monopole
  std::vector<int>   m_tagged;

  std::vector<double> m_monoEta;
  std::vector<double> m_monoPhi;
  std::vector<double> m_monoTime;
  std::vector<double> m_monoPt;
  std::vector<int> m_monoPID;

  static constexpr double s_EcalR = 1.29;
  static constexpr double s_EEz = 3.144;

  double m_dRcut;  
  bool m_tagEB;
  unsigned m_nClusters;

};


template <class T>
void GenMonoClusterTagger::tag(const unsigned nObs, const T * obs)
{
  assert( obs );
  const unsigned nMonopoles = m_monoPID.size();

  if ( nObs > m_tagged.size() ) {
    m_tagged.resize(nObs);
    m_dR.resize(nObs);
    m_matchTime.resize(nObs);
    m_monoMatch.resize(nObs);
  }
  
  for ( unsigned c=0; c != nObs; c++ ) {
    double minDR = DBL_MAX;
    double minTime = DBL_MAX;
    int minPID = 0;
    const T & cluster = obs[c];
    const double eta = cluster.eta();
    const double phi = cluster.phi();
    for ( unsigned m=0; m != nMonopoles; m++ ) {
      const double dR = reco::deltaR(m_monoEta[m],m_monoPhi[m],eta,phi);
      if ( dR < minDR ) {
	minDR = dR;
	minTime = m_monoTime[m];
	minPID = m_monoPID[m];
      }
    }
    m_dR[c] = minDR;
    m_matchTime[c] = minTime;
    m_monoMatch[c] = minPID;
    const bool tagged = minDR < m_dRcut;
    m_tagged[c] = tagged;
    m_nClusters++;
  } 

}



// ---------------------------------------------------------
// MonoEcalObs0 
class MonoEcalObs0 {

public:

  inline MonoEcalObs0(const edm::ParameterSet &ps) 
    :m_seedLength(ps.getParameter<unsigned>("StripSeedLength") )
    ,m_clustLength(ps.getParameter<unsigned>("ClusterLength") )
    ,m_threshold(ps.getParameter<double>("SeedThreshold") )
    //,m_calibName(ps.getParameter<std::string>("EnergyCalibrationName") )
    //,m_tCalibName(ps.getParameter<std::string>("TimeCalibrationName") )
    ,m_wsSize(50U)
    {
      m_seedFinder = StripSeedFinder(m_seedLength,m_clustLength,m_threshold,m_ecalMap.nCells());
      m_seedFinder.initialize();
     
      //loadHMatTables(); 

      m_workspace.resize(m_wsSize);

      double pars[3] = {1.,0.4,0.4};
      m_functor.setParameters(3,pars);
    } 

  inline virtual ~MonoEcalObs0()
    {
     // if ( m_ecalMap ) delete m_ecalMap;
    }


  // calculate observable method
  // The results per each cluster are returned to inside the vectors
  // pass as arguments.  The energy beta is the first vector argument, whilst
  // the time beta is the second vector argument.
  double calculate(const edm::EventSetup &,const edm::Event &, std::vector<double> *, std::vector<double> *);

  // accessor methods
  inline const StripSeedFinder & finder() const { return m_seedFinder; }
  inline const ClusterBuilder & clusterBuilder() const { return m_clusterBuilder; }
  inline const EBmap & ecalMap() const { return m_ecalMap; }

  // set the functor parameters
  inline void setClusterParameters(const unsigned N, const double * pars)
    {
      m_functor.setParameters(N,pars);
    }

  // set the functor
  inline void setFunctor(const EnergyFlowFunctor &functor)
    {
      m_functor = EnergyFlowFunctor(functor);
    }


private:

  // -- private member functions
  
  // return the expected energy in bin i 
  double eBarI(unsigned i);

  // find Beta ij
  double betaij(unsigned i, unsigned j);

  // calculate M ij
  double mij(unsigned i, unsigned j);

  // load H matrix lookup tables
  inline void loadHMatTables() {
    MonoEcalCalibReader reader;
    reader.readCalib(m_calibName,&m_hMatMap,&m_Eavg);
    reader.readCalib(m_tCalibName,&m_hTMatMap,&m_Tavg);
  }


  // -- private member data
  // seed length
  unsigned m_seedLength;
  unsigned m_clustLength;
  // seed threshold
  double m_threshold;
  

  // ecalMap
  EBmap m_ecalMap;

  // calibration file name
  std::string m_calibName;
  std::string m_tCalibName;

  // some workspace
  unsigned m_wsSize;
  std::vector<double> m_workspace;

  // the seed finder
  StripSeedFinder m_seedFinder;

  // the cluster builder
  ClusterBuilder m_clusterBuilder;

  // H matrix look up table map
  MIJType m_hMatMap;
  MIJType m_Eavg;
  MIJType m_hTMatMap;
  MIJType m_Tavg;

  // energy flow functor
  EnergyFlowFunctor m_functor;

};  


//
// Calibrator class for the Monopole Ecal observable
class MonoEcalObs0Calibrator {

  typedef std::map<ClustCategorizer,std::vector<std::vector<double> > > MIJNType;
  //typedef std::map<ClustCategorizer,std::vector<double> >	          MIJType;

public:

  inline MonoEcalObs0Calibrator(const edm::ParameterSet &ps)
    :m_seedLength(ps.getParameter<unsigned>("StripSeedLength"))
    ,m_clustLength(ps.getParameter<unsigned>("ClusterLength"))
    ,m_threshold(ps.getParameter<double>("SeedThreshold"))
    ,m_calibName(ps.getParameter<std::string>("EnergyCalibrationName")) 
    ,m_tCalibName(ps.getParameter<std::string>("TimeCalibrationName")) 
    ,m_wsSize(50U)
    {
      m_seedFinder = StripSeedFinder(m_seedLength,m_clustLength,m_threshold,m_ecalMap.nCells());
      m_seedFinder.initialize();

      m_workspace.resize(m_wsSize);

      double pars[3] = {1.,0.4,0.4};
      m_functor.setParameters(3,pars);
    }

  inline virtual ~MonoEcalObs0Calibrator() { }

  // calculate M_ij^n called for every event
  void calculateMijn(const edm::EventSetup &, const edm::Event &);
  void fillClust(const edm::EventSetup &, const edm::Event &);

  // computed H_ij at end of run
  void calculateHij();

  // dump calibration
  inline void dumpCalibration()
  { 
    MonoEcalCalibReader reader;
    reader.dumpCalib(m_calibName,m_hij,m_Eavg);
    reader.dumpCalib(m_tCalibName,m_hTij,m_Tavg);
  }


  // set the functor parameters
  inline void setClusterParameters(const unsigned N, const double * pars)
    {
      m_functor.setParameters(N,pars);
    }

  // set the functor
  inline void setFunctor(const EnergyFlowFunctor &functor)
    {
      m_functor = EnergyFlowFunctor(functor);
    }

private:

  // -- private member functions
  // compute M_ij at end of run in order to find H_ij
  void computeMij();

  // -- private member data
  MIJType m_hij;
  MIJType m_Mij;
  MIJNType m_Mijn;
  MIJType m_hTij;
  MIJType m_MTij;

  // seed length
  unsigned m_seedLength;
  unsigned m_clustLength;
  // seed threshold
  double m_threshold;
  // output calibration file name
  std::string m_calibName;
  std::string m_tCalibName;

  EBmap m_ecalMap;
  StripSeedFinder m_seedFinder;
  ClusterBuilder m_clusterBuilder;

  unsigned m_wsSize;
  std::vector<double> m_workspace;

  // energy flow 
  EnergyFlowFunctor m_functor;

  // output calibration files
  std::string m_hOutput;
  std::string m_tOutput;

  // test using averages rather than fits
  MIJNType  m_Eclusts;
  MIJNType  m_Tclusts;
  MIJType   m_Eavg;
  MIJType   m_Tavg;

};  // end calibrator class


}  // end Mono namespace

#endif
