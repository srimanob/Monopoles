// -*- C++ -*-
//
// Package:    MonoNtupleDumper
// Class:      MonoNtupleDumper
// 
/**\class MonoNtupleDumper MonoNtupleDumper.cc Monopoles/MonoNtupleDumper/src/MonoNtupleDumper.cc
   
   Description: [one line class summary]

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  Christopher Cowden
//         Created:  Tue Feb  7 16:21:08 CST 2012
// $Id: MonoNtupleDumper.cc,v 1.7 2013/06/13 21:45:27 cowden Exp $
//
// June 2017: Phat Srimanobhas: Modify for 80X framework
//
//
// July 2021 : Lin Shih:fix HepMC -> GenParticle
//      
// 2025: Matheus added pileUp information and new triggers for TRG Eff
//
// 2025: Thales included the Type-1 MET corrected via pat;
// 2025: Added L1 functions
//
// system include files
#include <vector>
#include <string>
#include <set>
#include <memory>
#include <algorithm>
#include <sstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// math
#include "DataFormats/Math/interface/deltaR.h"

// geometry
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"

// data formats
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/Common/interface/SortedCollection.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
//
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

// vertex
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

// Ecal includes
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
#include "RecoCaloTools/Navigation/interface/CaloRectangle.h"

// Hcal includes
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

// trigger includes
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

// Monopole analysis includes
//#include "Monopoles/MonoAlgorithms/interface/NPVHelper.h"
#include "Monopoles/MonoAlgorithms/interface/MonoDefs.h" 
#include "Monopoles/MonoAlgorithms/interface/MonoEcalObs0.h"
#include "Monopoles/MonoAlgorithms/interface/ClustCategorizer.h"
#include "Monopoles/MonoAlgorithms/interface/MonoTrackMatcher.h"
#include "Monopoles/MonoAlgorithms/interface/MonoGenTrackExtrapolator.h"
#include "Monopoles/TrackCombiner/interface/MplTracker.h"

//Track + Tracker
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"

//TrackerRecHit2D
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
#include "RecoTracker/DeDx/interface/DeDxTools.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"


// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TF1.h"
#include "TF2.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TFitResult.h"
#include "TVirtualFitter.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/ErrorFrameTransformer.h"
using namespace std;
using namespace edm;
typedef std::vector<Trajectory> TrajectoryCollection;

bool passFilter(const float eta,const float phi,const trigger::TriggerEvent& trigEvent,const std::string& filter,const float maxDR2);
int getPho175TrigCode(const float eta,const float phi,const trigger::TriggerEvent& trigEvent,const float maxDR2=0.2*0.2);
int getPho200TrigCode(const float eta,const float phi,const trigger::TriggerEvent& trigEvent,const float maxDR2=0.2*0.2);
bool evtPassesPho175L1(const trigger::TriggerEvent& trigEvent);
bool evtPassesPho200L1(const trigger::TriggerEvent& trigEvent);
bool evtPassesPFMET250L1(const trigger::TriggerEvent& trigEvent);
bool evtPassesPFMET300L1(const trigger::TriggerEvent& trigEvent);

#define DEBUG

//
// class declaration
//
//namespace edm { class ConsumesCollector; }

class MonoNtupleDumper : public edm::EDAnalyzer {
  
  typedef std::vector<reco::BasicCluster> BasicClusterCollection;
  typedef std::vector<reco::Photon> PhotonCollection;
  //typedef std::vector<reco::Electron> ElectronCollection;
  typedef std::vector<reco::GsfElectron> ElectronCollection;
  typedef reco::GsfElectron Electron;
  typedef reco::Photon Photon;
  
public:
  explicit MonoNtupleDumper(const edm::ParameterSet&);
  ~MonoNtupleDumper();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

  // clear tree variables
  void clear();
  
  //
  void rematch(); 
  //void fill_trigInfo(const edm::TriggerResults& triggerResults, const edm::TriggerNames& trigNames);

  // ----------member data ---------------------------
  
  std::string m_output;
  TFile *m_outputFile;
  
  // HLTConfigProvider
  HLTConfigProvider m_hltConfig;
  
  // Input tags
  edm::EDGetTokenT< edm::TriggerResults > m_hltResults;
  edm::EDGetTokenT< trigger::TriggerEvent > trigEvent;
//  edm::EDGetTokenT< edm::HepMCProduct > m_mcproduct;

  edm::EDGetTokenT< vector<reco::GenParticle> > m_genParticles;

  edm::EDGetTokenT< reco::VertexCollection > m_PVTag;
  edm::EDGetTokenT< EBRecHitCollection > m_TagEcalEB_RecHits;
  edm::EDGetTokenT< EERecHitCollection > m_TagEcalEE_RecHits;
  edm::EDGetTokenT< HBHERecHitCollection > m_TagHcalHBHE_RecHits;
  edm::EDGetTokenT< reco::TrackCollection > m_TrackTag;
  edm::EDGetTokenT< std::vector<Trajectory> > m_TrajectoryTag;
  edm::EDGetTokenT< TrajTrackAssociationCollection > m_assoMapTag;
  edm::EDGetTokenT< reco::PFJetCollection > m_Tag_Jets;
  edm::EDGetTokenT< PhotonCollection > m_Tag_Photons;
  edm::EDGetTokenT< ElectronCollection > m_Tag_Electrons;
  edm::EDGetTokenT< vector<reco::PFCandidate>> m_Tag_PF;
  edm::EDGetTokenT< std::vector<reco::PFMET> > m_Tag_MET;
//MET
  edm::EDGetTokenT< vector<reco::GenMET>> m_Tag_GenMET;
  edm::EDGetTokenT< vector<reco::CaloMET> > m_Tag_CaloMET;

  edm::EDGetTokenT< reco::BasicClusterCollection > m_Tag_bClusters;
  edm::EDGetTokenT< reco::BasicClusterCollection > m_Tag_cClusters;
  edm::EDGetTokenT< reco::BasicClusterCollection > m_Tag_combClusters;
  edm::EDGetTokenT< reco::BasicClusterCollection > m_Tag_eeClean;
  edm::EDGetTokenT< reco::BasicClusterCollection > m_Tag_eeUnclean;
  edm::EDGetTokenT< reco::BasicClusterCollection > m_Tag_eeComb;

  // Pileup 
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> m_puInfoToken;


  // PAT objects
  edm::EDGetTokenT< std::vector<pat::Jet>> m_Tag_PatJets;
  edm::EDGetTokenT< std::vector<pat::MET>> m_Tag_PatMETs;

  bool m_isData;

  // Monopole Ecal Observables
  Mono::MonoEcalObs0 m_ecalObs;
  
  edm::Handle<trigger::TriggerEvent> m_trigEventHandle;

  //NPV

  //Generator
  static const int MONO_PID = MONOID;

  //Tracker
  MplTracker *_Tracker;

  // TFileService
  //edm::Service<TFileService> m_fs;

  // map cluster category (lengthxwidth thing) to a histogram
  // showing the average energy or time  in each cell
  //TFileDirectory *m_avgDir;
  std::map<Mono::ClustCategorizer,TH2D *> m_clustEMap;
  std::map<Mono::ClustCategorizer,TH2D *> m_clustTMap;
  std::map<Mono::ClustCategorizer,unsigned> m_clustCatCount;

  std::vector<double> m_betas;
  std::vector<double> m_betaTs;

  //trigger infos
  bool passHLT_Photon175_ ; //2016
  bool passHLT_Photon200_ ; //2017, 2018
  bool passHLT_DoublePhoton60_ ; //2016
  bool passHLT_DoublePhoton70_ ; //2017, 2018
  bool passHLT_PFMET300_ ; //2016
  bool passHLT_PFMET170_HBHE_BeamHaloCleaned_; //2016
  bool passHLT_PFMET140_PFMHT140_IDTight_; //2017, 2018
  bool passHLT_PFMET250_HBHECleaned_; //2017, 2018
  bool passHLT_PFMET300_HBHECleaned_; //2017, 2018
  //
  bool passHLT_PFMET200_HBHE_BeamHaloCleaned_;
  bool passHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_;
  bool passHLT_CaloMET300_HBHECleaned_;

  // For Trigger Efficiency
  bool passHLT_IsoMu24_;
  bool passHLT_Ele27_WPTight_Gsf_;
  bool passHLT_PFJet500_;

  //
  TTree * m_tree;
  
  bool _ClustHitOutput;
  bool _EleJetPhoOutput;
  //bool _TrackHitOutput;

  // Event information
  unsigned m_run;
  unsigned m_lumi;
  unsigned m_event;
  unsigned m_NPV;
  bool m_isFirstEvent;
  unsigned m_NTrigs;
  std::vector<bool> m_trigResults;
  std::vector<std::string> m_trigNames;  
  bool m_evtPassesPho175L1;
  bool m_evtPassesPho200L1;
  bool m_evtPassesPFMET250L1;
  bool m_evtPassesPFMET300L1;
 
  // Ecal Observable information
  unsigned m_nClusters;
  std::vector<double> m_clust_E;
  std::vector<double> m_clust_eta;
  std::vector<double> m_clust_phi;
  std::vector<double> m_clust_L;
  std::vector<double> m_clust_W;
  std::vector<double> m_clust_N;
  std::vector<double> m_clust_sigEta;
  std::vector<double> m_clust_sigPhi;
  std::vector<double> m_clust_meanEta;
  std::vector<double> m_clust_meanPhi;
  std::vector<double> m_clust_skewEta;
  std::vector<double> m_clust_skewPhi;
  std::vector<double> m_clust_seedFrac;
  std::vector<double> m_clust_firstFrac;
  std::vector<double> m_clust_secondFrac;
  std::vector<double> m_clust_thirdFrac;
  std::vector<double> m_clust_phiStripFrac;
  std::vector<double> m_clust_matchDR;
  std::vector<double> m_clust_tagged;
  std::vector<double> m_clust_matchPID;
  std::vector<double> m_clust_matchTime;
  std::vector<double> m_clust_matchPt;
  std::vector<double> m_clust_hsE;
  std::vector<double> m_clust_hsTime;
  std::vector<int>    m_clust_hsInSeed;
  std::vector<int>    m_clust_hsWeird;
  std::vector<int>    m_clust_hsDiWeird;


  //int and bool for the new function . int for SC and bool is per event 


  // Treat these arrays as 3D arrays
  // There is space for 15 clusters of 100 total elements in each cluster
  // One must use m_nClusters, m_clust_L, and m_clust_W when unpacking
  // the cluster from the TTree.
  static const unsigned WS = 100;
  static const unsigned SS = 15*100;
  double m_clust_Ecells[1500]; 
  double m_clust_Tcells[1500]; 
  

  //June 25 2021 add e2x5 eLFTB...
  // energy in the 2x5 strip right of the max crystal (does not contain max crystal)
  // 2 crystals wide in eta, 5 wide in phi.
  // ---e2x5Max----
  // energy in a 2x5 strip containing the seed (max) crystal.
  // 2 crystals wide in eta, 5 wide in phi.
  // it is the maximum of either (1x5left + 1x5center) or (1x5right + 1x5center)
  // ---eLeft/eRight/eTop/eBottom---
  // energies in the crystal left, right, top, bottom w.r.t. to the most energetic crystal


  // Ecal hybrid clusters
  unsigned m_nClusterEgamma;
  std::vector<double> m_egClust_E;
  std::vector<double> m_egClust_size;
  std::vector<double> m_egClust_eta;
  std::vector<double> m_egClust_phi;
  std::vector<double> m_egClust_frac51;
  std::vector<double> m_egClust_frac15;
  std::vector<double> m_egClust_e55;
  std::vector<double> m_egClust_e2x5Right;
  std::vector<double> m_egClust_e2x5Left;
  std::vector<double> m_egClust_e2x5Top;
  std::vector<double> m_egClust_e2x5Bottom;
  std::vector<double> m_egClust_e2x5Max;
  std::vector<double> m_egClust_eLeft;
  std::vector<double> m_egClust_eRight;
  std::vector<double> m_egClust_eTop;
  std::vector<double> m_egClust_eBottom;    
  std::vector<double> m_egClust_eMax;
  std::vector<double> m_egClust_matchDR;
  std::vector<double> m_egClust_tagged;
  std::vector<double> m_egClust_matchPID;
  std::vector<double> m_egClust_hcalIso;
  std::vector<double> m_egClust_SwissCross;
  
  
  // Ecal hybrid clusters (cleaned collection)
  unsigned m_nCleanEgamma;
  std::vector<double> m_egClean_E;
  std::vector<double> m_egClean_size;
  std::vector<double> m_egClean_eta;
  std::vector<double> m_egClean_phi;
  std::vector<double> m_egClean_frac51;
  std::vector<double> m_egClean_frac15;
  std::vector<double> m_egClean_e55;
  std::vector<double> m_egClean_e2x5Right;
  std::vector<double> m_egClean_e2x5Left;
  std::vector<double> m_egClean_e2x5Top;
  std::vector<double> m_egClean_e2x5Bottom;
  std::vector<double> m_egClean_e2x5Max;
  std::vector<double> m_egClean_eLeft;
  std::vector<double> m_egClean_eRight;
  std::vector<double> m_egClean_eTop;
  std::vector<double> m_egClean_eBottom;
  std::vector<double> m_egClean_eMax;
  std::vector<double> m_egClean_matchDR;
  std::vector<double> m_egClean_tagged;
  std::vector<double> m_egClean_matchPID;
  std::vector<double> m_egClean_hcalIso;
  std::vector<double> m_egClean_SwissCross;
  
  // Ecal hybrid clusters (combined collection)
  unsigned m_nCombEgamma;
  std::vector<double> m_egComb_E;
  std::vector<double> m_egComb_size;
  std::vector<double> m_egComb_eta;
  std::vector<double> m_egComb_phi;
  std::vector<double> m_egComb_frac51;
  std::vector<double> m_egComb_frac15;
  std::vector<double> m_egComb_e55;
  std::vector<double> m_egComb_e99;
  std::vector<double> m_egComb_e2x5Right;
  std::vector<double> m_egComb_e2x5Left;
  std::vector<double> m_egComb_e2x5Top;
  std::vector<double> m_egComb_e2x5Bottom;
  std::vector<double> m_egComb_e2x5Max;
  std::vector<double> m_egComb_eLeft;
  std::vector<double> m_egComb_eRight;
  std::vector<double> m_egComb_eTop;
  std::vector<double> m_egComb_eBottom;
  std::vector<double> m_egComb_eMax;
  std::vector<double> m_egComb_e25Right;
  std::vector<double> m_egComb_e25Left;
  std::vector<double> m_egComb_e25Top;
  std::vector<double> m_egComb_e25Bottom;
  std::vector<double> m_egComb_e25Max;
  std::vector<double> m_egComb_matchDR;
  std::vector<double> m_egComb_tagged;
  std::vector<double> m_egComb_matchPID;
  std::vector<double> m_egComb_hcalIso;
  std::vector<double> m_egComb_SwissCross;
  
  // EE clusters (cleaned)
  unsigned m_nCleanEE;
  std::vector<double> m_eeClean_E;
  std::vector<double> m_eeClean_size;
  std::vector<double> m_eeClean_eta;
  std::vector<double> m_eeClean_phi;
  std::vector<double> m_eeClean_frac51;
  std::vector<double> m_eeClean_frac15;
  std::vector<double> m_eeClean_e55;
  std::vector<double> m_eeClean_e2x5Right;
  std::vector<double> m_eeClean_e2x5Left;
  std::vector<double> m_eeClean_e2x5Top;
  std::vector<double> m_eeClean_e2x5Bottom;
  std::vector<double> m_eeClean_e2x5Max;
  std::vector<double> m_eeClean_eLeft;
  std::vector<double> m_eeClean_eRight;
  std::vector<double> m_eeClean_eTop;
  std::vector<double> m_eeClean_eBottom;
  std::vector<double> m_eeClean_eMax;
  std::vector<double> m_eeClean_matchDR;
  std::vector<double> m_eeClean_tagged;
  std::vector<double> m_eeClean_matchPID;
  std::vector<double> m_eeClean_hcalIso;
  std::vector<double> m_eeClean_SwissCross;
  
  // EE clusters (uncleanOnly)
  unsigned m_nUncleanEE;
  std::vector<double> m_eeUnclean_E;
  std::vector<double> m_eeUnclean_size;
  std::vector<double> m_eeUnclean_eta;
  std::vector<double> m_eeUnclean_phi;
  std::vector<double> m_eeUnclean_frac51;
  std::vector<double> m_eeUnclean_frac15;
  std::vector<double> m_eeUnclean_e55;
  std::vector<double> m_eeUnclean_e2x5Right;
  std::vector<double> m_eeUnclean_e2x5Left;
  std::vector<double> m_eeUnclean_e2x5Top;
  std::vector<double> m_eeUnclean_e2x5Bottom;
  std::vector<double> m_eeUnclean_e2x5Max;
  std::vector<double> m_eeUnclean_eLeft;
  std::vector<double> m_eeUnclean_eRight;
  std::vector<double> m_eeUnclean_eTop;
  std::vector<double> m_eeUnclean_eBottom;
  std::vector<double> m_eeUnclean_eMax;
  std::vector<double> m_eeUnclean_matchDR;
  std::vector<double> m_eeUnclean_tagged;
  std::vector<double> m_eeUnclean_matchPID;
  std::vector<double> m_eeUnclean_hcalIso;
  std::vector<double> m_eeUnclean_SwissCross;
  
  // EE clusters (combined)
  unsigned m_nCombEE;
  std::vector<double> m_eeComb_E;
  std::vector<double> m_eeComb_size;
  std::vector<double> m_eeComb_eta;
  std::vector<double> m_eeComb_phi;
  std::vector<double> m_eeComb_frac51;
  std::vector<double> m_eeComb_frac15;
  std::vector<double> m_eeComb_e55;
  std::vector<double> m_eeComb_e99;
  std::vector<double> m_eeComb_e2x5Right;
  std::vector<double> m_eeComb_e2x5Left;
  std::vector<double> m_eeComb_e2x5Top;
  std::vector<double> m_eeComb_e2x5Bottom;
  std::vector<double> m_eeComb_e2x5Max;
  std::vector<double> m_eeComb_eLeft;
  std::vector<double> m_eeComb_eRight;
  std::vector<double> m_eeComb_eTop;
  std::vector<double> m_eeComb_eBottom;
  std::vector<double> m_eeComb_eMax;
  std::vector<double> m_eeComb_e25Left;
  std::vector<double> m_eeComb_e25Right;
  std::vector<double> m_eeComb_e25Top;
  std::vector<double> m_eeComb_e25Bottom;
  std::vector<double> m_eeComb_e25Max;
  std::vector<double> m_eeComb_matchDR;
  std::vector<double> m_eeComb_tagged;
  std::vector<double> m_eeComb_matchPID;
  std::vector<double> m_eeComb_hcalIso;
  std::vector<double> m_eeComb_SwissCross;
  
  // Ecal RecHits
  std::vector<double> m_ehit_eta;
  std::vector<double> m_ehit_phi;
  std::vector<double> m_ehit_time;
  std::vector<double> m_ehit_energy;
  std::vector<double> m_ehit_otEnergy;
  std::vector<double> m_ehit_flag;
  std::vector<double> m_ehit_kWeird;
  std::vector<double> m_ehit_kDiWeird;
  std::vector<double> m_ehit_chi2;

  std::vector<double> m_mono_ehit_eta;
  std::vector<double> m_mono_ehit_phi;
  std::vector<double> m_mono_ehit_time;
  std::vector<double> m_mono_ehit_energy;
  std::vector<double> m_mono_ehit_otEnergy;
  std::vector<double> m_mono_ehit_flag;
  std::vector<double> m_mono_ehit_kWeird;
  std::vector<double> m_mono_ehit_kDiWeird;


  std::vector<double> m_test_ehit_eta;
  std::vector<double> m_test_ehit_phi;
  std::vector<double> m_test_ehit_time;
  std::vector<double> m_test_ehit_energy;
  std::vector<double> m_test_ehit_otEnergy;
  std::vector<double> m_test_ehit_flag;
  std::vector<double> m_test_ehit_kWeird;
  std::vector<double> m_test_ehit_kDiWeird;

  
  // Jet information
  unsigned m_jet_N;
  std::vector<double> m_jet_E;
  std::vector<double> m_jet_p;
  std::vector<double> m_jet_pt;
  std::vector<double> m_jet_px;
  std::vector<double> m_jet_py;
  std::vector<double> m_jet_pz;
  std::vector<double> m_jet_eta;
  std::vector<double> m_jet_phi;
  std::vector<double> m_jet_matchDR;
  std::vector<int>  m_jet_tagged;
  std::vector<int>  m_jet_matchPID;
  
  // Photon information
  unsigned m_pho_N;
  std::vector<double> m_pho_E;
  std::vector<double> m_pho_p;
  std::vector<double> m_pho_pt;
  std::vector<double> m_pho_px;
  std::vector<double> m_pho_py;
  std::vector<double> m_pho_pz;
  std::vector<double> m_pho_eta;
  std::vector<double> m_pho_phi;
  std::vector<double> m_pho_matchDR;
  std::vector<int>  m_pho_tagged;
  std::vector<int>  m_pho_matchPID;
  
  // Electron information
  unsigned m_ele_N;
  std::vector<double> m_ele_E;
  std::vector<double> m_ele_p;
  std::vector<double> m_ele_pt;
  std::vector<double> m_ele_px;
  std::vector<double> m_ele_py;
  std::vector<double> m_ele_pz;
  std::vector<double> m_ele_eta;
  std::vector<double> m_ele_phi;
//  std::vector<int>    m_ele_isPF;
  std::vector<double> m_ele_matchDR;
  std::vector<int>  m_ele_tagged;
  std::vector<int>  m_ele_matchPID;
  
  // PFCandidate information
  //add by Lin 7 28 2021

  unsigned m_pf_N;
  std::vector<double> m_pf_E;
  std::vector<double> m_pf_p;
  std::vector<double> m_pf_pt;
  std::vector<double> m_pf_px;
  std::vector<double> m_pf_py;
  std::vector<double> m_pf_pz;
  std::vector<double> m_pf_eta;
  std::vector<double> m_pf_phi;
//  std::vector<reco::PFCandidate::ParticleType> m_pf_pdgId;
  std::vector<int> m_pf_pdgId;
  std::vector<int> m_pf_status;

  // Pileup
  double m_nTrueInteractions;
  
  // MET information
  double m_mpt;
  double m_mpPhi;

  double m_GenMpt;
  double m_GenMpx;
  double m_GenMpy;
  double m_GenMETPhi;

  double m_CaloMpt;
  double m_CaloMpx;
  double m_CaloMpy;
  double m_CaloMETPhi;

  // MET from PAT
  double m_pat_mpt;
  double m_pat_mpPhi;
  double m_pat_CaloMpt;
  double m_pat_CaloMpx;
  double m_pat_CaloMpy;
  double m_pat_CaloMETPhi;
  
  // Generator level Branches
  double m_mono_p;
  double m_mono_eta;
  double m_mono_phi;
  double m_mono_m;
  double m_mono_px;
  double m_mono_py;
  double m_mono_pz;
  //Lin added  July 2021
  double m_mono_pt;
  double m_mono_Et; 
  double m_mono_Et2; 
  double m_mono_E; 
  double m_mono_mt; 
  double m_mono_mt2; 
  int    m_mono_status;
  int    m_mono_pdgId; 
  size_t m_mono_Nmother;
  // int for monopole pass trigger 
  double m_monoExp_eta;
  double m_monoExp_phi;
  double m_monoExpEE_eta;
  double m_monoExpEE_phi;
  
  double m_amon_p;
  double m_amon_eta;
  double m_amon_phi;
  double m_amon_m;
  double m_amon_px;
  double m_amon_py;
  double m_amon_pz;
  double m_amon_pt;
  double m_amon_Et; 
  double m_amon_Et2; 
  double m_amon_E; 
  double m_amon_mt; 
  double m_amon_mt2; 
  int    m_amon_status;
  int    m_amon_pdgId; 
  size_t m_amon_Nmother;
  double m_amonExp_eta;
  double m_amonExp_phi;
  double m_amonExpEE_eta;
  double m_amonExpEE_phi; 

  
////////////////////////////////////////
  // Branches to hold the combined monopole objects
  unsigned m_nCandidates;
  std::vector<double> m_candDist;
  std::vector<double> m_candSubHits;
  std::vector<double> m_candSatSubHits;
  std::vector<double> m_canddEdXSig;
  std::vector<double> m_candTIso;
  std::vector<double> m_candSeedFrac;
  std::vector<double> m_candf15;
  std::vector<double> m_candE55;
  std::vector<double> m_candE99;
  std::vector<double> m_candSwissCross;
  std::vector<double> m_candHIso;
  std::vector<double> m_candXYPar0;
  std::vector<double> m_candXYPar1;
  std::vector<double> m_candXYPar2;
  std::vector<double> m_candRZPar0;
  std::vector<double> m_candRZPar1;
  std::vector<double> m_candRZPar2;
  std::vector<double> m_candEta;
  std::vector<double> m_candPhi;  
  std::vector<int>   m_candPho175TrigCode;
  std::vector<int>   m_candPho200TrigCode;
  //SC or monopole cand function for HLT 


  // Extra variables
  std::vector<double> m_cande2x5Right;
  std::vector<double> m_cande2x5Left;
  std::vector<double> m_cande2x5Top;
  std::vector<double> m_cande2x5Bottom;
  std::vector<double> m_cande2x5Max;



  std::vector<double> m_candeLeft;
  std::vector<double> m_candeRight;
  std::vector<double> m_candeTop;
  std::vector<double> m_candeBottom;
  std::vector<double> m_candeMax;
  std::vector<double> m_cande25Right;
  std::vector<double> m_cande25Left;




 

  



};

//
// constants, enums and typedefs
//

namespace chow {
  double mag ( double x, double y) {
    return sqrt( x*x + y*y );
  }
  double mag ( double x, double y, double z){
    return sqrt( x*x + y*y + z*z );
  }
}

//
// static data member definitions
//

//
// constructors and destructor
//
MonoNtupleDumper::MonoNtupleDumper(const edm::ParameterSet& iConfig)
  :m_output(iConfig.getParameter<std::string>("Output"))
  ,m_hltResults( consumes< edm::TriggerResults >( iConfig.getParameter< edm::InputTag >( "TriggerResults" ) ) )
  ,trigEvent(consumes<trigger::TriggerEvent>( iConfig.getParameter<edm::InputTag>("TriggerEvent")))
//  ,m_mcproduct( consumes< edm::HepMCProduct >( iConfig.getParameter< edm::InputTag >( "GeneratorTag" ) ) )
  ,m_genParticles( consumes <vector<reco::GenParticle>>( iConfig.getParameter< edm::InputTag >( "GeneratorTag" ) ) )
  ,m_PVTag( consumes< reco::VertexCollection >( iConfig.getParameter< edm::InputTag >( "PrimaryVertices" ) ) )
  ,m_TagEcalEB_RecHits( consumes< EBRecHitCollection >( iConfig.getParameter<edm::InputTag>("EcalEBRecHits") ) )
  ,m_TagEcalEE_RecHits( consumes< EERecHitCollection >( iConfig.getParameter<edm::InputTag>("EcalEERecHits") ) )
  ,m_TagHcalHBHE_RecHits( consumes< HBHERecHitCollection >( iConfig.getParameter<edm::InputTag>("HBHERecHits") ) )
  ,m_TrackTag( consumes< reco::TrackCollection >( iConfig.getParameter< edm::InputTag >( "TrackTag" ) ) )
  ,m_TrajectoryTag(consumes<std::vector<Trajectory> >(iConfig.getParameter<edm::InputTag>("TrackTag") ) )
  ,m_assoMapTag(consumes<TrajTrackAssociationCollection>(iConfig.getParameter<edm::InputTag>("TrackTag") ) )
  ,m_Tag_Jets( consumes< reco::PFJetCollection >( iConfig.getParameter<edm::InputTag>("JetTag") ) )
  ,m_Tag_Photons( consumes< PhotonCollection >( iConfig.getParameter<edm::InputTag>("PhotonTag") ) )
  ,m_Tag_Electrons( consumes< ElectronCollection >( iConfig.getParameter<edm::InputTag>("ElectronTag") ) )
  ,m_Tag_PF( consumes< vector<reco::PFCandidate> >( iConfig.getParameter<edm::InputTag>("PFTag") ) )
  ,m_Tag_MET( consumes< std::vector<reco::PFMET> >(iConfig.getParameter<edm::InputTag>("METTag") ) )
  ,m_Tag_GenMET( consumes< std::vector<reco::GenMET> >(iConfig.getParameter<edm::InputTag>("GenMETTag") ) )
  ,m_Tag_CaloMET( consumes< std::vector<reco::CaloMET> >(iConfig.getParameter<edm::InputTag>("CaloMETTag") ) )
  ,m_Tag_bClusters( consumes< reco::BasicClusterCollection >(iConfig.getParameter<edm::InputTag>("bcClusterTag") ) )
  ,m_Tag_cClusters( consumes< reco::BasicClusterCollection >(iConfig.getParameter<edm::InputTag>("ccClusterTag") ) )
  ,m_Tag_combClusters( consumes< reco::BasicClusterCollection >(iConfig.getParameter<edm::InputTag>("combClusterTag") ) )
  ,m_Tag_eeClean( consumes< reco::BasicClusterCollection >(iConfig.getParameter<edm::InputTag>("eeCleanTag") ) )
  ,m_Tag_eeUnclean( consumes< reco::BasicClusterCollection >(iConfig.getParameter<edm::InputTag>("eeUncleanTag") ) )
  ,m_Tag_eeComb( consumes< reco::BasicClusterCollection >(iConfig.getParameter<edm::InputTag>("eeCombTag") ) )
  ,m_puInfoToken( consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupInfoTag") ))
  ,m_Tag_PatJets( consumes< std::vector<pat::Jet> >(iConfig.getParameter<edm::InputTag>("PatJetTag")))
  ,m_Tag_PatMETs( consumes< std::vector<pat::MET> >(iConfig.getParameter<edm::InputTag>("PatMETTag")))
  ,m_isData(iConfig.getParameter<bool>("isData") )
  ,m_ecalObs(iConfig)
{
  m_isFirstEvent = true;
  _Tracker         = new MplTracker(iConfig,consumesCollector());
  _ClustHitOutput  = iConfig.getUntrackedParameter<bool>("ClustHitOutput", true);
  _EleJetPhoOutput = iConfig.getUntrackedParameter<bool>("EleJetPhoOutput", true);
}


MonoNtupleDumper::~MonoNtupleDumper()
{ 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void MonoNtupleDumper::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  clear();
  m_run = iEvent.id().run();
  m_lumi = iEvent.id().luminosityBlock();
  m_event = iEvent.id().event();

#ifdef DEBUG
  std::cout<<"Start ntuple production"<<std::endl;
  std::cout<<"PASS: EventSetup"<<std::endl;
#endif
 


  // Pileup - Implementation
  edm::Handle<std::vector<PileupSummaryInfo>> puInfoH;
  iEvent.getByToken(m_puInfoToken, puInfoH);
  
  m_nTrueInteractions = -1.0;
  if (puInfoH.isValid()) {
  for (auto & PVI : *puInfoH) {
      if (PVI.getBunchCrossing() == 0) {
        m_nTrueInteractions = PVI.getTrueNumInteractions();
        break;
       }
    }
 }
  
  ////////////////////////////////////
  // get a handle on the trigger results
  //void MonoNtupleDumper::set_trigInfo(TTree* tree){

  // tree->Branch("passHLT_Photon175" , &passHLT_Photon175_ , "passHLT_Photon175/O");
  // }

  //=============================================================
  //
  //             Trigger Info
  //    
  //=============================================================     
  
  edm::Handle < edm::TriggerResults > triggerResultsHandle;
  iEvent.getByToken(m_hltResults, triggerResultsHandle);
  const edm::TriggerResults triggerResults = *triggerResultsHandle.product();
  const edm::TriggerNames& trigNames = iEvent.triggerNames(triggerResults);
  passHLT_Photon175_ = 0;
  passHLT_Photon200_ = 0;
  passHLT_DoublePhoton60_ = 0;
  passHLT_DoublePhoton70_ = 0;
  passHLT_PFMET300_ = 0;
  passHLT_PFMET170_HBHE_BeamHaloCleaned_ = 0;
  passHLT_PFMET140_PFMHT140_IDTight_ = 0;
  passHLT_PFMET250_HBHECleaned_ = 0;
  passHLT_PFMET300_HBHECleaned_ = 0;
  passHLT_PFMET200_HBHE_BeamHaloCleaned_ = 0;
  passHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_ = 0;
  passHLT_CaloMET300_HBHECleaned_ = 0;

  // For trigger efficiency
  passHLT_IsoMu24_ = 0;
  passHLT_Ele27_WPTight_Gsf_ = 0;
  passHLT_PFJet500_ = 0;

  for (size_t i = 0; i < trigNames.size(); ++i) {
    const std::string &name = trigNames.triggerName(i);
    bool fired = triggerResults.accept(i);
    
    if(!fired) continue;
    passHLT_Photon175_ |= name.find("HLT_Photon175_v") != std::string::npos;
    passHLT_Photon200_ |= name.find("HLT_Photon200_v") != std::string::npos;
    passHLT_DoublePhoton60_ |= name.find("HLT_DoublePhoton60_v") != std::string::npos;
    passHLT_DoublePhoton70_ |= name.find("HLT_DoublePhoton70_v") != std::string::npos;
    passHLT_PFMET300_ |= name.find("HLT_PFMET300_v") != std::string::npos;
    passHLT_PFMET170_HBHE_BeamHaloCleaned_ |= name.find("HLT_PFMET170_HBHE_BeamHaloCleaned_v") != std::string::npos;
    passHLT_PFMET140_PFMHT140_IDTight_ |= name.find("HLT_PFMET140_PFMHT140_IDTight_v") != std::string::npos;
    passHLT_PFMET250_HBHECleaned_ |= name.find("HLT_PFMET250_HBHECleaned_v") != std::string::npos;
    passHLT_PFMET300_HBHECleaned_ |= name.find("HLT_PFMET300_HBHECleaned_v") != std::string::npos;
    passHLT_PFMET200_HBHE_BeamHaloCleaned_ |= name.find("HLT_PFMET200_HBHE_BeamHaloCleaned_v") != std::string::npos;
    passHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_ |= name.find("HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_v") != std::string::npos;
    passHLT_CaloMET300_HBHECleaned_ |= name.find("HLT_CaloMET300_HBHECleaned_v") != std::string::npos;
    
    passHLT_IsoMu24_ |= name.find("HLT_IsoMu24_v") != std::string::npos;
    passHLT_Ele27_WPTight_Gsf_ |= name.find("HLT_Ele27_WPTight_Gsf_v") != std::string::npos;
    passHLT_PFJet500_ |= name.find("HLT_PFJet500_v") != std::string::npos;
  }
#ifdef DEBUG
  cout<< "HLT_Photon175 = " <<passHLT_Photon175_<<endl;
  cout<< "HLT_Photon200 = " <<passHLT_Photon200_<<endl;
  cout<< "HLT_DoublePhoton60 = " <<passHLT_DoublePhoton60_<<endl;
  cout<< "HLT_DoublePhoton70 = " <<passHLT_DoublePhoton70_<<endl;
  cout<< "HLT_PFMET300 = " <<passHLT_PFMET300_<<endl;
  cout<< "HLT_PFMET170_HBHE_BeamHaloCleaned = " <<passHLT_PFMET170_HBHE_BeamHaloCleaned_<<endl;
  cout<< "HLT_PFMET140_PFMHT140_IDTight = " <<passHLT_PFMET140_PFMHT140_IDTight_<<endl;
  cout<< "HLT_PFMET250_HBHECleaned = " <<passHLT_PFMET250_HBHECleaned_<<endl;
  cout<< "HLT_PFMET300_HBHECleaned = " <<passHLT_PFMET300_HBHECleaned_<<endl;
  cout<< "HLT_PFMET200_HBHE_BeamHaloCleaned = " <<passHLT_PFMET200_HBHE_BeamHaloCleaned_<<endl;
  cout<< "HLT_PFMETTypeOne200_HBHE_BeamHaloCleaned = " <<passHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_<<endl;
  cout<< "HLT_CaloMET300_HBHECleaned = " <<passHLT_CaloMET300_HBHECleaned_<<endl;
  cout<< "HLT_IsoMu24 = " <<passHLT_IsoMu24_<<endl;
  cout<< "HLT_Ele27_WPTight_Gsf = " <<passHLT_Ele27_WPTight_Gsf_<<endl;
  cout<< "HLT_PFJet500 = " <<passHLT_PFJet500_<<endl;


#endif

  // get a handle on the trigger results                                                    
  Handle<edm::TriggerResults> HLTR;
  iEvent.getByToken(m_hltResults,HLTR);

  //get a list of trigger names                                                            
  //std::vector<std::string> hltNames = m_hltConfig.triggerNames();                         
  edm::TriggerNames hltNames(iEvent.triggerNames(*HLTR));
  const std::vector<std::string> & hltNameVec = hltNames.triggerNames();

  // fill trigger names and set the size of the trigger result vector                       
  if (m_isFirstEvent) {
    m_isFirstEvent=false;
    m_NTrigs = HLTR->size();
    m_trigResults.resize(m_NTrigs);
    m_trigNames.resize(m_NTrigs);
    for ( unsigned i=0; i != m_NTrigs; i++ )
      m_trigNames[i] = hltNameVec[i];
    //m_trigNames = hltNameVec;
  }
  
  //std::string pathName="HLT_Photon175_v";
  // bool passTrig=m_trigResults->accept(triggerNames.triggerIndex(pathName));         
 
  //Triggerevent handle 
  iEvent.getByToken(trigEvent, m_trigEventHandle);
  //  const trigger::TriggerEvent & trigEvent = *(handle_trigEvent.product());
  
  // cycle over trigger names                                                               
  const unsigned nNames = HLTR->size();
  //  std::string pathName="HLT_Photon175_v"; 
  for ( unsigned i=0; i != nNames; i++ ) {
    m_trigResults[i] = HLTR->accept(i);
    // m_trigResults = HLTR->accept(H);  
    
    // a very verbose debug statement!!                                                     
    //cout << m_trigNames[i] << " " << HLTR->accept(i) << endl;                             
    
    //cout << m_trigNames("HLT_Photon1756_v") << HLTR->accept(i) << endl;
    //std::string pathName="HLT_Photon175_v";
    //    cout << m_trigNames << " " << HLTR->accept(HLT_Photon175_v) << endl; 
  }
  m_evtPassesPho175L1 = evtPassesPho175L1(*m_trigEventHandle);
  m_evtPassesPho200L1 = evtPassesPho200L1(*m_trigEventHandle);
  m_evtPassesPFMET250L1 = evtPassesPFMET250L1(*m_trigEventHandle);
  m_evtPassesPFMET300L1 = evtPassesPFMET300L1(*m_trigEventHandle);
  
  /*
    Handle<edm::TriggerResults> triggerResultsHandle; //our trigger result object
    edm::InputTag trigResultsTag("TriggerResults","","HLT"); //make sure have correct process on MC
    //data process=HLT, MC depends, Spring11 is REDIGI311X
    iEvent.getByToken(m_hltResults,triggerResultsHandle);
    const edm::TriggerNames& trigNames = iEvent.triggerNames(triggerResults);   
    std::string pathName="HLT_Photon175_v";
    bool passTrig=triggerResults->accept(trigNames.triggerIndex(pathName)); 
  */
  
  /////////////////////////////////////
  // get NPV for this event
  //m_NPV = _NPVHelper->getNPV(iEvent,iSetup);
  m_NPV = 0;
  edm::Handle<reco::VertexCollection> handlePV;
  iEvent.getByToken(m_PVTag,handlePV);
  reco::VertexCollection::const_iterator pv = handlePV->begin();
  for ( ; pv != handlePV->end(); pv++ ) {
    if ( !pv->isFake() && pv->ndof() > 4.0 ) {
      ++m_NPV;
    }
  }
  
#ifdef DEBUG
  std::cout<<"Pass: NPV"<<std::endl;
#endif

  // execute observable calculations
  //(unused) double monoObs = m_ecalObs.calculate(iSetup,iEvent,&m_betas,&m_betaTs);
  //(unused) const Mono::EBmap & ebMap = m_ecalObs.ecalMap();

  //#ifdef DEBUG
  //std::cout<<"Pass: EBmap"<<std::endl;
  //#endif
  
  // limits in eta and phi of ebmap
  //(unused) const unsigned nEta = ebMap.nEta();
  //(unused) const unsigned nPhi = ebMap.nPhi();

  // initialize the MC monopole tagger to tag MC monopoels fo RECO objects
  //  Mono::GenMonoClusterTagger tagger(0.3);
  // tagger.initialize(iEvent,iSetup);

#ifdef DEBUG
  std::cout<<"Pass: Mono::GenMonoClusterTagger"<<std::endl;
#endif

  // get RecHit collection
  Handle<EBRecHitCollection > ecalRecHits;
  iEvent.getByToken(m_TagEcalEB_RecHits,ecalRecHits);
  assert( ecalRecHits->size() > 0 );

  // get EE RecHit collection
  Handle<EERecHitCollection > eeRecHits;
  iEvent.getByToken(m_TagEcalEE_RecHits,eeRecHits);
  assert( eeRecHits->size() > 0 );

  // get HB RecHit Collection
  Handle<HBHERecHitCollection> hbRecHits;
  iEvent.getByToken(m_TagHcalHBHE_RecHits,hbRecHits);
  assert( hbRecHits->size() > 0 );

  // get calo geometry and topology
  ESHandle<CaloGeometry> calo;
  iSetup.get<CaloGeometryRecord>().get(calo);
  const CaloGeometry *m_caloGeo = (const CaloGeometry*)calo.product();
  const CaloSubdetectorGeometry *geom = m_caloGeo->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);

  ESHandle<CaloTopology> topo;
  iSetup.get<CaloTopologyRecord>().get(topo);
  const CaloTopology * topology = (const CaloTopology*)topo.product();

  // get HE geometry and topology
  // get HB geometry and topology
  //SimpleCaloRecHitMetaCollection< HBHERecHitCollection > mhbrh(hbRecHits.product());
  EgammaHcalIsolation egIso(0.4,0.1,10.,10.,10.,10.,calo,*hbRecHits);

  /* Comment to test */
  // fill RecHit Branches
  EBRecHitCollection::const_iterator itHit = ecalRecHits->begin();
  for ( ; itHit != ecalRecHits->end(); itHit++ ) {

    EBDetId detId( (*itHit).id() );
    //const CaloCellGeometry *cell = geom->getGeometry( detId );
    //m_ehit_eta.push_back( cell->getPosition().eta() );
    //m_ehit_phi.push_back( cell->getPosition().phi() );

    m_ehit_eta.push_back( geom->getGeometry(detId)->getPosition().eta() );
    m_ehit_phi.push_back( geom->getGeometry(detId)->getPosition().phi() );
    m_ehit_energy.push_back( (*itHit).energy() );
    m_ehit_time.push_back( (*itHit).time() );
    // the outOfTimeEnergy method has been removed from the EcalRecHit class
    // in CMSSW_7.  I leave this commment here is a note/reminder this is something
    // I don't immediately know how to fix, but this analyzer is not used much
    // so it doesn't need to be fixed at the moment.
    //m_ehit_otEnergy.push_back( (*itHit).outOfTimeEnergy() );

    m_ehit_kWeird.push_back( (*itHit).checkFlag(EcalRecHit::kWeird) );
    m_ehit_kDiWeird.push_back( (*itHit).checkFlag(EcalRecHit::kDiWeird) );
    m_ehit_flag.push_back( (*itHit).recoFlag() );


    m_ehit_chi2.push_back( (*itHit).chi2() );

    //std::cout << "Chi2: " << (*itHit).chi2() << " RecoFlag: " << (*itHit).recoFlag() << std::endl;

    if ((*itHit).energy() > 4.0) {
      m_test_ehit_time.push_back( (*itHit).time() );
      m_test_ehit_energy.push_back( (*itHit).energy() );
      m_test_ehit_eta.push_back( geom->getGeometry(detId)->getPosition().eta() );
      m_test_ehit_phi.push_back( geom->getGeometry(detId)->getPosition().phi() );

      m_test_ehit_kWeird.push_back( (*itHit).checkFlag(EcalRecHit::kWeird) );
      m_test_ehit_kDiWeird.push_back( (*itHit).checkFlag(EcalRecHit::kDiWeird) );
      m_test_ehit_flag.push_back( (*itHit).recoFlag() );

    }


    if ((*itHit).energy() > 200.0) {
      m_mono_ehit_time.push_back( (*itHit).time() );
      m_mono_ehit_energy.push_back( (*itHit).energy() );
      m_mono_ehit_eta.push_back( geom->getGeometry(detId)->getPosition().eta() );
      m_mono_ehit_phi.push_back( geom->getGeometry(detId)->getPosition().phi() );

      m_mono_ehit_kWeird.push_back( (*itHit).checkFlag(EcalRecHit::kWeird) );
      m_mono_ehit_kDiWeird.push_back( (*itHit).checkFlag(EcalRecHit::kDiWeird) );
      m_mono_ehit_flag.push_back( (*itHit).recoFlag() );

    }



  }
  //*/








  // get BasicCluster Collection
  Handle<BasicClusterCollection> bClusters;
  iEvent.getByToken(m_Tag_bClusters,bClusters);
  const unsigned nbClusters = bClusters->size();

  EcalClusterTools ecalTool;
  EcalTools ecalTools;
  std::vector<int> exclFlags;
  std::vector<int> sevExcl;

  std::vector<const reco::CaloCluster *> ebUncleanClusters;

  //  tagger.clearTags();
  // if ( !m_isData && nbClusters ) tagger.tag(nbClusters,&(*bClusters)[0]);

  unsigned nClusterCount=0;
  for ( unsigned i=0; i != nbClusters; i++ ) {
  //  if ( (*bClusters)[i].energy() < 50. ) continue;

    ebUncleanClusters.push_back( &(*bClusters)[i] );

    nClusterCount++;
    m_egClust_E.push_back( (*bClusters)[i].energy() );
    m_egClust_size.push_back( (*bClusters)[i].size() );
    m_egClust_eta.push_back( (*bClusters)[i].eta() );
    m_egClust_phi.push_back( (*bClusters)[i].phi() );
    const float e55 = ecalTool.e5x5((*bClusters)[i],ecalRecHits.product(),topology);
    const float e51 = ecalTool.e5x1((*bClusters)[i],ecalRecHits.product(),topology);
    const float e15 = ecalTool.e1x5((*bClusters)[i],ecalRecHits.product(),topology);
    const float eMax = ecalTool.eMax((*bClusters)[i],ecalRecHits.product());
    const float e2x5Right  = ecalTool.e2x5Right((*bClusters)[i],ecalRecHits.product(),topology);
    const float e2x5Left   = ecalTool.e2x5Left((*bClusters)[i],ecalRecHits.product(),topology);
    const float e2x5Top    = ecalTool.e2x5Top((*bClusters)[i],ecalRecHits.product(),topology);
    const float e2x5Bottom = ecalTool.e2x5Bottom((*bClusters)[i],ecalRecHits.product(),topology);
    const float e2x5Max    = ecalTool.e2x5Max((*bClusters)[i],ecalRecHits.product(),topology);
    const float eLeft      = ecalTool.eLeft((*bClusters)[i],ecalRecHits.product(),topology);
    const float eRight     = ecalTool.eRight((*bClusters)[i],ecalRecHits.product(),topology);
    const float eTop       = ecalTool.eTop((*bClusters)[i],ecalRecHits.product(),topology);
    const float eBottom    = ecalTool.eBottom((*bClusters)[i],ecalRecHits.product(),topology);
    float SwissCross ;
    SwissCross = 1.-(eLeft+eRight+eTop+eBottom)/eMax;
    m_egClust_frac51.push_back( e51/e55 );
    m_egClust_frac15.push_back( e15/e55 );
    m_egClust_e55.push_back(e55);
    m_egClust_eMax.push_back(eMax/e55);
    m_egClust_hcalIso.push_back( egIso.getHcalESum((*bClusters)[i].position()) );
    m_egClust_SwissCross.push_back( SwissCross);
    if ( !m_isData ) {
      //      m_egClust_matchDR.push_back(tagger.matchDR()[i]);
      // m_egClust_tagged.push_back(tagger.tagResult()[i]);
      // m_egClust_matchPID.push_back(tagger.matchPID()[i]);
   
      // m_egClust_matchDR.push_back(tagger.matchDR()[i]);                                               
      // m_egClust_tagged.push_back(tagger.tagResult()[i]);                                              
      //      m_egCluㄕㄠst_matchPID.push_back(tagger.matchPID()[i]); 



    }
  }
  m_nClusterEgamma = nClusterCount;

  // get BasicCluster Collection (cleaned)
  Handle<BasicClusterCollection> cClusters;
  iEvent.getByToken(m_Tag_cClusters,cClusters);
  const unsigned ncClusters = cClusters->size();

  std::vector<const reco::CaloCluster *> ebCleanClusters;

  //  tagger.clearTags();
  // if ( !m_isData && ncClusters ) tagger.tag(ncClusters,&(*cClusters)[0]);

  nClusterCount=0;
  for ( unsigned i=0; i != ncClusters; i++ ) {
  //  if ( (*cClusters)[i].energy() < 50. ) continue;

    ebCleanClusters.push_back( &(*cClusters)[i] );

    nClusterCount++;
    m_egClean_E.push_back( (*cClusters)[i].energy() );
    m_egClean_size.push_back( (*cClusters)[i].size() );
    m_egClean_eta.push_back( (*cClusters)[i].eta() );
    m_egClean_phi.push_back( (*cClusters)[i].phi() );

    const float e55 = ecalTool.e5x5((*cClusters)[i],ecalRecHits.product(),topology);
    const float e51 = ecalTool.e5x1((*cClusters)[i],ecalRecHits.product(),topology);
    const float e15 = ecalTool.e1x5((*cClusters)[i],ecalRecHits.product(),topology);
    const float eMax = ecalTool.eMax((*cClusters)[i],ecalRecHits.product());
    const float e2x5Right  = ecalTool.e2x5Right((*cClusters)[i],ecalRecHits.product(),topology);
    const float e2x5Left   = ecalTool.e2x5Left((*cClusters)[i],ecalRecHits.product(),topology);
    const float e2x5Top    = ecalTool.e2x5Top((*cClusters)[i],ecalRecHits.product(),topology);
    const float e2x5Bottom = ecalTool.e2x5Bottom((*cClusters)[i],ecalRecHits.product(),topology);
    const float e2x5Max    = ecalTool.e2x5Max((*cClusters)[i],ecalRecHits.product(),topology);
    const float eLeft      = ecalTool.eLeft((*cClusters)[i],ecalRecHits.product(),topology);
    const float eRight     = ecalTool.eRight((*cClusters)[i],ecalRecHits.product(),topology);
    const float eTop       = ecalTool.eTop((*cClusters)[i],ecalRecHits.product(),topology);
    const float eBottom    = ecalTool.eBottom((*cClusters)[i],ecalRecHits.product(),topology);
    float SwissCross ;
    SwissCross = 1.-(eLeft+eRight+eTop+eBottom)/eMax;
    m_egClean_frac51.push_back( e51/e55 );
    m_egClean_frac15.push_back( e15/e55 );
    m_egClean_e55.push_back(e55);
    m_egClean_eMax.push_back(eMax/e55);
    m_egClean_hcalIso.push_back( egIso.getHcalESum((*cClusters)[i].position()) );
    m_egClean_SwissCross.push_back( SwissCross);

    if ( !m_isData ) {
      // m_egClean_matchDR.push_back(tagger.matchDR()[i]);
      // m_egClean_tagged.push_back(tagger.tagResult()[i]);
      // m_egClean_matchPID.push_back(tagger.matchPID()[i]);
    }
  }
  m_nCleanEgamma = nClusterCount;

  // get BasicCluster Collection (combined)
  Handle<reco::BasicClusterCollection> combClusters;
  iEvent.getByToken(m_Tag_combClusters,combClusters);
  const unsigned ncombClusters = combClusters->size();

  std::vector<const reco::CaloCluster *> ebClusters;
  //matching was commenting 
  // tagger.clearTags();//comment this if there is an error 
  // if ( !m_isData && ncombClusters ) tagger.tag(ncombClusters,&(*combClusters)[0]);//comment this is there is an error

  nClusterCount=0;
  for ( unsigned i=0; i != ncombClusters; i++ ) {
   // if ( (*combClusters)[i].energy() < 50. ) continue;

  // Step 1: Identify the seed crystal (maximum energy hit)
   DetId seedId = ecalTool.getMaximum((*combClusters)[i], ecalRecHits.product()).first;

   //std::cout << "seedId: " << seedId.rawId() << std::endl;

   // Step 2: Define a 9x9 rectangle around the seed
   int size = 4;
   CaloRectangleRange<DetId> rectangleRange(size, seedId, *topology);



    ebClusters.push_back( &(*combClusters)[i] );

    nClusterCount++;
    m_egComb_E.push_back( (*combClusters)[i].energy() );
    m_egComb_size.push_back( (*combClusters)[i].size() );
    m_egComb_eta.push_back( (*combClusters)[i].eta() );
    m_egComb_phi.push_back( (*combClusters)[i].phi() );

    const float e55 = ecalTool.e5x5((*combClusters)[i],ecalRecHits.product(),topology);
   // const float e99 = ecalTool.e9x9((*combClusters)[i],ecalRecHits.product(),topology);
   // float e99 = ecalTool.matrixEnergy((*combClusters)[i], ecalRecHits.product(),topology, 4, 4);
   // Step 4: Sum the energy within the 9x9 matrix
    float e99 = 0.0;
      for (const auto& detId : rectangleRange) {
        auto hit = ecalRecHits.product()->find(detId);
          if (hit != ecalRecHits.product()->end()) {
              e99 += hit->energy();
      }
    } 
    //std::cout << "e99: " << e99 << std::endl;
    //std::cout << "e55: " << e55 << std::endl;
    const float e51 = ecalTool.e5x1((*combClusters)[i],ecalRecHits.product(),topology);
    const float e15 = ecalTool.e1x5((*combClusters)[i],ecalRecHits.product(),topology);
    const float eMax = ecalTool.eMax((*combClusters)[i],ecalRecHits.product());
    const float e2x5Right  = ecalTool.e2x5Right((*combClusters)[i],ecalRecHits.product(),topology);
    const float e2x5Left   = ecalTool.e2x5Left((*combClusters)[i],ecalRecHits.product(),topology);
    const float e2x5Top    = ecalTool.e2x5Top((*combClusters)[i],ecalRecHits.product(),topology);
    const float e2x5Bottom = ecalTool.e2x5Bottom((*combClusters)[i],ecalRecHits.product(),topology);
    const float e2x5Max    = ecalTool.e2x5Max((*combClusters)[i],ecalRecHits.product(),topology);
    const float eLeft      = ecalTool.eLeft((*combClusters)[i],ecalRecHits.product(),topology);
    const float eRight     = ecalTool.eRight((*combClusters)[i],ecalRecHits.product(),topology);
    const float eTop       = ecalTool.eTop((*combClusters)[i],ecalRecHits.product(),topology);
    const float eBottom    = ecalTool.eBottom((*combClusters)[i],ecalRecHits.product(),topology);
    float SwissCross ;
    SwissCross = 1.-(eLeft+eRight+eTop+eBottom)/eMax;
    m_egComb_frac51.push_back( e51/e55 );
    m_egComb_frac15.push_back( e15/e55 );
    m_egComb_e55.push_back(e55);
    m_egComb_e99.push_back(e99);
    m_egComb_eMax.push_back(eMax/e55);
    m_egComb_e25Right.push_back(ecalTool.e2x5Right((*combClusters)[i],ecalRecHits.product(),topology));
    m_egComb_e25Left.push_back(ecalTool.e2x5Left((*combClusters)[i],ecalRecHits.product(),topology));
    m_egComb_e25Top.push_back(e2x5Top);
    m_egComb_e25Bottom.push_back(e2x5Bottom);
    m_egComb_e25Max.push_back(e2x5Max);
    m_egComb_eRight.push_back(eRight);
    m_egComb_eLeft.push_back(eLeft);
    m_egComb_eTop.push_back(eTop);
    m_egComb_eBottom.push_back(eBottom);

    m_egComb_hcalIso.push_back( egIso.getHcalESum((*combClusters)[i].position()) );
    m_egComb_SwissCross.push_back( SwissCross);

    if ( !m_isData ) {//comment the below three lines if error 
      //  m_egComb_matchDR.push_back(tagger.matchDR()[i]);
      // m_egComb_tagged.push_back(tagger.tagResult()[i]);
      //  m_egComb_matchPID.push_back(tagger.matchPID()[i]);
    }
  }
  m_nCombEgamma = nClusterCount;

  for (const auto& cluster : *combClusters) {
    // Loop over hits and fractions
    for (const auto& hitFractionPair : cluster.hitsAndFractions()) {
      const DetId& detid = hitFractionPair.first;
      float fraction = hitFractionPair.second;

      const EcalRecHit* recHit = nullptr;

      // Check if the hit is from the EB or EE
      if (detid.subdetId() == EcalBarrel) {
          auto it = ecalRecHits->find(detid);  // Iterator to the RecHit
          if (it != ecalRecHits->end()) {
              recHit = &(*it);  // Dereference the iterator and take the address
          }
      } else if (detid.subdetId() == EcalEndcap) {
          auto it = eeRecHits->find(detid);  // Iterator to the RecHit
          if (it != eeRecHits->end()) {
              recHit = &(*it);  // Dereference the iterator and take the address
          }
      }                                 
    }
  }
  
  //Handle<EBRecHitCollection > ecalRecHits;
  //Handle<EERecHitCollection > eeRecHits;




  //
  // ------------- EE clusters ----------------------
  // get basic clusters
  Handle<BasicClusterCollection> eeClean;
  iEvent.getByToken(m_Tag_eeClean,eeClean);
  const unsigned nEECleanClusters = eeClean->size();

  std::vector<const reco::CaloCluster *> eeCleanClusters;
  //comment the following three lines if there is an error 
  // tagger.clearTags();
  // Mono::GenMonoClusterTagger eetagger(0.3,false);
  // if ( !m_isData && nEECleanClusters ) eetagger.tag(nEECleanClusters,&(*eeClean)[0]);

  nClusterCount=0;
  for ( unsigned i=0; i != nEECleanClusters; i++ ) {
    //if ( (*eeClean)[i].energy() < 50. ) continue;

    eeCleanClusters.push_back( &(*eeClean)[i] );

    nClusterCount++;
    m_eeClean_E.push_back( (*eeClean)[i].energy() );
    m_eeClean_size.push_back( (*eeClean)[i].size() );
    m_eeClean_eta.push_back( (*eeClean)[i].eta() );
    m_eeClean_phi.push_back( (*eeClean)[i].phi() );

    const float e55 = ecalTool.e5x5((*eeClean)[i],eeRecHits.product(),topology);
    const float e51 = ecalTool.e5x1((*eeClean)[i],eeRecHits.product(),topology);
    const float e15 = ecalTool.e1x5((*eeClean)[i],eeRecHits.product(),topology);
    const float eMax = ecalTool.eMax((*eeClean)[i],eeRecHits.product());
    const float e2x5Right  = ecalTool.e2x5Right((*eeClean)[i],eeRecHits.product(),topology);
    const float e2x5Left   = ecalTool.e2x5Left((*eeClean)[i],eeRecHits.product(),topology);
    const float e2x5Top    = ecalTool.e2x5Top((*eeClean)[i],eeRecHits.product(),topology);
    const float e2x5Bottom = ecalTool.e2x5Bottom((*eeClean)[i],eeRecHits.product(),topology);
    const float e2x5Max    = ecalTool.e2x5Max((*eeClean)[i],eeRecHits.product(),topology);
    const float eLeft      = ecalTool.eLeft((*eeClean)[i],eeRecHits.product(),topology);
    const float eRight     = ecalTool.eRight((*eeClean)[i],eeRecHits.product(),topology);
    const float eTop       = ecalTool.eTop((*eeClean)[i],eeRecHits.product(),topology);
    const float eBottom    = ecalTool.eBottom((*eeClean)[i],eeRecHits.product(),topology);
    float SwissCross ;
    SwissCross = 1.-(eLeft+eRight+eTop+eBottom)/eMax;
    m_eeClean_frac51.push_back( e51/e55 );
    m_eeClean_frac15.push_back( e15/e55 );
    m_eeClean_e55.push_back(e55);
    m_eeClean_eMax.push_back(eMax/e55);
    m_eeClean_hcalIso.push_back( egIso.getHcalESum((*eeClean)[i].position()) );
    m_eeClean_SwissCross.push_back( SwissCross);

    if ( !m_isData ) {//comment if error 
      // m_eeClean_matchDR.push_back(eetagger.matchDR()[i]);
      // m_eeClean_tagged.push_back(eetagger.tagResult()[i]);
      // m_eeClean_matchPID.push_back(eetagger.matchPID()[i]);
    }
  }
  m_nCleanEE = nClusterCount;

  // unclean only clusters (EE)
  Handle<BasicClusterCollection> eeUnclean;
  iEvent.getByToken(m_Tag_eeUnclean,eeUnclean);
  const unsigned nEEUncleanClusters = eeUnclean->size();

  std::vector<const reco::CaloCluster *> eeUncleanClusters;
  //comment if error 
  // eetagger.clearTags();
  // if ( !m_isData && nEEUncleanClusters ) eetagger.tag(nEEUncleanClusters,&(*eeUnclean)[0]);

  nClusterCount=0;
  for ( unsigned i=0; i != nEEUncleanClusters; i++ ) {
   // if ( (*eeUnclean)[i].energy() < 50. ) continue;

    eeUncleanClusters.push_back( &(*eeUnclean)[i] );

    nClusterCount++;
    m_eeUnclean_E.push_back( (*eeUnclean)[i].energy() );
    m_eeUnclean_size.push_back( (*eeUnclean)[i].size() );
    m_eeUnclean_eta.push_back( (*eeUnclean)[i].eta() );
    m_eeUnclean_phi.push_back( (*eeUnclean)[i].phi() );

    const float e55 = ecalTool.e5x5((*eeUnclean)[i],eeRecHits.product(),topology);
    const float e51 = ecalTool.e5x1((*eeUnclean)[i],eeRecHits.product(),topology);
    const float e15 = ecalTool.e1x5((*eeUnclean)[i],eeRecHits.product(),topology);
    const float eMax = ecalTool.eMax((*eeUnclean)[i],eeRecHits.product());
    const float e2x5Right  = ecalTool.e2x5Right((*eeUnclean)[i],eeRecHits.product(),topology);
    const float e2x5Left   = ecalTool.e2x5Left((*eeUnclean)[i],eeRecHits.product(),topology);
    const float e2x5Top    = ecalTool.e2x5Top((*eeUnclean)[i],eeRecHits.product(),topology);
    const float e2x5Bottom = ecalTool.e2x5Bottom((*eeUnclean)[i],eeRecHits.product(),topology);
    const float e2x5Max    = ecalTool.e2x5Max((*eeUnclean)[i],eeRecHits.product(),topology);
    const float eLeft      = ecalTool.eLeft((*eeUnclean)[i],eeRecHits.product(),topology);
    const float eRight     = ecalTool.eRight((*eeUnclean)[i],eeRecHits.product(),topology);
    const float eTop       = ecalTool.eTop((*eeUnclean)[i],eeRecHits.product(),topology);
    const float eBottom    = ecalTool.eBottom((*eeUnclean)[i],eeRecHits.product(),topology);
    float SwissCross ;
    SwissCross = 1.-(eLeft+eRight+eTop+eBottom)/eMax;
    m_eeUnclean_frac51.push_back( e51/e55 );
    m_eeUnclean_frac15.push_back( e15/e55 );
    m_eeUnclean_e55.push_back(e55);
    m_eeUnclean_eMax.push_back(eMax/e55);
    m_eeUnclean_hcalIso.push_back( egIso.getHcalESum((*eeUnclean)[i].position()) );
    m_eeUnclean_SwissCross.push_back( SwissCross);

    if ( !m_isData ) {//commnet if error 
      // m_eeUnclean_matchDR.push_back(eetagger.matchDR()[i]);
      // m_eeUnclean_tagged.push_back(eetagger.tagResult()[i]);
      // m_eeUnclean_matchPID.push_back(eetagger.matchPID()[i]);
    }
  }
  m_nUncleanEE = nClusterCount;

  // combined EE clusters
  Handle<BasicClusterCollection> eeComb;
  iEvent.getByToken(m_Tag_eeComb,eeComb);
  const unsigned nEECombClusters = eeComb->size();

  std::vector<const reco::CaloCluster *> eeCombClusters;
  //commnet if erro r
  // eetagger.clearTags();
  // if ( !m_isData && nEECombClusters ) eetagger.tag(nEECombClusters,&(*eeComb)[0]);

  nClusterCount=0;
  for ( unsigned i=0; i != nEECombClusters; i++ ) {
    //if ( (*eeComb)[i].energy() < 50. ) continue;

    eeCombClusters.push_back( &(*eeComb)[i] );

    nClusterCount++;
    m_eeComb_E.push_back( (*eeComb)[i].energy() );
    m_eeComb_size.push_back( (*eeComb)[i].size() );
    m_eeComb_eta.push_back( (*eeComb)[i].eta() );
    m_eeComb_phi.push_back( (*eeComb)[i].phi() );

    // Step 1: Identify the seed crystal (maximum energy hit)
   DetId seedId = ecalTool.getMaximum((*eeComb)[i], eeRecHits.product()).first;

   //std::cout << "seedId: " << seedId.rawId() << std::endl;

   // Step 2: Define a 9x9 rectangle around the seed
   int size = 4;
   CaloRectangleRange<DetId> rectangleRange(size, seedId, *topology);

   // Step 4: Sum the energy within the 9x9 matrix
    float e99 = 0.0;
      for (const auto& detId : rectangleRange) {
        auto hit = ecalRecHits.product()->find(detId);
          if (hit != ecalRecHits.product()->end()) {
              e99 += hit->energy();
      }
    } 

    const float e55 = ecalTool.e5x5((*eeComb)[i],eeRecHits.product(),topology);
    const float e51 = ecalTool.e5x1((*eeComb)[i],eeRecHits.product(),topology);
    const float e15 = ecalTool.e1x5((*eeComb)[i],eeRecHits.product(),topology);
    const float eMax = ecalTool.eMax((*eeComb)[i],eeRecHits.product());
    const float e2x5Right  = ecalTool.e2x5Right((*eeComb)[i],eeRecHits.product(),topology);
    const float e2x5Left   = ecalTool.e2x5Left((*eeComb)[i],eeRecHits.product(),topology);
    const float e2x5Top    = ecalTool.e2x5Top((*eeComb)[i],eeRecHits.product(),topology);
    const float e2x5Bottom = ecalTool.e2x5Bottom((*eeComb)[i],eeRecHits.product(),topology);
    const float e2x5Max    = ecalTool.e2x5Max((*eeComb)[i],eeRecHits.product(),topology);
    const float eLeft      = ecalTool.eLeft((*eeComb)[i],eeRecHits.product(),topology);
    const float eRight     = ecalTool.eRight((*eeComb)[i],eeRecHits.product(),topology);
    const float eTop       = ecalTool.eTop((*eeComb)[i],eeRecHits.product(),topology);
    const float eBottom    = ecalTool.eBottom((*eeComb)[i],eeRecHits.product(),topology);
    float SwissCross ;
    SwissCross = 1.-(eLeft+eRight+eTop+eBottom)/eMax;
    m_eeComb_frac51.push_back( e51/e55 );
    m_eeComb_frac15.push_back( e15/e55 );
    m_eeComb_e55.push_back(e55);
    m_eeComb_e99.push_back(e99);
    m_eeComb_eMax.push_back(eMax/e55);
    m_eeComb_e25Right.push_back(ecalTool.e2x5Right((*eeComb)[i],eeRecHits.product(),topology));
    m_eeComb_e25Left.push_back(ecalTool.e2x5Left((*eeComb)[i],eeRecHits.product(),topology));
    m_eeComb_hcalIso.push_back( egIso.getHcalESum((*eeComb)[i].position()) );
    m_eeComb_e25Top.push_back(e2x5Top);
    m_eeComb_e25Bottom.push_back(e2x5Bottom);
    m_eeComb_e25Max.push_back(e2x5Max);

    m_eeComb_eRight.push_back(eRight);
    m_eeComb_eLeft.push_back(eLeft);
    m_eeComb_eTop.push_back(eTop);
    m_eeComb_eBottom.push_back(eBottom);

    m_eeComb_SwissCross.push_back( SwissCross);




    if ( !m_isData ) {
      //      m_eeComb_matchDR.push_back(eetagger.matchDR()[i]);
      // m_eeComb_tagged.push_back(eetagger.tagResult()[i]);
      // m_eeComb_matchPID.push_back(eetagger.matchPID()[i]);
    }
  }
  m_nCombEE = nClusterCount;
  
#ifdef DEBUG
  std::cout<<"Pass: Mono::GenMonoClusterTagger"<<std::endl;
#endif 

  ////////////////////////////////
  // Tracking analysis
  _Tracker->analyze(iEvent,iSetup);
  //  const Mono::MonoEcalCluster * clusters = clusterBuilder.clusters(); //comment if error
  //_Tracker->doMatch(m_nClusters,clusters,ebMap); //comment if error 
  _Tracker->doMatch(m_nCombEgamma,&ebClusters[0],fEBCombined);
  _Tracker->doMatch(m_nCleanEgamma,&ebCleanClusters[0],fEBClean);
  _Tracker->doMatch(m_nClusterEgamma,&ebUncleanClusters[0],fEBUnclean);
  _Tracker->doMatch(m_nCombEE,&eeCombClusters[0],fEECombined);
  _Tracker->doMatch(m_nCleanEE,&eeCleanClusters[0],fEEClean);
  _Tracker->doMatch(m_nUncleanEE,&eeUncleanClusters[0],fEEUnclean);
  
#ifdef DEBUG
  std::cout<<"Pass: Tracker"<<std::endl;
#endif

  ////////////////////////////////    
  // get jet collection
  Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(m_Tag_Jets,jets);
  const unsigned nJets = jets->size();

  // fill jet Branches
  //  tagger.clearTags();
  // if ( !m_isData && nJets ) tagger.tag(nJets,&(*jets)[0]);
  for ( unsigned i=0; i != nJets; i++ ) {

    const reco::PFJet & jet = (*jets)[i];

    m_jet_E.push_back( jet.energy() );
    m_jet_p.push_back( jet.p() );
    m_jet_pt.push_back( jet.pt() );
    m_jet_px.push_back( jet.px() );
    m_jet_py.push_back( jet.py() );
    m_jet_pz.push_back( jet.pz() );
    m_jet_eta.push_back( jet.eta() );
    m_jet_phi.push_back( jet.phi() );

    if ( !m_isData ) {
      //  m_jet_matchDR.push_back( tagger.matchDR()[i] );
      // m_jet_tagged.push_back(tagger.tagResult()[i]);
      // m_jet_matchPID.push_back(tagger.matchPID()[i]);
    }

  }
  m_jet_N = nJets;

  // get photon collection
  Handle<PhotonCollection> photons;
  iEvent.getByToken(m_Tag_Photons,photons);
  const unsigned nPhotons = photons->size();

  // fill photon Branches
  //  tagger.clearTags();
  // if ( !m_isData && nPhotons ) tagger.tag(nPhotons,&(*photons)[0]); 
  for ( unsigned i=0; i != nPhotons; i++ ) {
    
    const Photon & pho = (*photons)[i];

    m_pho_E.push_back( pho.energy() );
    m_pho_p.push_back( pho.p() );
    m_pho_pt.push_back( pho.pt() );
    m_pho_px.push_back( pho.px() );
    m_pho_py.push_back( pho.py() );
    m_pho_pz.push_back( pho.pz() );
    m_pho_eta.push_back( pho.eta() );
    m_pho_phi.push_back( pho.phi() );

    if ( !m_isData ) {
      //      m_pho_matchDR.push_back( tagger.matchDR()[i] );
      // m_pho_tagged.push_back( tagger.tagResult()[i] );
      // m_pho_matchPID.push_back( tagger.matchPID()[i] );
    }

  }
  m_pho_N = nPhotons;

  // get electron collection
  Handle<ElectronCollection> electrons;
  iEvent.getByToken(m_Tag_Electrons,electrons);
  const unsigned nElectrons = electrons->size();

  // fill electron Branches
  // tagger.clearTags();
  //  if ( !m_isData && nElectrons ) tagger.tag(nElectrons,&(*electrons)[0]);
  for (unsigned i=0; i != nElectrons; i++) {

    const Electron & ele = (*electrons)[i];

    m_ele_E.push_back( ele.energy() );
    m_ele_p.push_back( ele.p() );
    m_ele_pt.push_back( ele.pt() );
    m_ele_px.push_back( ele.px() );
    m_ele_py.push_back( ele.py() );
    m_ele_pz.push_back( ele.pz() );
    m_ele_eta.push_back( ele.eta() );
    m_ele_phi.push_back( ele.phi() );
    if ( !m_isData ) {
      //     m_ele_matchDR.push_back( tagger.matchDR()[i] );
      //  m_ele_tagged.push_back( tagger.tagResult()[i] );
      // m_ele_matchPID.push_back( tagger.matchPID()[i] );
    }

  }
  m_ele_N = nElectrons;

  // add by Lin 7 28 2021
  //get PFCandidtae 

  Handle< vector<reco::PFCandidate>> pfCandidate;
  iEvent.getByToken(m_Tag_PF,pfCandidate);
  const unsigned nPFCandidate = pfCandidate->size();
 /* for (unsigned  i=0; i != nPFCandidate; i++){

    const reco::PFCandidate &pf = (*pfCandidate)[i];
    m_pf_E.push_back( pf.energy() );
    m_pf_p.push_back( pf.p() );
    m_pf_pt.push_back( pf.pt() );
    m_pf_px.push_back( pf.px() );
    m_pf_py.push_back( pf.py() );
    m_pf_pz.push_back( pf.pz() );
    m_pf_eta.push_back( pf.eta() );
    m_pf_phi.push_back( pf.phi() );
    // reco::PFCandidate::ParticleType & type = (*pfCandidate)[i];
    m_pf_pdgId.push_back( pf.particleId());
    m_pf_status.push_back( pf.status());
  } */
  m_pf_N = nPFCandidate;

  // get MET collection
  Handle<std::vector<reco::PFMET> > met;
  iEvent.getByToken(m_Tag_MET,met);
  //assert (met[0].pt() > 150); // ALSO DO WITH PAT MET

  Handle<std::vector<reco::CaloMET> > Calomet;
  iEvent.getByToken(m_Tag_CaloMET,Calomet);
  
  // fill MET Branches
  m_mpt = (*met)[0].pt();
  m_mpPhi = (*met)[0].phi();
 


  m_CaloMpt = (*Calomet)[0].pt();
  m_CaloMpx = (*Calomet)[0].px();
  m_CaloMpy = (*Calomet)[0].py();
  m_CaloMETPhi = (*Calomet)[0].phi();

  Handle<std::vector<pat::MET> > patmet;
  iEvent.getByToken(m_Tag_PatMETs,patmet);
 
// Such filtering should be only applied when running on DATA 
 //if (patmet->at(0).pt() <= 150) {
  //   edm::LogInfo("LowMET") << "MET pt = " << patmet->at(0).pt() << " (below threshold)";
  //   return; // Skip this event as it doesn't meet the threshold
  //}  


  // fill MET Branches
  m_pat_mpt = (*patmet)[0].pt();
  m_pat_mpPhi = (*patmet)[0].phi();


  m_pat_CaloMpt = ((*patmet)[0]).caloMETPt();
  m_pat_CaloMpx = ((*patmet)[0]).caloMETP2().px;
  m_pat_CaloMpy = ((*patmet)[0]).caloMETP2().py;
  m_pat_CaloMETPhi = ((*patmet)[0]).caloMETPhi();

  // fill generator Branches
  if (!m_isData){
    Handle<std::vector<reco::GenMET> > Genmet;
    iEvent.getByToken(m_Tag_GenMET,Genmet);

    m_GenMpt = (*Genmet)[0].pt(); 
    m_GenMpx = (*Genmet)[0].px(); 
    m_GenMpy = (*Genmet)[0].py(); 
    m_GenMETPhi = (*Genmet)[0].phi(); 

    Mono::MonoGenTrackExtrapolator extrap; 

    //added by Lin July 21
    Handle< vector<reco::GenParticle> > genParticles;
    iEvent.getByToken(m_genParticles, genParticles);
    for(size_t i = 0; i < genParticles->size(); ++ i) {

      const reco::GenParticle & p = (*genParticles)[i];
  

      if(abs(p.pdgId())== MONO_PID &&  p.status()==1){
        if(p.pdgId()>0){  
          m_mono_p   = p.p();
          m_mono_eta = p.eta();
          m_mono_phi = p.phi();
          m_mono_m   = p.mass();
          m_mono_px  = p.px();
          m_mono_py  = p.py();
          m_mono_pz  = p.pz();
          m_mono_pt  = p.pt();
          m_mono_Et  = p.et();  
          m_mono_Et2 = p.et2();//transverse energy squared, use this for cut  
          m_mono_E   = p.energy();  
          m_mono_mt  = p.mt();//transverse mass
          m_mono_mt2 = p.mtSqr();//transverse mass squared
          m_mono_status   = p.status(); 
          m_mono_pdgId    = p.pdgId();
          #ifdef DEBUG
            std::cout<<"PDG ID : "<<p.pdgId()<<endl;
            std::cout<<"Status : "<<p.status()<<endl;
          #endif 
        }
        else{
          m_amon_p   = p.p();
          m_amon_eta = p.eta();
          m_amon_phi = p.phi();
          m_amon_m   = p.mass();
          m_amon_px  = p.px();
          m_amon_py  = p.py();
          m_amon_pz  = p.pz();
          m_amon_pt  = p.pt();
          m_amon_Et  = p.et();  
          m_amon_Et2 = p.et2();//transverse energy squared, use this for cut  
          m_amon_E   = p.energy();  
          m_amon_mt  = p.mt();//transverse mass
          m_amon_mt2 = p.mtSqr();//transverse mass squared
          m_amon_status   = p.status(); 
          m_amon_pdgId    = p.pdgId();
          #ifdef DEBUG
            std::cout<<"PDG ID : "<<p.pdgId()<<endl;
            std::cout<<"Status : "<<p.status()<<endl;
          #endif
        }
      }
    }
  }
  //////////////////////////////////////////////////////////
  // Perform track to cluster matching and form monopole candidates
  //
  rematch();

  // fill tree, must go last in this function
  m_tree->Fill();

}

// ------------ method called once each job just before starting event loop  ------------
void 
MonoNtupleDumper::beginJob()
{
  m_outputFile = new TFile(m_output.c_str(), "recreate");

  //m_avgDir = new TFileDirectory( m_fs->mkdir("avgClusterMaps"));

  m_tree = new TTree("monopoles","Monopole Variables");

//  gInterpreter->GenerateDictionary("vector<reco::PFcandidate::ParticleType>", "vector");
  m_tree->Branch("run",&m_run,"run/i");
  m_tree->Branch("lumiBlock",&m_lumi,"lumiBlock/i");
  m_tree->Branch("event",&m_event,"evnet/i");

  m_tree->Branch("NPV",&m_NPV,"NPV/i");

  // trigger information
  m_tree->Branch("trigResult",&m_trigResults);
  m_tree->Branch("trigNames",&m_trigNames);
  m_tree->Branch("passHLT_Photon175" , &passHLT_Photon175_ , "passHLT_Photon175/O");
  m_tree->Branch("passHLT_Photon200" , &passHLT_Photon200_ , "passHLT_Photon200/O");
  m_tree->Branch("passHLT_DoublePhoton60" , &passHLT_DoublePhoton60_ , "passHLT_DoublePhoton60/O");
  m_tree->Branch("passHLT_DoublePhoton70" , &passHLT_DoublePhoton70_ , "passHLT_DoublePhoton70/O");
  m_tree->Branch("passHLT_PFMET300" , &passHLT_PFMET300_ , "passHLT_PFMET300/O");
  m_tree->Branch("passHLT_PFMET170_HBHE_BeamHaloCleaned" , &passHLT_PFMET170_HBHE_BeamHaloCleaned_ , "passHLT_PFMET170_HBHE_BeamHaloCleaned/O");
  m_tree->Branch("passHLT_PFMET140_PFMHT140_IDTight" , &passHLT_PFMET140_PFMHT140_IDTight_ , "passHLT_PFMET140_PFMHT140_IDTight/O");
  m_tree->Branch("passHLT_PFMET250_HBHECleaned" , &passHLT_PFMET250_HBHECleaned_ , "passHLT_PFMET250_HBHECleaned/O");
  m_tree->Branch("passHLT_PFMET300_HBHECleaned" , &passHLT_PFMET300_HBHECleaned_ , "passHLT_PFMET300_HBHECleaned/O");
  m_tree->Branch("passHLT_PFMET200_HBHE_BeamHaloCleaned" , &passHLT_PFMET200_HBHE_BeamHaloCleaned_ , "passHLT_PFMET200_HBHE_BeamHaloCleaned/O");
  m_tree->Branch("passHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned" , &passHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned_ , "passHLT_PFMETTypeOne200_HBHE_BeamHaloCleaned/O");
  m_tree->Branch("passHLT_CaloMET300_HBHECleaned" , &passHLT_CaloMET300_HBHECleaned_ , "passHLT_CaloMET300_HBHECleaned/O"); 
  // trigger efficiency paths
  m_tree->Branch("passHLT_IsoMu24" , &passHLT_IsoMu24_ , "passHLT_IsoMu24/O");
  m_tree->Branch("passHLT_Ele27_WPTight_Gsf" , &passHLT_Ele27_WPTight_Gsf_ , "passHLT_Ele27_WPTight_Gsf/O");
  m_tree->Branch("passHLT_PFJet500" , &passHLT_PFJet500_ , "passHLT_PFJet500/O");

  //L1
  m_tree->Branch("evtPassesPho175L1" , &m_evtPassesPho175L1);
  m_tree->Branch("evtPassesPho200L1" , &m_evtPassesPho200L1);
  m_tree->Branch("evtPassesPFMET250L1" , &m_evtPassesPFMET250L1);
  m_tree->Branch("evtPassesPFMET300L1" , &m_evtPassesPFMET300L1);
  
  // combined candidates
  m_tree->Branch("cand_N",&m_nCandidates,"cand_N/i");
  m_tree->Branch("cand_dist",&m_candDist);
  m_tree->Branch("cand_SubHits",&m_candSubHits);
  m_tree->Branch("cand_SatSubHits",&m_candSatSubHits);
  m_tree->Branch("cand_dEdXSig",&m_canddEdXSig);
  m_tree->Branch("cand_TIso",&m_candTIso);
  m_tree->Branch("cand_f51",&m_candSeedFrac);
  m_tree->Branch("cand_f15",&m_candf15);
  m_tree->Branch("cand_e55",&m_candE55);
  m_tree->Branch("cand_e99",&m_candE99);
  m_tree->Branch("cand_SwissCross",&m_candSwissCross);
  m_tree->Branch("cand_HIso",&m_candHIso);
  m_tree->Branch("cand_XYPar0",&m_candXYPar0);
  m_tree->Branch("cand_XYPar1",&m_candXYPar1);
  m_tree->Branch("cand_XYPar2",&m_candXYPar2);
  m_tree->Branch("cand_RZPar0",&m_candRZPar0);
  m_tree->Branch("cand_RZPar1",&m_candRZPar1);
  m_tree->Branch("cand_RZPar2",&m_candRZPar2);
  m_tree->Branch("cand_eta",&m_candEta);
  m_tree->Branch("cand_phi",&m_candPhi);
  m_tree->Branch("cand_pho175TrigCode",&m_candPho175TrigCode);
  m_tree->Branch("cand_pho200TrigCode",&m_candPho200TrigCode);


  // Extra variables
  m_tree->Branch("cand_e2x5Right",&m_cande2x5Right);
  m_tree->Branch("cand_e2x5Left",&m_cande2x5Left);
  m_tree->Branch("cand_e2x5Top",&m_cande2x5Top);
  m_tree->Branch("cand_e2x5Bottom",&m_cande2x5Bottom);
  m_tree->Branch("cand_e2x5Max",&m_cande2x5Max);


  m_tree->Branch("cand_eRight",&m_candeRight);
  m_tree->Branch("cand_eLeft",&m_candeLeft);
  m_tree->Branch("cand_eTop",&m_candeTop);
  m_tree->Branch("cand_eBottom",&m_candeBottom);
  m_tree->Branch("cand_eMax",&m_candeMax);


 
//  m_tree->Branch("clust_N",&m_nClusters,"clust_N/i");
//  m_tree->Branch("clust_E",&m_clust_E);
//  m_tree->Branch("clust_eta",&m_clust_eta);
//  m_tree->Branch("clust_phi",&m_clust_phi);
//  m_tree->Branch("clust_L",&m_clust_L);
//  m_tree->Branch("clust_W",&m_clust_W);
//  m_tree->Branch("clust_sigEta",&m_clust_sigEta);
//  m_tree->Branch("clust_sigPhi",&m_clust_sigPhi);
//  m_tree->Branch("clust_skewEta",&m_clust_skewEta);
//  m_tree->Branch("clust_skewPhi",&m_clust_skewPhi);
//  m_tree->Branch("clust_seedFrac",&m_clust_seedFrac);
//  m_tree->Branch("clust_firstFrac",&m_clust_firstFrac);
//  m_tree->Branch("clust_secondFrac",&m_clust_secondFrac);
//  m_tree->Branch("clust_thirdFrac",&m_clust_thirdFrac);
//  m_tree->Branch("clust_matchDR",&m_clust_matchDR);
//  m_tree->Branch("clust_matchTime",&m_clust_matchTime);
//  m_tree->Branch("clust_matchPt",&m_clust_matchPt);
//  m_tree->Branch("clust_matchPID",&m_clust_matchPID);
//  m_tree->Branch("clust_tagged",&m_clust_tagged);
//  m_tree->Branch("clust_hsE",&m_clust_hsE);
//  m_tree->Branch("clust_hsTime",&m_clust_hsTime);
//  m_tree->Branch("clust_hsInSeed",&m_clust_hsInSeed);
//  m_tree->Branch("clust_hsWeird",&m_clust_hsWeird);
//  m_tree->Branch("clust_hsDiWeird",&m_clust_hsDiWeird);
//   if(_ClustHitOutput){
//   m_tree->Branch("clust_Ecells",&m_clust_Ecells,"clust_Ecells[1500]/D");
//   m_tree->Branch("clust_Tcells",&m_clust_Tcells,"clust_Tcells[1500]/D");
//  }
   m_tree->Branch("egClust_N",&m_nCleanEgamma,"egClust_N/i");
   m_tree->Branch("egClust_E",&m_egClust_E);
   m_tree->Branch("egClust_size",&m_egClust_size);
   m_tree->Branch("egClust_eta",&m_egClust_eta);
   m_tree->Branch("egClust_phi",&m_egClust_phi);
   m_tree->Branch("egClust_frac51",&m_egClust_frac51);
   m_tree->Branch("egClust_frac15",&m_egClust_frac15);
   m_tree->Branch("egClust_e55",&m_egClust_e55);
   m_tree->Branch("egClust_eMax",&m_egClust_eMax);
   m_tree->Branch("egClust_matchDR",&m_egClust_matchDR);
   m_tree->Branch("egClust_matchPID",&m_egClust_matchPID);
   m_tree->Branch("egClust_tagged",&m_egClust_tagged);
   m_tree->Branch("egClust_hcalIso",&m_egClust_hcalIso);
   m_tree->Branch("egClust_SwissCross",&m_egClust_SwissCross); 

   m_tree->Branch("egClean_N",&m_nClusterEgamma,"egClean_N/i");
   m_tree->Branch("egClean_E",&m_egClean_E);
   m_tree->Branch("egClean_size",&m_egClean_size);
   m_tree->Branch("egClean_eta",&m_egClean_eta);
   m_tree->Branch("egClean_phi",&m_egClean_phi);
   m_tree->Branch("egClean_frac51",&m_egClean_frac51);
   m_tree->Branch("egClean_frac15",&m_egClean_frac15);
   m_tree->Branch("egClean_e55",&m_egClean_e55);
   m_tree->Branch("egClean_eMax",&m_egClean_eMax);
   m_tree->Branch("egClean_matchDR",&m_egClean_matchDR);
   m_tree->Branch("egClean_matchPID",&m_egClean_matchPID);
   m_tree->Branch("egClean_tagged",&m_egClean_tagged);
   m_tree->Branch("egClean_hcalIso",&m_egClean_hcalIso);
   m_tree->Branch("egClean_SwissCross",&m_egClean_SwissCross); 
   
  m_tree->Branch("egComb_N",&m_nCombEgamma,"egComb_N/i");
  m_tree->Branch("egComb_E",&m_egComb_E);
  m_tree->Branch("egComb_size",&m_egComb_size);
  m_tree->Branch("egComb_eta",&m_egComb_eta);
  m_tree->Branch("egComb_phi",&m_egComb_phi);
  m_tree->Branch("egComb_frac51",&m_egComb_frac51);
  m_tree->Branch("egComb_frac15",&m_egComb_frac15);
  m_tree->Branch("egComb_e55",&m_egComb_e55);
  m_tree->Branch("egComb_e99",&m_egComb_e99);
  m_tree->Branch("egComb_eRight",&m_egComb_eRight);
  m_tree->Branch("egComb_eLeft",&m_egComb_eLeft);
  m_tree->Branch("egComb_eTop",&m_egComb_eTop);
  m_tree->Branch("egComb_eBottom",&m_egComb_eBottom);
  m_tree->Branch("egComb_eMax",&m_egComb_eMax);
  m_tree->Branch("egComb_e25Right",&m_egComb_e25Right);
  m_tree->Branch("egComb_e25Left",&m_egComb_e25Left);
  m_tree->Branch("egComb_e25Top",&m_egComb_e25Top);
  m_tree->Branch("egComb_e25Bottom",&m_egComb_e25Bottom);
  m_tree->Branch("egComb_e25Max",&m_egComb_e25Max);
  m_tree->Branch("egComb_matchDR",&m_egComb_matchDR);
  m_tree->Branch("egComb_matchPID",&m_egComb_matchPID);
  m_tree->Branch("egComb_tagged",&m_egComb_tagged);
  m_tree->Branch("egComb_hcalIso",&m_egComb_hcalIso);
   m_tree->Branch("egComb_SwissCross",&m_egComb_SwissCross); 
  
   m_tree->Branch("eeClean_N",&m_nCleanEE,"eeClean_N/i");
   m_tree->Branch("eeClean_E",&m_eeClean_E);
   m_tree->Branch("eeClean_size",&m_eeClean_size);
   m_tree->Branch("eeClean_eta",&m_eeClean_eta);
   m_tree->Branch("eeClean_phi",&m_eeClean_phi);
   m_tree->Branch("eeClean_frac51",&m_eeClean_frac51);
   m_tree->Branch("eeClean_frac15",&m_eeClean_frac15);
   m_tree->Branch("eeClean_eMax",&m_eeClean_eMax);
   m_tree->Branch("eeClean_e55",&m_eeClean_e55);
   m_tree->Branch("eeClean_matchDR",&m_eeClean_matchDR);
   m_tree->Branch("eeClean_matchPID",&m_eeClean_matchPID);
   m_tree->Branch("eeClean_tagged",&m_eeClean_tagged);
   m_tree->Branch("eeClean_hcalIso",&m_eeClean_hcalIso);
   m_tree->Branch("eeClean_SwissCross",&m_eeClean_SwissCross); 

   m_tree->Branch("eeUnclean_N",&m_nUncleanEE,"eeUnclean_N/i");
   m_tree->Branch("eeUnclean_E",&m_eeUnclean_E);
   m_tree->Branch("eeUnclean_size",&m_eeUnclean_size);
   m_tree->Branch("eeUnclean_eta",&m_eeUnclean_eta);
   m_tree->Branch("eeUnclean_phi",&m_eeUnclean_phi);
   m_tree->Branch("eeUnclean_frac51",&m_eeUnclean_frac51);
   m_tree->Branch("eeUnclean_frac15",&m_eeUnclean_frac15);
   m_tree->Branch("eeUnclean_eMax",&m_eeUnclean_eMax);
   m_tree->Branch("eeUnclean_e55",&m_eeUnclean_e55);
   m_tree->Branch("eeUnclean_matchDR",&m_eeUnclean_matchDR);
   m_tree->Branch("eeUnclean_matchPID",&m_eeUnclean_matchPID);
   m_tree->Branch("eeUnclean_tagged",&m_eeUnclean_tagged);
   m_tree->Branch("eeUnclean_hcalIso",&m_eeUnclean_hcalIso);
   m_tree->Branch("eeUnclean_SwissCross",&m_eeUnclean_SwissCross); 
  
  m_tree->Branch("eeComb_N",&m_nCombEE,"eeComb_N/i");
  m_tree->Branch("eeComb_E",&m_eeComb_E);
  m_tree->Branch("eeComb_size",&m_eeComb_size);
  m_tree->Branch("eeComb_eta",&m_eeComb_eta);
  m_tree->Branch("eeComb_phi",&m_eeComb_phi);
  m_tree->Branch("eeComb_frac51",&m_eeComb_frac51);
  m_tree->Branch("eeComb_frac15",&m_eeComb_frac15);
  m_tree->Branch("eeComb_eMax",&m_eeComb_eMax);
  m_tree->Branch("eeComb_eRight",&m_eeComb_eRight);
  m_tree->Branch("eeComb_eLeft",&m_eeComb_eLeft);
  m_tree->Branch("eeComb_eTop",&m_eeComb_eTop);
  m_tree->Branch("eeComb_eBottom",&m_eeComb_eBottom);
  m_tree->Branch("eeComb_e55",&m_eeComb_e55);
  m_tree->Branch("eeComb_e99",&m_eeComb_e99);
  m_tree->Branch("eeComb_e25Left",&m_eeComb_e25Left);
  m_tree->Branch("eeComb_e25Right",&m_eeComb_e25Right);
  m_tree->Branch("eeComb_e25Top",&m_eeComb_e25Top);
  m_tree->Branch("eeComb_e25Bottom",&m_eeComb_e25Bottom);
  m_tree->Branch("eeComb_e25Max",&m_eeComb_e25Max);
  m_tree->Branch("eeComb_matchDR",&m_eeComb_matchDR);
  m_tree->Branch("eeComb_matchPID",&m_eeComb_matchPID);
  m_tree->Branch("eeComb_tagged",&m_eeComb_tagged);
  m_tree->Branch("eeComb_hcalIso",&m_eeComb_hcalIso);
   m_tree->Branch("eeComb_SwissCross",&m_eeComb_SwissCross); 
  
   if(_ClustHitOutput){
   m_tree->Branch("ehit_eta",&m_ehit_eta);
   m_tree->Branch("ehit_phi",&m_ehit_phi);
   m_tree->Branch("ehit_time",&m_ehit_time);
   m_tree->Branch("ehit_E",&m_ehit_energy);
   m_tree->Branch("ehit_kWeird",&m_ehit_kWeird);
   m_tree->Branch("ehit_kDiWeird",&m_ehit_kDiWeird);
   m_tree->Branch("ehit_flag",&m_ehit_flag);
   m_tree->Branch("ehit_chi2",&m_ehit_chi2);
   

   m_tree->Branch("mono_ehit_eta",&m_mono_ehit_eta);
   m_tree->Branch("mono_ehit_phi",&m_mono_ehit_phi);
   m_tree->Branch("mono_ehit_time",&m_mono_ehit_time);
   m_tree->Branch("mono_ehit_energy",&m_mono_ehit_energy);
   m_tree->Branch("mono_ehit_kWeird",&m_mono_ehit_kWeird);
   m_tree->Branch("mono_ehit_kDiWeird",&m_mono_ehit_kDiWeird);
   m_tree->Branch("mono_ehit_flag",&m_mono_ehit_flag);


   m_tree->Branch("test_ehit_eta",&m_test_ehit_eta);
   m_tree->Branch("test_ehit_phi",&m_test_ehit_phi);
   m_tree->Branch("test_ehit_time",&m_test_ehit_time);
   m_tree->Branch("test_ehit_energy",&m_test_ehit_energy);
   m_tree->Branch("test_ehit_kWeird",&m_test_ehit_kWeird);
   m_tree->Branch("test_ehit_kDiWeird",&m_test_ehit_kDiWeird);
   m_tree->Branch("test_ehit_flag",&m_test_ehit_flag);



   }

  //_Tracker->beginJob(m_tree);
  
   if(_EleJetPhoOutput){
  m_tree->Branch("jet_N",&m_jet_N,"jet_N/i");
  m_tree->Branch("jet_E",&m_jet_E);
  m_tree->Branch("jet_p",&m_jet_p);
  m_tree->Branch("jet_pt",&m_jet_pt);
  m_tree->Branch("jet_px",&m_jet_px);
  m_tree->Branch("jet_py",&m_jet_py);
  m_tree->Branch("jet_pz",&m_jet_pz);
  m_tree->Branch("jet_eta",&m_jet_eta);
  m_tree->Branch("jet_phi",&m_jet_phi);
//  m_tree->Branch("jet_matchDR",&m_jet_matchDR);
//  m_tree->Branch("jet_tagged",&m_jet_tagged);
//  m_tree->Branch("jet_matchPID",&m_jet_matchPID);

  m_tree->Branch("pho_N",&m_pho_N,"pho_N/i");
  m_tree->Branch("pho_E",&m_pho_E);
  m_tree->Branch("pho_p",&m_pho_p);
  m_tree->Branch("pho_pt",&m_pho_pt);
  m_tree->Branch("pho_px",&m_pho_px);
  m_tree->Branch("pho_py",&m_pho_py);
  m_tree->Branch("pho_pz",&m_pho_pz);
  m_tree->Branch("pho_eta",&m_pho_eta);
  m_tree->Branch("pho_phi",&m_pho_phi);
//  m_tree->Branch("pho_matchDR",&m_pho_matchDR);
//  m_tree->Branch("pho_tagged",&m_pho_tagged);
//  m_tree->Branch("pho_matchPID",&m_pho_matchPID);

  m_tree->Branch("ele_N",&m_ele_N,"ele_N/i");
  m_tree->Branch("ele_E",&m_ele_E);
  m_tree->Branch("ele_p",&m_ele_p);
  m_tree->Branch("ele_pt",&m_ele_pt);
  m_tree->Branch("ele_px",&m_ele_px);
  m_tree->Branch("ele_py",&m_ele_py);
  m_tree->Branch("ele_pz",&m_ele_pz);
  m_tree->Branch("ele_eta",&m_ele_eta);
  m_tree->Branch("ele_phi",&m_ele_phi);
//  m_tree->Branch("ele_matchDR",&m_ele_matchDR);
//  m_tree->Branch("ele_tagged",&m_ele_tagged);
//  m_tree->Branch("ele_matchPID",&m_ele_matchPID);

  //add by Lin 7 28 2021

  m_tree->Branch("pf_N",&m_pf_N,"pf_N/i");
  m_tree->Branch("pf_E",&m_pf_E);
  m_tree->Branch("pf_p",&m_pf_p);
  m_tree->Branch("pf_pt",&m_pf_pt);
  m_tree->Branch("pf_px",&m_pf_px);
  m_tree->Branch("pf_py",&m_pf_py);
  m_tree->Branch("pf_pz",&m_pf_pz);
  m_tree->Branch("pf_eta",&m_pf_eta);
  m_tree->Branch("pf_phi",&m_pf_phi);
  m_tree->Branch("pf_pdgId",&m_pf_pdgId);
  m_tree->Branch("pf_status",&m_pf_status);
  }

  m_tree->Branch("mpt_pt",&m_mpt,"mpt_pt/D");
  m_tree->Branch("mpt_phi",&m_mpPhi,"mpt_phi/D");

  m_tree->Branch("GenMET_pt",&m_GenMpt,"GenMET_pt/D");
  m_tree->Branch("GenMET_px",&m_GenMpx,"GenMET_px/D");
  m_tree->Branch("GenMET_py",&m_GenMpy,"GenMET_py/D");
  m_tree->Branch("GenMET_phi",&m_GenMETPhi,"GenMET_phi/D");
  m_tree->Branch("CaloMET_pt",&m_CaloMpt,"CaloMET_pt/D");
  m_tree->Branch("CaloMET_px",&m_CaloMpx,"CaloMET_px/D");
  m_tree->Branch("CaloMET_py",&m_CaloMpy,"CaloMET_py/D");
  m_tree->Branch("CaloMET_phi",&m_CaloMETPhi,"CaloMET_phi/D");

  m_tree->Branch("PAT_mpt_pt",     &m_pat_mpt,  "PAT_mpt_pt/D");//PAT
  m_tree->Branch("PAT_mpt_phi",    &m_pat_mpPhi,"PAT_mpt_phi/D");//PAT
  m_tree->Branch("PAT_CaloMET_pt", &m_pat_CaloMpt,    "PAT_CaloMET_pt/D");
  m_tree->Branch("PAT_CaloMET_px", &m_pat_CaloMpx,    "PAT_CaloMET_px/D");
  m_tree->Branch("PAT_CaloMET_py", &m_pat_CaloMpy,    "PAT_CaloMET_py/D");
  m_tree->Branch("PAT_CaloMET_phi",&m_pat_CaloMETPhi,"PAT_CaloMET_phi/D");

 
  if ( !m_isData ) {
    m_tree->Branch("mono_eta", &m_mono_eta, "mono_eta/D");
    m_tree->Branch("mono_phi", &m_mono_phi, "mono_phi/D");
    m_tree->Branch("mono_m", &m_mono_m, "mono_m/D");
    m_tree->Branch("mono_mt",&m_mono_mt,"mono_mt/D");
    m_tree->Branch("mono_mt2",&m_mono_mt2,"mono_mt2/D");
    m_tree->Branch("mono_px",&m_mono_px, "mono_px/D");
    m_tree->Branch("mono_py",&m_mono_py, "mono_py/D");
    m_tree->Branch("mono_pz",&m_mono_pz, "mono_pz/D");
    m_tree->Branch("mono_pt",&m_mono_pt,"mono_pt/D");
    m_tree->Branch("mono_p", &m_mono_p, "mono_p/D");
    m_tree->Branch("mono_Et",&m_mono_Et,"mono_Et/D");
    m_tree->Branch("mono_Et2",&m_mono_Et2,"mono_Et2/D");
    m_tree->Branch("mono_E",&m_mono_E,"mono_E/D");
    m_tree->Branch("mono_status",&m_mono_status,"mono_status/D");
    m_tree->Branch("mono_pdgId",&m_mono_pdgId,"mono_pdgId/D");
    
    
    m_tree->Branch("amon_eta", &m_amon_eta, "amon_eta/D");
    m_tree->Branch("amon_phi", &m_amon_phi, "amon_phi/D");
    m_tree->Branch("amon_m", &m_amon_m, "amon_m/D");
    m_tree->Branch("amon_mt",&m_amon_mt,"amon_mt/D");
    m_tree->Branch("amon_mt2",&m_amon_mt2,"amon_mt2/D");
    m_tree->Branch("amon_px",&m_amon_px, "amon_px/D");
    m_tree->Branch("amon_py",&m_amon_py, "amon_py/D");
    m_tree->Branch("amon_pz",&m_amon_pz, "amon_pz/D");
    m_tree->Branch("amon_pt",&m_amon_pt,"amon_pt/D");
    m_tree->Branch("amon_p", &m_amon_p, "amon_p/D");
    m_tree->Branch("amon_Et",&m_amon_Et,"amon_Et/D");
    m_tree->Branch("amon_Et2",&m_amon_Et2,"amon_Et2/D");
    m_tree->Branch("amon_E",&m_amon_E,"amon_E/D");
    m_tree->Branch("amon_status",&m_amon_status,"amon_status/D");
    m_tree->Branch("amon_pdgId",&m_amon_pdgId,"amon_pdgId/D");

  } 
}

// ------------ method called once each job just after ending the event loop  ------------
void 
MonoNtupleDumper::endJob() 
{
  m_outputFile->cd();
  m_tree->Write();

  for(std::map<Mono::ClustCategorizer,TH2D*>::iterator i = m_clustEMap.begin(); i != m_clustEMap.end(); i++){
    i->second->Write();
  }
  for(std::map<Mono::ClustCategorizer,TH2D*>::iterator i = m_clustTMap.begin(); i != m_clustTMap.end(); i++){
    i->second->Write();
  }
  
  //_Tracker->endJob();

  m_outputFile->Close();
}

// ------------ method called when starting to processes a run  ------------
void 
MonoNtupleDumper::beginRun(edm::Run const&, edm::EventSetup const&)
{
}



void MonoNtupleDumper::clear()
{ 
  m_run = 0;
  m_lumi = 0;
  m_event = 0;
  
  m_NPV = 0;
  
  m_nClusters = 0;
  m_clust_E.clear();
  m_clust_eta.clear();
  m_clust_phi.clear();
  m_clust_L.clear();
  m_clust_W.clear();
  m_clust_sigEta.clear();
  m_clust_sigPhi.clear();
  m_clust_skewEta.clear();
  m_clust_skewPhi.clear();
  m_clust_seedFrac.clear();
  m_clust_firstFrac.clear();
  m_clust_secondFrac.clear();
  m_clust_thirdFrac.clear();
  m_clust_matchDR.clear();
  m_clust_matchTime.clear();
  m_clust_matchPt.clear();
  m_clust_matchPID.clear();
  m_clust_tagged.clear();
  m_clust_hsE.clear();
  m_clust_hsTime.clear();
  m_clust_hsInSeed.clear();
  m_clust_hsWeird.clear();
  m_clust_hsDiWeird.clear();
  for ( unsigned i=0; i != SS; i++ ){
    m_clust_Ecells[i] = -999.;
    m_clust_Tcells[i] = -999.;
  }
  
  m_nClusterEgamma = 0;
  m_egClust_E.clear();
  m_egClust_size.clear();
  m_egClust_eta.clear();
  m_egClust_phi.clear();
  m_egClust_frac51.clear();
  m_egClust_frac15.clear();
  m_egClust_e55.clear();
  m_egClust_e2x5Right.clear();
  m_egClust_e2x5Left.clear();
  m_egClust_e2x5Top.clear();
  m_egClust_e2x5Bottom.clear();
  m_egClust_e2x5Max.clear();
  m_egClust_eLeft.clear();
  m_egClust_eRight.clear();
  m_egClust_eTop.clear();
  m_egClust_eBottom.clear();    
  m_egClust_eMax.clear();
  m_egClust_matchDR.clear();
  m_egClust_matchPID.clear();
  m_egClust_tagged.clear();
  m_egClust_hcalIso.clear();
  m_egClust_SwissCross.clear();
  
  m_nCleanEgamma = 0;
  m_egClean_E.clear();
  m_egClean_size.clear();
  m_egClean_eta.clear();
  m_egClean_phi.clear();
  m_egClean_frac51.clear();
  m_egClean_frac15.clear();
  m_egClean_e55.clear();
  m_egClean_e2x5Right.clear();
  m_egClean_e2x5Left.clear();
  m_egClean_e2x5Top.clear();
  m_egClean_e2x5Bottom.clear();
  m_egClean_e2x5Max.clear();
  m_egClean_eLeft.clear();
  m_egClean_eRight.clear();
  m_egClean_eTop.clear();
  m_egClean_eBottom.clear();    
  m_egClean_eMax.clear();
  m_egClean_matchDR.clear();
  m_egClean_matchPID.clear();
  m_egClean_tagged.clear();
  m_egClean_hcalIso.clear();
  m_egClean_SwissCross.clear();
  
  m_nCombEgamma = 0;
  m_egComb_E.clear();
  m_egComb_size.clear();
  m_egComb_eta.clear();
  m_egComb_phi.clear();
  m_egComb_frac51.clear();
  m_egComb_frac15.clear();
  m_egComb_e55.clear();
  m_egComb_e99.clear();
  m_egComb_e2x5Right.clear();
  m_egComb_e2x5Left.clear();
  m_egComb_e2x5Top.clear();
  m_egComb_e2x5Bottom.clear();
  m_egComb_e2x5Max.clear();
  m_egComb_eLeft.clear();
  m_egComb_eRight.clear();
  m_egComb_eTop.clear();
  m_egComb_eBottom.clear();    
  m_egComb_eMax.clear();
  m_egComb_e25Right.clear();
  m_egComb_e25Left.clear();
  m_egComb_e25Top.clear();
  m_egComb_e25Bottom.clear();
  m_egComb_e25Max.clear();
  m_egComb_matchDR.clear();
  m_egComb_matchPID.clear();
  m_egComb_tagged.clear();
  m_egComb_hcalIso.clear();
  m_egComb_SwissCross.clear();
  

  m_nCleanEE = 0;
  m_eeClean_E.clear();
  m_eeClean_size.clear();
  m_eeClean_eta.clear();
  m_eeClean_phi.clear();
  m_eeClean_frac51.clear();
  m_eeClean_frac15.clear();
  m_eeClean_eMax.clear();
  m_eeClean_e55.clear();
  m_eeClean_e2x5Right.clear();
  m_eeClean_e2x5Left.clear();
  m_eeClean_e2x5Top.clear();
  m_eeClean_e2x5Bottom.clear();
  m_eeClean_e2x5Max.clear();
  m_eeClean_eLeft.clear();
  m_eeClean_eRight.clear();
  m_eeClean_eTop.clear();
  m_eeClean_eBottom.clear();    
  m_eeClean_matchDR.clear();
  m_eeClean_matchPID.clear();
  m_eeClean_tagged.clear();
  m_eeClean_hcalIso.clear();
  m_eeClean_SwissCross.clear();
  
  m_nUncleanEE = 0;
  m_eeUnclean_E.clear();
  m_eeUnclean_size.clear();
  m_eeUnclean_eta.clear();
  m_eeUnclean_phi.clear();
  m_eeUnclean_frac51.clear();
  m_eeUnclean_frac15.clear();
  m_eeUnclean_eMax.clear();
  m_eeUnclean_e55.clear();
  m_eeUnclean_e2x5Right.clear();
  m_eeUnclean_e2x5Left.clear();
  m_eeUnclean_e2x5Top.clear();
  m_eeUnclean_e2x5Bottom.clear();
  m_eeUnclean_e2x5Max.clear();
  m_eeUnclean_eLeft.clear();
  m_eeUnclean_eRight.clear();
  m_eeUnclean_eTop.clear();
  m_eeUnclean_eBottom.clear();    
  m_eeUnclean_matchDR.clear();
  m_eeUnclean_matchPID.clear();
  m_eeUnclean_tagged.clear();
  m_eeUnclean_hcalIso.clear();
  m_eeUnclean_SwissCross.clear();
  
  m_nCombEE = 0;
  m_eeComb_E.clear();
  m_eeComb_size.clear();
  m_eeComb_eta.clear();
  m_eeComb_phi.clear();
  m_eeComb_frac51.clear();
  m_eeComb_frac15.clear();
  m_eeComb_eMax.clear();
  m_eeComb_e55.clear();
  m_eeComb_e99.clear();
  m_eeComb_e2x5Right.clear();
  m_eeComb_e2x5Left.clear();
  m_eeComb_e2x5Top.clear();
  m_eeComb_e2x5Bottom.clear();
  m_eeComb_e2x5Max.clear();
  m_eeComb_eLeft.clear();
  m_eeComb_eRight.clear();
  m_eeComb_eTop.clear();
  m_eeComb_eBottom.clear();    
  m_eeComb_e25Left.clear();
  m_eeComb_e25Right.clear();
  m_eeComb_e25Top.clear();
  m_eeComb_e25Bottom.clear();
  m_eeComb_e25Max.clear();
  m_eeComb_matchDR.clear();
  m_eeComb_matchPID.clear();
  m_eeComb_tagged.clear();
  m_eeComb_hcalIso.clear();
  m_eeComb_SwissCross.clear();
  
  // Ecal RecHits
  m_ehit_eta.clear();
  m_ehit_phi.clear();
  m_ehit_time.clear();
  m_ehit_energy.clear();
  m_ehit_otEnergy.clear();
  m_ehit_kWeird.clear();
  m_ehit_kDiWeird.clear();
  m_ehit_flag.clear();
  m_ehit_chi2.clear();

  m_test_ehit_eta.clear();
  m_test_ehit_phi.clear();
  m_test_ehit_time.clear();
  m_test_ehit_energy.clear();
  m_test_ehit_otEnergy.clear();
  m_test_ehit_kWeird.clear();
  m_test_ehit_kDiWeird.clear();
  m_test_ehit_flag.clear();


  m_mono_ehit_eta.clear();
  m_mono_ehit_phi.clear();
  m_mono_ehit_time.clear();
  m_mono_ehit_energy.clear();
  m_mono_ehit_otEnergy.clear();
  m_mono_ehit_kWeird.clear();
  m_mono_ehit_kDiWeird.clear();
  m_mono_ehit_flag.clear();
  
  // Jet information
  m_jet_N = 0;
  m_jet_E.clear();
  m_jet_p.clear();
  m_jet_pt.clear();
  m_jet_px.clear();
  m_jet_py.clear();
  m_jet_pz.clear();
  m_jet_eta.clear();
  m_jet_phi.clear();
  m_jet_matchDR.clear();
  m_jet_tagged.clear();
  m_jet_matchPID.clear();
  
  // Photon information
  m_pho_N = 0;
  m_pho_E.clear();
  m_pho_p.clear();
  m_pho_pt.clear();
  m_pho_px.clear();
  m_pho_py.clear();
  m_pho_pz.clear();
  m_pho_eta.clear();
  m_pho_phi.clear();
  
  // Electron information
  m_ele_N = 0;
  m_ele_E.clear();
  m_ele_p.clear();
  m_ele_pt.clear();
  m_ele_px.clear();
  m_ele_py.clear();
  m_ele_pz.clear();
  m_ele_eta.clear();
  m_ele_phi.clear();
  //m_ele_isPF.clear();

  //PF information

  m_pf_N = 0;
  m_pf_E.clear();
  m_pf_p.clear();
  m_pf_pt.clear();
  m_pf_px.clear();
  m_pf_py.clear();
  m_pf_pz.clear();
  m_pf_eta.clear();
  m_pf_phi.clear();
  m_pf_pdgId.clear();
  m_pf_status.clear();

  m_mpt = 0.;
  m_mpPhi = 0.;
 m_GenMpt=0.;
 m_GenMpx=0.;
 m_GenMpy=0.;
 m_GenMETPhi=0.;
               
 m_CaloMpt=0.;
 m_CaloMpx=0.;
 m_CaloMpy=0.;
 m_CaloMETPhi=0.;

  m_pat_mpt=0.;
  m_pat_mpPhi=0.;
  m_pat_CaloMpt=0.;
  m_pat_CaloMpx=0.;
  m_pat_CaloMpy=0.;
  m_pat_CaloMETPhi=0.;
 
  m_mono_p = 0.;
  m_mono_eta = 0.;
  m_mono_phi = 0.;
  m_mono_m = 0.;
  m_mono_px = 0.;
  m_mono_py = 0.;
  m_mono_pz = 0.;
  m_mono_pt = 0.;
  m_mono_Et = 0.; 
  m_mono_Et2 = 0.; 
  m_mono_E = 0.; 
  m_mono_mt = 0.; 
  m_mono_mt2 = 0.; 
  m_mono_status = 0;
  m_mono_pdgId = 0; 
  m_monoExp_eta = 0.;
  m_monoExp_phi = 0.;

  m_amon_p = 0.;
  m_amon_eta = 0.;
  m_amon_phi = 0.;
  m_amon_m = 0.;
  m_amon_px = 0.;
  m_amon_py = 0.;
  m_amon_pz = 0.;
  m_amon_pt = 0.;
  m_amon_Et = 0.; 
  m_amon_Et2 = 0.; 
  m_amon_E = 0.; 
  m_amon_mt = 0.; 
  m_amon_mt2 = 0.; 
  m_amon_status = 0;
  m_amon_pdgId = 0; 
  m_amonExp_eta = 0.;
  m_amonExp_phi = 0.;

  // combined candidates
  m_nCandidates = 0U;
  m_candDist.clear();
  m_candSubHits.clear();
  m_candSatSubHits.clear();
  m_canddEdXSig.clear();
  m_candTIso.clear();
  m_candSeedFrac.clear();
  m_candf15.clear();
  m_candE55.clear();
  m_candE99.clear();
  m_cande2x5Right.clear();
  m_cande2x5Left.clear();
  m_cande2x5Top.clear();
  m_cande2x5Max.clear();
  m_cande2x5Bottom.clear();

  m_candeLeft.clear();
  m_candeRight.clear();
  m_candeTop.clear();
  m_candeBottom.clear();
  m_candeMax.clear();

  m_candSwissCross.clear();
  m_candHIso.clear();
  m_candXYPar0.clear();
  m_candXYPar1.clear();
  m_candXYPar2.clear();
  m_candRZPar0.clear();
  m_candRZPar1.clear();
  m_candRZPar2.clear();
  m_candEta.clear();
  m_candPhi.clear();
  m_candPho175TrigCode.clear();
  m_candPho200TrigCode.clear();
}

void MonoNtupleDumper::rematch()
{
  const double EcalR = 129.0;

  // Get a handle on the monopole tracks
  std::vector<Mono::MonoTrack> tracks;
  _Tracker->getTracks(tracks);

  const std::vector<int> & subHits = _Tracker->getSubHits();
  const std::vector<int> & satSubHits = _Tracker->getSatSubHits();
  const std::vector<float> & tIso = _Tracker->getTIso();
  const std::vector<float> & nDof = _Tracker->getNdof();

  const unsigned nTracks = tracks.size();

  for(unsigned i=0; i < nTracks; i++){

    // if Ndof <= 3 continue
    if ( nDof[i] <= 3 ) continue;

    const Mono::MonoTrack & tr = tracks[i];

    // Calculate expected Eta, Phi for ECAL cluster
    float ThisZ = tr.rzp0() + EcalR*tr.rzp1() + EcalR*EcalR*tr.rzp2();
    float ThisEta = TMath::ASinH(ThisZ / EcalR);
    float ThisPhi = tr.xyp1() - TMath::ASin((EcalR*EcalR-tr.xyp0()*tr.xyp2())/(2*EcalR*(tr.xyp2()-tr.xyp0())));

    float MinDR = 999;
    int BestEBCluster=-1, BestEECluster=-1;
    const unsigned nEBclusters = m_nCombEgamma;
    for(unsigned j=0; j < nEBclusters; j++){
      float DEta = ThisEta-m_egComb_eta[j];
      float DPhi = ThisPhi-m_egComb_phi[j];
      while(DPhi < -3.14159) DPhi += 2*3.14159;
      while(DPhi >  3.14159) DPhi -= 2*3.14159;

      float ThisDR = sqrt(pow(DEta,2) + pow(DPhi,2));
      if(ThisDR < MinDR){
        MinDR = ThisDR;
        BestEBCluster = j;
      }
    }

    const unsigned nEEclusters = m_nCombEE;
    for(unsigned j=0; j < nEEclusters; j++){
      float DEta = ThisEta-m_eeComb_eta[j];
      float DPhi = ThisPhi-m_eeComb_phi[j];
      while(DPhi < -3.14159) DPhi += 2*3.14159;
      while(DPhi >  3.14159) DPhi -= 2*3.14159;

      float ThisDR = sqrt(pow(DEta,2) + pow(DPhi,2));
      if(ThisDR < MinDR){
        MinDR = ThisDR;
        BestEECluster = j;
        BestEBCluster = -1;
      }
      //mai added this recenlty

      //MatchGenLeptons(ntuple_.get_lep2_eta(), ntuple_.get_lep2_phi(), mu) : make_pair(-999., -999.);   

      //      getPho175TrigCode(DEta, DPhi,ThisDR);
    }

    // fill place holders of matches
    int matchEB = BestEBCluster;
    int matchEE = BestEECluster;
    double distEB = BestEBCluster >= 0 ? MinDR : 999;
    double distEE = BestEECluster >= 0 ? MinDR : 999;

    // continue to next track if match = -1
    if ( matchEB == -1 && matchEE == -1 ) continue;

    // calculate dE/dX significance
    const double dEdXSig = sqrt(-TMath::Log(TMath::BinomialI(0.07, subHits[i], satSubHits[i])));
 
    /////////////////////////////////////////
    // assign values to Branches
    m_nCandidates++;
    m_candSubHits.push_back( subHits[i] );
    m_candSatSubHits.push_back( satSubHits[i] );
    m_canddEdXSig.push_back( dEdXSig );
    m_candTIso.push_back( tIso[i] );
    m_candXYPar0.push_back( tr.xyp0() );
    m_candXYPar1.push_back( tr.xyp1() );
    m_candXYPar2.push_back( tr.xyp2() );
    m_candRZPar0.push_back( tr.rzp0() );
    m_candRZPar1.push_back( tr.rzp1() );
    m_candRZPar2.push_back( tr.rzp2() );

    if ( distEB < distEE ) {
      m_candDist.push_back( distEB );
      
      m_candSeedFrac.push_back( m_egComb_frac51[matchEB] );
      m_candf15.push_back( m_egComb_frac15[matchEB] );
      m_candE55.push_back( m_egComb_e55[matchEB] );
      m_candE99.push_back( m_egComb_e99[matchEB] );
      m_candSwissCross.push_back( m_egComb_SwissCross[matchEB]);
      m_candHIso.push_back( m_egComb_hcalIso[matchEB] );
      m_candEta.push_back( m_egComb_eta[matchEB] );
      m_candPhi.push_back( m_egComb_phi[matchEB] );
      m_candPho175TrigCode.push_back ( getPho175TrigCode( m_egComb_eta[matchEB] , m_egComb_phi[matchEB] , *m_trigEventHandle ) );
      m_candPho200TrigCode.push_back ( getPho200TrigCode( m_egComb_eta[matchEB] , m_egComb_phi[matchEB] , *m_trigEventHandle ) );
    
      m_cande2x5Right.push_back( m_egComb_e25Right[matchEB] );
      m_cande2x5Left.push_back( m_egComb_e25Left[matchEB]);
      m_cande2x5Top.push_back( m_egComb_e25Top[matchEB] );
      m_cande2x5Bottom.push_back( m_egComb_e25Bottom[matchEB] );
      m_cande2x5Max.push_back( m_egComb_e25Max[matchEB] );

      m_candeRight.push_back( m_egComb_eRight[matchEB] );
      m_candeLeft.push_back( m_egComb_eLeft[matchEB]);
      m_candeTop.push_back( m_egComb_eTop[matchEB] );
      m_candeBottom.push_back( m_egComb_eBottom[matchEB] );
      m_candeMax.push_back( m_egComb_eMax[matchEB] );



    } else {
      m_candDist.push_back( distEE );

      m_candSeedFrac.push_back( m_eeComb_frac51[matchEE] );
      m_candf15.push_back( m_eeComb_frac15[matchEE] );
      m_candE55.push_back( m_eeComb_e55[matchEE] );
      m_candE99.push_back( m_eeComb_e99[matchEE] );
      m_candSwissCross.push_back( m_eeComb_SwissCross[matchEE]);
      m_candHIso.push_back( m_eeComb_hcalIso[matchEE] );
      m_candEta.push_back( m_eeComb_eta[matchEE] );
      m_candPhi.push_back( m_eeComb_phi[matchEE] );
      m_candPho175TrigCode.push_back (getPho175TrigCode(m_eeComb_eta[matchEE] , m_eeComb_phi[matchEE] , *m_trigEventHandle ) );
      m_candPho200TrigCode.push_back (getPho200TrigCode(m_eeComb_eta[matchEE] , m_eeComb_phi[matchEE] , *m_trigEventHandle ) );

      m_cande2x5Right.push_back( m_eeComb_e25Right[matchEE] );
      m_cande2x5Left.push_back( m_eeComb_e25Left[matchEE] );
      m_cande2x5Top.push_back( m_eeComb_e25Top[matchEE] );
      m_cande2x5Bottom.push_back( m_eeComb_e25Bottom[matchEE] );
      m_cande2x5Max.push_back( m_eeComb_e25Max[matchEE] );


      m_candeRight.push_back( m_eeComb_eRight[matchEB] );
      m_candeLeft.push_back( m_eeComb_eLeft[matchEE]);
      m_candeTop.push_back( m_eeComb_eTop[matchEE] );
      m_candeBottom.push_back( m_eeComb_eBottom[matchEE] );
      m_candeMax.push_back( m_eeComb_eMax[matchEE] );

    }

 }
}



// ------------ method called when ending the processing of a run  ------------
void MonoNtupleDumper::endRun(edm::Run const&, edm::EventSetup const&)
{
  //m_ecalCalib.calculateHij();
  //m_ecalCalib.dumpCalibration();

  std::map<Mono::ClustCategorizer,unsigned>::iterator counts = m_clustCatCount.begin();
  std::map<Mono::ClustCategorizer,unsigned>::iterator end = m_clustCatCount.end();
  for ( ; counts != end; counts++ ) {
    const Mono::ClustCategorizer & cat = counts->first;
    unsigned & count = counts->second;
    m_clustEMap[cat]->Scale(1.0/count);  
    m_clustTMap[cat]->Scale(1.0/count);
  }

}

// ------------ method called when starting to processes a luminosity block  ------------
void 
MonoNtupleDumper::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
MonoNtupleDumper::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MonoNtupleDumper::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MonoNtupleDumper);

bool passFilter(const float eta,const float phi,const trigger::TriggerEvent& trigEvent,const std::string& filter,const float maxDR2)
{
  edm::InputTag filterTag(filter,"",trigEvent.usedProcessName());
  trigger::size_type filterIndex = trigEvent.filterIndex(filterTag); 
  if(filterIndex<trigEvent.sizeFilters()){ //check that filter is in triggerEvent
    const trigger::Keys& trigKeys = trigEvent.filterKeys(filterIndex); 
    const trigger::TriggerObjectCollection & trigObjColl(trigEvent.getObjects());
    for(trigger::Keys::const_iterator keyIt=trigKeys.begin();keyIt!=trigKeys.end();++keyIt){ 
      const trigger::TriggerObject& obj = trigObjColl[*keyIt];
      if(reco::deltaR2(eta,phi,obj.eta(),obj.phi())<maxDR2) return true;
    }//end loop over keys
  }//end check that filter is valid and in trigEvent
  return false;
}

//
//fragment.HLTPhoton175Sequence = cms.Sequence( fragment.HLTDoFullUnpackingEgammaEcalSequence + fragment.HLTPFClusteringForEgamma + fragment.hltEgammaCandidates + fragment.hltEGL1SingleEGNonIsoOrWithJetAndTauFilter + fragment.hltEG175EtFilter + fragment.HLTDoLocalHcalSequence + fragment.HLTFastJetForEgamma + fragment.hltEgammaHoverE + fragment.hltEG175HEFilter )
//fragment.HLTPhoton200Sequence = cms.Sequence( fragment.HLTDoFullUnpackingEgammaEcalSequence + fragment.HLTPFClusteringForEgamma + fragment.hltEgammaCandidates + fragment.hltEGL1SingleEGNonIsoOrWithJetAndTauFilter + fragment.hltEG200EtFilter + fragment.HLTDoLocalHcalSequence + fragment.HLTFastJetForEgamma + fragment.hltEgammaHoverE + fragment.hltEG200HEFilter )
//=0 failed
//=1 passed L1
//=2 reco supercluster & matched to L1
//=3 passed Et
//=4 passed H/E (and thus the trigger)
int getPho175TrigCode(const float eta,const float phi,const trigger::TriggerEvent& trigEvent,const float maxDR2)
{
  //the filters of pho175 in order
  int retCode=0;
  const std::vector<std::string> hltFilters = {"hltL1sSingleEGNonIsoOrWithJetAndTau",
                 "hltEGL1SingleEGNonIsoOrWithJetAndTauFilter",
                 "hltEG175EtFilter",
                 "hltEG175HEFilter"};
  for(size_t filterNr=0;filterNr<hltFilters.size();filterNr++){
    if(passFilter(eta,phi,trigEvent,hltFilters[filterNr],maxDR2)) retCode=filterNr+1;
  }
  return retCode;
}

int getPho200TrigCode(const float eta,const float phi,const trigger::TriggerEvent& trigEvent,const float maxDR2)
{
  //the filters of pho175 in order
  int retCode=0;
  const std::vector<std::string> hltFilters = {"hltL1sSingleEGNonIsoOrWithJetAndTau",
                                               "hltEGL1SingleEGNonIsoOrWithJetAndTauFilter",
                                               "hltEG200EtFilter",
                                               "hltEG200HEFilter"};
  for(size_t filterNr=0;filterNr<hltFilters.size();filterNr++){
    if(passFilter(eta,phi,trigEvent,hltFilters[filterNr],maxDR2)) retCode=filterNr+1;
  }
  return retCode;
}

bool evtPassesPho175L1(const trigger::TriggerEvent& trigEvent)
{
  //the filters of pho175 in order
  edm::InputTag filterTag("hltL1sSingleEGNonIsoOrWithJetAndTau","",trigEvent.usedProcessName());
  trigger::size_type filterIndex = trigEvent.filterIndex(filterTag); 
  if(filterIndex<trigEvent.sizeFilters()){
    return trigEvent.filterKeys(filterIndex).size()>=1;//at least 1 object to pass
  } else return false;
}

bool evtPassesPho200L1(const trigger::TriggerEvent& trigEvent)
{
  // Filters for pho200 
  edm::InputTag filterTag("hltL1sSingleEGNonIsoOrWithJetAndTau","",trigEvent.usedProcessName());
  trigger::size_type filterIndex = trigEvent.filterIndex(filterTag); 
  if(filterIndex<trigEvent.sizeFilters()){
    return trigEvent.filterKeys(filterIndex).size()>=1;//at least 1 object to pass
  } else return false;                    
}

bool evtPassesPFMET250L1(const trigger::TriggerEvent& trigEvent)
{
  edm::InputTag filterTag("hltL1sAllETMHFSeeds","",trigEvent.usedProcessName());
  //Filter for 2017: hltL1sAllETMHadSeeds
  trigger::size_type filterIndex = trigEvent.filterIndex(filterTag); 
  if(filterIndex<trigEvent.sizeFilters()){
    return trigEvent.filterKeys(filterIndex).size()>=1;//at least 1 object to pass
  } else return false;                    
}


bool evtPassesPFMET300L1(const trigger::TriggerEvent& trigEvent)
{
  edm::InputTag filterTag("hltL1sETM60IorETM70IorETM80IorETM90IorETM100IorETM120","",trigEvent.usedProcessName());  
  //hltL1sAllETMHadSeeds
  trigger::size_type filterIndex = trigEvent.filterIndex(filterTag); 
  if(filterIndex<trigEvent.sizeFilters()){
    return trigEvent.filterKeys(filterIndex).size()>=1;//at least 1 object to pass
  } else return false;                    
}

