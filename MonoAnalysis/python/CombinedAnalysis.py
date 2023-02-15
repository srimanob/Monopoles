import FWCore.ParameterSet.Config as cms
 
process = cms.Process("Mpl")

### standard MessageLoggerConfiguration
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'DEBUG'
process.MessageLogger.categories.append('EcalCleaning')
process.MessageLogger.cerr.DEBUG = cms.untracked.PSet( limit=cms.untracked.int32(-1) )

### Standard Configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
#process.load('Configuration.StandardSequences.GeometryDB_cff')


## Fitter-smoother: loosen outlier rejection as for first data-taking with LHC "collisions"
process.KFFittingSmootherWithOutliersRejectionAndRK.BreakTrajWith2ConsecutiveMissing = False
process.KFFittingSmootherWithOutliersRejectionAndRK.EstimateCut = 1000



### Conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'MC_53_V8::All'

### Track refitter specific stuff
process.load("RecoTracker.TrackProducer.TrackRefitters_cff") #the correct one


process.maxEvents = cms.untracked.PSet(
     input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#	'file:/data1/swilbur/Pythia_534_monopole_1000GeV_RECO_RAW.root'
	'file:Pythia_536p1_monopole_1700GeV_DY_lhe_QUICK_RAW_RECO.root'
    )
    #,skipEvents = cms.untracked.uint32(21)
)

### Construct combined (clean and uncleanOnly Ecal clusters)
process.load("RecoEcal.EgammaClusterProducers.uncleanSCRecovery_cfi")
process.uncleanEERecovered = cms.EDProducer('UncleanSCRecoveryProducer',

            # input collections:
            cleanBcCollection = cms.InputTag('multi5x5SuperClusters','multi5x5EndcapBasicClusters'),
            cleanScCollection = cms.InputTag('multi5x5SuperClusters','multi5x5EndcapSuperClusters'),
                                    
            uncleanBcCollection = cms.InputTag('multi5x5SuperClusters','uncleanOnlyMulti5x5EndcapBasicClusters'),
            uncleanScCollection = cms.InputTag('multi5x5SuperClusters','uncleanOnlyMulti5x5EndcapSuperClusters'),
            # names of collections to be produced:
            bcCollection = cms.string('uncleanEndcapBasicClusters'),
            scCollection = cms.string('uncleanEndcapSuperClusters'),

            )

import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process.Monopoler = cms.EDAnalyzer('MonoNtupleDumper'
  ,isData = cms.bool(False)
  ,Output = cms.string("Combined.test.root")
  ,EcalEBRecHits = cms.InputTag("ecalRecHit","EcalRecHitsEB") 
  ,EcalEERecHits = cms.InputTag("ecalRecHit","EcalRecHitsEE") 
  ,HBHERecHits = cms.InputTag("hbhereco","")
  ,JetTag = cms.InputTag("ak5PFJets","")
  ,ElectronTag = cms.InputTag("gsfElectrons","")
  ,PhotonTag = cms.InputTag("photons","")
  ,METTag = cms.InputTag("pfMet","")
  ,StripSeedLength = cms.uint32(3)
  ,ClusterLength = cms.uint32(5)
  ,SeedThreshold = cms.double(50.)
  ,TrackSource=cms.string("TrackRefitter")
  ,TrackChi2Cut=cms.untracked.double(7.5)
  ,TrackPtCut=cms.untracked.double(3.0)
  ,TrackDeDxCut=cms.untracked.double(0)
  ,TrackDefaultError=cms.untracked.double(0.05)
  ,TrackErrorFudge=cms.untracked.double(0.02)
  ,TrackHitOutput=cms.untracked.bool(True)
)


process.ecalCombine_step = cms.Path(process.uncleanSCRecovered)
process.ecalCombineEE_step = cms.Path(process.uncleanEERecovered)
process.refit_step = cms.Path(process.TrackRefitter)
process.mpl_step = cms.Path(process.Monopoler)

process.options = cms.untracked.PSet(     wantSummary = cms.untracked.bool(True) )


process.p1 = cms.Schedule(process.ecalCombine_step,process.ecalCombineEE_step,process.refit_step,
			process.mpl_step
)
#process.outpath = cms.EndPath(process.TRACKS)
