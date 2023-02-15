import FWCore.ParameterSet.Config as cms

process = cms.Process("Mpl")

### standard MessageLoggerConfiguration
process.load("FWCore.MessageService.MessageLogger_cfi")

### Standard Configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

## Fitter-smoother: loosen outlier rejection as for first data-taking with LHC "collisions"
process.KFFittingSmootherWithOutliersRejectionAndRK.BreakTrajWith2ConsecutiveMissing = False
process.KFFittingSmootherWithOutliersRejectionAndRK.EstimateCut = 1000

### Conditions
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v17', '')

### Track refitter specific stuff
#from RecoTracker.TrackProducer.TrackRefitters_cff import *
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi")

### unclean EE
process.uncleanEERecovered = cms.EDProducer(
    'UncleanSCRecoveryProducer',
    # input collections
    cleanBcCollection = cms.InputTag('multi5x5SuperClusters','multi5x5EndcapBasicClusters'),
    cleanScCollection = cms.InputTag('multi5x5SuperClusters','multi5x5EndcapSuperClusters'),
    
    uncleanBcCollection = cms.InputTag('multi5x5SuperClusters','uncleanOnlyMulti5x5EndcapBasicClusters'),
    uncleanScCollection = cms.InputTag('multi5x5SuperClusters','uncleanOnlyMulti5x5EndcapSuperClusters'),
    # names of collections to be produced:
    bcCollection = cms.string('uncleanEndcapBasicClusters'),
    scCollection = cms.string('uncleanEndcapSuperClusters'),
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
 	'root://xrootd-cms.infn.it///store/relval/CMSSW_10_6_1_patch1/SinglePhoton/USER/EXOMONOPOLE-106X_dataRun2_v17_RelVal_2017B-v1/20000/EBA96788-2C1E-794F-8613-A75B55E1E4F3.root',
        'root://xrootd-cms.infn.it///store/relval/CMSSW_10_6_1_patch1/SinglePhoton/USER/EXOMONOPOLE-106X_dataRun2_v17_RelVal_2017B-v1/20000/C8185FD8-6847-F940-8C9C-F0769982DBC2.root',
        'root://xrootd-cms.infn.it///store/relval/CMSSW_10_6_1_patch1/SinglePhoton/USER/EXOMONOPOLE-106X_dataRun2_v17_RelVal_2017B-v1/20000/C5C1AB63-183D-A54C-A07C-588DB2B86BF6.root'
    ),
    duplicateCheckMode = cms.untracked.string('checkEachRealDataFile') 
)

### Construct combined (clean and uncleanOnly Ecal clusters)
process.load("RecoEcal.EgammaClusterProducers.uncleanSCRecovery_cfi")

import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes

process.Monopoler = cms.EDAnalyzer(
    'MonoNtupleDumper'
    ,isData = cms.bool(True)
    ,Output = cms.string("MonoNtuple_Data.root")
    ,TriggerResults = cms.InputTag("TriggerResults","","HLT")
    ,TriggerEvent = cms.InputTag("hltTriggerSummaryAOD","","HLT")
    ,GeneratorTag = cms.InputTag("generatorSmeared","")
    ,PrimaryVertices = cms.InputTag("offlinePrimaryVertices","")                                 
    ,EcalEBRecHits = cms.InputTag("ecalRecHit","EcalRecHitsEB") 
    ,EcalEERecHits = cms.InputTag("ecalRecHit","EcalRecHitsEE") 
    ,HBHERecHits = cms.InputTag("hbhereco","")
    ,JetTag = cms.InputTag("ak4PFJets","")
    ,ElectronTag = cms.InputTag("gedGsfElectrons","")
    ,PhotonTag = cms.InputTag("photons","")
    ,METTag = cms.InputTag("pfMet","")
    ,bcClusterTag = cms.InputTag("hybridSuperClusters","uncleanOnlyHybridBarrelBasicClusters") 
    ,ccClusterTag = cms.InputTag("hybridSuperClusters","hybridBarrelBasicClusters")
    ,combClusterTag = cms.InputTag("uncleanSCRecovered","uncleanHybridBarrelBasicClusters") 
    ,eeCleanTag = cms.InputTag("multi5x5SuperClusters","multi5x5EndcapBasicClusters")
    ,eeUncleanTag = cms.InputTag("multi5x5SuperClusters","uncleanOnlyMulti5x5EndcapBasicClusters") 
    ,eeCombTag = cms.InputTag("uncleanEERecovered","uncleanEndcapBasicClusters")
    ,StripSeedLength = cms.uint32(3)
    ,ClusterLength = cms.uint32(5)
    ,SeedThreshold = cms.double(50.)
    ,TrackTag=cms.InputTag("TrackRefitter")                                 
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
#process.refit_step = cms.Path(process.TrackRefitter)
process.refit_step = cms.Path(process.MeasurementTrackerEvent * process.TrackRefitter)
process.mpl_step = cms.Path(process.Monopoler)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.suppressWarning.append("Monopoler")

process.p1 = cms.Schedule(
    process.ecalCombine_step
    ,process.ecalCombineEE_step
    ,process.refit_step
    ,process.mpl_step
)
#process.outpath = cms.EndPath(process.TRACKS)

