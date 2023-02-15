import FWCore.ParameterSet.Config as cms
 
process = cms.Process("Mpl")

### standard MessageLoggerConfiguration
process.load("FWCore.MessageService.MessageLogger_cfi")

### Standard Configurations
process.load("Configuration.StandardSequences.Services_cff")
process.load('Configuration/StandardSequences/GeometryIdeal_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
 

## Fitter-smoother: loosen outlier rejection as for first data-taking with LHC "collisions"
process.KFFittingSmootherWithOutliersRejectionAndRK.BreakTrajWith2ConsecutiveMissing = False
process.KFFittingSmootherWithOutliersRejectionAndRK.EstimateCut = 1000



### Conditions
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = 'GR09_P_V6::All'
process.GlobalTag.globaltag = 'MC_53_V8::All'

### Track refitter specific stuff
process.load("RecoTracker.TrackProducer.TrackRefitters_cff") #the correct one
#process.load("RecoTracker.TrackProducer.RefitterWithMaterial_cff") #the one for backward compatibility


process.maxEvents = cms.untracked.PSet(
     input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
### tracks from collisions                            
#    fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/user/c/chiochia/09_beam_commissioning/BSCskim_123151_Express.root') 
### reco mpl mc
    fileNames = cms.untracked.vstring('file:combiner.1000.root')
) 

#process.TRACKS = cms.OutputModule("PoolOutputModule",
#                                outputCommands = cms.untracked.vstring('drop *_*_*_*', 
#                                                                       'keep recoTracks_*_*_*',
#                                                                       'keep recoTrackExtras_*_*_*',
#                                                                       'keep TrackingRecHitsOwned_*_*_*'),
#
#                                fileName = cms.untracked.string('refitting.root')
#                                )

process.Combiner = cms.EDAnalyzer("TrackCombinerReco",
				Source=cms.string("generalTracks"),
				Output=cms.string("CombinedTracks.reco.1000.root"),
				Chi2Cut=cms.untracked.double(7.5),
                                PtCut=cms.untracked.double(3.0),
                                DeDxCut=cms.untracked.double(0),
				DefaultError=cms.untracked.double(0.05),
				ErrorFudge=cms.untracked.double(0.02)
)

process.refit_step = cms.Path(process.TrackRefitter)
process.combine_step = cms.Path(process.Combiner)

process.options = cms.untracked.PSet(     wantSummary = cms.untracked.bool(True) )


process.p1 = cms.Schedule(process.refit_step,
                      process.combine_step
                      #process.TrackRefitterP5
                      #process.TrackRefitterBHM
)
#process.outpath = cms.EndPath(process.TRACKS)
