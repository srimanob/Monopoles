import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )



#from EXODisplacedPhotonSkim2012D import *
#process.source = source
process.source = cms.Source("PoolSource",
    #skipEvents = cms.untracked.uint32(2003),
    # replace 'myfile.root' with the source file you want to use
    fileNames = cms.untracked.vstring(
     'file:Pythia_536p1_monopole_1000GeV_DY_lhe_RAW_RECO.root'
     #['file:../CMSSW_5_2_6_patch1/data/Pythia_526_monopole_1000GeV_RECO_RAW_1.root',
     #'file:../CMSSW_5_2_6_patch1/data/Pythia_526_monopole_1000GeV_RECO_RAW_2.root']
    ),
    duplicateCheckMode = cms.untracked.string('checkEachRealDataFile') 
)

import FWCore.PythonUtilities.LumiList as LumiList
import FWCore.ParameterSet.Types as CfgTypes
# YOU DON"T WANT the FOLLOWING LINES WITH CRAB USAGE!!!!!!
#process.souce.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#myLumis = LumiList.LumiList(filename='jsonFile.json').getCMSSWString().split(',')
#process.source.lumisToProcess.extend(myLumis)

process.TFileService = cms.Service("TFileService"
  #, fileName = cms.string('singlePhotonAnalysis_EXODisp_2012D.root')
  , fileName = cms.string('testDump.root')
  , closeFileFast = cms.untracked.bool( True )
)



process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Other statements
#process.GlobalTag.globaltag = 'MC_44_V7::All'
process.GlobalTag.globaltag = 'MC_52_V11::All'



process.demo = cms.EDAnalyzer('MonoNtupleDumper'
  ,isData = cms.bool(False)
  ,EcalEBRecHits = cms.InputTag("ecalRecHit","EcalRecHitsEB") 
  ,ElectronTag = cms.InputTag("gsfElectrons","")
  ,StripSeedLength = cms.uint32(3)
  ,ClusterLength = cms.uint32(5)
  ,SeedThreshold = cms.double(50.)
)


process.p = cms.Path(process.demo)
