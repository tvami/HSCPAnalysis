import os
import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018


process = cms.Process("PIXEL",Run2_2018)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.categories.append('HLTrigReport')
#process.MessageLogger.categories.append('L1GtTrigReport')
process.options = cms.untracked.PSet( 
    #SkipEvent = cms.untracked.vstring('ProductNotFound'), wantSummary = cms.untracked.bool(True) 
)


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_upgrade2018_realistic_v11', '')

# -- Conditions
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

# -- Input files
process.source = cms.Source(
    "PoolSource",
    # replace with your files
    #lastRun = cms.untracked.uint32(64789),
    #timetype = cms.string('runnumber'),
    #firstRun = cms.untracked.uint32(64108),
    #interval = cms.uint32(1),
    fileNames = cms.untracked.vstring(
        #"file:DY_13TeV_M1000_Q6.root",
	#"file:3HSCP_Gluino_Mass1800_RAWtoRECO_TrackingOnly_SignalOnly.root",
        "file:HSCP_Gluino_Mass1800_RECO.root",
    )
    )
# -- number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

#---------------------- Refitter -----------------------
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.TrackRefitter.src = 'generalTracks'
process.TrackRefitter.TrajectoryInEvent = True
#process.TrackRefitter.NavigationSchool = ""
#process.MeasurementTrackerEvent.stripClusterProducer = "SiStripClustersFromRawFacility"

#-------------------------------------------------------


process.PixelTree = cms.EDAnalyzer(
    "PixelTree",
    verbose                      = cms.untracked.int32(0),
    rootFileName                 = cms.untracked.string("4HSCP_Gluino_Mass1800BasedpixelTreeWithSimInfo.root"),
    associateRecoTracks = cms.bool(False),
    associateStrip = cms.bool(False),
    associatePixel = cms.bool(True),
    #RecHitProducer = cms.string('siStripMatchedRecHits'),
    pixelSimLinkSrc = cms.InputTag("simSiPixelDigis"),
    stripSimLinkSrc = cms.InputTag("simSiStripDigis"),
    ROUList = cms.vstring(
        'TrackerHitsPixelBarrelLowTof', 
        'TrackerHitsPixelBarrelHighTof', 
        'TrackerHitsPixelEndcapLowTof', 
        'TrackerHitsPixelEndcapHighTof'),
    #type                         = cms.untracked.string(getDataset(process.source.fileNames[0])),
    globalTag                    = process.GlobalTag.globaltag,
    dumpAllEvents                = cms.untracked.int32(0),
    PrimaryVertexCollectionLabel = cms.untracked.InputTag('offlinePrimaryVertices'),
    muonCollectionLabel          = cms.untracked.InputTag('muons'),
    trajectoryInputLabel         = cms.untracked.InputTag('TrackRefitter::PIXEL'),
    trackCollectionLabel         = cms.untracked.InputTag('generalTracks'),
    #pixelClusterLabel            = cms.untracked.InputTag('siPixelClusters'),
    pixelRecHitLabel             = cms.untracked.InputTag('siPixelRecHits'),
    HLTProcessName               = cms.untracked.string('RECO'),
    L1GTReadoutRecordLabel       = cms.untracked.InputTag('gtDigis'),
    hltL1GtObjectMap             = cms.untracked.InputTag('hltL1GtObjectMap'),
    HLTResultsLabel              = cms.untracked.InputTag('TriggerResults::RECO'),
    accessSimHitInfo             = cms.untracked.bool(True),
    )

# -- Path
process.PixelTree_step = cms.Path(
    process.PixelTree
    )

process.Refitter_step = cms.Path(
  process.MeasurementTrackerEvent*
  process.TrackRefitter*
  process.siPixelRecHits
)

process.schedule = cms.Schedule(
  process.Refitter_step,
  process.PixelTree_step
)

process.options.numberOfThreads=cms.untracked.uint32(8)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

