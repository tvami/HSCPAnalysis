import os
import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run2_2017_cff import Run2_2017


process = cms.Process("PIXEL",Run2_2017)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.threshold = 'INFO'
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.categories.append('HLTrigReport')
#process.MessageLogger.categories.append('L1GtTrigReport')
process.options = cms.untracked.PSet( 
    SkipEvent = cms.untracked.vstring('ProductNotFound'), wantSummary = cms.untracked.bool(True) 
)


process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v27', '')

# -- Conditions
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")

# -- Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring("root://cms-xrd-global.cern.ch//store/data/Run2017B/SingleMuon/ALCARECO/SiPixelCalSingleMuon-ForPixelALCARECO_UL2017-v1/70000/A8F0BEC6-C0A5-2040-BAC6-A63C8860B838.root")
    )
# -- number of events
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )

#---------------------- Refitter -----------------------
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cff")
process.load("RecoTracker.TrackProducer.TrackRefitters_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("RecoTracker.MeasurementDet.MeasurementTrackerEventProducer_cfi") 

process.MeasurementTrackerEvent.pixelClusterProducer = 'ALCARECOSiPixelCalSingleMuon'
process.MeasurementTrackerEvent.stripClusterProducer = 'ALCARECOSiPixelCalSingleMuon'
process.MeasurementTrackerEvent.inactivePixelDetectorLabels = cms.VInputTag()
process.MeasurementTrackerEvent.inactiveStripDetectorLabels = cms.VInputTag()

process.TrackRefitter.src = 'ALCARECOSiPixelCalSingleMuon'
process.TrackRefitter.TrajectoryInEvent = True


#-------------------------------------------------------


process.PixelTree = cms.EDAnalyzer(
    "PixelTree",
    verbose                      = cms.untracked.int32(0),
    rootFileName                 = cms.untracked.string("PixelTree.root"),
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
    trackCollectionLabel         = cms.untracked.InputTag('ALCARECOSiPixelCalSingleMuon'),
    #pixelClusterLabel            = cms.untracked.InputTag('siPixelClusters'),
    pixelRecHitLabel             = cms.untracked.InputTag('siPixelRecHits'),
    HLTProcessName               = cms.untracked.string('HLT'),
    L1GTReadoutRecordLabel       = cms.untracked.InputTag('gtDigis'),
    hltL1GtObjectMap             = cms.untracked.InputTag('hltL1GtObjectMap'),
    accessSimHitInfo             = cms.untracked.bool(False),
    )

# -- Path
process.TrackRefitter_step = cms.Path(
  process.offlineBeamSpot*
  process.MeasurementTrackerEvent*
  process.TrackRefitter
)

process.PixelTree_step = cms.Path(
  process.PixelTree
)

process.schedule = cms.Schedule(
  process.TrackRefitter_step,
  process.PixelTree_step
)

#process.options.numberOfThreads=cms.untracked.uint32(8)
#process.options.numberOfStreams=cms.untracked.uint32(0)
#process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)

