#
# rechits are not persisent anymore, so one should run one of the CPEs
# on clusters ot do the track fitting. 11/08 d.k.
#
import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

process = cms.Process('HSCPAnalyzer',Run2_2018)

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('CalibTracker.SiPixelESProducers.SiPixelTemplateDBObjectESProducer_cfi')
process.load('CalibTracker.SiPixelESProducers.SiPixel2DTemplateDBObjectESProducer_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

process.MessageLogger = cms.Service("MessageLogger",
    debugModules = cms.untracked.vstring('HSCPAnalyzer'),
    destinations = cms.untracked.vstring('cout'),
#    destinations = cms.untracked.vstring("log","cout"),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('ERROR')
    )
#    log = cms.untracked.PSet(
#        threshold = cms.untracked.string('DEBUG')
#    )
)
process.source = cms.Source("EmptySource",
    firstRun = cms.untracked.uint32(1),
    #firstRun = cms.untracked.uint32(317664),
)


from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v27', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '110X_mcRun3_2021_realistic_v6', '')
#process.GlobalTag.DumpStat = cms.untracked.bool(True)


# read rechits
process.analysis = cms.EDAnalyzer("HSCPStudy",
    Verbosity = cms.untracked.int32(0),
    rootFileName = cms.untracked.string("file:/afs/cern.ch/work/t/tvami/public/HSCP/HSCP/CMSSW_10_6_1_patch1/src/0RECOChain/HSCPSingalNoPU/4HSCP_Gluino_Mass800BasedpixelTreeWithSimInfo.root")
)

process.p = cms.Path(process.analysis)

# test the DB object, works
#process.load("RecoLocalTracker.SiPixelRecHits.PixelCPEESProducers_cff")
##process.load("RecoLocalTracker.SiPixelRecHits.SiPixelRecHits_cfi")
##process.load("CalibTracker.SiPixelESProducers.SiPixelFakeTemplateDBObjectESSource_cfi"
##process.load("CalibTracker.SiPixelESProducers.SiPixelFakeCPEGenericErrorParmESSource_cfi"
#process.test = cms.EDAnalyzer("CPEAccessTester",
##    PixelCPE = cms.string('PixelCPEGeneric'),
#    PixelCPE = cms.string('PixelCPETemplateReco'),
#)
#process.p = cms.Path(process.test)





