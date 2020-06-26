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
    debugModules = cms.untracked.vstring('PCP'),
    destinations = cms.untracked.vstring('cout'),
                                    #    destinations = cms.untracked.vstring("log","cout"),
    cout = cms.untracked.PSet(
        threshold = cms.untracked.string('ERROR')
    )
#    log = cms.untracked.PSet(
#        threshold = cms.untracked.string('DEBUG')
#    )
                                    )
process.source = cms.Source(
    "EmptySource",
    #firstRun = cms.untracked.uint32(320064),
    #firstRun = cms.untracked.uint32(320804),
    firstRun = cms.untracked.uint32(320804),    
    )


from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '110X_upgrade2018_realistic_v9', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v27', '')
# Ultra Legacy
#process.GlobalTag = GlobalTag(process.GlobalTag, '110X_dataRun2_v12', '')

#process.GlobalTag = GlobalTag(process.GlobalTag, '110X_mcRun3_2021_realistic_v6', '')
#process.GlobalTag.DumpStat = cms.untracked.bool(True)

# read rechits
process.analysis = cms.EDAnalyzer(
    "LorentzAngleAnalyser",
    #ListFile = cms.untracked.string("./320804.list"),
    ListFile = cms.untracked.string("./HSCP_MC_2018.list"),    
    #RootOutfile = cms.untracked.string("./320804.root"),
    RootOutfile = cms.untracked.string("./HSCP_MC_2018.root"),    
    Verbosity = cms.untracked.int32(0),
    year = cms.untracked.int32(2018),
    NBinsDrift = cms.untracked.int32(285),
    MinDrift = cms.untracked.double(-1000), 
    MaxDrift = cms.untracked.double(1000),
    NBinsDepth = cms.untracked.int32(50),
    MinDepth = cms.untracked.double(0), 
    MaxDepth = cms.untracked.double(285.),
    MinFit = cms.untracked.double(5), 
    MaxFit = cms.untracked.double(280.),    
    #MinFit = cms.untracked.double(50), 
    #MaxFit = cms.untracked.double(250.),    
    InstLumiMin = cms.untracked.double(0.), 
    InstLumiMax  = cms.untracked.double(100.),    
    InstLumiFile  = cms.untracked.string("./run_ls_instlumi_pileup_2018.txt"),    
    PtCut = cms.untracked.double(3.),
    ClusterSizeYCut = cms.untracked.double(0.1),
    #ClusterSizeUpperYCut = cms.untracked.double(11.5),
    ClusterSizeUpperYCut = cms.untracked.double(10000),    
    #ClusterSizeXCut = cms.untracked.double(5),
    ClusterSizeXCut = cms.untracked.double(10000),
    ResidualCut = cms.untracked.double(0.1), 
    NormChi2Cut = cms.untracked.double(100), 
    #ClusterChargeCut = cms.untracked.double(50),
    ClusterChargeCut = cms.untracked.double(1000),
    ProbabilityCut= cms.untracked.double(0.01),    
    LargePix = cms.untracked.int32(1),   
    #AnalysisTypeString = cms.untracked.string("WithBowing")
    ApplyBowingCorrection = cms.untracked.int32(1)
   )

process.p = cms.Path(process.analysis)






