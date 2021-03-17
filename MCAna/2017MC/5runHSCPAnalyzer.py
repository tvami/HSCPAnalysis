#
# rechits are not persisent anymore, so one should run one of the CPEs
# on clusters ot do the track fitting. 11/08 d.k.
#
import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run2_2017_cff import Run2_2017

process = cms.Process('HSCPAnalyzer',Run2_2017)

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
    
import FWCore.ParameterSet.VarParsing as opts
opt = opts.VarParsing ('analysis')
opt.register('runNumber',           1,
	     opts.VarParsing.multiplicity.singleton, opts.VarParsing.varType.int,
	     'Using a predetermined run number')

opt.parseArguments()
runNumber = opt.runNumber
process.source = cms.Source("EmptySource",
    firstRun = cms.untracked.uint32(runNumber),
)

from Configuration.AlCa.GlobalTag import GlobalTag
if runNumber >1:
    process.GlobalTag = GlobalTag(process.GlobalTag, '106X_dataRun2_v27', '')
else:
    process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mc2017_realistic_v7', '')

#process.GlobalTag.DumpStat = cms.untracked.bool(True)

process.TFileService = cms.Service('TFileService',
                                 fileName = cms.string('HSCPAna_Run'+str(runNumber)+'.root'),
                                 )

if (runNumber >1):
    FileName = cms.untracked.string("file:/afs/cern.ch/work/t/tvami/public/MonitoringTrees/CMSSW_11_0_0_pre9/src/2018ALCARECOBasedPixelTree/pixelTree_Run"+str(runNumber)+".root")
else:
    #FileName = cms.untracked.string("file:4HSCP_2017_Gluino_Mass1800_RECOwIncThreshold_PixelTrees.root")
    #FileName = cms.untracked.string("file:/eos/cms/store/caf/user/tvami/HSCP/4HSCP_2017_Gluino_Mass1800_FullRECO_PixelTrees.root")
    FileName = cms.untracked.string("file:4HSCP_Gluino_Mass1800BasedpixelTreeWithSimInfo.root")
        
# read rechits
process.analysis = cms.EDAnalyzer("HSCPStudy",
    Verbosity = cms.untracked.int32(0),
    rootFileName = FileName
)

process.p = cms.Path(process.analysis)

