
# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --filein file:HSCP_Gluino_Mass1800_GENSIM.root --fileout file:HSCP_Gluino_Mass1800_RECO.root --eventcontent RECOSIM --datatier GEN-SIM-RECO --conditions 106X_mc2017_realistic_v7 --step DIGI,DIGI2RAW,RAW2DIGI,L1,RECO --nThreads 8 --geometry DB:Extended --era Run2_2017
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run2_2017_cff import Run2_2017

process = cms.Process('RECO',Run2_2017)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(80)
)

# Input source
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring('root://cmseos.fnal.gov//store/user/tvami/HSCP/HSCPgluino_M_1800/HSCP_Gluino_Mass1800_GENSIM/200414_222359/0000/HSCP_Gluino_Mass1800_GENSIM_119.root'),
    inputCommands = cms.untracked.vstring(
        'keep *', 
        'drop *_genParticles_*_*', 
        'drop *_genParticlesForJets_*_*', 
        'drop *_kt4GenJets_*_*', 
        'drop *_kt6GenJets_*_*', 
        'drop *_iterativeCone5GenJets_*_*', 
        'drop *_ak4GenJets_*_*', 
        'drop *_ak7GenJets_*_*', 
        'drop *_ak8GenJets_*_*', 
        'drop *_ak4GenJetsNoNu_*_*', 
        'drop *_ak8GenJetsNoNu_*_*', 
        'drop *_genCandidatesForMET_*_*', 
        'drop *_genParticlesForMETAllVisible_*_*', 
        'drop *_genMetCalo_*_*', 
        'drop *_genMetCaloAndNonPrompt_*_*', 
        'drop *_genMetTrue_*_*', 
        'drop *_genMetIC5GenJs_*_*'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('GEN-SIM-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:HSCP_Gluino_Mass1800_RECO.root'),
    outputCommands = cms.untracked.vstring(
        'keep *', 
        'drop *_genParticlesForJets_*_*', 
        'drop *_kt4GenJets_*_*', 
        'drop *_kt6GenJets_*_*', 
        'drop *_iterativeCone5GenJets_*_*', 
        'drop *_ak4GenJets_*_*', 
        'drop *_ak7GenJets_*_*', 
        'drop *_ak8GenJets_*_*', 
        'drop *_ak4GenJetsNoNu_*_*', 
        'drop *_ak8GenJetsNoNu_*_*', 
        'drop *_genCandidatesForMET_*_*', 
        'drop *_genParticlesForMETAllVisible_*_*', 
        'drop *_genMetCalo_*_*', 
        'drop *_genMetCaloAndNonPrompt_*_*', 
        'drop *_genMetTrue_*_*', 
        'drop *_genMetIC5GenJs_*_*',
        'drop *CSC*_*_*_*',
        'drop *CTPPS*_*_*_*',
        'drop *alo*_*_*_*',
        'drop *_*alo*_*_*',
        'drop *DT*_*_*_*',
        'drop *EB*_*_*_*',
        'drop *EE*_*_*_*',
        'drop *cal*_*_*_*',
        'drop *mix*_*_*_*',
        'drop *_*mix_*_*',
        'drop *_*cal_*_*', #is :FEDRawDataCollection_hcalRawDatauHTR__RECO.obj or FEDRawDataCollection_caloLayer1RawFed1356__RECO.obj or :ESDigiCollection_simEcalUnsuppressedDigis__RECO.obj there?
        'drop *_*ecal*_*_*',
        'drop *GEM*_*_*_*',
        'drop *HF*_*_*_*',
        'drop *Hcal*_*_*_*',
        'drop *RPC*_*_*_*',
        'drop *Muon*_*_*_*',
        'drop *_*_*Muon*_*',
        'drop *_*_*GEM*_*',
        'drop *_*gem*_*_*', #is FEDRawDataCollection_gemPacker__RECO there?
        'drop *Totem*_*_*_*',
        'drop *_*_*Totem*_*',
        'drop *_*_*PLT*_*',
        'drop *_*_*Timer*_*',
        'drop *_*_*BCM*_*',
        'drop *_*_*BHM*_*',
        'drop *_*_*BSC*_*',
        'drop *l1t*_*_*_*',
        'drop *HBHE*_*_*_*',
        'drop *_*_*Hcal*_*',
        'drop *_*_*hltIter*_*',
        'drop *CompositeCandidates*_*_*_*',
        'drop *astor*_*_*_*',
        'drop *_*astor*_*_*', # is FEDRawDataCollection_castorRawData__RECO.obj there?
        'drop *METs*_*_*_*',
        'drop *PFTaus*_*_*_*',
        'drop *PFJets*_*_*_*',
        'drop *TriggerFilterObject*_*_*_*',
        'drop *_*gctDigis*_*_*',
    ),
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mc2017_realistic_v7', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.digi2raw_step = cms.Path(process.DigiToRaw)
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.reconstruction_step = cms.Path(process.reconstruction_trackingOnly)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RECOSIMoutput_step = cms.EndPath(process.RECOSIMoutput)

#process.DigiToRawTask_trackerOnly = cms.Task(process.siPixelRawData, process.SiStripDigiToRaw, process.rawDataCollector)
#process.DigiToRaw_trackerOnly = cms.Sequence(process.DigiToRawTask_trackerOnly)
#process.digi2raw_step = cms.Path(process.DigiToRaw_trackerOnly)

#process.RawToDigiTask_trackerOnly = cms.Task(process.siPixelDigis,process.siStripDigis)
#process.RawToDigi_trackerOnly = cms.Sequence(process.RawToDigiTask_trackerOnly)
#process.raw2digi_step = cms.Path(process.RawToDigi_trackerOnly)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.digi2raw_step,process.raw2digi_step,process.L1simulation_step,process.reconstruction_step,process.endjob_step,process.RECOSIMoutput_step)
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(8)
process.options.numberOfStreams=cms.untracked.uint32(0)
process.options.numberOfConcurrentLuminosityBlocks=cms.untracked.uint32(1)


# Customisation from command line

#Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
