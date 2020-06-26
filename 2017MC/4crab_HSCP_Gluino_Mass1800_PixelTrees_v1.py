import os                                                                                                              
import glob                                                                                                            

from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'PrivateHSCP_2017_Gluino_Mass1800_PixelTrees_v1'

config.section_('JobType')
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '4pixelTree-mc.py'
config.JobType.outputFiles = ['4HSCP_Gluino_Mass1800BasedpixelTreeWithSimInfo.root']
config.JobType.disableAutomaticOutputCollection = True
config.JobType.maxMemoryMB = 6000
config.JobType.numCores = 8

config.section_('Data')
config.Data.inputDBS = 'phys03'
config.Data.inputDataset = '/HSCPgluino_M_1800/tvami-crab_PrivateHSCP_2017_Gluino_Mass1800_DIGI2RECO_NoPU_v1-1297ab9d779455c655abd3d1857271b2/USER'
config.Data.outLFNDirBase = '/store/user/tvami/HSCP/'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 10
config.Data.publication = True

config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'

