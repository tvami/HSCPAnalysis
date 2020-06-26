import os                                                                                                              
import glob                                                                                                            

from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.transferLogs = True
config.General.requestName = 'HSCP_Gluino_Mass1800_GENSIM_2018_13TeV_v2' #can be anything

config.section_('JobType')
config.JobType.allowUndistributedCMSSW = True
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = '1HSCP_Gluino_Mass1800_GENSIM_2018.py'
config.JobType.inputFiles = []
config.JobType.outputFiles = ['HSCP_Gluino_Mass1800_GENSIM.root']
config.JobType.disableAutomaticOutputCollection = True
#config.JobType.priority = -1
config.JobType.maxMemoryMB = 6000
config.JobType.numCores = 8
#config.JobType.eventsPerLumi = 125

config.section_('Data')
config.Data.outLFNDirBase = '/store/user/tvami/HSCP/'
config.Data.outputPrimaryDataset = 'HSCPgluino_M_1800'
config.Data.outputDatasetTag = config.General.requestName
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 1000
config.Data.totalUnits = 1000000
config.Data.ignoreLocality = True
config.Data.publication = True

config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'
config.Site.whitelist   = ['T3_US_Baylor', 'T2_US_Caltech', 'T3_US_Colorado', 'T3_US_FIT', 'T1_US_FNAL', 'T3_US_FNALLPC', 'T2_US_Florida', 'T3_US_JHU', 'T3_US_Kansas', 'T2_US_MIT', 'T3_US_NotreDame', 'T2_US_Nebraska', 'T3_US_OSU', 'T2_US_Purdue', 'T3_US_Rice', 'T3_US_Rutgers', 'T3_US_MIT', 'T3_US_NERSC', 'T3_US_PSC', 'T3_US_SDSC', 'T3_US_OSG', 'T3_US_TACC', 'T3_US_TAMU', 'T3_US_TTU', 'T3_US_UCR', 'T2_US_UCSD', 'T3_US_UMD', 'T3_US_UMiss', 'T3_US_PuertoRico', 'T3_US_VC3_NotreDame', 'T2_US_Vanderbilt', 'T2_US_Wisconsin'] # All US sites

