# Usage: 
# voms-proxy-init -rfc -voms cms -valid 192:00
# python 2getListsRunDir.py -l crab_ALCARECO_RunsDir.txt 
# where crab_ALCARECO_RunsDir.txt is coming from 
# dpns-ls /dpm/kfki.hu/home/cms/phedex/store/user/tvami/PixelTreesRuns/SingleMuon/ 

import argparse
import ROOT
import os
from tqdm import tqdm

def extractRunNumber():
  os.system("dpns-ls /dpm/kfki.hu/home/cms/phedex/store/user/tvami/PixelTreesRuns/SingleMuon/ > 1crab_ALCARECO_RunsDir.txt")
#  with open("1crab_ALCARECO_RunsDir.txt", "r") as file:
  with open("1crab_ALCARECO_RunsDir2018Only.txt", "r") as file:
    for f_name in tqdm(file):
      rline = f_name.rstrip().split()
      dir1 = rline[0]
      lsCommandforDates = "dpns-ls /dpm/kfki.hu/home/cms/phedex/store/user/tvami/PixelTreesRuns/SingleMuon/" + dir1+" > 2crab_ALCARECO_DatesDir.txt"
      os.system(lsCommandforDates)
      with open("2crab_ALCARECO_DatesDir.txt", "r") as fileDatesDir:
        for f_nameFileDatesDir in (fileDatesDir):
          fileDatesDirLine = f_nameFileDatesDir.rstrip().split()
          dir2 = fileDatesDirLine[0]
          lsCommandForFilesDir = "dpns-ls /dpm/kfki.hu/home/cms/phedex/store/user/tvami/PixelTreesRuns/SingleMuon/" + dir1+ "/"+ dir2 + "/0000/" + " >> 3crab_ALCARECO_FilesDir.txt"
          os.system(lsCommandForFilesDir)
          with open("3crab_ALCARECO_FilesDir.txt", "r") as fileFilesDir:
	    for f_name in (fileFilesDir):
              rl = f_name.rstrip().split()
              dir3 = rl[0]
              print("root://cms-xrd-global.cern.ch//store/user/tvami/PixelTreesRuns/SingleMuon/" + dir1 + "/"+ dir2+"/0000/" + dir3)
          os.system("rm 3crab_ALCARECO_FilesDir.txt")
      os.system("rm 2crab_ALCARECO_DatesDir.txt")
  os.system("rm 1crab_ALCARECO_RunsDir.txt")
		
if __name__ == '__main__':

  ROOT.gROOT.SetBatch(True)
  extractRunNumber()

