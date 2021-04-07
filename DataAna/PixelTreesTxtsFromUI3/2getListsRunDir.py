# Usage: 
# voms-proxy-init -rfc -voms cms -valid 192:00
# python 2getListsRunDir.py -l crab_ALCARECO_RunsDir.txt 
# where crab_ALCARECO_RunsDir.txt is coming from 
# dpns-ls /dpm/kfki.hu/home/cms/phedex/store/user/tvami/PixelTreesRuns/SingleMuon/ 

import argparse
import ROOT
import os
from tqdm import tqdm

def extractRunNumber(inlist):
	with open(inlist) as file:
		for f_name in tqdm(file):
			rline = f_name.rstrip().split()
			dir1 = rline[0]
			lsCommand = "dpns-ls /dpm/kfki.hu/home/cms/phedex/store/user/tvami/PixelTreesRuns/SingleMuon/" + dir1  +" >> 2crab_ALCARECO_DatesDir.txt"
			os.system(lsCommand)
		
if __name__ == '__main__':
    
    ROOT.gROOT.SetBatch(True)

    parser = argparse.ArgumentParser(description='Run number extractor')
    parser.add_argument('-l', action="store", dest = 'inlist', default = 'HSCP_file.list')
    
    p = parser.parse_args()
    extractRunNumber(p.inlist)

