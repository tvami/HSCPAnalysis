# Usage: 
# voms-proxy-init -rfc -voms cms -valid 192:00
# python 3getListsDatesDir.py -l crab_ALCARECO_RunsDir.txt -l2 crab_ALCARECO_DatesDir.txt 

import argparse
import ROOT
import os
from tqdm import tqdm

def extractRunNumber(inlist, inlist2):
	with open(inlist) as file1, open (inlist2) as file2:
		for f_name1, f_name2 in tqdm(zip(file1,file2)):
			rline = f_name1.rstrip().split()
                        rline2 = f_name2.rstrip().split()
			dir1 = rline[0]
                        dir2 = rline2[0]
			lsCommand = "dpns-ls /dpm/kfki.hu/home/cms/phedex/store/user/tvami/PixelTreesRuns/SingleMuon/" + dir1  + "/"+ dir2 + "/0000/" + " >> 3crab_ALCARECO_FilesDir.txt"
                        os.system(lsCommand)
		
if __name__ == '__main__':
    
    ROOT.gROOT.SetBatch(True)

    parser = argparse.ArgumentParser(description='Run number extractor')
    parser.add_argument('-l', action="store", dest = 'inlist', default = 'crab_ALCARECO_RunsDir.txt')
    parser.add_argument('-l2', action="store", dest = 'inlist2', default = 'crab_ALCARECO_DatesDir.txt')
    
    p = parser.parse_args()
    extractRunNumber(p.inlist, p.inlist2)

