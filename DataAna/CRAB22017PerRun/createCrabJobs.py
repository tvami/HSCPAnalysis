# Usage: 
# voms-proxy-init -rfc -voms cms -valid 192:00
# python3 extractRunNumber.py -l HSCP_file2.list

import argparse
import ROOT
import os
from tqdm import tqdm

def extractRunNumber(inlist):
	with open(inlist) as file:
		os.system("eval `scramv1 runtime -csh`")
		os.system("voms-proxy-init -rfc -voms cms -valid 192:00")
		for f_name in tqdm(file):
			rline = f_name.rstrip().split()
			rn = rline[0]
			newConfigName = "4crab_ALCARECO_2018A_Run"+rn+"_PixelTrees_v1.py"
			copyCommand = "cp 4crab_ALCARECO_2018A_RunBela_PixelTrees_v1.py " +newConfigName
			os.system(copyCommand)
			sedCommand = "sed -i 's/Bela/"+rn+"/g'i " + newConfigName
			os.system(sedCommand) 
			crabCommand = "crab submit -c " + newConfigName
			os.system(crabCommand)
		
if __name__ == '__main__':
    
    ROOT.gROOT.SetBatch(True)

    parser = argparse.ArgumentParser(description='Run number extractor')
    parser.add_argument('-l', action="store", dest = 'inlist', default = 'HSCP_file.list')
    
    p = parser.parse_args()
    extractRunNumber(p.inlist)

