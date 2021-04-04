# Usage: 
# voms-proxy-init -rfc -voms cms -valid 192:00
# python3 extractRunNumber.py -l HSCP_file2.list

import argparse
import ROOT

def extractRunNumber(inlist):
	with open(inlist) as file:
		for f_name in file:
			#print(f_name.rstrip())
			#filename = "root://cms-xrd-global.cern.ch/"+f_name.rstrip()
			#f=ROOT.TFile.Open("root://cms-xrd-global.cern.ch//store/user/tvami/PixelTrees/SingleMuon/crab_ALCARECO_2018A_PixelTrees_v1/200604_040318/0000/PixelTree_99.root")
			f=ROOT.TFile.Open(f_name.rstrip())
			if f.IsZombie():
				continue
			tree=f.Get("pixelTree")
			runHisto=ROOT.TH1F("runHisto","", 3300000,0,3300000 )
			tree.Draw("run>>runHisto")
			FirstNonZeroBin = runHisto.FindFirstBinAbove()
			runNumber = int(round(FirstNonZeroBin-1))
			print("inputFileName="+(f_name.rstrip())+" runNumber=" +str(runNumber))
			
		
if __name__ == '__main__':
    
    ROOT.gROOT.SetBatch(True)

    parser = argparse.ArgumentParser(description='Run number extractor')
    parser.add_argument('-l', action="store", dest = 'inlist', default = 'HSCP_file.list')
    
    p = parser.parse_args()
    extractRunNumber(p.inlist)

