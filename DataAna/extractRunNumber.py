import argparse
import ROOT
from array import array
import re 


def extractRunNumber(inlist):
	with open(inlist) as file:
		for f_name in file:
			#print(f_name.rstrip())
			f=ROOT.TFile(f_name.rstrip())
			tree=f.Get("pixelTree")
			runHisto=ROOT.TH1F("runHisto","", 50,0,5000000 )
			tree.Draw("run>>runHisto")
			mean = runHisto.GetMean()
			runNumber = round(mean)
			print("inputFileName="+(f_name.rstrip())+" runNumber=" +str(runNumber))
			
		
if __name__ == '__main__':
    
    ROOT.gROOT.SetBatch(True)
 #-------------------------------------------------------------------------
 #  PARSE COMMAND LINE ARGUMENT
 #------------------------------------------------------------------------- 
    parser = argparse.ArgumentParser(description='Charge MPV extractor')
    parser.add_argument('-l', action="store", dest = 'inlist', default = 'HSCP_file.list')
    
    p = parser.parse_args()
    extractRunNumber(p.inlist)