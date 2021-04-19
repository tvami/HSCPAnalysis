# Usage
# python listAllRootFiles.py > inputForHistoryPlots.txt



import os
from tqdm import tqdm

#directory_to_check = "/Users/tav/Documents/1Research/projects/PixelOffline/TrackerPublicPlots/AddCMS-DP-2013_014/tracker-plots/content/CMS-DP-2013_014/"
#directory_to_check = "/eos/cms/store/caf/user/tvami/HSCP/PixelTrees/Histos_2018BC_April11/"
directory_to_check = "/eos/cms/store/caf/user/tvami/HSCP/PixelTrees/HistosOn2021-04-10_2017Files/"

# Get all the subdirectories of directory_to_check recursively and store them in a list:
directories = [os.path.abspath(x[0]) for x in os.walk(directory_to_check)]
directories.remove(os.path.abspath(directory_to_check)) # If you don't want your main directory included

for i in tqdm(directories):
      os.chdir(i)         # Change working Directory
      RootFileHere = False
      for fname in os.listdir('.'):
             if fname.endswith('.root') : RootFileHere = True
      if (RootFileHere) : 
             RunInd = i.find("_Run")
             rn = i[RunInd+4:RunInd+10]
             print(str(rn)+" " + str(i) + "/PixelTree_BasedHistos.root")