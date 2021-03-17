#!/bin/bash
echo   "UnixTime-JobStart: "$(date +%s)

arch=slc7_amd64_gcc820
rel=CMSSW_10_6_11_patch1
sandbox=$(ls *.tar.bz2)
arguments=${@:1}
filelist=$1
cwd=$(pwd)

echo -e "------------------- START --------------------"
printf "Start time: "; TZ=CET /bin/date
printf "Job is running on node: "; /bin/hostname
printf "Job running as user: "; /usr/bin/id
printf "Job is running in directory: "; /bin/pwd -P
echo
echo -e "---------------- Environments ----------------"

echo -e "\n[0] source /cvmfs/cms.cern.ch/cmsset_default.sh"
source /cvmfs/cms.cern.ch/cmsset_default.sh

echo -e "\n[1] export SCRAM_ARCH=$arch"
export SCRAM_ARCH=$arch

echo -e "\n[2] scramv1 project CMSSW $rel"
scramv1 project CMSSW $rel

echo -e "\n[3] tar -xf $sandbox -C $rel/src/"
tar -xf $sandbox -C $rel/src/

echo -e "\n[4] cd $rel/src/"
cd $rel/src/

echo -e "\n[5] cmsenv"
eval `scramv1 runtime -sh`

echo -e "\n[5] scram b -j"
scram b -j

echo -e "\n------------------ Analyzer ------------------"

echo -e "\n[6] cd DataAna"
cd DataAna


echo   "UnixTime-AnalyzerStart: "$(date +%s)

echo -e "\n[8] cmsRun 5runHSCPAnalyzer.py $arguments 2>&1"
cmsRun 5runHSCPAnalyzer.py $arguments 2>&1

echo -e "\n[11] ls -ltr"
ls -ltr

echo -e "$arguments" > arguments.temp
era=$(cut -d '/' -f10 arguments.temp)
datum=$(cut -d '/' -f11 arguments.temp)


(cut -d '.' -f3 arguments.temp) > filename.temp
filename=$(cut -d '/' -f11 filename.temp)

echo $era
echo $datum
echo $filename

echo -e "\n[12] ls -ltr"
ls -ltr

echo -e "Copy root file to CAF"
xrdcp -f "output.root" "root://eoscms.cern.ch//eos/cms/store/caf/user/tvami/HSCP/PixelTrees/${era}/${datum}/${filename}BasedHistos.root" 


echo   "UnixTime-AnalyzerEnd: "$(date +%s)

echo -e "\n------------------ Cleanup -------------------"

echo -e "\n[9] cd ../../.."
cd ../../..

echo -e "\n[10] rm -r $rel $sandbox"
rm -r $rel $sandbox

echo -e "\n[11] ls -ltr"
ls -ltr

# delete output file if too small
# This prevents late, parallel and failing condor jobs overwriting the otherwise good output
if [ -f output.root ]; then 
   if [ $(ls -l output.root | awk '{ print $5 }' ) -lt 1000 ]; then 
       echo -e "\n[12] rm output.root"
       rm output.root
   fi
fi

echo -e "\n"
echo -e "-------------------- END ---------------------\n"
echo   "UnixTime-JobEnd: "$(date +%s)