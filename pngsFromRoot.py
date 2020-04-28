import ROOT, sys, os, time, re
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog fileName.root Type (HSCP/Norm)")
(opt,args) = parser.parse_args()

save = 1

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(1)

fileName = sys.argv[1]
Type = sys.argv[2]

# Only remove stuff for before matching
clean_dirname = {
    "Counts_vs_" : "",
    "Bins"       : "",
    "RazorLep"   : "Razor",
    "Razor1vl"   : "Razor",
    "Razorll"    : "Razor",
    "RazorNoPho" : "Razor",
    "JetAK8"     : "",
    }
replace_in_dirname = {
    #"METFine" : "dir:Fine_binning/MET",
    #"MRFine"  : "dir:Fine_binning/MRFine",
    #"R2Fine"  : "dir:Fine_binning/R2Fine",
    "ROCCurve_vs_"        : "",
    "EleFromW_"           : "dir:ROC_Curves/EleFromW",
    "MuFromW_"            : "dir:ROC_Curves/MuFromW",
    "EleFromTop_"         : "dir:ROC_Curves/EleFromTop",
    "EleNuFromTop_"       : "dir:ROC_Curves/EleNuFromTop",
    "MuFromTop_"          : "dir:ROC_Curves/MuFromTop",
    "MuNuFromTop_"        : "dir:ROC_Curves/MuNuFromTop",
    "EleFromHardProcess_" : "dir:ROC_Curves/EleFromHardProcess",
    "MuFromHardProcess_"  : "dir:ROC_Curves/MuFromHardProcess",
    "EleLepTop_"          : "dir:ROC_Curves/EleLepTop",
    "EleNuFromLepTop_"    : "dir:ROC_Curves/EleNuFromLepTop",
    "MuLepTop_"           : "dir:ROC_Curves/MuLepTop",
    "MuNuFromLepTop_"     : "dir:ROC_Curves/MuNuFromLepTop",
    "HadTop_"             : "dir:ROC_Curves/HadTop",
    "HadW_"               : "dir:ROC_Curves/HadW",
    "HadZ_"               : "dir:ROC_Curves/HadZ",
    "HadH_"               : "dir:ROC_Curves/HadH",
    }
# Modify plotname
replace_in_plotname = {
    # Boost SRs
    "StackPlot_" : "",
    "StackPlotSignal" : "",
    "StackPlotTopSignal" : "",
    "StackPlotVSignal" : "",
    "StackPlotHSignal" : "",
    "_JetHTMET" : "",
    "JetHTMET_" : "",
    "_BlindData" : "",
    "_Ratio" : "",
    "_0_" : "dir:SR_0Lepton",
    "_1_" : "dir:SR_1Lepton",
    # HLT
    "HadronicMeasurements" : "dir:HLT_Hadronic",
    "LeptonicMeasurements" : "dir:HLT_Leptonic",
    "_Baseline" : "",
    "_1Ele"     : "",
    "_1Muon"    : "",
    "Bins"      : "",
    "CR_Fake"             : "dir:Fake_rates",
    "CR_"                 : "dir:Control_regions",
    "SR_"                 : "dir:Signal_regions",
    "Val_"                : "dir:Validation_regions",
    }

# list of: plot directories, match strings
plot_these = [
    [
       #[".*Razor.*", ".*HT.*", ".*MET.*", ".*R2.*", ".*MR.*"],
       #[".*Razor.*", ".*R2.*", ".*MR.*"],
        [".*MRR2.*", ".*HT.*", ".*MET.*", ".*DeltaPhi.*", ".*NJet.*"],
        ["StackPlot_JetHTMET_Pre.*_Ratio",
            "StackPlot_JetHTMET_CR.*_Ratio",
            "StackPlot_JetHTMET_Val.*_Ratio",
            #"StackPlot.*Signal_BlindData_SR.*",
        ]
    ],
    # Preselection plots
    [
        ["JetPt.*", ".*MegaJet.*"],
        ["StackPlot_JetHTMET_Pre.*_Ratio",
            "StackPlot_JetHTMET_CR.*_Ratio",
            "StackPlot_JetHTMET_Val.*_Ratio",
            #"StackPlot.*Signal_BlindData_SR.*",
        ]
    ],
]
for plot_comb in plot_these:
    selected_dirs = plot_comb[0]
    selected_plots = plot_comb[1]
    
    # compile patterns
    dir_patterns = []
    for pattern in selected_dirs:
        dir_patterns.append(re.compile(pattern))
    plot_patterns = []
    for pattern in selected_plots:
        plot_patterns.append(re.compile(pattern))
    
    if (fileName!=None):
        print fileName
        input_file = fileName
        output_dir = "Plots/"

        if save and not os.path.exists(output_dir): os.mkdir(output_dir)
        
        f = ROOT.TFile.Open(input_file)
        dirs = []
        for i in range(0, f.GetListOfKeys().GetEntries()):
            # Remove/modify unnecessary stuff from the name of the plot that was required by SmartHistos to ditinguish plots
            dirname = f.GetListOfKeys().At(i).GetName()
            clean_dirname_to_match = dirname
            for part_to_replace in clean_dirname.keys():
                if part_to_replace in clean_dirname_to_match:
                    clean_dirname_to_match = clean_dirname_to_match.replace(part_to_replace, clean_dirname[part_to_replace])
            # Find the plot name to match
            match_dir = True #False
            for pattern in dir_patterns:
                if re.match(pattern, clean_dirname_to_match):
                    match_dir = True
                    break
            if match_dir:
                #print "Matched directory: "+dirname
                curr_dir = f.GetDirectory(dirname)
                for i in range(0, curr_dir.GetListOfKeys().GetEntries()):
                    # Match the plot of interest
                    keyname = curr_dir.GetListOfKeys().At(i).GetName()
                    match_plot = True #False
                    for pattern in plot_patterns:
                        if re.match(pattern, keyname):
                            match_plot = True
                            break
                    if match_plot:
                        # The plot should be TCanvas
                        obj = f.Get(dirname+"/"+keyname)
                        if (True):
                        #if obj.InheritsFrom("TCanvas"):
                            #can = obj
                            can = ROOT.TCanvas(dirname+"/"+keyname)
                            clean_keyname = keyname
                            subdirname = ""
                            for part_to_replace in replace_in_plotname.keys():
                                if part_to_replace in clean_keyname:
                                    if "dir:" in replace_in_plotname[part_to_replace]:
                                        subdirname += replace_in_plotname[part_to_replace].split(":")[1] + "/"
                                        clean_keyname = clean_keyname.replace(part_to_replace, "")
                                    else:
                                        clean_keyname = clean_keyname.replace(part_to_replace, replace_in_plotname[part_to_replace])
                            plotname = clean_dirname_to_match
                            for part_to_replace in replace_in_dirname.keys():
                                if part_to_replace in plotname:
                                    if "dir:" in replace_in_dirname[part_to_replace]:
                                        subdirname += replace_in_dirname[part_to_replace].split(":")[1] + "/"
                                        plotname = plotname.replace(part_to_replace, "")
                                    else:
                                        plotname = plotname.replace(part_to_replace, replace_in_dirname[part_to_replace])
                            newname = subdirname + plotname
                            if clean_keyname != "":
                                newname += "_"+ clean_keyname
                            # Save plot
                            if save:
                                name = output_dir + newname + ".png"
                                if (Type=="HSCP"):
                                    name =  name.replace("analysis","SignalOnly")
                                elif (Type=="Norm"):
                                    name =  name.replace("analysis","NormalTracks")
                                else:
                                    print("Please specify if Type is HSCP or Norm")
                                if not os.path.exists(os.path.dirname(name)): os.makedirs(os.path.dirname(name))
                                obj.Draw()
                                can.SaveAs(name)
                                #can.SaveAs(name.replace(".png",".pdf"))
                                #can.SaveAs(name.replace(".png",".C"))
                                can.Close()
                            else:
                                print keyname+"   "+newname
