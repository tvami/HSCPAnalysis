import ROOT, sys, os, time, re
from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog date_tag")
(opt,args) = parser.parse_args()

#name = "2019_09_15"
#name = "2019_09_29_part1_replot"
#name = "2019_11_04"
#name = "2019_11_06"
#name = "2020_01_20"
#name = "2020_01_23"
#name = "2020_01_25"
#name = "2020_01_29"
#name = "2020_02_03"
#name = "2020_02_05"
#name = "2020_03_02"
#name = "2020_03_11"
#name = "2020_03_21"
#name = "2020_04_10"
name = "2020_05_21-2"

save = 1

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetOptStat(1)
ROOT.gStyle.SetPalette(1)

if len(args)==0: args=[name]

# Only remove stuff for before matching
clean_dirname = {
    "Counts_vs_" : "",
    "Bins"       : "",
    #"RazorLep"   : "Razor",
    #"Razor1vl"   : "Razor",
    #"Razorll"    : "Razor",
    #"RazorNoPho" : "Razor",
    #"JetAK8"     : "",
    }
replace_in_dirname = {
    #"METFine" : "dir:Fine_binning/MET",
    #"MRFine"  : "dir:Fine_binning/MRFine",
    #"R2Fine"  : "dir:Fine_binning/R2Fine",
    "ROCCurve_vs_"        : "",
    "EleFromW_"           : "dir:ROC_Curves/EleFromW",
    "MuFromW_"            : "dir:ROC_Curves/MuFromW",
    #"EleFromTop_"         : "dir:ROC_Curves/EleFromTop",
    #"EleNuFromTop_"       : "dir:ROC_Curves/EleNuFromTop",
    #"MuFromTop_"          : "dir:ROC_Curves/MuFromTop",
    #"MuNuFromTop_"        : "dir:ROC_Curves/MuNuFromTop",
    #"EleFromHardProcess_" : "dir:ROC_Curves/EleFromHardProcess",
    #"MuFromHardProcess_"  : "dir:ROC_Curves/MuFromHardProcess",
    #"EleLepTop_"          : "dir:ROC_Curves/EleLepTop",
    #"EleNuFromLepTop_"    : "dir:ROC_Curves/EleNuFromLepTop",
    #"MuLepTop_"           : "dir:ROC_Curves/MuLepTop",
    #"MuNuFromLepTop_"     : "dir:ROC_Curves/MuNuFromLepTop",
    #"HadTop_"             : "dir:ROC_Curves/HadTop",
    #"HadW_"               : "dir:ROC_Curves/HadW",
    #"HadZ_"               : "dir:ROC_Curves/HadZ",
    #"HadH_"               : "dir:ROC_Curves/HadH",
    }
# Modify plotname
replace_in_plotname = {
    # Boost SRs
    "CR_Fake"             : "dir:Plots_"+name,
    "CR_"                 : "dir:Plots_"+name,
    "SR_"                 : "dir:Plots_"+name,
    "Val_"                : "dir:Plots_"+name,
    }

# list of: plot directories, match strings
plot_these = [
    # Stack plots
    [
        [".*Eta.*",".*PtBins.*"],
        [".*Gen*"],
    ]
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
    
    for  name in args:
        input_file = "Plotter_out_2020_05_21.root" 
        output_dir = "./"

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
            match_dir = False
            for pattern in dir_patterns:
                if re.match(pattern, clean_dirname_to_match):
                    #print(clean_dirname_to_match)
                    match_dir = True
                    break
            if match_dir:
                #print "Matched directory: "+dirname
                curr_dir = f.GetDirectory(dirname)
                for i in range(0, curr_dir.GetListOfKeys().GetEntries()):
                    # Match the plot of interest
                    keyname = curr_dir.GetListOfKeys().At(i).GetName()
                    #if re.match("StackPlot", keyname):
                    #match_plot = True 
                    #for pattern in plot_patterns:
                        #if (re.match("StackPlot", keyname) or re.match("HLT", keyname)):
                            #match_plot = False
                            #break
                            
                    #if re.match("StackPlot", keyname):
                        #print(keyname)
                        #continue
                    if re.match("HadronicMeasurements", keyname):
                        print("Skipped:",keyname)
                        continue
                    else:
                        #print(keyname)
                        # The plot should be TCanvas
                        obj = f.Get(dirname+"/"+keyname)
                        if obj.InheritsFrom("TCanvas"):
                            can = obj
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
                                if not os.path.exists(os.path.dirname(name)): os.makedirs(os.path.dirname(name))
                                can.Draw()
                                #print(name)
                                can.SaveAs(name)
                                #can.SaveAs(name.replace(".png",".pdf"))
                                #can.SaveAs(name.replace(".png",".C"))
                                can.Close()
                            else:
                                print keyname+"   "+newname
