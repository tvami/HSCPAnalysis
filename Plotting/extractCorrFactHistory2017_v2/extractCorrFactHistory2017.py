import argparse
import ROOT
from array import array
import re 
import tdrstyle

def make_filemap(infile):
    file_map = {}
    with open (infile) as ifp:
        for line in ifp:
            rline = line.rstrip().split()
            if len(rline) < 2:
                continue
            rn = int(rline[0])
            file_name = rline[1]
            file_map[rn] = file_name
            
    return file_map

def draw(file_map, h_name, h_name2, layer,  outfile):
    x = array('f')
    ex = array('f')
    y = array('f')
    ey = array('f')
    y2 = array('f')
    ey2 = array('f')
    correctionFactor = array('f')
    eCorrectionFactor = array('f')

    for rn, f_name in file_map.iteritems():
        f_in = ROOT.TFile(f_name)
        if f_in.IsZombie(): continue
        
        h1 = ROOT.gDirectory.Get(h_name)
        if not h1 : continue
        
        corrFactor = h1.GetMean()
        if (corrFactor == 0.0) : continue
        corrFactorUnc = h1.GetRMS()
        
#        print(str(rn)+" "+str(corrFactor)+" "+str(corrFactorUnc))
        
        x.append(float(rn))
        ex.append(0.)
        y.append(corrFactor)
        ey.append(corrFactorUnc)


    g = ROOT.TGraphErrors(len(x), x, y, ex, ey)
    
    #g.SetMarkerStyle(24)
    g.SetMarkerStyle(8)
    g.SetMarkerSize(.5)
    g.SetMarkerColor(2) #red
    g.SetLineStyle(1)
    g.SetLineWidth(1504)
    g.SetFillStyle(3005)
    g.SetFillColor(17)
    g.SetLineColor(17)
    
    mg =  ROOT.TMultiGraph()
    mg.Add(g,"p");
    title = "Charge correction factor on " + layer + " vs run numbers"
    mg.SetTitle(title+";Run numbers;Charge correction factors")
    mg.GetXaxis().SetMaxDigits(6)
    mg.SetMinimum(0.0)
    mg.SetMaximum(2.5)
    
    c = ROOT.TCanvas('c', 'c',1600,800)

    mg.Draw('AC')    
        
    
    #myline313041R=ROOT.TLine(313041,0,313041,40000)
    #myline314881R=ROOT.TLine(314881,0,314881,40000)
    myline290543R=ROOT.TLine(316758,0,316758,40000)
    myline317475R=ROOT.TLine(317475,0,317475,40000)
    #myline317485R=ROOT.TLine(317485,0,317485,40000)
    #myline317527R=ROOT.TLine(317527,0,317527,40000)
    #myline317661R=ROOT.TLine(317661,0,317661,40000)
    myline317664R=ROOT.TLine(317664,0,317664,40000)
    myline318227R=ROOT.TLine(318227,0,318227,40000)
    myline320377R=ROOT.TLine(320377,0,320377,40000)
    myline321831R=ROOT.TLine(321831,0,321831,40000)
    myline322510R=ROOT.TLine(322510,0,322510,40000)
    myline322603R=ROOT.TLine(322603,0,322603,40000)
    myline323232R=ROOT.TLine(323232,0,323232,40000)
    myline324245R=ROOT.TLine(324245,0,324245,40000)
    
    #Gain changes
    myline318227RG=ROOT.TLine(318227,0,318227,40000)
    myline319937RG=ROOT.TLine(319937,0,319937,40000)
    myline319943RG=ROOT.TLine(319943,0,319943,40000)
    myline320377RG=ROOT.TLine(320377,0,320377,40000)
    myline323232RG=ROOT.TLine(323232,0,323232,40000)
    myline326851RG=ROOT.TLine(326851,0,326851,40000)
    
    #myline314881R.SetLineColor(3)
    myline316758R.SetLineColor(3)
    myline317475R.SetLineColor(3)
    #myline317485R.SetLineColor(3)
    #myline317527R.SetLineColor(3)
    #myline317661R.SetLineColor(3)
    myline317664R.SetLineColor(3)
    myline318227R.SetLineColor(3)
    myline320377R.SetLineColor(3)
    myline321831R.SetLineColor(3)
    myline322510R.SetLineColor(3)
    myline322603R.SetLineColor(3)
    myline323232R.SetLineColor(3)
    myline324245R.SetLineColor(3)
    
    #gain
    myline318227RG.SetLineColor(7)
    myline319937RG.SetLineColor(7)
    myline319943RG.SetLineColor(7)
    myline320377RG.SetLineColor(7)
    myline323232RG.SetLineColor(7)
    myline326851RG.SetLineColor(7)
    
    #myline314881R.SetLineStyle(3)
    myline316758R.SetLineStyle(3)
    myline317475R.SetLineStyle(3)
    #myline317485R.SetLineStyle(3)
    #myline317527R.SetLineStyle(3)
    #myline317661R.SetLineStyle(3)
    myline317664R.SetLineStyle(3)
    myline318227R.SetLineStyle(3)
    myline320377R.SetLineStyle(3)
    myline321831R.SetLineStyle(3)
    myline322510R.SetLineStyle(3)
    myline322603R.SetLineStyle(3)
    myline323232R.SetLineStyle(3)
    myline324245R.SetLineStyle(3)
    
    #gain
    myline318227RG.SetLineStyle(3)
    myline319937RG.SetLineStyle(3)
    myline319943RG.SetLineStyle(3)
    myline320377RG.SetLineStyle(3)
    myline323232RG.SetLineStyle(3)
    myline326851RG.SetLineStyle(3)
    
    #myline313041R.Draw()
    #myline314881R.Draw()
    myline316758R.Draw()
    myline317475R.Draw()
    #myline317485R.Draw()
    #myline317527R.Draw()
    #myline317661R.Draw()
    myline317664R.Draw()
    myline318227R.Draw()
    myline320377R.Draw()
    myline321831R.Draw()
    myline322510R.Draw()
    myline322603R.Draw()
    myline323232R.Draw()
    myline324245R.Draw()
    
    #gain
    myline318227RG.Draw()
    myline319937RG.Draw()
    myline319943RG.Draw()
    myline320377RG.Draw()
    myline323232RG.Draw()
    myline326851RG.Draw()
    
    tex2 = ROOT.TLatex(0.1,0.96,"CMS");
    #tex2 = ROOT.TLatex(0.20,0.94,"CMS");#if there is 10^x
    tex2.SetNDC();
    tex2.SetTextFont(61);
    tex2.SetTextSize(0.0375);
    tex2.SetLineWidth(2);
    tex2.Draw("SAME");
    
    tex3 = ROOT.TLatex(0.14,0.96,"Internal 2017");
    #tex3 = ROOT.TLatex(0.28,0.94,"Work in Progress 2017"); #if there is 10^x
    tex3.SetNDC();
    tex3.SetTextFont(52);
    tex3.SetTextSize(0.0285);
    tex3.SetLineWidth(2);
    tex3.Draw("SAME");
    
    c.SaveAs(outfile)
    
if __name__ == '__main__':
    
    ROOT.gROOT.SetBatch(True)
    tdrstyle.setTDRStyle()

 #-------------------------------------------------------------------------
 #  PARSE COMMAND LINE ARGUMENT
 #------------------------------------------------------------------------- 
    
    parser = argparse.ArgumentParser(description='Charge MPV extractor')
    parser.add_argument('-l', action="store", dest = 'inlist', default = 'HSCP_file.list')
    parser.add_argument('-n', action="store", dest = 'h_name', default = 'analysis/h701_n1')
    parser.add_argument('-n2', action="store", dest = 'h_name2', default = 'analysis/h701_n1c')
    parser.add_argument('-lay', action="store", dest = 'layer', default = 'Layer 1')
    parser.add_argument('-o', action="store", dest = 'outfile', default = 'tmp.pdf')

    p = parser.parse_args()

    file_map = make_filemap(p.inlist)
    #draw(file_map, p.h_name, p.h_name2, p.outfile)

#    draw(file_map, "analysis/corrFactL1", "analysis/h701_n1c", "Layer 1", "L1ChargeCorrFactors.png")
#    draw(file_map, "analysis/corrFactL2", "analysis/h701_n2c", "Layer 2", "L2ChargeCorrFactors.png")
#    draw(file_map, "analysis/corrFactL3", "analysis/h701_n3c", "Layer 3", "L3ChargeCorrFactors.png")
#    draw(file_map, "analysis/corrFactL4", "analysis/h701_n4c", "Layer 4", "L4ChargeCorrFactors.png")
#    draw(file_map, "analysis/corrFactR1", "analysis/h702_nr1c", "Ring 1", "R1ChargeCorrFactors.png")
#    draw(file_map, "analysis/corrFactR2", "analysis/h702_nr2c", "Ring 2", "R2ChargeCorrFactors.png") 
