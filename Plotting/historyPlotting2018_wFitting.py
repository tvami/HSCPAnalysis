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
        #print(rn, f_name)
        f_in = ROOT.TFile(f_name)
        if f_in.IsZombie(): continue
        
        h1 = ROOT.gDirectory.Get(h_name)
        if not h1 :
            #print("h1 does not exists for the file")
            #print(rn, f_name)
            continue
        h2 = ROOT.gDirectory.Get(h_name2)
        if not h2 :
            #print("h2 does not exists for the file")
            #print(rn, f_name)
            continue
        
        fun_name = ""
        fun_name2 = ""
        o_list = h1.GetListOfFunctions()
        for o in o_list:
            o_name = o.GetName()
            if re.search('f[0-9]*', o_name):
                fun_name = o_name
        
        o_list2 = h2.GetListOfFunctions()
        for o in o_list2:
            o_name2 = o.GetName()
            if re.search('f[0-9]*', o_name2):
                fun_name2 = o_name2
                
        ff = h1.GetFunction(fun_name)
        ff2 = h2.GetFunction(fun_name2)
        if not ff : continue
        if not ff2: continue
        
        p2 = ff.GetParameter(2)
        p2_2 = ff2.GetParameter(2)
        p2_err = ff.GetParError(2)
        p2_err_2 = ff2.GetParError(2)
        #print p2, p2_err
        #print p2_2, p2_err_2
        
        if p2 <= 5000: continue
        if p2_2 <= 5000 : continue
        if p2 > 1000000 : continue
        if p2_2 > 1000000 : continue
        if p2_err > 1000 :
            #print("p2_err is :",p2_err)
            continue
        if p2_err_2 > 1000 :
            #print("p2_err_2 is :",p2_err_2)
            continue
        #print (rn, p2,p2_err)
        #print (rn, p2_2,p2_err_2)
        #print (rn, p2_2/p2,(p2_2/p2)*(p2_err/p2+p2_err_2/p2_2))

        x.append(float(rn))
        ex.append(0.)
        y.append(p2)
        y2.append(p2_2)
        ey.append(p2_err)
        ey2.append(p2_err_2)
        correctionFactor.append(p2_2/p2)
        eCorrectionFactor.append((p2_2/p2)*(p2_err/p2+p2_err_2/p2_2))


    g2 = ROOT.TGraphErrors(len(x), x, y2, ex, ey2)
    
    g2.SetMarkerStyle(20)
    g2.SetMarkerColor(4) #blue
    g2.SetMarkerSize(.5)
    
    mg =  ROOT.TMultiGraph()
    mg.Add(g2,"p");
    title = "MPV of charges on " + layer + " vs run numbers"
    mg.SetTitle(title+";Run numbers;MPV of charges [electrons]")
    mg.GetXaxis().SetMaxDigits(6)
    mg.SetMinimum(10)
    mg.SetMaximum(40000)
    
    c = ROOT.TCanvas('c', 'c',1600,800)
    #c.Divide(1,2)
    
    # Upper histogram plot is pad1
    #pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    #pad1.SetBottomMargin(0)  # joins upper and lower plot
    #pad1.SetGridx()
    #pad1.Draw()
    # Lower ratio plot is pad2
    #c.cd()  # returns to main canvas before defining pad2
    #pad2 = ROOT.TPad("pad2", "pad2", 0, 0.05, 1, 0.3)
    #pad2.SetTopMargin(0.001)  # joins upper and lower plot
    #pad2.SetBottomMargin(0.2)
    #pad2.Draw("SAME")
    #pad1.cd()
    
    #c.cd(1)
    mg.Draw('AP')


    #linFitIOVBefore0 = ROOT.TF1("linFitIOVBefore0","[0]*x+[1]",313041,314880)
    #linFitIOVBefore0.SetParameter(1,25000.)
    #linFitIOVBefore0.SetParameter(0,-10.)
    #mg.Fit("linFitIOVBefore0","RQW")
    #print("correctionOnTheCorrFactorIOVBefore0  = "+str(-linFitIOVBefore0.GetParameter(0)/(linFitIOVBefore0.GetParameter(0)*313041+linFitIOVBefore0.GetParameter(1))))

    linFitIOV0 = ROOT.TF1("linFitIOV0","[0]*x+[1]",314881,316757)
    linFitIOV0.SetParameter(1,25000.)
    linFitIOV0.SetParameter(0,-10.)
    mg.Fit("linFitIOV0","RQW")
    print("correctionOnTheCorrFactorIOV0  = "+str(-linFitIOV0.GetParameter(0)/(linFitIOV0.GetParameter(0)*314881+linFitIOV0.GetParameter(1))))
    
    linFitIOV1 = ROOT.TF1("linFitIOV1","[0]*x+[1]",316758,317475)
    linFitIOV1.SetParameter(1,25000.)
    linFitIOV1.SetParameter(0,-10.)
    mg.Fit("linFitIOV1","RQW")
    print("correctionOnTheCorrFactorIOV1  = "+str(-linFitIOV1.GetParameter(0)/(linFitIOV1.GetParameter(0)*316758+linFitIOV1.GetParameter(1))))
    
    
    linFitIOV2 = ROOT.TF1("linFitIOV2","[0]*x+[1]",317475,317664)
    linFitIOV2.SetParameter(1,25000.)
    linFitIOV2.SetParameter(0,-10.)
    mg.Fit("linFitIOV2","RQW")
    print("correctionOnTheCorrFactorIOV2  = "+str(-linFitIOV2.GetParameter(0)/(linFitIOV2.GetParameter(0)*317475+linFitIOV2.GetParameter(1))))
    
    linFitIOV3 = ROOT.TF1("linFitIOV3","[0]*x+[1]",317664,318227)
    linFitIOV3.SetParameter(1,25000.)
    linFitIOV3.SetParameter(0,-10.)
    mg.Fit("linFitIOV3","RQW")
    print("correctionOnTheCorrFactorIOV3  = "+str(-linFitIOV3.GetParameter(0)/(linFitIOV3.GetParameter(0)*317664+linFitIOV3.GetParameter(1))))
    
    linFitIOV4 = ROOT.TF1("linFitIOV4","[0]*x+[1]",318227,320377)
    linFitIOV4.SetParameter(1,25000.)
    linFitIOV4.SetParameter(0,-10.)
    mg.Fit("linFitIOV4","RQW")
    print("correctionOnTheCorrFactorIOV4  = "+str(-linFitIOV4.GetParameter(0)/(linFitIOV4.GetParameter(0)*318227+linFitIOV4.GetParameter(1))))

    linFitIOV5 = ROOT.TF1("linFitIOV5","[0]*x+[1]",320377,321831)
    linFitIOV5.SetParameter(1,25000.)
    linFitIOV5.SetParameter(0,-10.)
    mg.Fit("linFitIOV5","RQW")
    print("correctionOnTheCorrFactorIOV5  = "+str(-linFitIOV5.GetParameter(0)/(linFitIOV5.GetParameter(0)*320377+linFitIOV5.GetParameter(1))))
    
    linFitIOV6 = ROOT.TF1("linFitIOV6","[0]*x+[1]",321831,322510)
    linFitIOV6.SetParameter(1,25000.)
    linFitIOV6.SetParameter(0,-10.)
    mg.Fit("linFitIOV6","RQW")
    print("correctionOnTheCorrFactorIOV6  = "+str(-linFitIOV6.GetParameter(0)/(linFitIOV6.GetParameter(0)*321831+linFitIOV6.GetParameter(1))))
    
    linFitIOV7 = ROOT.TF1("linFitIOV7","[0]*x+[1]",322510,322603)
    linFitIOV7.SetParameter(1,25000.)
    linFitIOV7.SetParameter(0,-10.)
    mg.Fit("linFitIOV7","RQW")
    print("correctionOnTheCorrFactorIOV7  = "+str(-linFitIOV7.GetParameter(0)/(linFitIOV7.GetParameter(0)*322510+linFitIOV7.GetParameter(1))))
    
    linFitIOV8 = ROOT.TF1("linFitIOV8","[0]*x+[1]",322603,323232)
    linFitIOV8.SetParameter(1,25000.)
    linFitIOV8.SetParameter(0,-10.)
    mg.Fit("linFitIOV8","RQW")
    print("correctionOnTheCorrFactorIOV8  = "+str(-linFitIOV8.GetParameter(0)/(linFitIOV8.GetParameter(0)*322603+linFitIOV8.GetParameter(1))))
    
    linFitIOV9 = ROOT.TF1("linFitIOV9","[0]*x+[1]",323232,324245)
    linFitIOV9.SetParameter(1,25000.)
    linFitIOV9.SetParameter(0,-10.)
    mg.Fit("linFitIOV9","RQW")
    print("correctionOnTheCorrFactorIOV9  = "+str(-linFitIOV9.GetParameter(0)/(linFitIOV9.GetParameter(0)*323232+linFitIOV9.GetParameter(1))))
    
    linFitIOV10 = ROOT.TF1("linFitIOV10","[0]*x+[1]",324246,32500)
    linFitIOV10.SetParameter(1,25000.)
    linFitIOV10.SetParameter(0,-10.)
    mg.Fit("linFitIOV10","RQW")
    print("correctionOnTheCorrFactorIOV10  = "+str(-linFitIOV10.GetParameter(0)/(linFitIOV10.GetParameter(0)*324246+linFitIOV10.GetParameter(1))))
    
    leg = ROOT.TLegend(.63,.12,.87,.33)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.AddEntry(g2,"With correction from templates","LP")
    leg.Draw()
   
    #linFitIOVBefore0.Draw("SAME")
    linFitIOV0.Draw("SAME") 
    linFitIOV1.Draw("SAME")
    linFitIOV2.Draw("SAME")
    linFitIOV3.Draw("SAME")
    linFitIOV4.Draw("SAME")
    linFitIOV5.Draw("SAME")
    linFitIOV6.Draw("SAME")
    linFitIOV7.Draw("SAME")
    linFitIOV8.Draw("SAME")
    linFitIOV9.Draw("SAME")
    linFitIOV10.Draw("SAME")
    
    #myline313041R=ROOT.TLine(313041,0,313041,40000)
    #myline314881R=ROOT.TLine(314881,0,314881,40000)
    myline316758R=ROOT.TLine(316758,0,316758,40000)
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
    
    tex3 = ROOT.TLatex(0.14,0.96,"Internal 2018");
    #tex3 = ROOT.TLatex(0.28,0.94,"Work in Progress 2018"); #if there is 10^x
    tex3.SetNDC();
    tex3.SetTextFont(52);
    tex3.SetTextSize(0.0285);
    tex3.SetLineWidth(2);
    tex3.Draw("SAME");
    
    #pad2.cd()
    #c.cd(2)
    #factor = ROOT.TGraphErrors(len(x), x, correctionFactor, ex, eCorrectionFactor)
    #factor.GetXaxis().SetMaxDigits(6)
    #factor.SetTitle(";Run numbers;Corr factor")
    #factor.Draw("AP")
    
    #c.cd(0)
    
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
    draw(file_map, "analysis/h701_n1", "analysis/h701_n1c", "Layer 1", "L1ChargeHistory.png")
    draw(file_map, "analysis/h701_n2", "analysis/h701_n2c", "Layer 2", "L2ChargeHistory.png")
    draw(file_map, "analysis/h701_n3", "analysis/h701_n3c", "Layer 3", "L3ChargeHistory.png")
    draw(file_map, "analysis/h701_n4", "analysis/h701_n4c", "Layer 4", "L4ChargeHistory.png")
    draw(file_map, "analysis/h702_nr1", "analysis/h702_nr1c", "Ring 1", "R1ChargeHistory.png")
    draw(file_map, "analysis/h702_nr2", "analysis/h702_nr2c", "Ring 2", "R2ChargeHistory.png")
        
