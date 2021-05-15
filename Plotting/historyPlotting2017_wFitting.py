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
    mg.Draw('AP')
    
    linFitIOV1 = ROOT.TF1("linFitIOV1","[0]*x+[1]",290543,297281)
    linFitIOV1.SetParameter(1,25000.)
    linFitIOV1.SetParameter(0,-10.)
    mg.Fit("linFitIOV1","RQW")
    print("correctionOnTheCorrFactorIOV1  = "+str(-linFitIOV1.GetParameter(0)/(linFitIOV1.GetParameter(0)*290543+linFitIOV1.GetParameter(1))))
    
    
    linFitIOV2 = ROOT.TF1("linFitIOV2","[0]*x+[1]",297281,298653)
    linFitIOV2.SetParameter(1,25000.)
    linFitIOV2.SetParameter(0,-10.)
    mg.Fit("linFitIOV2","RQW")
    print("correctionOnTheCorrFactorIOV2  = "+str(-linFitIOV2.GetParameter(0)/(linFitIOV2.GetParameter(0)*297281+linFitIOV2.GetParameter(1))))
    
    linFitIOV3 = ROOT.TF1("linFitIOV3","[0]*x+[1]",298653,299443)
    linFitIOV3.SetParameter(1,25000.)
    linFitIOV3.SetParameter(0,-10.)
    mg.Fit("linFitIOV3","RQW")
    print("correctionOnTheCorrFactorIOV3  = "+str(-linFitIOV3.GetParameter(0)/(linFitIOV3.GetParameter(0)*298653+linFitIOV3.GetParameter(1))))
    
    linFitIOV4 = ROOT.TF1("linFitIOV4","[0]*x+[1]",299443,300389)
    linFitIOV4.SetParameter(1,25000.)
    linFitIOV4.SetParameter(0,-10.)
    mg.Fit("linFitIOV4","RQW")
    print("correctionOnTheCorrFactorIOV4  = "+str(-linFitIOV4.GetParameter(0)/(linFitIOV4.GetParameter(0)*299443+linFitIOV4.GetParameter(1))))

    linFitIOV5 = ROOT.TF1("linFitIOV5","[0]*x+[1]",300389,301046)
    linFitIOV5.SetParameter(1,25000.)
    linFitIOV5.SetParameter(0,-10.)
    mg.Fit("linFitIOV5","RQW")
    print("correctionOnTheCorrFactorIOV5  = "+str(-linFitIOV5.GetParameter(0)/(linFitIOV5.GetParameter(0)*300389+linFitIOV5.GetParameter(1))))
    
    linFitIOV6 = ROOT.TF1("linFitIOV6","[0]*x+[1]",301046,302131)
    linFitIOV6.SetParameter(1,25000.)
    linFitIOV6.SetParameter(0,-10.)
    mg.Fit("linFitIOV6","RQW")
    print("correctionOnTheCorrFactorIOV6  = "+str(-linFitIOV6.GetParameter(0)/(linFitIOV6.GetParameter(0)*301046+linFitIOV6.GetParameter(1))))
    
    linFitIOV7 = ROOT.TF1("linFitIOV7","[0]*x+[1]",302131,303790)
    linFitIOV7.SetParameter(1,25000.)
    linFitIOV7.SetParameter(0,-10.)
    mg.Fit("linFitIOV7","RQW")
    print("correctionOnTheCorrFactorIOV7  = "+str(-linFitIOV7.GetParameter(0)/(linFitIOV7.GetParameter(0)*302131+linFitIOV7.GetParameter(1))))
    
    linFitIOV8 = ROOT.TF1("linFitIOV8","[0]*x+[1]",303790,303998)
    linFitIOV8.SetParameter(1,25000.)
    linFitIOV8.SetParameter(0,-10.)
    mg.Fit("linFitIOV8","RQW")
    print("correctionOnTheCorrFactorIOV8  = "+str(-linFitIOV8.GetParameter(0)/(linFitIOV8.GetParameter(0)*303790+linFitIOV8.GetParameter(1))))
    
    linFitIOV9 = ROOT.TF1("linFitIOV9","[0]*x+[1]",303998,304911)
    linFitIOV9.SetParameter(1,25000.)
    linFitIOV9.SetParameter(0,-10.)
    mg.Fit("linFitIOV9","RQW")
    print("correctionOnTheCorrFactorIOV9  = "+str(-linFitIOV9.GetParameter(0)/(linFitIOV9.GetParameter(0)*303998+linFitIOV9.GetParameter(1))))
    
    linFitIOV10 = ROOT.TF1("linFitIOV10","[0]*x+[1]",304911,306211)
    linFitIOV10.SetParameter(1,25000.)
    linFitIOV10.SetParameter(0,-10.)
    mg.Fit("linFitIOV10","RQW")
    print("correctionOnTheCorrFactorIOV10  = "+str(-linFitIOV10.GetParameter(0)/(linFitIOV10.GetParameter(0)*304911+linFitIOV10.GetParameter(1))))
    


    
    
    
    myline300389R=ROOT.TLine(300389,0,300389,40000)
    myline301046R=ROOT.TLine(301046,0,301046,40000)
    myline302131R=ROOT.TLine(302131,0,302131,40000)
    myline303790R=ROOT.TLine(303790,0,303790,40000)
    myline303998R=ROOT.TLine(303998,0,303998,40000)
    myline304911R=ROOT.TLine(304911,0,304911,40000)
    
    leg = ROOT.TLegend(.63,.12,.87,.33)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.AddEntry(g2,"With correction from templates","LP")
    leg.Draw()
    
    linFitIOV1.Draw("SAME")
    linFitIOV2.Draw("SAME")
    linFitIOV3.Draw("SAME")
    linFitIOV4.Draw("SAME")
    linFitIOV5.Draw("SAME")
    linFitIOV6.Draw("SAME")
    linFitIOV7.Draw("SAME")
    linFitIOV8.Draw("SAME")
    linFitIOV9.Draw("SAME")
                
    
    
    
       #myline313041R=ROOT.TLine(313041,0,313041,40000)
    #myline314881R=ROOT.TLine(314881,0,314881,40000)
    myline290543R=ROOT.TLine(290543,0,290543,40000)
    myline297281R=ROOT.TLine(297281,0,297281,40000)
    myline298653R=ROOT.TLine(298653,0,298653,40000)
    myline299443R=ROOT.TLine(299443,0,299443,40000)
    myline300389R=ROOT.TLine(300389,0,300389,40000)
    myline301046R=ROOT.TLine(301046,0,301046,40000)
    myline302131R=ROOT.TLine(302131,0,302131,40000)
    myline303790R=ROOT.TLine(303790,0,303790,40000)
    myline303998R=ROOT.TLine(303998,0,303998,40000)
    myline304911R=ROOT.TLine(304911,0,304911,40000)
    
    #Gain changes
    myline299443RG=ROOT.TLine(299443,0,299443,40000)
    myline319937RG=ROOT.TLine(319937,0,319937,40000)
    myline319943RG=ROOT.TLine(319943,0,319943,40000)
    myline300389RG=ROOT.TLine(300389,0,300389,40000)
    myline303998RG=ROOT.TLine(303998,0,303998,40000)
    myline326851RG=ROOT.TLine(326851,0,326851,40000)
    
    #myline314881R.SetLineColor(3)
    myline290543R.SetLineColor(3)
    myline297281R.SetLineColor(3)
    #myline317485R.SetLineColor(3)
    #myline317527R.SetLineColor(3)
    #myline317661R.SetLineColor(3)
    myline298653R.SetLineColor(3)
    myline299443R.SetLineColor(3)
    myline300389R.SetLineColor(3)
    myline301046R.SetLineColor(3)
    myline302131R.SetLineColor(3)
    myline303790R.SetLineColor(3)
    myline303998R.SetLineColor(3)
    myline304911R.SetLineColor(3)
    
    #gain
    myline299443RG.SetLineColor(7)
    myline319937RG.SetLineColor(7)
    myline319943RG.SetLineColor(7)
    myline300389RG.SetLineColor(7)
    myline303998RG.SetLineColor(7)
    myline326851RG.SetLineColor(7)
    
    #myline314881R.SetLineStyle(3)
    myline290543R.SetLineStyle(3)
    myline297281R.SetLineStyle(3)
    #myline317485R.SetLineStyle(3)
    #myline317527R.SetLineStyle(3)
    #myline317661R.SetLineStyle(3)
    myline298653R.SetLineStyle(3)
    myline299443R.SetLineStyle(3)
    myline300389R.SetLineStyle(3)
    myline301046R.SetLineStyle(3)
    myline302131R.SetLineStyle(3)
    myline303790R.SetLineStyle(3)
    myline303998R.SetLineStyle(3)
    myline304911R.SetLineStyle(3)
    
    #gain
    myline299443RG.SetLineStyle(3)
    myline319937RG.SetLineStyle(3)
    myline319943RG.SetLineStyle(3)
    myline300389RG.SetLineStyle(3)
    myline303998RG.SetLineStyle(3)
    myline326851RG.SetLineStyle(3)
    
    #myline313041R.Draw()
    #myline314881R.Draw()
    myline290543R.Draw()
    myline297281R.Draw()
    #myline317485R.Draw()
    #myline317527R.Draw()
    #myline317661R.Draw()
    myline298653R.Draw()
    myline299443R.Draw()
    myline300389R.Draw()
    myline301046R.Draw()
    myline302131R.Draw()
    myline303790R.Draw()
    myline303998R.Draw()
    myline304911R.Draw()
    
    #gain
#    myline299443RG.Draw()
#    myline319937RG.Draw()
#    myline319943RG.Draw()
#    myline300389RG.Draw()
#    myline303998RG.Draw()
#    myline326851RG.Draw()
    
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
