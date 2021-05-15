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
        
        runnumber = float(rn)

        correctionOnTheCorrFactorIOV1  = 3.04610088475e-05
        correctionOnTheCorrFactorIOV2  = -3.73547306594e-05
        correctionOnTheCorrFactorIOV3  = 0.00020407567075
        correctionOnTheCorrFactorIOV4  = 3.08601351651e-05
        correctionOnTheCorrFactorIOV5  = 1.37114905558e-05
        correctionOnTheCorrFactorIOV6  = 6.42378603724e-07
        correctionOnTheCorrFactorIOV7  = 5.12817895229e-05
        correctionOnTheCorrFactorIOV8  = -0.000154040990698
        correctionOnTheCorrFactorIOV9  = 6.74145838168e-05
        correctionOnTheCorrFactorIOV10  = 4.04053326176e-05

        correctedCorrFactor = 0.0
        x.append(runnumber)
        ex.append(0.)
        if   (runnumber>=290543 and runnumber<297281) :
            correctedCorrFactor = (runnumber-290543)*correctionOnTheCorrFactorIOV1*corrFactor+corrFactor
        elif (runnumber>=297281 and runnumber<298653) :
            correctedCorrFactor = (runnumber-297281)*correctionOnTheCorrFactorIOV2*corrFactor+corrFactor
        elif (runnumber>=298653 and runnumber<299443) :
            correctedCorrFactor = (runnumber-298653)*correctionOnTheCorrFactorIOV3*corrFactor+corrFactor
        elif (runnumber>=299443 and runnumber<300389) :
            correctedCorrFactor = (runnumber-299443)*correctionOnTheCorrFactorIOV4*corrFactor+corrFactor
        elif (runnumber>=300389 and runnumber<301046) :
            correctedCorrFactor = (runnumber-300389)*correctionOnTheCorrFactorIOV5*corrFactor+corrFactor
        elif (runnumber>=301046 and runnumber<302131) :
            correctedCorrFactor = (runnumber-301046)*correctionOnTheCorrFactorIOV6*corrFactor+corrFactor
        elif (runnumber>=302131 and runnumber<303790) :
            correctedCorrFactor = (runnumber-302131)*correctionOnTheCorrFactorIOV7*corrFactor+corrFactor
        elif (runnumber>=303790 and runnumber<303998) :
            correctedCorrFactor = (runnumber-303790)*correctionOnTheCorrFactorIOV8*corrFactor+corrFactor
        elif (runnumber>=303998 and runnumber<304911) :
            correctedCorrFactor = (runnumber-303998)*correctionOnTheCorrFactorIOV9*corrFactor+corrFactor
        elif (runnumber>=304911) :
            correctedCorrFactor = (runnumber-304911)*correctionOnTheCorrFactorIOV10*corrFactor+corrFactor
            
        y.append(correctedCorrFactor)
        ey.append(corrFactorUnc)
        print(str(rn)+" "+str(correctedCorrFactor)+" "+str(corrFactorUnc))


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

#    draw(file_map, "analysis/corrFactL2", "analysis/h701_n2c", "Layer 2", "L2ChargeCorrFactors.png")
#    draw(file_map, "analysis/corrFactL3", "analysis/h701_n3c", "Layer 3", "L3ChargeCorrFactors.png")
    draw(file_map, "analysis/corrFactL4", "analysis/h701_n4c", "Layer 4", "L4ChargeCorrFactors.png")
#    draw(file_map, "analysis/corrFactR1", "analysis/h702_nr1c", "Ring 1", "R1ChargeCorrFactors.png")
#    draw(file_map, "analysis/corrFactR2", "analysis/h702_nr2c", "Ring 2", "R2ChargeCorrFactors.png") 
