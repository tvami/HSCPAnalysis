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
        
        p2 = ff.GetParameter(2)
        p2_2 = ff2.GetParameter(2)
        p2_err = ff.GetParError(2)
        p2_err_2 = ff2.GetParError(2)
        #print p2, p2_err
        #print p2_2, p2_err_2
        if p2 >  10E20: continue
        if p2_2 >  10E20 : continue
        
        if p2 <= 0: continue
        if p2_2 <= 0 : continue
        if p2 > 1000000 : continue
        if p2_2 > 1000000 : continue
        if p2_err > 10000 :
            #print("p2_err is :",p2_err)
            continue
        if p2_err_2 > 10000 :
            #print("p2_err_2 is :",p2_err_2)
            continue
        #print (rn, p2,p2_err)
        #print (rn, p2_2,p2_err_2)
        x.append(float(rn))
        ex.append(0.)
        y.append(p2)
        y2.append(p2_2)
        ey.append(p2_err)
        ey2.append(p2_err_2)


    g = ROOT.TGraphErrors(len(x), x, y, ex, ey)
    g2 = ROOT.TGraphErrors(len(x), x, y2, ex, ey2)
    
    g.GetYaxis().SetLimits( 0, 40000)
    g2.GetYaxis().SetLimits( 0, 40000)

    g.SetMarkerStyle(20)
    g.SetMarkerColor(2)
    g2.SetMarkerStyle(20)
    g2.SetMarkerColor(4)
    
    mg =  ROOT.TMultiGraph()
    mg.Add(g,"p");
    mg.Add(g2,"p");
    title = "MPV of charges on " + layer + " vs run numbers"
    mg.SetTitle(title+";Run numbers;MPV of charges [electrons]")
    mg.GetYaxis().SetLimits( 0, 40000)
    
    c = ROOT.TCanvas('c', 'c')

    mg.Draw('AP')
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
        
