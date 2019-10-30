//! PGM to read Urs' pixel tree and produce summary files.
//! Update to version r23 of the pixel tree
//! Make run specific output files
//! Add hit position to angle output info
//! Require barrel layer one on all tracks
//! Write out hit position and angles
//! Loop over several input files
//! Update to version r28 (not including simhit info)
//! Modify to handle multiple vertices
//! Fix vertex selection bug (3-Dec-2011)
//! Keep 3 layers for each subdetector, add layer info to output file
//! Test 2011 end of year templates
//! Test 2012 June templates
//! Test 2016 summer templates
//! Do different FB templates
//! Do full ring and layer templates
//! Read more general file names
//! Update to 1st 2017 IOV
//! Update for new L2 voltages and improved FPix LAs
//! Update for new L2 voltages and new. new other voltages.
//! Charge scale and threhsold changes
#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include "boost/multi_array.hpp"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <sys/time.h>
#include "SiPixelTemplate.cc"
static int theVerboseLevel = {0};
#include "VVIObjF.cc"
#include "SiPixelTemplateReco.cc"

#include "pixelTree.h"

#include <cmath>
#include <cstdlib>

#include "TROOT.h"
#include "TH1.h"
#include "TH2.h" 
#include "TProfile.h" 
#include "TVector3.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TPostScript.h"


#define PVMAX 100 
#define MUMAX 100 
#define TRACKMAX 10000 
#define CLPERTRACKMAX 20 
#define CLUSTERMAX 42000
#define DGPERCLMAX 100  
#define TKPERCLMAX 100  
#define DIGIMAX 200000
// 0 -- Normal tracks, 1 -- HSCP only
#define HSCPONLY 0


using namespace std;

  struct Pixel {
    float gx;
    float gy;
    float gz;
  };


int rocId(int col, int row) {
  int rocRow = row/80;
  int rocCol = col/52;
  int rocId = rocCol + rocRow*8;
  return rocId;
}

main(int argc, char **argv) {

  const bool SELECT_BX = false;
  const bool SELECT_LUMI = false;
  const bool SELECT_ORBIT = false;
  const bool SELECT_ON_TRACKS  = true;
  const bool SELECT_OFF_TRACKS = false;

  // Select a module 
  const bool SELECT_MODULE = false;
  const int selectLayer=1, selectLadder=-1, selectModule=-3;

//   const bool use_MUONs = false;
//   const bool use_HF = false;
//   const bool use_Tracks = false;
//   const bool use_Clus = false;
//   const bool use_ClusSim = false;
//   const bool use_Digis = true;

  const bool PRINT = false;
  int debug = 0;

  int maxClus = 30000;
  int maxDigi = 150000;
  int maxTrk =  7000;
  int maxVtx =  100;
  int maxDigiPerClu = 100;
  
  // Local variables 
  bool ydouble[TYSIZE], xdouble[TXSIZE];
  static float qscale, qscaleB, qscaleF, pcut, tkpcut, probQ, xs, ys, probQonTrack, probQonTrackTerm;
  static float probQonTrackWMulti = 1 ;
  static float xhit, yhit, xrec, yrec, sigmax, sigmay, probx, proby, cotalpha, cotbeta, locBx, locBz, xoff, yoff, xtemp, ytemp;  
  static int sfile, nrun, external, sizex, sizey, layer, llayer, module, ladder, offladder, side, disk, blade, onblade, panel, lowpt;
  static int tladp1[4], qlad[4]={3, 7, 11, 16};
  static int lumiMin = 100000, lumiMax = 0;
  static vector<int> nbin(5,0);
  int i, j, ierr, qbin, ID;
  bool bpix;
  double log10probXY, log10probQ, log10probXYQ, logprobQonTrackWMulti, qclust, qnorm, proba, probXY, probXYQ, dx, dy, TkP, xhmod, yhmod;
  static int iy, ix, ngood, nbad, speed, IDF1, ring;	
  static int bptmp[4][8], fptmp[2][3][4], bpnew;
  static int acctrk[4][11], alltrk[4][11];
  static char infile[150], header[80], outfile[150], outfile0[80], outfile1[80], outfile2[80];
  float cluster[TXSIZE][TYSIZE];
  int mrow = TXSIZE, mcol = TYSIZE;    
  float xpitch, ypitch;
  int factorial(int);

  TFile hFile( "histo.root", "RECREATE" );

  int nx=120;  
  gStyle->SetOptFit(101);
  gStyle->SetHistLineWidth(2);
  static vector<TH1F*> hp(58);
     
  hp[0] = new TH1F("h201","number of vertices",60,-0.5,59.5);
  hp[1] = new TH1F("h202","vertex x",20,-1.,1.);      
  hp[2] = new TH1F("h203","vertex y",20,-1.,1.);      
  hp[3] = new TH1F("h204","vertex z",80,-8.,8.);     
  hp[4] = new TH1F("h205","chi2/dof",50,0.,10.);      
  hp[5] = new TH1F("h206","ProbXY",50,0.,1.);      
  hp[6] = new TH1F("h207","ProbQ",50,0.,1.);
  hp[7] = new TH1F("h208","ProbXYQ",50,0.,1.);      
  hp[8] = new TH1F("h209","log10(ProbXY)",nx,-12.,0.);      
  hp[9] = new TH1F("h210","log10(ProbQ)",nx,-12.,0.);     
  hp[10] = new TH1F("h211","log10(ProbXYQ)",nx,-12.,0.);      
  hp[11] = new TH1F("h212","log10(ProbXY) (ProbQ>0.01)",nx,-12.,0.);      
  hp[12] = new TH1F("h701","Normalized Cluster Charge (BPix)",nx,0.,60000.);
  hp[13] = new TH1F("h702","Normalized Cluster Charge (FPix)",nx,0.,60000.);
  hp[14] = new TH1F("h703","Cluster Charge (ProbXY > 0.01)",nx,0.,200000.);      
  hp[15] = new TH1F("h704","Cluster Charge (ProbQ>0.01)",nx,0.,200000.);
  hp[16] = new TH1F("h705","Cluster Charge (ProbXYQ>0.01)",nx,0.,200000.);      
  hp[17] = new TH1F("h106","dx_temp (bpix); #Deltax (#mum)",nx,-50.,50.);
  hp[18] = new TH1F("h117","dy_temp (bpix); #Deltay (#mum)",nx,-50.,50.);
  hp[19] = new TH1F("h118","Track pT; pT (GeV)",nx,0.,12.);
  hp[20] = new TH1F("h119","Track momentum; p (GeV)",nx,0.,12.);
  hp[21] = new TH1F("h120","dof; ndof",20,-0.5,19.5);
  hp[22] = new TH1F("h121","fraction RH w/ good probs; frac",50,0.,1.);
  hp[23] = new TH1F("h122","Vertex NDOF; NDOF",50,-0.5,49.5);
  hp[24] = new TH1F("h123","Vertex Chi2/NDof; Chi2/Dof",40,0.,20.);
  hp[25] = new TH1F("h124","TkVx-PVx",50,-0.5,0.5);      
  hp[26] = new TH1F("h125","TkVy-PVy",50,-0.5,0.5);      
  hp[27] = new TH1F("h126","TkVz-PVz",100,-1.,1.);     
  hp[28] = new TH1F("h127","Track Quality",17,-1.5,15.5);     
  hp[29] = new TH1F("h213","Number of tracks",200,-0.5,999.5);
  hp[30] = new TH1F("h301","ProbQ on tracks (w/ multipication)",50,0.,1.);
  hp[31] = new TH1F("h302","log(probQ) on tracks (w/ multipication)",nx,-12.,0.);     
  hp[32] = new TH1F("h303","ProbQ on tracks (w/ combine)",50,0.,1.);
  hp[33] = new TH1F("h304","FPix yhit (mod pixel); y (#mum)",160,-80,80);
  hp[34] = new TH1F("h218","After cuts: Track pT; pT (GeV)",nx,0.,12.);
  hp[35] = new TH1F("h219","After cuts: Track momentum; p (GeV)",nx,0.,12.);
  hp[36] = new TH1F("h220","After cuts: BPix cotalpha; cot(#alpha)",100,-1.0,1.0);
  hp[37] = new TH1F("h221","After cuts: BPix cotbeta; cot(#beta)",200,-10.,10.);
  hp[38] = new TH1F("h222","After cuts: FPix cotalpha; cot(#alpha)",50,0.0,1.0);
  hp[39] = new TH1F("h223","After cuts: FPix |cotbeta|; |cot(#beta)|",50,0.0,1.0);
  hp[40] = new TH1F("h107","dx_temp (fpix); #Deltax (#mum)",nx,-100.,100.);
  hp[41] = new TH1F("h138","dy_temp (fpix); #Deltay (#mum)",nx,-100.,100.);
  hp[42] = new TH1F("h801","dx_temp (BP L1); #Deltax (#mum)",nx,-100.,100.);
  hp[43] = new TH1F("h802","dy_temp (BP L1); #Deltay (#mum)",nx,-100.,100.);
  hp[44] = new TH1F("h803","dx_temp (BP L2); #Deltax (#mum)",nx,-100.,100.);
  hp[45] = new TH1F("h804","dy_temp (BP L2); #Deltay (#mum)",nx,-100.,100.);
  hp[46] = new TH1F("h805","dx_temp (BP L3); #Deltax (#mum)",nx,-100.,100.);
  hp[47] = new TH1F("h806","dy_temp (BP L3); #Deltay (#mum)",nx,-100.,100.);
  hp[48] = new TH1F("h807","dx_temp (BP L4); #Deltax (#mum)",nx,-100.,100.);
  hp[49] = new TH1F("h808","dy_temp (BP L4); #Deltay (#mum)",nx,-100.,100.);
  hp[50] = new TH1F("h809","dx_temp (FP R2P1); #Deltax (#mum)",nx,-100.,100.);
  hp[51] = new TH1F("h810","dy_temp (FP R2P1); #Deltay (#mum)",nx,-100.,100.);
  hp[52] = new TH1F("h811","dx_temp (FP R1P1); #Deltax (#mum)",nx,-100.,100.);
  hp[53] = new TH1F("h812","dy_temp (FP R1P1); #Deltay (#mum)",nx,-100.,100.);
  hp[54] = new TH1F("h813","dx_temp (FP R1P2); #Deltax (#mum)",nx,-100.,100.);
  hp[55] = new TH1F("h814","dy_temp (FP R1P2); #Deltay (#mum)",nx,-100.,100.);
  hp[56] = new TH1F("h815","dx_temp (FP R2P2); #Deltax (#mum)",nx,-100.,100.);
  hp[57] = new TH1F("h816","dy_temp (FP R2P2); #Deltay (#mum)",nx,-100.,100.);
  
  // Set style for the the histograms  
  for(i=0; i<58; ++i) {
    hp[i]->SetLineColor(2);
    hp[i]->SetFillColor(38);
    hp[i]->SetMinimum(0.0);
  }
  
  TChain chain("pixelTree");
  chain.Add("pixelTreeWithSimInfo.root");
//  chain.Add("4pixelTree_Signal.root"); 
//   chain.Add("4pixelTree.root"); 
  // Template reading
//   printf("enter probability cut, momentum cut, lowpt (1 = yes), bpix scale factor (pixelav/data), fpix scale factor \n");

//   scanf("%f %f %d %f %f", &pcut, &tkpcut, &lowpt, &qscaleB, &qscaleF);
  pcut = 0.01, tkpcut = 3.0, lowpt = 1, qscaleB=1., qscaleF=1.;

  printf("probability cut = %f, momentum cut = %f, lowpt = %d, bpix scale factor = %f, fpix scale factor = %f \n", pcut, tkpcut, lowpt, qscaleB, qscaleF);
  
  for(int lay=0; lay<4; ++lay) {
      tladp1[lay] = 5*qlad[lay]+1;
  }
  
   bool newmodule[4][64][8];
   for(int lay=0; lay<4; ++lay) {
      for(int lad=0; lad<64; ++lad) {
         for(int mod=0; mod<8; ++mod) {
            newmodule[lay][lad][mod] = false;
         }
      }
   }
  
  // Initialize template store
  // Q: where are these comming from?, [model][layer][IOV]
  int bpl1 = 70; 	int bpl2 = 1513; int bpl3 = 1721; int bpl4 = 7731; int fpr2p1 = 1841; int fpr1p1 = 1951; int fpr1p2 = 1961; int fpr2p2 = 1871;
  
  //Q: what's the second index here? --> module number
  bptmp[0][0] = bpl1; bptmp[0][1] = bpl1; bptmp[0][2] = bpl1; bptmp[0][3] = bpl1; bptmp[0][4] = bpl1; bptmp[0][5] = bpl1; bptmp[0][6] = bpl1; bptmp[0][7] = bpl1; bptmp[1][0] = bpl2; bptmp[1][1] = bpl2; bptmp[1][2] = bpl2; bptmp[1][3] = bpl2; bptmp[1][4] = bpl2; bptmp[1][5] = bpl2; bptmp[1][6] = bpl2; bptmp[1][7] = bpl2; bptmp[2][0] = bpl3; bptmp[2][1] = bpl3; bptmp[2][2] = bpl3; bptmp[2][3] = bpl3; bptmp[2][4] = bpl3; bptmp[2][5] = bpl3; bptmp[2][6] = bpl3; bptmp[2][7] = bpl3; bptmp[3][0] = bpl4; bptmp[3][1] = bpl4; bptmp[3][2] = bpl4; bptmp[3][3] = bpl4; bptmp[3][4] = bpl4; bptmp[3][5] = bpl4; bptmp[3][6] = bpl4; bptmp[3][7] = bpl4;
  
  bpnew = bpl1;

// fptmp[side][disk][rindex] ;rindex = 0 [R2,P1], 1 [R1,P1], 2 [R1,P2], 3 [R2,P2]  
  fptmp[0][0][0] = fpr2p1; fptmp[0][1][0] = fpr2p1; fptmp[0][2][0] = fpr2p1; fptmp[1][0][0] = fpr2p1; fptmp[1][1][0] = fpr2p1; fptmp[1][2][0] = fpr2p1; fptmp[0][0][1] = fpr1p1; fptmp[0][1][1] = fpr1p1; fptmp[0][2][1] = fpr1p1; fptmp[1][0][1] = fpr1p1; fptmp[1][1][1] = fpr1p1; fptmp[1][2][1] = fpr1p1; fptmp[0][0][2] = fpr1p2; fptmp[0][1][2] = fpr1p2; fptmp[0][2][2] = fpr1p2; fptmp[1][0][2] = fpr1p2; fptmp[1][1][2] = fpr1p2; fptmp[1][2][2] = fpr1p2; fptmp[0][0][3] = fpr2p2; fptmp[0][1][3] = fpr2p2; fptmp[0][2][3] = fpr2p2; fptmp[1][0][3] = fpr2p2; fptmp[1][1][3] = fpr2p2; fptmp[1][2][3] = fpr2p2; 
  
  std::vector< SiPixelTemplateStore > thePixelTemp_;
  SiPixelTemplate templ(thePixelTemp_);
  
  
  
  ID = bpl1;      templ.pushfile(ID,thePixelTemp_);
  templ.interpolate(ID, 0.f, 0.f, 1.f, 1.f); // Q: interpolate to what exactly?
  xpitch = templ.xsize();
  ypitch = templ.ysize();
  printf("\n xpitch/ypitch = %f/%f \n\n", xpitch, ypitch);
  ID = bpl2;     templ.pushfile(ID,thePixelTemp_);  
  ID = bpl3;      templ.pushfile(ID,thePixelTemp_);
  ID = bpl4;      templ.pushfile(ID,thePixelTemp_);
  ID = fpr2p1;    templ.pushfile(ID,thePixelTemp_);
  ID = fpr1p1;    templ.pushfile(ID,thePixelTemp_);
  ID = fpr1p2;    templ.pushfile(ID,thePixelTemp_);
  ID = fpr2p2;    templ.pushfile(ID,thePixelTemp_);

  // Other inits
  speed = -2;
  ngood = 0; nbad = 0;
    
  float cotaminb = 100.f, cotamaxb = -100.f, cotaminf = 100.f, cotamaxf = -100.f;
  float cotbminb = 100.f, cotbmaxb = -100.f, cotbminf = 100.f, cotbmaxf = -100.f;
  float cotabminf = 100.f, cotabmaxf = 0.f;
    
  for(i=0; i<4; ++i) {
      for(j=0; j<11; ++j) {
          acctrk[i][j] = 0; alltrk[i][j] = 0;
    }
  }
  //Q: acctrk? --> accepted track on a given layer with a given pt cut

  // Set branch adresses in memory outside of this file (PixelTree.cc and PixelTree.h)
  pixelTree *ana=new pixelTree(&chain);

  // Loop over events
  Long64_t nentries=chain.GetEntries(); 
  Long64_t nbytes = 0, nb = 0;

  //nentries = 1;

  cout << "Running over " << nentries << " events" << endl;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
//         if (jentry!=96) continue; 

    // Set current entry. If -2 Entry does not exist
    Long64_t ientry = chain.LoadTree(jentry);
    if (ientry < 0) {
      cout<<"ientry "<<ientry<<"does not exist!"<<endl;
      break;
    }


    // Read all branches of entry. Return total number of bytes read.
//     cout<< "jentry " << jentry<<" and ientry "<<ientry<<" ";
    nb = chain.GetEntry(jentry);  // very very slow!
//     cout << "number of bytes read: " <<nb<<" ";
    nbytes += nb;
    
    int event = ana->event;
    int bx = ana->bx;
    int orbit = ana->orbit;
    int lumiSect = ana->lumiblock; //lumiblock in Morris code
    int run = ana->run;

    if(jentry%1000 == 0) cout<<" jentry: "<< jentry <<" event: "<<event<<" run: "<<run<<" lumiSect: "<<lumiSect<<endl;
    cout<<"jentry: "<<jentry <<" event: "<<event<<" run: "<<run<<" lumiSect: "<<lumiSect<<" bx: "<<bx<<endl;
//     continue;

    const int SELECT_EVENT = 379330015;
//     if(event!=SELECT_EVENT) continue;

    if(SELECT_ORBIT) {
      const int SELECT_MIN=0, SELECT_MAX=99900000;  //18ns, run123596
      if(orbit>=SELECT_MIN && orbit<=SELECT_MAX ) {
        cout<<" found orbit "<<lumiSect<<endl;
      } else {
        continue;
      }
    } // SELECT_ORBIT

    if(SELECT_LUMI) {
      const int SELECT_MIN=0, SELECT_MAX=999999;  // dummy
      //const int SELECT = 0, SELECT_MIN=313, SELECT_MAX=318;  //  132477 wbc156
      if(lumiSect>=SELECT_MIN && lumiSect<=SELECT_MAX ) {
        bool found=true;
	//cout<<" found lumi "<<lumiSect<<endl;
      } else {
	//if(found) break;  // assumes same lumiSect events are together 
	//else continue;
        continue;
      }
    } // SELECT_LUMI
     
    int nclus = ana->ClN;
    int ndigis = ana->DgN;
    int iTkN = ana->TkN;
    int iPvN = ana->PvN;
//     cout<< " ClN " << nclus << " TkN " << iTkN <<" DgN "<<ndigis<<" pvs "<<iPvN<<endl;

    if( nclus>maxClus) {
      cout<<" Too many clus, extend array from "<<maxClus<<" to "<<nclus<<endl;
      maxClus=nclus;
    }

    if( ndigis>maxDigi) {
      cout<<" Too many digis, extend array from "<<maxDigi<<" to "<<ndigis<<endl;
      maxDigi=ndigis;
    }
     
    if (iPvN > PVMAX) {
      cout<<" Too many vtxs, extend array from "<<PVMAX<<" to "<<iPvN<<endl;
      iPvN = PVMAX;
      continue;
    }
     
    if( iTkN>maxTrk) {
      cout<<" Too many tracks, extend array from "<<maxTrk<<" to "<<iTkN<<endl;
       maxTrk=iTkN;
    }

    hp[0]->Fill(ana->PvN);  if(ana->PvN < 1) continue;  vector<bool> vtxok(iPvN,false);
    for(i = 0; i<iPvN; ++i) {
      hp[1]->Fill(ana->PvX[i]);
      hp[2]->Fill(ana->PvY[i]);
      hp[3]->Fill(ana->PvZ[i]);
      hp[23]->Fill(ana->PvNdof[i]);
      hp[24]->Fill((ana->PvChi2[i]/ana->PvNdof[i]));

      if(ana->PvIsFake[i] == 1) continue;
      if(fabsf(ana->PvX[i]) > 0.3) continue;
      if(fabsf(ana->PvY[i]) > 0.3) continue;
      if(fabsf(ana->PvZ[i]) > 6.) continue;
      if(ana->PvNdof[i] < 5.) continue;
      vtxok[i] = true;	
    }
    
    hp[29]->Fill(iTkN);

//     if(pvs>0) {
//        if(pvs>maxVtx) {
// 	 
// 	 pvs=maxVtx;
// 	 continue; // skip event 
//        }
//        for(int n=0; n<pvs;++n) {
//          //cout<<ana->PvIsFake[n]<<" "<<ana->PvChi2[n]<<" "<<ana->PvX[n]<<" "<<ana->PvY[n]
//          //  <<" "<<ana->PvZ[n]<<endl;
//          if(ana->PvIsFake[n]==1) continue; // skip fakes
//          float pvz = ana->PvZ[n];
//          hpvz->Fill(pvz);
//          if(pvz>-20. && pvz<20.) {
//            float pvr = sqrt(ana->PvX[n]*ana->PvX[n] + ana->PvY[n]*ana->PvY[n]); 
//            hpvr->Fill(pvr);
//            if(pvr<0.6) pvsTrue++;
//          }
//        }
//      }

//      if(PRINT) {
//        cout << "new event: " << jentry << "/  "<<event <<" "<<nb<<" "<<nbytes<<" "
// 	    <<bx<<" "<<orbit<<" "<<lumiSect<<endl;
//        //cout<<hex<<l1t<<" "<<l1ta0<<" "<<l1ta1<<" "<<l1ta2<<" "<<l1ta3<<dec<<endl;
//        cout<< " ClN: " << nclus << " TkN: " << ana->TkN <<" DgN "<<ndigis<<endl;
//      }
//      int maxC = min(nclus,maxClus);

// Loop over all tracks in the event
//iTkN = 10;      
    for(int iTk = 0; iTk<iTkN; ++iTk) {
      // Require at least 3 pixel clusters on the track   
      int iTkClN = ana->TkClN[iTk];
      int TkClN = iTkClN; if(TkClN < 2) continue;
        
      // Require at least one to be in barrel layer 1 (to suppress conversions/interactions)
      bool layer1 = false; 
        
      for(int iTkCl = 0; iTkCl < TkClN; ++iTkCl) {
        // iCl is the index in the cluster arrays        
        int iCl = ana->TkClI[iTk][iTkCl];
          
        // See if any cluster in the first barrel layer    
// printf("track cluster %d in layer %d, module %d, disk %d \n", iTkCl, fClLayer[iCl], fClModule[iCl], fClDisk[iCl]);
          
        if(ana->ClLayer[iCl] == 1) layer1 = true;
          
      }
        
      if(!layer1) continue;
        
      // Add track selection requirements
      float fTkChi2 = ana->TkChi2[iTk];
      int iTkNdof = ana->TkNdof[iTk];
      hp[4]->Fill((fTkChi2/iTkNdof));
      hp[21]->Fill(iTkNdof);

      double dvx = (double)(ana->TkVx[iTk]-(ana->PvX[0]));
      double dvy = (double)(ana->TkVy[iTk]-(ana->PvY[0]));
      double dvz = (double)(ana->TkVz[iTk]-(ana->PvZ[0]));
      double dv3 = sqrt(dvx*dvx+dvy*dvy+dvz*dvz);
      bool vtxgood = vtxok[0];
      if(ana->PvN > 1) {
        for(i=1; i<iPvN; ++i) {
          double dx = (double)(ana->TkVx[iTk]-ana->PvX[i]);
          double dy = (double)(ana->TkVy[iTk]-ana->PvY[i]);
          double dz = (double)(ana->TkVz[iTk]-ana->PvZ[i]);
          double d3 = sqrt(dx*dx+dy*dy+dz*dz);
          if(d3 >= dv3) continue;
          dvx = dx;
          dvy = dy;
          dvz = dz;
          dv3 = d3;
          vtxgood = vtxok[i];
        }
      }
      if(!vtxgood) continue;              
      hp[25]->Fill(dvx);
      hp[26]->Fill(dvy);
      hp[27]->Fill(dvz);
      hp[28]->Fill(ana->TkQuality[iTk]);
      
      // Make tight vertex cuts to suppress conversions and secondard interactions
        
      if(fabs(dvx) > 0.10) continue;
      if(fabs(dvy) > 0.10) continue;
      if(fabs(dvz) > 0.20) continue;
      hp[19]->Fill(ana->TkPt[iTk]);
      TkP = (double)ana->TkPt[iTk]/sin((double)ana->TkTheta[iTk]);
      hp[20]->Fill(TkP);
      if(ana->TkChi2[iTk]/ana->TkNdof[iTk] > 2.) continue;
      if(lowpt != 1) {
//         int high_purity = 0x4;
//         if((fTkQuality[iTk] & high_purity) == 0x0) continue;
        if(ana->TkNdof[iTk] < 5.) continue;
        if(ana->TkPt[iTk] < 0.1) continue;
      } else {
         if(ana->TkNdof[iTk] < 1.) continue;                    
         if(ana->TkPt[iTk] < 0.05) continue;
      }
        
      // Calculate the track momentum      
      if(TkP < tkpcut) continue;
                
      // p,pt after cuts
      hp[34]->Fill(ana->TkPt[iTk]);
      hp[35]->Fill(TkP);
        
      int iTkPt = (int)(ana->TkPt[iTk]/0.1f);
      if(iTkPt > 10) iTkPt = 10;
        
      // Add initializations here
      bool large_pix = false;
      probQonTrackWMulti = 1 ;
        
      // Loop over pixel clusters
      int nprobQOnTrack = 0;
      for(int iTkCl = 0; iTkCl < TkClN; ++iTkCl) {
          
        cotalpha = TMath::Tan(TMath::Pi()/2. - ana->TkAlpha[iTk][iTkCl]); 
        cotbeta = TMath::Tan(TMath::Pi()/2. - ana->TkBeta[iTk][iTkCl]);
        if(ana->ClDisk[ana->TkClI[iTk][iTkCl]] < -10) {
          if(cotalpha < cotaminb) cotaminb = cotalpha;
          if(cotalpha > cotamaxb) cotamaxb = cotalpha;
          if(cotbeta < cotbminb) cotbminb = cotbeta;
          if(cotbeta > cotbmaxb) cotbmaxb = cotbeta;
          hp[36]->Fill(cotalpha);
          hp[37]->Fill(cotbeta);
        } else {
          if(cotalpha < cotaminf) cotaminf = cotalpha;
          if(cotalpha > cotamaxf) cotamaxf = cotalpha;
          if(cotbeta < cotbminf) cotbminf = cotbeta;
          if(cotbeta > cotbmaxf) cotbmaxf = cotbeta;
          if(fabsf(cotbeta) < cotabminf) cotabminf = fabsf(cotbeta);
          if(fabsf(cotbeta) > cotabmaxf) cotabmaxf = fabsf(cotbeta);
          hp[38]->Fill(cotalpha);
          hp[39]->Fill(fabs(cotbeta));
        }
          
        //  If track angle is larger than the nominal acceptance, skip it        
        if(fabsf(cotbeta) > 6.0) continue;
          
          // iCl is the index in the cluster arrays        
        int iCl = ana->TkClI[iTk][iTkCl];
        sizex = ana->ClSizeX[iCl];
        sizey = ana->ClSizeY[iCl];
        if(sizex < 0 || sizey < 0) continue;
           
        //Make sure that the cluster can fit into the container        
        if(sizex > (TXSIZE-2)) continue;
        if(sizey > (TYSIZE-2)) continue;
          
        // Make sure that there aren't bad pixels
        if(ana->ClRhHasBadPixels[iCl] != 0) continue;
        
        //ClCharge is given in ke
        qclust = 1000.*ana->ClCharge[iCl];
          
        // Decide if this is a BPix or an FPix cluster
        bpix = true;
        if(ana->ClDisk[iCl] < -10) {
          qscale = qscaleB;
          layer = ana->ClLayer[iCl];
          module = ana->ClModule[iCl];
          ladder = ana->ClLadder[iCl];
          if(ladder < 0) {
            offladder = qlad[layer-1] + abs(ladder);
          } else {
            int qladp1 = qlad[layer-1] + 1;
            if(ladder < qladp1) {
              offladder = qladp1 - ladder;
            } else {
              offladder = tladp1[layer-1] - ladder;
            }
          }

 //       printf("layer/module/ladder/offladder = %d/%d/%d/%d \n", layer, module, ladder, offladder);
          if(module < 0) {ring = module+4;} else {ring = module+3;}
          ID = bptmp[layer-1][ring];
//        printf("layer/module/ring/ladder/ID = %d/%d/%d/%d/%d \n", layer, module, ring, ladder, ID);
          llayer = layer;
          if(newmodule[layer-1][offladder-1][ring]){
            ID = bpnew;
            llayer = layer+4;
          }

       } else {
                  bpix = false;
                  qscale = qscaleF;
//                  printf("disk/ring/panel/blade = %d/%d/%d/%d \n", fClDisk[iCl], fClRing[iCl], fClPanel[iCl], fClBlade[iCl]);
                  disk = ana->ClDisk[iCl];
                  if(disk < 0) {side = 1;} else {side = 2;}
                  disk = abs(disk);
                  ring = abs(ana->ClRing[iCl]);
                  panel = abs(ana->ClPanel[iCl]);
                  onblade = ana->ClBlade[iCl];
                  if(ring == 2) {
                     if(panel == 1) {layer = 1;} else {layer = 4;}
                  } else {
                     if(panel == 1) {layer = 2;} else {layer = 3;}
                  }
                  llayer = layer;
                  ID = fptmp[side-1][disk-1][layer-1];
//  Convert the online blade number to an offline blade number
                  if(ring == 1) {
                     if(onblade < 0) {
                        blade = 6 - onblade;
                     } else {
                        if(onblade < 7) {blade = 7 - onblade;} else {blade = 29 - onblade;}
                     }
                   } else {
                     if(onblade < 0) {
                        blade = 31 - onblade;
             } else {
               if(onblade < 10) {blade = 32 - onblade;} else {blade = 66 - onblade;}
             }
           }

                  
//                  cout << "ID/fClDisk/fClPanel/fClPlaquette/layer = " << ID << "/" << fClDisk[iCl] << "/" << fClPanel[iCl] << "/" << fClPlaquette[iCl] << "/" << layer << endl;
//                  cout << "side/disk/fClBlade/ring/blade = " << side << "/" << disk << "/" << fClBlade[iCl] << "/" << ring << "/" << blade << endl << endl;
         }

                    
          // Look for clusters that traverse minimal material (layer 1, flipped modules)        
          //        if(ana->ClLayer[iCl] != 1 || fClFlipped[iCl] != 1) continue;
          
          for(j=0; j<TXSIZE; ++j) {for(i=0; i<TYSIZE; ++i) {cluster[j][i] = 0.;} }        
          for(j=0; j<TXSIZE; ++j) {xdouble[j] = false;}  
          for(i=0; i<TYSIZE; ++i) {ydouble[i] = false;}
          large_pix = false;
          
          // Loop over pixels in cluster
          
          int ClDgN = ana->ClDgN[iCl];
          if(ClDgN <= 0) continue;
          int minPixelRow = 161;
          int minPixelCol = 417;
          for(int iClDg = 0; iClDg < ClDgN; ++iClDg) {
            
            // iDg is the index in the digi arrays        
            int iDg = ana->ClDgI[iCl][iClDg];
            int row = ana->DgRow[iDg];
            int col = ana->DgCol[iDg];
            
            //  Flag edge or double pixels           
            if(row == 0 || row == 79 || row == 80 || row == 159) large_pix = true;
            if(col % 52 == 0 || col % 52 == 51) large_pix = true;
            
            //  Find lower left corner pixel and its coordinates
            if(row < minPixelRow) {
              minPixelRow = row; 
              xoff = 10000.*ana->DgLx[iDg] + xpitch/2.;
              if(row == 79 || row == 80) xoff += xpitch/2.;
            }
            if(col < minPixelCol) {
              minPixelCol = col;
              yoff = 10000.*ana->DgLy[iDg] + ypitch/2.;
              if(col % 52 == 0 || col % 52 == 51) yoff += ypitch/2.;
            }
          }
          
          // Discard clusters with double-size or edge pixels from this analysis
          if(large_pix) continue;
          
          // Now fill the cluster buffer with charges
          for(int iClDg = 0; iClDg < ClDgN; ++iClDg) {
            
            // iDg is the index in the digi arrays        
            int iDg = ana->ClDgI[iCl][iClDg];
            int row = ana->DgRow[iDg];
            int col = ana->DgCol[iDg];
            ix = row - minPixelRow;
            if(ix >= TXSIZE) continue;
            iy = col - minPixelCol;
            if(iy >= TYSIZE) continue;
            
            // Set pixel charges          
            cluster[ix][iy] = qscale*1000.*ana->DgCharge[iDg];
            
            // and double pixel flags          
            if (row == 79 || row == 80){
              xdouble[ix] = true;
            }
            if (col % 52 == 0 || col % 52 == 51 ){
              ydouble[iy] = true;
            }
          }

          xhit = ana->ClRhLx[iCl]*10000.;
          yhit = ana->ClRhLy[iCl]*10000.;
          xhmod = (double)xhit - floor((double)xhit/8100.)*8100.;
          yhmod = (double)yhit - floor((double)yhit/8100.)*8100.;
          
// remove hits with reconstructed positions near edges to avoid bias
          
          if((xhmod < 400.) || (xhmod > 7700.)) continue;
          if(yhmod < 600. || yhmod > 7500.) continue;
#if HSCPONLY == 1
          // check if there are HSCPs on the track, exit if not
          if (ana->ClSimHitPID[iCl][0]<1000000||ana->ClSimHitPID[iCl][0]>10000000) continue;
#else
          // no cut if the mixure/normal track are studied
#endif
          
          
          // Do the template analysis on the cluster 
          SiPixelTemplateReco::ClusMatrix clusterPayload{&cluster[0][0], xdouble, ydouble, mrow,mcol};
          locBx = 1.;
          if(cotbeta < 0.) locBx = -1.;
          locBz = locBx;
          if(cotalpha < 0.) locBz = -locBx;
          // qbin 0-4 bins in charge, speed 
          ierr = PixelTempReco1D(ID, cotalpha, cotbeta, locBz, locBx, clusterPayload, templ, yrec, sigmay, proby, xrec, sigmax, probx, qbin, speed, probQ);  
          if(ierr != 0) {
            ++nbad;
            if(nbad < 50) {printf("ID %d reco of cotalpha/cotbeta = %f/%f failed with error %d \n", ID, cotalpha, cotbeta, ierr);}
          } else {
            xtemp = xoff + xrec;
            ytemp = yoff + yrec;
            dx = xtemp - xhit;
            dy = ytemp - yhit;
            
// Check resolution and weights 
//        if(qbin > 3) {printf(" qbin = %d \n", qbin);}
            if (debug>2) cout << "probQ in cluster " << iTkCl << " is " << probQ <<endl;
            proba = probx*proby;
            probXY = proba*(1.-log(proba));
            proba = probx*proby*probQ;
            probXYQ = proba*(1.-log(proba)+0.5*log(proba)*log(proba));
            log10probQ = log10((double)probQ);
            log10probXY = log10(probXY);
            log10probXYQ = log10(probXYQ);
            hp[5]->Fill(probXY);
            hp[6]->Fill(probQ);
            
            
            hp[7]->Fill(probXYQ);
            hp[8]->Fill(log10probXY);
            hp[9]->Fill(log10probQ);
            hp[10]->Fill(log10probXYQ);
            if(probQ > 0.01) hp[11]->Fill(log10probXY);
            qnorm = qclust/sqrt((double)(1.+cotbeta*cotbeta+cotalpha*cotalpha)); 
            if(bpix) {hp[12]->Fill(qnorm);} else {hp[13]->Fill(qnorm);};
            if(probXY > 0.01) {hp[14]->Fill(qclust);}
            if(probQ > 0.01) {hp[15]->Fill(qclust);}
            if(probXYQ > 0.01) {hp[16]->Fill(qclust);}
            if(bpix) {
              hp[17]->Fill(dx);
              hp[18]->Fill(dy);
              int ihist = 42 + (layer-1)*2;
              hp[ihist]->Fill(dx);
              hp[ihist+1]->Fill(dy);
              ++alltrk[layer-1][iTkPt];
              if(fabs(cotalpha) <=0.3f) {++acctrk[layer-1][iTkPt];}
            } else {
              hp[40]->Fill(dx);
              hp[41]->Fill(dy);
              int ihist = 50 + (llayer-1)*2;
              hp[ihist]->Fill(dx);
              hp[ihist+1]->Fill(dy);
            }
            if(probXY < pcut) continue;
            if (probQ!=0) nprobQOnTrack++;
            probQonTrackWMulti *= probQ;
//             cout << "probQonTrackWMulti " << probQonTrackWMulti <<endl;
//             if(probQonTrackWMulti > 0.00000001) cout << "log probQonTrackWMulti " << log(probQonTrackWMulti) << endl;
          } // end if there is no error from the template analysis
        } // end of loop over pixel clusters
        if(nprobQOnTrack < 2 || nprobQOnTrack!=TkClN) continue;
	if (debug>2) cout << "nprobQOnTrack: " << nprobQOnTrack << endl;
        if (debug>1) cout << "probQonTrackWMulti outside the cluster loop is " << probQonTrackWMulti <<endl;
        hp[30]->Fill(probQonTrackWMulti);
        logprobQonTrackWMulti = log(probQonTrackWMulti);
        if(probQ > 0.00001) hp[31]->Fill(logprobQonTrackWMulti);
        probQonTrackTerm = 0;
        for(int iTkCl = 0; iTkCl < nprobQOnTrack; ++iTkCl) {
            probQonTrackTerm += ((pow(-logprobQonTrackWMulti,iTkCl))/(factorial(iTkCl)));
//             cout << "For cluster " << iTkCl << " pow(-logprobQonTrackWMulti,iTkCl) is " << pow(-logprobQonTrackWMulti,iTkCl) << endl;
            if (debug>1) cout << "For cluster " << iTkCl << " probQonTrackTerm is " << probQonTrackTerm << endl;
        }
          
        probQonTrack = probQonTrackWMulti*probQonTrackTerm;
        if (debug>0) cout << " probQonTrack " << probQonTrack << endl;
        if (probQonTrack>0.98) cout << "probQonTrack is more then 0.98 in event " << event << endl;;
        hp[32]->Fill(probQonTrack);
      } // end of loop over tracks
   } //end cycle for events

  cout << "max number of track/cluster/digis = " << maxTrk << "/" << maxClus << "/" << maxDigi << endl;
  cout << "number of all track/on-track-cluster/digis = " << ana->TkN << "/" << ana->ClN << "/" << ana->DgN << endl;
  cout << "minimum LumiBlock = " << lumiMin << ", maximum lumiBlock = " << lumiMax << endl;
  cout << " number of trigger events = " << ngood << ", template failures = " << nbad << endl;
  cout << " BPix: min/max cot(alpha) = " << cotaminb << "/" << cotamaxb << ", min/max cot(beta) = " << cotbminb << "/" << cotbmaxb << endl;
  cout << " FPix: min/max cot(alpha) = " << cotaminf << "/" << cotamaxf << ", min/max cot(beta) = " << cotbminf << "/" << cotbmaxf << endl;
  cout << " FPix: min/max abs(cot(beta)) = " << cotabminf << "/" << cotabmaxf << endl;
   
   for(i=0; i<4; ++i) {
      // Q: what is frac exactly? 
      float frac[11];
      for(j=0; j<11; ++j) {
         frac[j] = 0.f;
         if(alltrk[i][j] > 10) frac[j] = ((float)acctrk[i][j])/((float)alltrk[i][j]);
       }
       cout << "layer " << i << ": " << frac[0] << " " << frac[1] << " " << frac[2] << " " 
       << frac[3] << " " << frac[4] << " " << frac[5] << " " << frac[6] << " " << frac[7] << " "
       << frac[8] << " " << frac[9] << " " << frac[10] << endl;
    }
    
 
/*
 * Histograms plotting
 */
//   for(i=0; i<5; ++i) {hp[i]->Fit("gaus"); hp[i+10]->Fit("gaus");}
//    for(i=42; i<58; ++i) {hp[i]->Fit("gaus");} // this was the final on Oct 24, 2019

  
//  Create an output filename for this run 
#if HSCPONLY == 0
   sprintf(outfile0,"1DReco_NormalTracks.pdf[");
   sprintf(outfile1,"1DReco_NormalTracks.pdf");
   sprintf(outfile2,"1DReco_NormalTracks.pdf]");
#elif HSCPONLY == 1
   sprintf(outfile0,"1DReco_SignalOnly.pdf[");
   sprintf(outfile1,"1DReco_SignalOnly.pdf");
   sprintf(outfile2,"1DReco_SignalOnly.pdf]");
#endif
   TCanvas c1("c1", header);
   c1.SetFillStyle(4000);
   c1.Print(outfile0);
   for(i=0; i<58; ++i) {
     hp[i]->Draw();
     c1.Print(outfile1);
   }
   c1.Print(outfile2);
  int limit = nentries/100;

  
   hFile.Write();
  //gDirectory->pwd();
  //hFile.ls();
   hFile.Close();

   return 0;
} // end of main()

int factorial(int n) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
