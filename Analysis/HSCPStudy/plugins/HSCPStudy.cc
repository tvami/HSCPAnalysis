 // -*- C++ -*-
//
// Package:    Analysis/HSCPStudy
// Class:      HSCPStudy
//
/**\class HSCPStudy HSCPStudy.cc Analysis/HSCPStudy/plugins/HSCPStudy.cc

    Description: Code for reading external trees (like PixelTrees) and run template reco on the clusters. Then take some low level quantities to look for Heavy Stable Charged Particles (HSCP). 

    Implementation:
    Original code taken from Morris Swartz
    Implementing as EDAnalyzer by Tamas Vami (Tav)
*/
//
// Original Author:  Tamas Almos Vami
//         Created:  Mon, 10 Feb 2020 21:06:13 GMT
//
//


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

#include "pixelTree.h"

#include <cmath>
#include <cstdlib>

#include "TROOT.h"
#include "Math/Vavilov.h"
#include "Math/VavilovAccurate.h"
#include "Math/SpecFuncMathCore.h"
#include "Math/SpecFuncMathMore.h"

#include "TStyle.h"
#include "TH1.h"
#include "TF1.h"
#include "TH2.h" 
#include "TProfile.h" 
#include "TVector3.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TPostScript.h"

// user include files
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include <FWCore/Framework/interface/EventSetup.h>
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// // Geometry services
// #include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
// #include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"

#include "RecoLocalTracker/SiPixelRecHits/src/PixelCPEBase.cc"
#include "RecoLocalTracker/SiPixelRecHits/interface/PixelCPEBase.h"

// #include "CondFormats/SiPixelTransient/src/SiPixelTemplate.cc"
// #include "RecoLocalTracker/SiPixelRecHits/src/SiPixelTemplate.cc"
#include "RecoLocalTracker/SiPixelRecHits/src/SiPixelTemplateReco.cc"

// #include "CondFormats/SiPixelTransient/src/SiPixelTemplate2D.cc"
// #include "RecoLocalTracker/SiPixelRecHits/src/SiPixelTemplate2D.cc"
#include "RecoLocalTracker/SiPixelRecHits/src/SiPixelTemplateReco2D.cc"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CalibTracker/Records/interface/SiPixelTemplateDBObjectESProducerRcd.h"
#include "CalibTracker/Records/interface/SiPixel2DTemplateDBObjectESProducerRcd.h"

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
#define TwoDTempAna
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using namespace std;
using namespace edm;
using namespace ROOT;
using namespace Math;

class HSCPStudy : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit HSCPStudy(const edm::ParameterSet&);
      ~HSCPStudy();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      int             verbosity;
      double          ptCutLabel;
      bool            removePixelLayer1Label;
      std::string     fileName; 
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HSCPStudy::HSCPStudy(const edm::ParameterSet& iConfig):
   verbosity(iConfig.getUntrackedParameter<int>("Verbosity", 0)),
   ptCutLabel(iConfig.getUntrackedParameter<double>("ptCut", 3.0)),
   removePixelLayer1Label(iConfig.getUntrackedParameter<bool>("RemovePixelLayer1", false)),
   fileName(iConfig.getUntrackedParameter<string>("rootFileName", string("pixelTree.root")))
{
}

HSCPStudy::~HSCPStudy()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HSCPStudy::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::cout << "Opening input file " << fileName.c_str() << std::endl;
  TChain chain("pixelTree");
  chain.Add(fileName.c_str());

  int maxClus = 30000;
  int maxDigi = 150000;
  int maxTrk =  7000;
//   int maxVtx =  100;
//   int maxDigiPerClu = 100;
//   
  // Local variables 
  bool ydouble[TYSIZE], xdouble[TXSIZE];
  static float qscale, qscaleB, qscaleF, pcut, tkpcut, probQ, /*xs, ys,*/ probQonTrack, probQonTrackTerm, probXYonTrack, probXYonTrackTerm, dEdxEstimator;
  static float probQonTrackWMulti, probXYonTrackWMulti, corrFactor;
  static float xhit, yhit, xrec, yrec, sigmax, sigmay, probx, proby, cotalpha, cotbeta, locBx, locBz, xoff, yoff, xtemp, ytemp;  
  static int /*sfile, nrun, external,*/ size, sizex, sizey, layer, llayer, /*module, ladder, offladder, side,*/ disk, /*blade, onblade,*/ panel, lowpt;
//   static int tladp1[4], qlad[4]={3, 7, 11, 16};
  static int lumiMin = 100000, lumiMax = 0;
  static vector<int> nbin(5,0);
  int i, j, ierr, ierr2, qbin;
  bool bpix;
   
  double log10probXY, log10probQ, log10probXYQ, logprobQonTrackWMulti, logprobXYonTrackWMulti, qclust, qnorm, qnormcorr, proba, probXY, probXYQ, dx, dy, TkP, xhmod, yhmod;
  static int iy, ix, ngood, nbad, speed, /*IDF1,*/ ring;	
  static char /*infile[150],*/ header[80], /*outfile[150],*/ outfile0[80], outfile1[80], outfile2[80];
  float cluster[TXSIZE][TYSIZE];
  int mrow = TXSIZE, mcol = TYSIZE;
  float xpitch, ypitch, logTkP;
  int factorial(int);
  const std::string currentDate();
    
  bool removePixelLayer1 = removePixelLayer1Label;
  

#ifdef TwoDTempAna
  // For the 2D templates
  int qbin2D;
  static float probQonTrackWMulti2D;
  float cluster2d[TXSIZE][TYSIZE];
  
  static int edgeflagy = 0, edgeflagx = 0;
  static float xrec2D, yrec2D, sigmax2D, sigmay2D, probQ2D, probQonTrack2D, probXYonTrack2D, probXYonTrackWMulti2D, probQonTrackTerm2D, probXYonTrackTerm2D, deltay, probXY2D;
  double /*log10probXY2D, log10probQ2D, log10probXYQ2D,*/ logprobQonTrackWMulti2D, logprobXYonTrackWMulti2D, probXYQ2D;
#endif
  
  int nx=120;  
  gStyle->SetOptFit(101);
  gStyle->SetHistLineWidth(2);
  static vector<TH1F*> hp(85);
  static vector<TH2F*> hp2(5);
  edm::Service<TFileService> fs;
     
  hp[0] = fs->make<TH1F>("numVertex","number of vertices",60,-0.5,59.5);
  hp[1] = fs->make<TH1F>("vertexSizeX","vertex x",20,-1.,1.);
  hp[2] = fs->make<TH1F>("vertexSizeY","vertex y",20,-1.,1.);
  hp[3] = fs->make<TH1F>("corrNormCluChargeBPix","Corrected normalized Cluster Charge (BPix);Cluster charge;Entries (1/bin)",nx,0.,120000.);
  hp[4] = fs->make<TH1F>("corrNormCluChargeFPix","Corrected normalized Cluster Charge (FPix);Cluster charge;Entries (1/bin)",nx,0.,120000.);
  hp[5] = fs->make<TH1F>("probXY","ProbXY;Shape probability;Entries (1/bin)",50,0.,1.);
  hp[6] = fs->make<TH1F>("probQ","ProbQ;Charge probability;Entries (1/bin)",50,0.,1.);
  hp[7] = fs->make<TH1F>("probXYQ","ProbXYQ",50,0.,1.);
  hp[8] = fs->make<TH1F>("log10ProbXY","log10(ProbXY)",nx,-12.,0.);
  hp[9] = fs->make<TH1F>("log10ProbQ","log10(ProbQ)",nx,-12.,0.);
  hp[10] = fs->make<TH1F>("probXYQ","log10(ProbXYQ)",nx,-12.,0.);
  hp[11] = fs->make<TH1F>("log10ProbXY","log10(ProbXY) (ProbQ>0.01)",nx,-12.,0.);
  hp[12] = fs->make<TH1F>("corrNormCluChargeBPix","Normalized Cluster Charge (BPix);Cluster charge;Entries (1/bin)",nx,0.,120000.);
  hp[13] = fs->make<TH1F>("corrNormCluChargeFPix","Normalized Cluster Charge (FPix);Cluster charge;Entries (1/bin)",nx,0.,120000.);
  hp[14] = fs->make<TH1F>("cluChargeProbQHigh","Cluster Charge (ProbQ > 0.98);Cluster charge;Entries (1/bin)",nx,0.,200000.);
  hp[15] = fs->make<TH1F>("cluSizeProbQHigh","Cluster Size  (ProbQ > 0.98);Cluster size;Entries (1/bin)",20,0.,20.);
  hp[16] = fs->make<TH1F>("cluSizeXProbQHigh","Cluster Size X  (ProbQ > 0.98);Cluster size X;Entries (1/bin)",20,0.,20.);
  hp[17] = fs->make<TH1F>("cluSizeYProbQHigh","Cluster Size Y  (ProbQ > 0.98);Cluster size Y;Entries (1/bin)",20,0.,20.);
  hp[18] = fs->make<TH1F>("hNotFilled14","dy_temp (bpix); #Deltay (#mum)",nx,-50.,50.);
  hp[19] = fs->make<TH1F>("pt","Track tranverse momentum; pT (GeV)",nx,0.,12.);
  hp[20] = fs->make<TH1F>("pfull","Track full momentum; p (GeV)",nx,0.,12.);
  hp[21] = fs->make<TH1F>("multiProbQ","ProbQ on tracks (w/ multiplication);Multiplied on-track charge probability;Entries (1/bin)",100,0.,1.);
  hp[22] = fs->make<TH1F>("combProbXY","ProbXY on tracks (w/ combine);Combined on-track shape probability;Entries (1/bin)",100,0.,1.);
  hp[23] = fs->make<TH1F>("combProbQ","ProbQ on tracks (w/ combine);Combined on-track charge probability;Entries (1/bin)",100,0.,1.);
  hp[24] = fs->make<TH1F>("probQ2D","ProbQ2D;2D charge probability;Entries (1/bin)",50,0.,1.);
  hp[25] = fs->make<TH1F>("multiProbQ2D","ProbQ2D on tracks (w/ multiplication);Multiplied on-track 2D charge probability;Entries (1/bin)",100,0.,1.);
  hp[26] = fs->make<TH1F>("logprobQ2D","ln(probQ2D) on tracks (w/ multiplication);Log of the multiplied on-track 2D charge probability;Entries (1/bin)",nx,-12.,0.);
  hp[27] = fs->make<TH1F>("combProbQ2D","ProbQ2D on tracks (w/ combine);Combined on-track 2D charge probability;Entries (1/bin)",100,0.,1.);
  hp[28] = fs->make<TH1F>("probXY2D","ProbXY2D",50,0.,1.);
  hp[29] = fs->make<TH1F>("probXYQ2D","ProbXY2D x ProbQ2D",50,0.,1.);
  
  hp[30] = fs->make<TH1F>("logprobXY2D","ln(probXY2D) on tracks (w/ multiplication);Log of the multiplied on-track 2D shape probability;Entries (1/bin)",nx,-12.,0.);
  hp[31] = fs->make<TH1F>("combProbXY2D","ProbXY2D on tracks (w/ combine);Combined on-track 2D shape probability;Entries (1/bin)",nx,0.,1.);
  hp[32] = fs->make<TH1F>("hNotFilled16","ProbQ on tracks (w/ combine);Combined on-track charge probability;Entries (1/bin)",50,0.,1.);
  hp[33] = fs->make<TH1F>("hTrackQual","Track Quality",17,-1.5,15.5);
  hp[34] = fs->make<TH1F>("hPtAfterCuts","After cuts: Track pT; pT (GeV)",150,0.,1500.);
  hp[35] = fs->make<TH1F>("hNumTracks","Number of tracks",200,-0.5,999.5);
  
  hp[36] = fs->make<TH1F>("hNotFilled20","After cuts: BPix cotalpha; cot(#alpha)",100,-1.0,1.0);
  hp[37] = fs->make<TH1F>("hNotFilled21","After cuts: BPix cotbeta; cot(#beta)",200,-10.,10.);
  hp[38] = fs->make<TH1F>("hNotFilled22","After cuts: FPix cotalpha; cot(#alpha)",50,0.0,1.0);
  hp[39] = fs->make<TH1F>("hNotFilled23","After cuts: FPix |cotbeta|; |cot(#beta)|",50,0.0,1.0);
  
  hp[40] = fs->make<TH1F>("h107","dx_temp (fpix); #Deltax (#mum)",nx,-100.,100.);
  hp[41] = fs->make<TH1F>("h138","dy_temp (fpix); #Deltay (#mum)",nx,-100.,100.);
  hp[42] = fs->make<TH1F>("h801","dx_temp (BP L1); #Deltax (#mum)",nx,-100.,100.);
  hp[43] = fs->make<TH1F>("h802","dy_temp (BP L1); #Deltay (#mum)",nx,-100.,100.);
  hp[44] = fs->make<TH1F>("h803","dx_temp (BP L2); #Deltax (#mum)",nx,-100.,100.);
  hp[45] = fs->make<TH1F>("h804","dy_temp (BP L2); #Deltay (#mum)",nx,-100.,100.);
  hp[46] = fs->make<TH1F>("h805","dx_temp (BP L3); #Deltax (#mum)",nx,-100.,100.);
  hp[47] = fs->make<TH1F>("h806","dy_temp (BP L3); #Deltay (#mum)",nx,-100.,100.);
  hp[48] = fs->make<TH1F>("h807","dx_temp (BP L4); #Deltax (#mum)",nx,-100.,100.);
  hp[49] = fs->make<TH1F>("h808","dy_temp (BP L4); #Deltay (#mum)",nx,-100.,100.);
  hp[50] = fs->make<TH1F>("h809","dx_temp (FP R2P1); #Deltax (#mum)",nx,-100.,100.);
  hp[51] = fs->make<TH1F>("h810","dy_temp (FP R2P1); #Deltay (#mum)",nx,-100.,100.);
  hp[52] = fs->make<TH1F>("h811","dx_temp (FP R1P1); #Deltax (#mum)",nx,-100.,100.);
  hp[53] = fs->make<TH1F>("h812","dy_temp (FP R1P1); #Deltay (#mum)",nx,-100.,100.);
  hp[54] = fs->make<TH1F>("h813","dx_temp (FP R1P2); #Deltax (#mum)",nx,-100.,100.);
  hp[55] = fs->make<TH1F>("h814","dy_temp (FP R1P2); #Deltay (#mum)",nx,-100.,100.);
  hp[56] = fs->make<TH1F>("h815","dx_temp (FP R2P2); #Deltax (#mum)",nx,-100.,100.);
  hp[57] = fs->make<TH1F>("h816","dy_temp (FP R2P2); #Deltay (#mum)",nx,-100.,100.);
  hp[58] = fs->make<TH1F>("h701_n","Normalized Cluster Charge (BPix);Cluster charge;Entries (1/bin)",nx,0.,120000.);
  hp[59] = fs->make<TH1F>("h702_n","Normalized Cluster Charge (FPix);Cluster charge;Entries (1/bin)",nx,0.,120000.);
  hp[60] = fs->make<TH1F>("h701_n1","Normalized Cluster Charge (BPix L1);Cluster charge;Entries (1/bin)",nx,0.,120000.);
  hp[61] = fs->make<TH1F>("h701_n1c","Corr Normalized Cluster Charge (BPix L1);Cluster charge;Entries (1/bin)",nx,0.,120000.);
  hp[62] = fs->make<TH1F>("h701_n2","Normalized Cluster Charge (BPix L2);Cluster charge;Entries (1/bin)",nx,0.,120000.);
  hp[63] = fs->make<TH1F>("h701_n2c","Corr Normalized Cluster Charge (BPix L2);Cluster charge;Entries (1/bin)",nx,0.,120000.);
  hp[64] = fs->make<TH1F>("h701_n3","Normalized Cluster Charge (BPix L3);Cluster charge;Entries (1/bin)",nx,0.,120000.);
  hp[65] = fs->make<TH1F>("h701_n3c","Corr Normalized Cluster Charge (BPix L3);Cluster charge;Entries (1/bin)",nx,0.,120000.);
  hp[66] = fs->make<TH1F>("h701_n4","Normalized Cluster Charge (BPix L4);Cluster charge;Entries (1/bin)",nx,0.,120000.);
  hp[67] = fs->make<TH1F>("h701_n4c","Corr Normalized Cluster Charge (BPix L4);Cluster charge;Entries (1/bin)",nx,0.,120000.);
  hp[68] = fs->make<TH1F>("h702_nr1","Normalized Cluster Charge (FPix R1);Cluster charge;Entries (1/bin)",nx,0.,120000.);
  hp[69] = fs->make<TH1F>("h702_nr1c","Corr Normalized Cluster Charge (FPix R1);Cluster charge;Entries (1/bin)",nx,0.,120000.);
  hp[70] = fs->make<TH1F>("h702_nr2","Normalized Cluster Charge (FPix R2);Cluster charge;Entries (1/bin)",nx,0.,120000.);
  hp[71] = fs->make<TH1F>("h702_nr2c","Corr Normalized Cluster Charge (FPix R2);Cluster charge;Entries (1/bin)",nx,0.,120000.);
  hp[72] = fs->make<TH1F>("hdEdxEstimator","dEdxEstimator; dEdx (MeV/cm);Entries (1/bin)",nx,0.,20.);
  hp[73] = fs->make<TH1F>("clustSize","Cluster size;Cluster size;Entries (1/bin)",20,0.,20.);
  hp[74] = fs->make<TH1F>("clustSizeX","Cluster size in x;Cluster size;Entries (1/bin)",20,0.,20.);
  hp[75] = fs->make<TH1F>("clustSizeY","Cluster size in y;Cluster size;Entries (1/bin)",20,0.,20.);
  hp[76] = fs->make<TH1F>("clustSize","Cluster size (after cuts);Cluster size;Entries (1/bin)",20,0.,20.);
  hp[77] = fs->make<TH1F>("clustSizeX","Cluster size in x (after cuts);Cluster size;Entries (1/bin)",20,0.,20.);
  hp[78] = fs->make<TH1F>("clustSizeY","Cluster size in y (after cuts);Cluster size;Entries (1/bin)",20,0.,20.);

  hp[79] = fs->make<TH1F>("corrFactL1","Correction factor on Layer 1",100,0,5);
  hp[80] = fs->make<TH1F>("corrFactL2","Correction factor on Layer 2",100,0,5);
  hp[81] = fs->make<TH1F>("corrFactL3","Correction factor on Layer 3",100,0,5);
  hp[82] = fs->make<TH1F>("corrFactL4","Correction factor on Layer 4",100,0,5);
  hp[83] = fs->make<TH1F>("corrFactR1","Correction factor on Ring 1",100,0,5);
  hp[84] = fs->make<TH1F>("corrFactR2","Correction factor on Ring 2",100,0,5);

#if HSCPONLY == 1
  hp2[0] = fs->make<TH2F>("h2dEdxEstimatorVsP","dEdx vs track momentum; p (GeV);dEdx (MeV/cm)",350,0.,3500.,100,0.0,6.0);
  hp2[1] = fs->make<TH2F>("h2dEdxEstimatorVsLogP","dEdx vs log track momentum; ln(p);dEdx (MeV/cm)",100,0.,10.,100,0.0,6.0);
#else
  hp2[0] = fs->make<TH2F>("h2dEdxEstimatorVsP","dEdx vs track momentum; p (GeV);dEdx (MeV/cm)",350,0.,3500.,100,0.0,6.0);
  hp2[1] = fs->make<TH2F>("h2dEdxEstimatorVsLogP","dEdx vs log track momentum; ln(p) (GeV);dEdx (MeV/cm)",100,0.,10.,100,0.0,6.0);
#endif
  hp2[2] = fs->make<TH2F>("h2ProbQvsProbXY","ProbQ vs ProbXY; Combined on-track charge probability;Combined on-track shape probability",50,0.,1.,50,0.,1.);
  hp2[3] = fs->make<TH2F>("h2ProbQvsProbQ2D","ProbQ vs ProbQ2D; Combined on-track charge probability;Combined on-track 2D charge probability",50,0.,1.,50,0.,1.);
  hp2[4] = fs->make<TH2F>("h2ProbQ2DvsProbXY2D","ProbQ2D vs ProbXY2D; Combined on-track 2D charge probability;Combined on-track 2D shape probability",50,0.,1.,50,0.,1.);

  // Set style for the the histograms  
  for(unsigned int i=0; i<hp.size(); ++i) {
    hp[i]->SetLineColor(2);
    hp[i]->SetFillColor(38);
     if (i==23 || i == 27) {
      hp[i]->SetMinimum(1);
     } else {
      hp[i]->SetMinimum(0.0);
    }
  }
  
  for(unsigned int i=0; i<hp2.size(); ++i) {
    hp2[i]->SetLineColor(2);
    hp2[i]->SetFillColor(38);
     if (i==0 || i == 27) {
      hp2[i]->SetMinimum(0.0);
     } else {
      hp2[i]->SetMinimum(0.0);
    }
  }
  
  // Template input info
//   printf("enter probability cut, momentum cut, lowpt (1 = yes), bpix scale factor (pixelav/data), fpix scale factor \n");
//   scanf("%f %f %d %f %f", &pcut, &tkpcut, &lowpt, &qscaleB, &qscaleF);
//   pcut = 0.01, tkpcut = 10., lowpt = 1, qscaleB=1., qscaleF=1.;
//   pcut = 0.5, tkpcut = 5., lowpt = 0, qscaleB=1., qscaleF=1.;
//     pcut = 0.01, tkpcut = 3, lowpt = 0, qscaleB=1., qscaleF=1.;
     pcut = 0.01, lowpt = 0, qscaleB=1., qscaleF=1.;
    tkpcut = ptCutLabel;
//

  printf("probability cut = %f, momentum cut = %f, lowpt = %d, bpix scale factor = %f, fpix scale factor = %f \n", pcut, tkpcut, lowpt, qscaleB, qscaleF);
  
//     for(int lay=0; lay<4; ++lay) {
//       tladp1[lay] = 5*qlad[lay]+1;
//   }
  
//    bool newmodule[4][64][8];
//    for(int lay=0; lay<4; ++lay) {
//       for(int lad=0; lad<64; ++lad) {
//          for(int mod=0; mod<8; ++mod) {
//             newmodule[lay][lad][mod] = false;
//          }
//       }
//    }
   
// Initialize 1D templates
  const SiPixelTemplateDBObject* templateDBobject_;
  edm::ESHandle<SiPixelTemplateDBObject> templateDBobject;
  iSetup.get<SiPixelTemplateDBObjectESProducerRcd>().get(templateDBobject);
  templateDBobject_ = templateDBobject.product();
  std::vector< SiPixelTemplateStore > thePixelTemp_;
  SiPixelTemplate templ(thePixelTemp_);
  
  if (verbosity>2) cout << " ---------  PixelCPETemplateReco: Loading templates from database (DB) --------- " << endl;

  if (!SiPixelTemplate::pushfile(*templateDBobject_, thePixelTemp_))
      cout << "\nERROR: Templates not filled correctly. Check the sqlite file. Using SiPixelTemplateDBObject version "
        << (*templateDBobject_).version() << "\n\n";
          
// Initialize 2D templates
  const SiPixel2DTemplateDBObject* templateDBobject2D_;
  edm::ESHandle<SiPixel2DTemplateDBObject> templateDBobject2D;
  iSetup.get<SiPixel2DTemplateDBObjectRcd>().get("numerator", templateDBobject2D);
  templateDBobject2D_ = templateDBobject2D.product();
  std::vector< SiPixelTemplateStore2D > thePixelTemp2D_;
  SiPixelTemplate2D templ2D(thePixelTemp2D_);
  
  if (!SiPixelTemplate2D::pushfile( *templateDBobject2D_, thePixelTemp2D_))
      cout << "\nERROR: 2D Templates not filled correctly. Check the sqlite file. Using SiPixel2DTemplateDBObject version "
        << (*templateDBobject2D_).version() << "\n\n";

  // Other inits
  xpitch = 100;
  ypitch = 150;
  speed = -2;
  ngood = 0; nbad = 0;
    
  float cotaminb = 100.f, cotamaxb = -100.f, cotaminf = 100.f, cotamaxf = -100.f;
  float cotbminb = 100.f, cotbmaxb = -100.f, cotbminf = 100.f, cotbmaxf = -100.f;
  float cotabminf = 100.f, cotabmaxf = 0.f;
  
  // Set branch adresses in memory outside of this file (PixelTree.cc and PixelTree.h)
  pixelTree *ana=new pixelTree(&chain);
  int runNumber;
  // Loop over events
  Long64_t nentries=chain.GetEntries(); 
  Long64_t nbytes = 0, nb = 0;

//   nentries = 1; //DEBUG

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
//     if (event!=223718312) continue; // this event has CombinedProb<0.99
//     int bx = ana->bx;
//     int orbit = ana->orbit;
    int lumiSect = ana->lumiblock; //lumiblock in Morris code
    int run = ana->run;
    runNumber = run;

    if(jentry%1000 == 0) cout<<" jentry: "<< jentry <<" event: "<<event<<" run: "<<run<<" lumiSect: "<<lumiSect<<endl;
	cout<<"jentry: "<<jentry <<" event: "<<event<<" run: "<<run<<" lumiSect: "<<lumiSect<<endl;
     
    int nclus = ana->ClN;
    int ndigis = ana->DgN;
    int iTkN = ana->TkN;
    int iPvN = ana->PvN;
    if (verbosity>1) cout<< " ClN " << nclus << " TkN " << iTkN <<" DgN "<<ndigis<<" pvs "<<iPvN<<endl;

    hp[0]->Fill(ana->PvN);  if(ana->PvN < 1) continue;  vector<bool> vtxok(iPvN,false);
    for(i = 0; i<iPvN; ++i) {
      hp[1]->Fill(ana->PvX[i]);
      hp[2]->Fill(ana->PvY[i]);

      if(ana->PvIsFake[i] == 1) continue;
      if(fabsf(ana->PvX[i]) > 0.3) continue;
      if(fabsf(ana->PvY[i]) > 0.3) continue;
      if(fabsf(ana->PvZ[i]) > 6.) continue;
      if(ana->PvNdof[i] < 5.) continue;
      vtxok[i] = true;	
    }
    
    hp[35]->Fill(iTkN);
       
// Loop over all tracks in the event
//iTkN = 10;      
    if (verbosity>2) cout << "-------------------------------------------------------" << endl << "----------- Loop over tracks -----------" <<endl;  
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
      hp[33]->Fill(ana->TkQuality[iTk]);
      
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
            
        
      // Add initializations here
      bool large_pix = false;
      probQonTrackWMulti = 1;
      probXYonTrackWMulti = 1;
      float squareSumCharge = 0.0;
      int nprobQOnTrack = 0;
      
#ifdef TwoDTempAna
      probQonTrackWMulti2D = 1;
      probXYonTrackWMulti2D = 1;
      int nprobQOnTrack2D = 0;
#endif
      
        
      // Loop over pixel clusters
      if (verbosity>2) cout << "----------- Loop over clusters on track -----------" <<endl;  
      for(int iTkCl = 0; iTkCl < TkClN; ++iTkCl) {
        
        cotalpha = TMath::Tan(TMath::Pi()/2. - ana->TkAlpha[iTk][iTkCl]); 
        cotbeta = TMath::Tan(TMath::Pi()/2. - ana->TkBeta[iTk][iTkCl]);
        if(ana->ClDisk[ana->TkClI[iTk][iTkCl]] < -10) {
          if(cotalpha < cotaminb) cotaminb = cotalpha;
          if(cotalpha > cotamaxb) cotamaxb = cotalpha;
          if(cotbeta < cotbminb) cotbminb = cotbeta;
          if(cotbeta > cotbmaxb) cotbmaxb = cotbeta;
//          hp[36]->Fill(cotalpha);
//          hp[37]->Fill(cotbeta);
        } else {
          if(cotalpha < cotaminf) cotaminf = cotalpha;
          if(cotalpha > cotamaxf) cotamaxf = cotalpha;
          if(cotbeta < cotbminf) cotbminf = cotbeta;
          if(cotbeta > cotbmaxf) cotbmaxf = cotbeta;
          if(fabsf(cotbeta) < cotabminf) cotabminf = fabsf(cotbeta);
          if(fabsf(cotbeta) > cotabmaxf) cotabmaxf = fabsf(cotbeta);
//          hp[38]->Fill(cotalpha);
//          hp[39]->Fill(fabs(cotbeta));
        }
          
        //  If track angle is larger than the nominal acceptance, skip it        
        if(fabsf(cotbeta) > 6.0) continue;
          
          // iCl is the index in the cluster arrays        
        int iCl = ana->TkClI[iTk][iTkCl];
        size = ana->ClSize[iCl];
        sizex = ana->ClSizeX[iCl];
        sizey = ana->ClSizeY[iCl];
           
        //Make sure that the cluster can fit into the container        
        if(sizex > (TXSIZE-2)) continue;
        if(sizey > (TYSIZE-2)) continue;
        hp[73]->Fill(size);
        hp[74]->Fill(sizex);
        hp[75]->Fill(sizey);
          
        // Make sure that there aren't bad pixels
        if(ana->ClRhHasBadPixels[iCl] != 0) continue;
        
        //ClCharge is given in ke
        qclust = 1000.*ana->ClCharge[iCl];
          
        // Decide if this is a BPix or an FPix cluster
        bpix = true;
        if(ana->ClDisk[iCl] < -10) {
          qscale = qscaleB;
          layer = ana->ClLayer[iCl];
          }

        else {
          bpix = false;
          qscale = qscaleF;
//                  printf("disk/ring/panel/blade = %d/%d/%d/%d \n", fClDisk[iCl], fClRing[iCl], fClPanel[iCl], fClBlade[iCl]);
          disk = ana->ClDisk[iCl];
//                   if(disk < 0) {side = 1;} else {side = 2;}
          disk = abs(disk);
          ring = abs(ana->ClRing[iCl]);
          panel = abs(ana->ClPanel[iCl]);
//                   onblade = ana->ClBlade[iCl];
         if(ring == 2) {
           if(panel == 1) {layer = -1;} else {layer = -4;}
         } else {
           if(panel == 1) {layer = -2;} else {layer = -3;}
           }
           llayer = layer;
        }

                    
          // Look for clusters that traverse minimal material (layer 1, flipped modules)
          //        if(ana->ClLayer[iCl] != 1 || fClFlipped[iCl] != 1) continue;
          
          for(j=0; j<TXSIZE; ++j) {for(i=0; i<TYSIZE; ++i) {cluster[j][i] = 0.;} } 
#ifdef TwoDTempAna
          for(j=0; j<TXSIZE; ++j) {for(i=0; i<TYSIZE; ++i) {cluster2d[j][i] = 0.;} }  
#endif
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
#ifdef TwoDTempAna
            cluster2d[ix][iy] = qscale*1000.*ana->DgCharge[iDg];
#endif
            
            // and double pixel flags          
            if (row == 79 || row == 80){ xdouble[ix] = true; }
            if (col % 52 == 0 || col % 52 == 51 ){ ydouble[iy] = true; }
          }

          xhit = ana->ClRhLx[iCl]*10000.;
          yhit = ana->ClRhLy[iCl]*10000.;
          xhmod = (double)xhit - floor((double)xhit/8100.)*8100.;
          yhmod = (double)yhit - floor((double)yhit/8100.)*8100.;
          
          // remove hits with reconstructed positions near edges to avoid bias
          if((xhmod < 400.) || (xhmod > 7700.)) continue;
          if(yhmod < 600. || yhmod > 7500.) continue;

          // check if there are HSCPs on the track, exit if not == check if there are SM particles, exit if yes
//           if ((ana->ClSimHitPID[iCl][0]<1000000||ana->ClSimHitPID[iCl][0]>10000000)) continue;
          int PID0=abs(ana->ClSimHitPID[iCl][0]);
          int PID1=abs(ana->ClSimHitPID[iCl][1]);
          int PID2=abs(ana->ClSimHitPID[iCl][2]);
          int PID3=abs(ana->ClSimHitPID[iCl][3]);
          int PID4=abs(ana->ClSimHitPID[iCl][4]);
          int PID5=abs(ana->ClSimHitPID[iCl][5]);
          int PID6=abs(ana->ClSimHitPID[iCl][6]);
          int PID7=abs(ana->ClSimHitPID[iCl][7]);
          int PID8=abs(ana->ClSimHitPID[iCl][8]);
          int PID9=abs(ana->ClSimHitPID[iCl][9]);
          bool PID0isHSCP = (PID0==17||PID0==1000993||PID0==1009213||PID0==1009313||PID0==1009323||PID0==1009113||PID0==1009223||PID0==1009333||PID0==1091114||PID0==1092114||PID0==1092214||PID0==1092224||PID0==1093114||PID0==1093214||PID0==1093224||PID0==1093314||PID0==1093324||PID0==1093334);
          bool PID1isHSCP = (PID1==17||PID1==1000993||PID1==1009213||PID1==1009313||PID1==1009323||PID1==1009113||PID1==1009223||PID1==1009333||PID1==1091114||PID1==1092114||PID1==1092214||PID1==1092224||PID1==1093114||PID1==1093214||PID1==1093224||PID1==1093314||PID1==1093324||PID1==1093334);
          bool PID2isHSCP = (PID2==17||PID2==1000993||PID2==1009213||PID2==1009313||PID2==1009323||PID2==1009113||PID2==1009223||PID2==1009333||PID2==1091114||PID2==1092114||PID2==1092214||PID2==1092224||PID2==1093114||PID2==1093214||PID2==1093224||PID2==1093314||PID2==1093324||PID2==1093334);
          bool PID3isHSCP = (PID3==17||PID3==1000993||PID3==1009213||PID3==1009313||PID3==1009323||PID3==1009113||PID3==1009223||PID3==1009333||PID3==1091114||PID3==1092114||PID3==1092214||PID3==1092224||PID3==1093114||PID3==1093214||PID3==1093224||PID3==1093314||PID3==1093324||PID3==1093334);
          bool PID4isHSCP = (PID4==17||PID4==1000993||PID4==1009213||PID4==1009313||PID4==1009323||PID4==1009113||PID4==1009223||PID4==1009333||PID4==1091114||PID4==1092114||PID4==1092214||PID4==1092224||PID4==1093114||PID4==1093214||PID4==1093224||PID4==1093314||PID4==1093324||PID4==1093334);
          bool PID5isHSCP = (PID5==17||PID5==1000993||PID5==1009213||PID5==1009313||PID5==1009323||PID5==1009113||PID5==1009223||PID5==1009333||PID5==1091114||PID5==1092114||PID5==1092214||PID5==1092224||PID5==1093114||PID5==1093214||PID5==1093224||PID5==1093314||PID5==1093324||PID5==1093334);
          bool PID6isHSCP = (PID6==17||PID6==1000993||PID6==1009213||PID6==1009313||PID6==1009323||PID6==1009113||PID6==1009223||PID6==1009333||PID6==1091114||PID6==1092114||PID6==1092214||PID6==1092224||PID6==1093114||PID6==1093214||PID6==1093224||PID6==1093314||PID6==1093324||PID6==1093334);
          bool PID7isHSCP = (PID7==17||PID7==1000993||PID7==1009213||PID7==1009313||PID7==1009323||PID7==1009113||PID7==1009223||PID7==1009333||PID7==1091114||PID7==1092114||PID7==1092214||PID7==1092224||PID7==1093114||PID7==1093214||PID7==1093224||PID7==1093314||PID7==1093324||PID7==1093334);
          bool PID8isHSCP = (PID8==17||PID8==1000993||PID8==1009213||PID8==1009313||PID8==1009323||PID8==1009113||PID8==1009223||PID8==1009333||PID8==1091114||PID8==1092114||PID8==1092214||PID8==1092224||PID8==1093114||PID8==1093214||PID8==1093224||PID8==1093314||PID8==1093324||PID8==1093334);
          bool PID9isHSCP = (PID9==17||PID9==1000993||PID9==1009213||PID9==1009313||PID9==1009323||PID9==1009113||PID9==1009223||PID9==1009333||PID9==1091114||PID9==1092114||PID9==1092214||PID9==1092224||PID9==1093114||PID9==1093214||PID9==1093224||PID9==1093314||PID9==1093324||PID9==1093334);
#if HSCPONLY == 1
          if (!(PID0isHSCP|| PID1isHSCP || PID2isHSCP || PID3isHSCP || PID4isHSCP || PID5isHSCP || PID6isHSCP || PID7isHSCP || PID8isHSCP || PID9isHSCP)
              ) continue;
#else
           if ((PID0isHSCP|| PID1isHSCP || PID2isHSCP || PID3isHSCP || PID4isHSCP || PID5isHSCP || PID6isHSCP || PID7isHSCP || PID8isHSCP || PID9isHSCP)
              ) continue;
          // no cut if the mixure/normal track are studied
#endif
          // p,pt after cuts
          hp[34]->Fill(ana->TkPt[iTk]);
//           cout << "ana->TkPt[iTk]: " << ana->TkPt[iTk] << endl;
          
          qnorm = qclust/sqrt((double)(1.+cotbeta*cotbeta+cotalpha*cotalpha));
          if(qnorm < 10000.) continue;
          
          // 1D templat analysis
          SiPixelTemplateReco::ClusMatrix clusterPayload{&cluster[0][0], xdouble, ydouble, mrow,mcol};
          locBx = 1.;
          if(cotbeta < 0.) locBx = -1.;
          locBz = locBx;
          if(cotalpha < 0.) locBz = -locBx;
          
          int TemplID1 = -9999;
          TemplID1 = templateDBobject_->getTemplateID(ana->ClDetId[iCl]);
          templ.interpolate(TemplID1, 0.f, 0.f, 1.f, 1.f);

          // Running the actualy 1D Template Reco
          ierr = PixelTempReco1D(TemplID1,
                                 cotalpha, 
                                 cotbeta, 
                                 locBz, 
                                 locBx, 
                                 clusterPayload, 
                                 templ, 
                                 yrec, 
                                 sigmay, 
                                 proby, 
                                 xrec, 
                                 sigmax, 
                                 probx, 
                                 qbin, 
                                 speed, 
                                 probQ);  
          if(ierr != 0) {
            ++nbad;
            if(nbad < 50) {printf("ID %d reco of cotalpha/cotbeta = %f/%f failed with error %d \n", TemplID1, cotalpha, cotbeta, ierr);}
            continue;
          } else {
//             if (verbosity>2) cout << "----------- 1D template analysis on a cluster -----------" <<endl;  
            xtemp = xoff + xrec;
            ytemp = yoff + yrec;
            dx = xtemp - xhit;
            dy = ytemp - yhit;
            
            // Check resolution and weights 
//        if(qbin > 3) {printf(" qbin = %d \n", qbin);}
            if (verbosity>2) {
                if (ana->ClLayer[iTkCl]!=-9999) {
                cout << "probQ in cluster " << iTkCl <<  " on Layer="  <<  ana->ClLayer[iTkCl] << " is " << probQ << endl;
                } else if (ana->ClDisk[iTkCl]!=-9999) {
                cout << "probQ in cluster " << iTkCl << " on Disk=" << ana->ClDisk[iTkCl] << " is " << probQ << endl;
                } else  {
                 cout << "Something is wrong, the cluster is on Layer=" <<  ana->ClLayer[iTkCl] << " AND on Disk="<< ana->ClDisk[iTkCl] << endl;
                }
            }
            
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


            if(bpix) {
                hp[12]->Fill(qnorm);
                if (layer==1) hp[60]->Fill(qnorm);
                if (layer==2) hp[62]->Fill(qnorm);
                if (layer==3) hp[64]->Fill(qnorm);
                if (layer==4) hp[66]->Fill(qnorm);
            } else {
                hp[13]->Fill(qnorm);
                if (ring==1) hp[68]->Fill(qnorm);
                if (ring==2) hp[70]->Fill(qnorm);
            };
            qnormcorr = (qnorm*templ.qscale())/templ.r_qMeas_qTrue();
            corrFactor = (templ.qscale())/templ.r_qMeas_qTrue();
            if(bpix) {
                hp[3]->Fill(qnormcorr);
                if (layer==1) hp[61]->Fill(qnormcorr);
                if (layer==2) hp[63]->Fill(qnormcorr);
                if (layer==3) hp[65]->Fill(qnormcorr);
                if (layer==4) hp[67]->Fill(qnormcorr);
                
                if (layer==1) hp[79]->Fill(corrFactor);
                if (layer==2) hp[80]->Fill(corrFactor);
                if (layer==3) hp[81]->Fill(corrFactor);
                if (layer==4) hp[82]->Fill(corrFactor);
            } else {
                hp[4]->Fill(qnormcorr);
                if (ring==1) hp[69]->Fill(qnormcorr);
                if (ring==2) hp[71]->Fill(qnormcorr);
                
                if (ring==1) hp[83]->Fill(corrFactor);
                if (ring==2) hp[84]->Fill(corrFactor);
            };
            
            
            squareSumCharge += pow((3.61e-05*qnormcorr),-2); //3.61e-06 multi to have the correct order of magnitudes
            if (verbosity>4) std::cout << "squareSumCharge: " << squareSumCharge << std::endl;

            if(bpix) {
              hp[17]->Fill(dx);
              hp[18]->Fill(dy);
              int ihist = 42 + (layer-1)*2;
              hp[ihist]->Fill(dx);
              hp[ihist+1]->Fill(dy);
            } else {
              hp[40]->Fill(dx);
              hp[41]->Fill(dy);
              int ihist = 50 + (llayer-1)*2;
              hp[ihist]->Fill(dx);
              hp[ihist+1]->Fill(dy);
            }
            if(probXY < pcut) continue; //TODO
            if(qnorm < 10000.) continue;
            hp[76]->Fill(size);
            hp[77]->Fill(sizex);
            hp[78]->Fill(sizey);

//             if (probQ<0.01) cout << "ProbQ < 0.01 for PID0: " <<  PID0 << endl;
//             if (probQ<0.01 && PID0==211) continue; // no pions
//             if (probQ<0.01 && PID0==11) continue; //no electrons
            if(removePixelLayer1 && ana->ClLayer[iTkCl]!=1) {
                if (verbosity>0) cout << "Layer 1 is not counted in 1D reco calc" << endl;
                if (probQ>0.98) {
                    hp[14]->Fill(qnorm);
                    hp[15]->Fill(size);
                    hp[16]->Fill(sizex);
                    hp[17]->Fill(sizey);
                    cout << "ProbQ is " << probQ << " and probXY is " << probXY  << " on Layer " << ana->ClLayer[iTkCl] << " or Disk " << ana->ClDisk[iTkCl]  << endl;
                    cout << " size is " << size << " and sizeX is " << sizex << " and sizeY " << sizey  << endl;
                    cout << " qnorm is " << qnorm << " and cotalpha " << cotalpha << " and cotbeta " << cotbeta << endl;
                    cout << " ----------------------------------------------------------------------------- " << endl;
                }
                if (probQ!=0 && probQ!=1 ) nprobQOnTrack++;
                if (probQ!=0 && probQ!=1) probQonTrackWMulti *= probQ; // \alpha_n in formula
                if (probXY!=0 && probXY!=1) probXYonTrackWMulti *= probXY; // \alpha_n in formula
            } else {
                if (probQ!=0 && probQ!=1) nprobQOnTrack++;
                if (probQ!=0 && probQ!=1 ) probQonTrackWMulti *= probQ; // \alpha_n in formula
                if (probXY!=0 && probXY!=1) probXYonTrackWMulti *= probXY; // \alpha_n in formula
            }
            
//             cout << "probQonTrackWMulti " << probQonTrackWMulti <<endl;
//             if(probQonTrackWMulti > 0.00000001) cout << "log probQonTrackWMulti " << log(probQonTrackWMulti) << endl;
          } // end if there is no error from the template analysis 1D
#ifdef TwoDTempAna
          // 2D template analysis
          int TemplID2 = -9999;
          TemplID2 = templateDBobject2D_->getTemplateID(ana->ClDetId[iCl]);
          templ2D.interpolate(TemplID2, 0.f, 0.f, 1.f, 1.f); // Q: interpolate to what exactly?

          SiPixelTemplateReco2D::ClusMatrix clusterPayload2d{&cluster2d[0][0], xdouble, ydouble, mrow,mcol};
          int npixels;
          ierr2 = PixelTempReco2D(TemplID2,
                                 cotalpha, 
                                 cotbeta, 
                                 locBz, 
                                 locBx, 
                                 edgeflagy, 
                                 edgeflagx, 
                                 clusterPayload2d, 
                                 templ2D, 
                                 yrec2D, 
                                 sigmay2D, 
                                 xrec2D, 
                                 sigmax2D, 
                                 probXY2D, 
                                 probQ2D, 
                                 qbin2D, 
                                 deltay,
                                 npixels); 
          if(ierr2 != 0) {
            ++nbad;
            if(nbad < 50) {printf("2D Template reco with ID %d of cotalpha/cotbeta = %f/%f failed with error %d \n", TemplID2, cotalpha, cotbeta, ierr2);}
          } else {
            if (verbosity>2) cout << "----------- 2D template analysis on a cluster -----------" <<endl;
            //xtemp = xoff + xrec2D;
            //ytemp = yoff + yrec2D;
            //dx = xrec2D - xhit;
            //dy = yrec2D - yhit;

            // Check resolution and weights 
//          if(qbin > 3) {printf(" qbin = %d \n", qbin);}
            if (verbosity>2) {
                if (ana->ClLayer[iTkCl]!=-9999) {
                cout << "probQ in cluster " << iTkCl <<  " on Layer="  <<  ana->ClLayer[iTkCl] << " is " << probQ << endl;
                } else if (ana->ClDisk[iTkCl]!=-9999) {
                cout << "probQ in cluster " << iTkCl << " on Disk=" << ana->ClDisk[iTkCl] << " is " << probQ << endl;
                } else  {
                 cout << "Something is wrong, the cluster is on Layer=" <<  ana->ClLayer[iTkCl] << " AND on Disk="<< ana->ClDisk[iTkCl] << endl;
                }
            }
            proba = probXY2D*probQ2D;
            probXYQ2D = proba*(1.-log(proba)+0.5*log(proba)*log(proba));
//             log10probQ2D = log10((double)probQ2D);
//             log10probXY2D = log10(probXY2D);
//             log10probXYQ2D = log10(probXYQ2D);

            hp[24]->Fill(probQ2D);
            hp[28]->Fill(probXY2D);
            hp[29]->Fill(probXYQ2D);
            
            if(probXY2D < pcut) continue; //TODO
            if(qnorm < 10000.) continue;
              
            if(removePixelLayer1) {
                if (verbosity>0) cout << "Layer 1 is not counted in 2D reco calc" << endl;
                if (probQ2D!=0 && probQ2D!=1 && ana->ClLayer[iTkCl]!=1) nprobQOnTrack2D++;
                if (probQ2D!=0 && probQ2D!=1 && ana->ClLayer[iTkCl]!=1) probQonTrackWMulti2D *= probQ2D; // \alpha_n in formula
                if (probXY2D!=0 && probXY2D!=1 && ana->ClLayer[iTkCl]!=1) probXYonTrackWMulti2D *= probXY2D; // \alpha_n in formula
            } else {
                if (probQ2D!=0 && probQ2D!=1) nprobQOnTrack2D++;
                if (probQ2D!=0 && probQ2D!=1) probQonTrackWMulti2D *= probQ2D;
                if (probXY2D!=0 && probXY2D!=1) probXYonTrackWMulti2D *= probXY2D; // \alpha_n in formula
              }
//             cout << "probQonTrackWMulti2D " << probQonTrackWMulti2D <<endl;
//             if(probQonTrackWMulti2D > 0.00000001) cout << "log probQonTrackWMulti2D " << log(probQonTrackWMulti2D) << endl;
          } // end if there is no error from the template analysis 2D
#endif
        } // end of loop over pixel clusters
        
        if (verbosity>2) cout <<endl << "----------- ON TRACK INFO -----------" <<endl;
//         if (removePixelLayer1) TkClN=TkClN-1;
        if(nprobQOnTrack!=TkClN) {
            if (verbosity>2) cout << "Error: nprobQOnTrack!=TkClN" << endl;
            continue; //TODO
        }
        if(nprobQOnTrack < 2) {
            if (verbosity>2) cout << "Error: nprobQOnTrack==TkClN but nprobQOnTrack < 2" << endl;
            continue;
        }
	    if (verbosity>2) cout << "nprobQOnTrack: " << nprobQOnTrack << endl;
        if (verbosity>1) cout << "probQonTrackWMulti outside the cluster loop is " << probQonTrackWMulti <<endl;
        hp[21]->Fill(probQonTrackWMulti);
        logprobQonTrackWMulti = log(probQonTrackWMulti);
        logprobXYonTrackWMulti = log(probXYonTrackWMulti);
        if(probQ > 0.00001) hp[22]->Fill(logprobQonTrackWMulti);
        probQonTrackTerm = 0;
        probXYonTrackTerm = 0;
        for(int iTkCl = 0; iTkCl < nprobQOnTrack; ++iTkCl) {
            probQonTrackTerm += ((pow(-logprobQonTrackWMulti,iTkCl))/(factorial(iTkCl)));
            probXYonTrackTerm += ((pow(-logprobXYonTrackWMulti,iTkCl))/(factorial(iTkCl)));
//             cout << "For cluster " << iTkCl << " pow(-logprobQonTrackWMulti,iTkCl) is " << pow(-logprobQonTrackWMulti,iTkCl) << endl;
            if (verbosity>1) cout << "For cluster " << iTkCl << " the probQonTrackTerm is " << probQonTrackTerm << " the probXYonTrackTerm is " << probXYonTrackTerm <<  endl;
        }
          
         probQonTrack = probQonTrackWMulti*probQonTrackTerm;
         probXYonTrack = probXYonTrackWMulti*probXYonTrackTerm;
        
        //if(probXYonTrack<0.98) continue;
        dEdxEstimator = pow((squareSumCharge/TkClN),-0.5); 
        if (verbosity>4) std::cout << "TkClN: " << TkClN << " dEdxEstimator: " << dEdxEstimator << std::endl;
        hp[72]->Fill(dEdxEstimator);
        
        logTkP = -9999.9;
        if (TkP!=0) logTkP = log(TkP);
        hp2[0]->Fill(TkP,dEdxEstimator);
        //if (logTkP!=-9999.9) cout << TkP << " and log of it" <<  logTkP << endl;
        if (logTkP!=-9999.9) hp2[1]->Fill(logTkP,dEdxEstimator);

        if (verbosity>0) cout << "   probQonTrack " << probQonTrack << endl;
        if (verbosity>3) if (probQonTrack>0.98) cout << "probQonTrack is more then 0.98 in event " << event << endl;;
        hp[23]->Fill(probQonTrack);
//        hp[32]->Fill(probQonTrack);
        hp[22]->Fill(probXYonTrack);
        hp2[2]->Fill(probQonTrack,probXYonTrack);
        
#ifdef TwoDTempAna
        if (verbosity>2) cout << "nprobQOnTrack2D: " << nprobQOnTrack2D << endl;
        if (verbosity>1) cout << "probQonTrackWMulti2D outside the cluster loop is " << probQonTrackWMulti2D <<endl;
        hp[25]->Fill(probQonTrackWMulti2D);
        logprobQonTrackWMulti2D = log(probQonTrackWMulti2D);
        logprobXYonTrackWMulti2D = log(probXYonTrackWMulti2D);
        if(probQ2D > 0.00001) hp[26]->Fill(logprobQonTrackWMulti2D);
        if(probXY2D > 0.00001) hp[30]->Fill(logprobXYonTrackWMulti2D);
        probQonTrackTerm2D = 0;
        probXYonTrackTerm2D = 0;
        for(int iTkCl = 0; iTkCl < nprobQOnTrack2D; ++iTkCl) {
            probQonTrackTerm2D += ((pow(-logprobQonTrackWMulti2D,iTkCl))/(factorial(iTkCl)));
            probXYonTrackTerm2D += ((pow(-logprobXYonTrackWMulti2D,iTkCl))/(factorial(iTkCl)));
//             cout << "For cluster " << iTkCl << " pow(-logprobQonTrackWMulti2D,iTkCl) is " << pow(-logprobQonTrackWMulti2D,iTkCl) << endl;
            if (verbosity>1) cout << "For cluster " << iTkCl  << " the probQonTrackTerm2D is " << probQonTrackTerm2D << endl;
            if (verbosity>1) cout << "For cluster " << iTkCl  << " the probXYonTrackTerm2D is " << probXYonTrackTerm2D << endl;
        }
          
        probQonTrack2D = probQonTrackWMulti2D*probQonTrackTerm2D;
        probXYonTrack2D  = probXYonTrackWMulti2D*probXYonTrackTerm2D;
        if (verbosity>0) cout << "   probQonTrack2D " << probQonTrack2D << endl;
        if (verbosity>0) cout << "   probXYonTrack2D " << probXYonTrack2D << endl;
        if (verbosity>3) if (probQonTrack2D>0.98) cout << "probQonTrack2D is more then 0.98 in event " << event << endl;;
        hp[27]->Fill(probQonTrack2D);
        hp[31]->Fill(probXYonTrack2D);
        hp2[3]->Fill(probQonTrack,probQonTrack2D);
        hp2[4]->Fill(probQonTrack2D,probXYonTrack2D);
#endif
      } // end of loop over tracks
   } //end cycle for events

  cout << "max number of track/cluster/digis = " << maxTrk << "/" << maxClus << "/" << maxDigi << endl;
  cout << "number of all track/on-track-cluster/digis = " << ana->TkN << "/" << ana->ClN << "/" << ana->DgN << endl;
  cout << "minimum LumiBlock = " << lumiMin << ", maximum lumiBlock = " << lumiMax << endl;
  cout << " number of trigger events = " << ngood << ", template failures = " << nbad << endl;
  cout << " BPix: min/max cot(alpha) = " << cotaminb << "/" << cotamaxb << ", min/max cot(beta) = " << cotbminb << "/" << cotbmaxb << endl;
  cout << " FPix: min/max cot(alpha) = " << cotaminf << "/" << cotamaxf << ", min/max cot(beta) = " << cotbminf << "/" << cotbmaxf << endl;
  cout << " FPix: min/max abs(cot(beta)) = " << cotabminf << "/" << cotabmaxf << endl;
       
/*
 * Histograms fitting
 */
//   for(i=0; i<5; ++i) {hp[i]->Fit("gaus"); hp[i+10]->Fit("gaus");}
//    for(i=42; i<58; ++i) {hp[i]->Fit("gaus");} // this was the final on Oct 24, 2019
// https://github.com/cms-sw/cmssw/blob/47c3d4c0989c1db4e371d8c888e1dbd66b0bc136/RecoLocalTracker/SiPixelRecHits/src/SiPixelTemplateReco.cc#L1188


// There are two Vavilov implementations, one from ROOT the other from vdt

struct Vavilov_Func { 
   Vavilov_Func() {}

   double operator() (const double *x, const double *p) { 
      double kappa = p[0]; 
      double beta2 = p[1];
      // float xvav = (qtotal-mpv)/sigmaQ;
      return p[4]*( pdf.Pdf( (x[0]-p[2])/p[3], kappa,beta2) );
   }

   ROOT::Math::VavilovAccurate pdf; 
};

struct Vavilov_FuncVVIObjF { 
   Vavilov_FuncVVIObjF() {}

   double operator() (const double *x, const double *p) { 
      double kappa = p[0]; 
      double beta2 = p[1];
      
//       templ.vavilov_pars(p[2], p[3]), kappa); // issue: templ is unknown at this point
      VVIObjF vvidist(kappa, beta2, 0);
      // float xvav = (qtotal-mpv)/sigmaQ;
      return (vvidist.fcn((x[0]-p[2])/p[3])*p[4] );
   }
 
};
//    TMinuit minuit(5);
//    minuit.SetFCN(Vavilov_FuncVVIObjF);
   bool fitVVf = true;
if (fitVVf) {
   static vector<TF1*> f(72);
//    static vector<Vavilov_FuncVVIObjF*> VVf(72);
   static vector<Vavilov_Func*> VVf(72);

   for(i=60; i<72; ++i) {
//         VVf[i] = new Vavilov_FuncVVIObjF();
       VVf[i] = new Vavilov_Func();
 
        
        std::string s = "f"+std::to_string(i);
        const char *name = s.c_str();
//         f[i] = new TF1(name,VVf[i], 0,90000,5,"Vavilov_FuncVVIObjF");
        f[i] = new TF1(name,VVf[i], 0,90000,5,"Vavilov_Func");
        f[i]->SetLineColor(416);
        
//         f[i]->SetParameters(kappa,1.,hp[i]->GetMean(),hp[i]->GetRMS(),hp[i]->GetEntries());
        f[i]->SetParameters(
            0.000001*hp[i]->GetMean(),   // p[0]: kappa [0.01:10]
            1.0,                         // p[1]: beta^2
            0.9*hp[i]->GetMean(),        // p[2]: mpv
            0.4*hp[i]->GetRMS(),         // p[3]: sigmaQ
            0.5*hp[i]->GetEntries()      // p[4]: norm
        );


        f[i]->SetParLimits(0, 0.01, 10.0);
        f[i]->SetParLimits(1, 0.9, 1.);
        f[i]->SetParLimits(2, 5000.0, 100000.0);
        f[i]->SetParLimits(3, 0.0, 10000.0);
        f[i]->SetParLimits(4, 0.0, 100000.0);
        hp[i]->Fit(f[i]);
   }
}
  
//  Create an output filename for this run    
#if HSCPONLY == 0
   sprintf(outfile0,"0NormalTracks.pdf[");
   sprintf(outfile1,"0NormalTracks.pdf");
   sprintf(outfile2,"0NormalTracks.pdf]");
#elif HSCPONLY == 1
   sprintf(outfile0,"0SignalOnly.pdf[");
   sprintf(outfile1,"0SignalOnly.pdf");
   sprintf(outfile2,"0SignalOnly.pdf]");
#elif HSCPONLY == 2
   sprintf(outfile0,"0Both.pdf[");
   sprintf(outfile1,"0Both.pdf");
   sprintf(outfile2,"0Both.pdf]");
#endif
   TCanvas c1("c2", header);
   c1.SetFillStyle(4000);
   c1.Print(outfile0);
   
   bool saveToEOS=false;
   
   if (saveToEOS) {
   string webDir ="https://tvami.web.cern.ch/tvami/projects/HSCP/"+currentDate()+"/"+to_string(runNumber)+"/";
   
   string mkDir ="mkdir /eos/user/t/tvami/www/projects/HSCP/"+currentDate();
   system(mkDir.c_str());
   string mkDir2 ="mkdir /eos/user/t/tvami/www/projects/HSCP/"+currentDate()+"/"+to_string(runNumber);
   system(mkDir2.c_str());
   string cpIndex ="cp /eos/user/t/tvami/www/projects/HSCP/index.php /eos/user/t/tvami/www/projects/HSCP/"+currentDate()+"/";
   system(cpIndex.c_str());
   string cpHt ="cp /eos/user/t/tvami/www/projects/HSCP/.htaccess /eos/user/t/tvami/www/projects/HSCP/"+currentDate()+"/";
   system(cpHt.c_str());
   
   string cpIndex2 ="cp /eos/user/t/tvami/www/projects/HSCP/index.php /eos/user/t/tvami/www/projects/HSCP/"+currentDate()+"/"+to_string(runNumber);
   system(cpIndex2.c_str());
   string cpHt2 ="cp /eos/user/t/tvami/www/projects/HSCP/.htaccess /eos/user/t/tvami/www/projects/HSCP/"+currentDate()+"/"+to_string(runNumber);
   system(cpHt2.c_str());
   
   for(unsigned int i=0; i<hp.size(); ++i) {
#if HSCPONLY == 0
     string name = "/eos/user/t/tvami/www/projects/HSCP/"+currentDate()+"/"+to_string(runNumber)+"/NormalTracks_hp"+to_string(i)+".png";
#elif HSCPONLY == 1
     string name = "/eos/user/t/tvami/www/projects/HSCP/"+currentDate()+"/"+to_string(runNumber)+"/SignalOnly_hp"+to_string(i)+".png";
#elif HSCPONLY == 2
     string name = "/eos/user/t/tvami/www/projects/HSCP/"+currentDate()+"/"+to_string(runNumber)+"/Both_hp"+to_string(i)+".png";
#endif
     if (i==23 || i == 27) {
       c1.SetLogy(1);
       hp[i]->Draw();
       hp[i]->Write();
       c1.Print(outfile1);
       c1.SaveAs(name.c_str());
     } else {
       c1.SetLogy(0);
       hp[i]->Draw();
       hp[i]->Write();
       c1.Print(outfile1);
       c1.SaveAs(name.c_str());
     }
   }
   
   for(unsigned int i=0; i<hp2.size(); ++i) {
#if HSCPONLY == 0
     string name = "/eos/user/t/tvami/www/projects/HSCP/"+currentDate()+"/"+to_string(runNumber)+"/NormalTracks_hp_2D_"+to_string(i)+".png";
#elif HSCPONLY == 1
     string name = "/eos/user/t/tvami/www/projects/HSCP/"+currentDate()+"/"+to_string(runNumber)+"/SignalOnly_hp_2D_"+to_string(i)+".png";
#elif HSCPONLY == 2
     string name = "/eos/user/t/tvami/www/projects/HSCP/"+currentDate()+"/"+to_string(runNumber)+"/Both_hp_2D_"+to_string(i)+".png";
#endif
   if (i==0) {
     c1.SetLogy(0);
     c1.SetLogx(1);
     hp2[i]->Draw("COLZ");
     hp2[i]->Write("COLZ");
     c1.Print(outfile1);
     c1.SaveAs(name.c_str());
     } else  if (i==2){
     c1.SetLogz(1);
     hp2[i]->Draw("COLZ");
     hp2[i]->Write("COLZ");
     c1.Print(outfile1);
     c1.SaveAs(name.c_str());
     } else {
     c1.SetLogy(0);
     c1.SetLogx(0);
     hp2[i]->Draw("COLZ");
     hp2[i]->Write("COLZ");
     c1.Print(outfile1);
     c1.SaveAs(name.c_str());
     }
   }
   
   c1.Print(outfile2);
   
   
#if HSCPONLY == 0
     string cpRootFile ="cp HSCPAna_Run"+to_string(runNumber)+".root /eos/user/t/tvami/www/projects/HSCP/"+currentDate()+"/"+to_string(runNumber)+"/HSCPAna_NormalTracks_Run"+to_string(runNumber)+".root";
#elif HSCPONLY == 1
     string cpRootFile ="cp HSCPAna_Run"+to_string(runNumber)+".root /eos/user/t/tvami/www/projects/HSCP/"+currentDate()+"/"+to_string(runNumber)+"/HSCPAna_SignalOnly_Run"+to_string(runNumber)+".root";
#endif
   std::cout << cpRootFile << std::endl;
   std::cout << "Plots can be seen in " << webDir << std::endl;
   }

}


// ------------ method called once each job just before starting event loop  ------------
void
HSCPStudy::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void
HSCPStudy::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HSCPStudy::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

int factorial(int n) {
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

const std::string currentDate() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d", &tstruct);

    return buf;
}




//define this as a plug-in
DEFINE_FWK_MODULE(HSCPStudy);
