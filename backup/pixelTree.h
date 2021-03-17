//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Oct 28 19:28:12 2019 by ROOT version 6.10/09
// from TTree pixelTree/pixelTree
// found on file: pixelTree-SimHitInfo.root
//////////////////////////////////////////////////////////

#ifndef pixelTree_h
#define pixelTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class pixelTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run;
   Int_t           lumiblock;
   Int_t           event;
   Int_t           bx;
   Int_t           orbit;
   Float_t         bz;
   UInt_t          tlo;
   UInt_t          thi;
   Float_t         Lumi;
   Float_t         LumiInt;
   UInt_t          fed1;
   UInt_t          fed2;
   UInt_t          l1t;
   UInt_t          l1ta[4];
   UInt_t          l1tt[4];
   UInt_t          hlta[8];
   Bool_t          ttA[64];
   Bool_t          l1A[128];
   Bool_t          hlA[1024];
   UInt_t          hlt;
   Int_t           PvN;
   Float_t         PvX[100];   //[PvN]
   Float_t         PvY[100];   //[PvN]
   Float_t         PvZ[100];   //[PvN]
   Float_t         PvXe[100];   //[PvN]
   Float_t         PvYe[100];   //[PvN]
   Float_t         PvZe[100];   //[PvN]
   Float_t         PvChi2[100];   //[PvN]
   Float_t         PvNdof[100];   //[PvN]
   Int_t           PvIsFake[100];   //[PvN]
   Int_t           MuN;
   Int_t           MuType[10];   //[MuN]
   Int_t           MuTkI[10];   //[MuN]
   Float_t         MuPt[10];   //[MuN]
   Float_t         MuTheta[10];   //[MuN]
   Float_t         MuPhi[10];   //[MuN]
   Float_t         MuT[10];   //[MuN]
   Float_t         MuTcorr[10];   //[MuN]
   Float_t         MuTerr[10];   //[MuN]
   Float_t         MuTmean;
   Int_t           MuTrigger;
   Float_t         HfEplus;
   Float_t         HfEminus;
   Float_t         HfTplus;
   Float_t         HfTminus;
   Int_t           TkN;
   Int_t           TkQuality[4169];   //[TkN]
   Int_t           TkCharge[4169];   //[TkN]
   Float_t         TkChi2[4169];   //[TkN]
   Float_t         TkNdof[4169];   //[TkN]
   Float_t         TkPt[4169];   //[TkN]
   Float_t         TkTheta[4169];   //[TkN]
   Float_t         TkEta[4169];   //[TkN]
   Float_t         TkPhi[4169];   //[TkN]
   Float_t         TkD0[4169];   //[TkN]
   Float_t         TkDz[4169];   //[TkN]
   Float_t         TkVx[4169];   //[TkN]
   Float_t         TkVy[4169];   //[TkN]
   Float_t         TkVz[4169];   //[TkN]
   Float_t         TkAlpha[4169][20];   //[TkN]
   Float_t         TkBeta[4169][20];   //[TkN]
   Float_t         TkResX[4169][20];   //[TkN]
   Float_t         TkResY[4169][20];   //[TkN]
   Float_t         TkResXe[4169][20];   //[TkN]
   Float_t         TkResYe[4169][20];   //[TkN]
   Float_t         TkRes2X[4169][20];   //[TkN]
   Float_t         TkRes2Xe[4169][20];   //[TkN]
   Int_t           TkClN[4169];   //[TkN]
   Int_t           TkClI[4169][20];   //[TkN]
   Int_t           TkType[4169];   //[TkN]
   Int_t           TkNHits[4169];   //[TkN]
   Int_t           TkLHits[4169];   //[TkN]
   Int_t           TkLHitsI[4169];   //[TkN]
   Int_t           TkLHitsO[4169];   //[TkN]
   Float_t         TkNHitFr[4169];   //[TkN]
   Int_t           TkMuI[4169];   //[TkN]
   Int_t           ClN;
   Float_t         ClRow[31107];   //[ClN]
   Float_t         ClCol[31107];   //[ClN]
   Float_t         ClLx[31107];   //[ClN]
   Float_t         ClLxE[31107];   //[ClN]
   Float_t         ClLy[31107];   //[ClN]
   Float_t         ClLyE[31107];   //[ClN]
   Float_t         ClGx[31107];   //[ClN]
   Float_t         ClGy[31107];   //[ClN]
   Float_t         ClGz[31107];   //[ClN]
   Int_t           ClSize[31107];   //[ClN]
   Int_t           ClSizeX[31107];   //[ClN]
   Int_t           ClSizeY[31107];   //[ClN]
   Int_t           ClFlipped[31107];   //[ClN]
   Int_t           ClLayer[31107];   //[ClN]
   Int_t           ClLadder[31107];   //[ClN]
   Int_t           ClModule[31107];   //[ClN]
   Int_t           ClDisk[31107];   //[ClN]
   Int_t           ClBlade[31107];   //[ClN]
   Int_t           ClPanel[31107];   //[ClN]
   Int_t           ClPlaquette[31107];   //[ClN]
   Int_t           ClRing[31107];   //[ClN]
   Int_t           ClDetId[31107];   //[ClN]
   Int_t           ClFedId[31107];   //[ClN]
   Int_t           ClROC[31107];   //[ClN]
   Int_t           ClChannel[31107];   //[ClN]
   Float_t         ClCharge[31107];   //[ClN]
   Float_t         ClChargeCorr[31107];   //[ClN]
   Int_t           ClType[31107];   //[ClN]
   Float_t         ClRhLx[31107];   //[ClN]
   Float_t         ClRhLxE[31107];   //[ClN]
   Float_t         ClRhLy[31107];   //[ClN]
   Float_t         ClRhLyE[31107];   //[ClN]
   Float_t         ClRhGx[31107];   //[ClN]
   Float_t         ClRhGy[31107];   //[ClN]
   Float_t         ClRhGz[31107];   //[ClN]
   Float_t         ClRhProb[31107];   //[ClN]
   Float_t         ClRhProbX[31107];   //[ClN]
   Float_t         ClRhProbY[31107];   //[ClN]
   UInt_t          ClRhQualWord[31107];   //[ClN]
   Int_t           ClRhqBin[31107];   //[ClN]
   Int_t           ClRhSpansTwoROCs[31107];   //[ClN]
   Int_t           ClRhIsOnEdge[31107];   //[ClN]
   Int_t           ClRhHasBadPixels[31107];   //[ClN]
   Int_t           ClTkN[31107];   //[ClN]
   Int_t           ClTkI[31107][100];   //[ClN]
   Int_t           ClDgN[31107];   //[ClN]
   Int_t           ClDgI[31107][100];   //[ClN]
   Int_t           ClSimHitN[31107];   //[ClN]
   Int_t           ClSimHitPID[31107][10];   //[ClN]
   Int_t           ClSimHitPRC[31107][10];   //[ClN]
   Int_t           ClSimHitTrkID[31107][10];   //[ClN]
   Float_t         ClSimHitLx[31107][10];   //[ClN]
   Float_t         ClSimHitLy[31107][10];   //[ClN]
   Float_t         ClSimHitThe[31107][10];   //[ClN]
   Float_t         ClSimHitPhi[31107][10];   //[ClN]
   Float_t         ClSimHitMom[31107][10];   //[ClN]
   Int_t           ClSimTrN[31107];   //[ClN]
   Int_t           ClSimTrID[31107][10];   //[ClN]
   Float_t         ClSimTrFr[31107][10];   //[ClN]
   Int_t           ClSimTrID2[31107][10];   //[ClN]
   Int_t           ClSimTrType[31107][10];   //[ClN]
   Int_t           ClSimTrQ[31107][10];   //[ClN]
   Float_t         ClSimTrPx[31107][10];   //[ClN]
   Float_t         ClSimTrPy[31107][10];   //[ClN]
   Float_t         ClSimTrPz[31107][10];   //[ClN]
   Float_t         ClSimTrEn[31107][10];   //[ClN]
   Float_t         ClSimTrEta[31107][10];   //[ClN]
   Float_t         ClSimTrPhi[31107][10];   //[ClN]
   Float_t         ClSimTrPt[31107][10];   //[ClN]
   Float_t         ClSimTrVx[31107][10];   //[ClN]
   Float_t         ClSimTrVy[31107][10];   //[ClN]
   Float_t         ClSimTrVz[31107][10];   //[ClN]
   Int_t           DgN;
   Int_t           DgRow[118776];   //[DgN]
   Int_t           DgCol[118776];   //[DgN]
   Int_t           DgDetid[118776];   //[DgN]
   Int_t           DgRoc[118776];   //[DgN]
   Int_t           DgRocR[118776];   //[DgN]
   Int_t           DgRocC[118776];   //[DgN]
   Float_t         DgLx[118776];   //[DgN]
   Float_t         DgLy[118776];   //[DgN]
   Float_t         DgGx[118776];   //[DgN]
   Float_t         DgGy[118776];   //[DgN]
   Float_t         DgGz[118776];   //[DgN]
   Float_t         DgAdc[118776];   //[DgN]
   Float_t         DgCharge[118776];   //[DgN]
   Int_t           DgClI[118776];   //[DgN]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_lumiblock;   //!
   TBranch        *b_event;   //!
   TBranch        *b_bx;   //!
   TBranch        *b_orbit;   //!
   TBranch        *b_bz;   //!
   TBranch        *b_tlo;   //!
   TBranch        *b_thi;   //!
   TBranch        *b_fLumi;   //!
   TBranch        *b_fLumiInt;   //!
   TBranch        *b_fed1;   //!
   TBranch        *b_fed2;   //!
   TBranch        *b_l1t;   //!
   TBranch        *b_l1ta;   //!
   TBranch        *b_l1tt;   //!
   TBranch        *b_hlta;   //!
   TBranch        *b_ttA;   //!
   TBranch        *b_l1A;   //!
   TBranch        *b_hlA;   //!
   TBranch        *b_hlt;   //!
   TBranch        *b_PvN;   //!
   TBranch        *b_PvX;   //!
   TBranch        *b_PvY;   //!
   TBranch        *b_PvZ;   //!
   TBranch        *b_PvXe;   //!
   TBranch        *b_PvYe;   //!
   TBranch        *b_PvZe;   //!
   TBranch        *b_PvChi2;   //!
   TBranch        *b_PvNdof;   //!
   TBranch        *b_PvIsFake;   //!
   TBranch        *b_MuN;   //!
   TBranch        *b_MuType;   //!
   TBranch        *b_MuTkI;   //!
   TBranch        *b_MuPt;   //!
   TBranch        *b_MuTheta;   //!
   TBranch        *b_MuPhi;   //!
   TBranch        *b_MuT;   //!
   TBranch        *b_MuTcorr;   //!
   TBranch        *b_MuTerr;   //!
   TBranch        *b_MuTmean;   //!
   TBranch        *b_MuTrigger;   //!
   TBranch        *b_HfEplus;   //!
   TBranch        *b_HfEminus;   //!
   TBranch        *b_HfTplus;   //!
   TBranch        *b_HfTminus;   //!
   TBranch        *b_TkN;   //!
   TBranch        *b_TkQuality;   //!
   TBranch        *b_TkCharge;   //!
   TBranch        *b_TkChi2;   //!
   TBranch        *b_TkNdof;   //!
   TBranch        *b_TkPt;   //!
   TBranch        *b_TkTheta;   //!
   TBranch        *b_TkEta;   //!
   TBranch        *b_TkPhi;   //!
   TBranch        *b_TkD0;   //!
   TBranch        *b_TkDz;   //!
   TBranch        *b_TkVx;   //!
   TBranch        *b_TkVy;   //!
   TBranch        *b_TkVz;   //!
   TBranch        *b_TkAlpha;   //!
   TBranch        *b_TkBeta;   //!
   TBranch        *b_TkResX;   //!
   TBranch        *b_TkResY;   //!
   TBranch        *b_TkResXe;   //!
   TBranch        *b_TkResYe;   //!
   TBranch        *b_TkRes2X;   //!
   TBranch        *b_TkRes2Xe;   //!
   TBranch        *b_TkClN;   //!
   TBranch        *b_TkClI;   //!
   TBranch        *b_TkType;   //!
   TBranch        *b_TkNHits;   //!
   TBranch        *b_TkLHits;   //!
   TBranch        *b_TkLHitsI;   //!
   TBranch        *b_TkLHitsO;   //!
   TBranch        *b_TkNHitFr;   //!
   TBranch        *b_TkMuI;   //!
   TBranch        *b_ClN;   //!
   TBranch        *b_ClRow;   //!
   TBranch        *b_ClCol;   //!
   TBranch        *b_ClLx;   //!
   TBranch        *b_ClLxE;   //!
   TBranch        *b_ClLy;   //!
   TBranch        *b_ClLyE;   //!
   TBranch        *b_ClGx;   //!
   TBranch        *b_ClGy;   //!
   TBranch        *b_ClGz;   //!
   TBranch        *b_ClSize;   //!
   TBranch        *b_ClSizeX;   //!
   TBranch        *b_ClSizeY;   //!
   TBranch        *b_ClFlipped;   //!
   TBranch        *b_ClLayer;   //!
   TBranch        *b_ClLadder;   //!
   TBranch        *b_ClModule;   //!
   TBranch        *b_ClDisk;   //!
   TBranch        *b_ClBlade;   //!
   TBranch        *b_ClPanel;   //!
   TBranch        *b_ClPlaquette;   //!
   TBranch        *b_ClRing;   //!
   TBranch        *b_ClDetId;   //!
   TBranch        *b_ClFedId;   //!
   TBranch        *b_ClROC;   //!
   TBranch        *b_ClChannel;   //!
   TBranch        *b_ClCharge;   //!
   TBranch        *b_ClChargeCorr;   //!
   TBranch        *b_ClType;   //!
   TBranch        *b_ClRhLx;   //!
   TBranch        *b_ClRhLxE;   //!
   TBranch        *b_ClRhLy;   //!
   TBranch        *b_ClRhLyE;   //!
   TBranch        *b_ClRhGx;   //!
   TBranch        *b_ClRhGy;   //!
   TBranch        *b_ClRhGz;   //!
   TBranch        *b_ClRhProb;   //!
   TBranch        *b_ClRhProbX;   //!
   TBranch        *b_ClRhProbY;   //!
   TBranch        *b_ClRhQualWord;   //!
   TBranch        *b_ClRhqBin;   //!
   TBranch        *b_ClRhSpansTwoROCs;   //!
   TBranch        *b_ClRhIsOnEdge;   //!
   TBranch        *b_ClRhHasBadPixels;   //!
   TBranch        *b_ClTkN;   //!
   TBranch        *b_ClTkI;   //!
   TBranch        *b_ClDgN;   //!
   TBranch        *b_ClDgI;   //!
   TBranch        *b_ClSimHitN;   //!
   TBranch        *b_ClSimHitPID;   //!
   TBranch        *b_ClSimHitPRC;   //!
   TBranch        *b_ClSimHitTrkID;   //!
   TBranch        *b_ClSimHitLx;   //!
   TBranch        *b_ClSimHitLy;   //!
   TBranch        *b_ClSimHitThe;   //!
   TBranch        *b_ClSimHitPhi;   //!
   TBranch        *b_ClSimHitMom;   //!
   TBranch        *b_ClSimTrN;   //!
   TBranch        *b_ClSimTrID;   //!
   TBranch        *b_ClSimTrFr;   //!
   TBranch        *b_ClSimTrID2;   //!
   TBranch        *b_ClSimTrType;   //!
   TBranch        *b_ClSimTrQ;   //!
   TBranch        *b_ClSimTrPx;   //!
   TBranch        *b_ClSimTrPy;   //!
   TBranch        *b_ClSimTrPz;   //!
   TBranch        *b_ClSimTrEn;   //!
   TBranch        *b_ClSimTrEta;   //!
   TBranch        *b_ClSimTrPhi;   //!
   TBranch        *b_ClSimTrPt;   //!
   TBranch        *b_ClSimTrVx;   //!
   TBranch        *b_ClSimTrVy;   //!
   TBranch        *b_ClSimTrVz;   //!
   TBranch        *b_DgN;   //!
   TBranch        *b_DgRow;   //!
   TBranch        *b_DgCol;   //!
   TBranch        *b_DgDetid;   //!
   TBranch        *b_DgRoc;   //!
   TBranch        *b_DgRocR;   //!
   TBranch        *b_DgRocC;   //!
   TBranch        *b_DgLx;   //!
   TBranch        *b_DgLy;   //!
   TBranch        *b_DgGx;   //!
   TBranch        *b_DgGy;   //!
   TBranch        *b_DgGz;   //!
   TBranch        *b_DgAdc;   //!
   TBranch        *b_DgCharge;   //!
   TBranch        *b_DgClI;   //!

   pixelTree(TTree *tree=0);
   virtual ~pixelTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef pixelTree_cxx
pixelTree::pixelTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("pixelTree-SimHitInfo.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("pixelTree-SimHitInfo.root");
      }
      f->GetObject("pixelTree",tree);

   }
   Init(tree);
}

pixelTree::~pixelTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t pixelTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t pixelTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void pixelTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("lumiblock", &lumiblock, &b_lumiblock);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("bx", &bx, &b_bx);
   fChain->SetBranchAddress("orbit", &orbit, &b_orbit);
   fChain->SetBranchAddress("bz", &bz, &b_bz);
   fChain->SetBranchAddress("tlo", &tlo, &b_tlo);
   fChain->SetBranchAddress("thi", &thi, &b_thi);
   fChain->SetBranchAddress("Lumi", &Lumi, &b_fLumi);
   fChain->SetBranchAddress("LumiInt", &LumiInt, &b_fLumiInt);
   fChain->SetBranchAddress("fed1", &fed1, &b_fed1);
   fChain->SetBranchAddress("fed2", &fed2, &b_fed2);
   fChain->SetBranchAddress("l1t", &l1t, &b_l1t);
   fChain->SetBranchAddress("l1ta", l1ta, &b_l1ta);
   fChain->SetBranchAddress("l1tt", l1tt, &b_l1tt);
   fChain->SetBranchAddress("hlta", hlta, &b_hlta);
   fChain->SetBranchAddress("ttA", ttA, &b_ttA);
   fChain->SetBranchAddress("l1A", l1A, &b_l1A);
   fChain->SetBranchAddress("hlA", hlA, &b_hlA);
   fChain->SetBranchAddress("hlt", &hlt, &b_hlt);
   fChain->SetBranchAddress("PvN", &PvN, &b_PvN);
   fChain->SetBranchAddress("PvX", PvX, &b_PvX);
   fChain->SetBranchAddress("PvY", PvY, &b_PvY);
   fChain->SetBranchAddress("PvZ", PvZ, &b_PvZ);
   fChain->SetBranchAddress("PvXe", PvXe, &b_PvXe);
   fChain->SetBranchAddress("PvYe", PvYe, &b_PvYe);
   fChain->SetBranchAddress("PvZe", PvZe, &b_PvZe);
   fChain->SetBranchAddress("PvChi2", PvChi2, &b_PvChi2);
   fChain->SetBranchAddress("PvNdof", PvNdof, &b_PvNdof);
   fChain->SetBranchAddress("PvIsFake", PvIsFake, &b_PvIsFake);
   fChain->SetBranchAddress("MuN", &MuN, &b_MuN);
   fChain->SetBranchAddress("MuType", MuType, &b_MuType);
   fChain->SetBranchAddress("MuTkI", MuTkI, &b_MuTkI);
   fChain->SetBranchAddress("MuPt", MuPt, &b_MuPt);
   fChain->SetBranchAddress("MuTheta", MuTheta, &b_MuTheta);
   fChain->SetBranchAddress("MuPhi", MuPhi, &b_MuPhi);
   fChain->SetBranchAddress("MuT", MuT, &b_MuT);
   fChain->SetBranchAddress("MuTcorr", MuTcorr, &b_MuTcorr);
   fChain->SetBranchAddress("MuTerr", MuTerr, &b_MuTerr);
   fChain->SetBranchAddress("MuTmean", &MuTmean, &b_MuTmean);
   fChain->SetBranchAddress("MuTrigger", &MuTrigger, &b_MuTrigger);
   fChain->SetBranchAddress("HfEplus", &HfEplus, &b_HfEplus);
   fChain->SetBranchAddress("HfEminus", &HfEminus, &b_HfEminus);
   fChain->SetBranchAddress("HfTplus", &HfTplus, &b_HfTplus);
   fChain->SetBranchAddress("HfTminus", &HfTminus, &b_HfTminus);
   fChain->SetBranchAddress("TkN", &TkN, &b_TkN);
   fChain->SetBranchAddress("TkQuality", TkQuality, &b_TkQuality);
   fChain->SetBranchAddress("TkCharge", TkCharge, &b_TkCharge);
   fChain->SetBranchAddress("TkChi2", TkChi2, &b_TkChi2);
   fChain->SetBranchAddress("TkNdof", TkNdof, &b_TkNdof);
   fChain->SetBranchAddress("TkPt", TkPt, &b_TkPt);
   fChain->SetBranchAddress("TkTheta", TkTheta, &b_TkTheta);
   fChain->SetBranchAddress("TkEta", TkEta, &b_TkEta);
   fChain->SetBranchAddress("TkPhi", TkPhi, &b_TkPhi);
   fChain->SetBranchAddress("TkD0", TkD0, &b_TkD0);
   fChain->SetBranchAddress("TkDz", TkDz, &b_TkDz);
   fChain->SetBranchAddress("TkVx", TkVx, &b_TkVx);
   fChain->SetBranchAddress("TkVy", TkVy, &b_TkVy);
   fChain->SetBranchAddress("TkVz", TkVz, &b_TkVz);
   fChain->SetBranchAddress("TkAlpha", TkAlpha, &b_TkAlpha);
   fChain->SetBranchAddress("TkBeta", TkBeta, &b_TkBeta);
   fChain->SetBranchAddress("TkResX", TkResX, &b_TkResX);
   fChain->SetBranchAddress("TkResY", TkResY, &b_TkResY);
   fChain->SetBranchAddress("TkResXe", TkResXe, &b_TkResXe);
   fChain->SetBranchAddress("TkResYe", TkResYe, &b_TkResYe);
   fChain->SetBranchAddress("TkRes2X", TkRes2X, &b_TkRes2X);
   fChain->SetBranchAddress("TkRes2Xe", TkRes2Xe, &b_TkRes2Xe);
   fChain->SetBranchAddress("TkClN", TkClN, &b_TkClN);
   fChain->SetBranchAddress("TkClI", TkClI, &b_TkClI);
   fChain->SetBranchAddress("TkType", TkType, &b_TkType);
   fChain->SetBranchAddress("TkNHits", TkNHits, &b_TkNHits);
   fChain->SetBranchAddress("TkLHits", TkLHits, &b_TkLHits);
   fChain->SetBranchAddress("TkLHitsI", TkLHitsI, &b_TkLHitsI);
   fChain->SetBranchAddress("TkLHitsO", TkLHitsO, &b_TkLHitsO);
   fChain->SetBranchAddress("TkNHitFr", TkNHitFr, &b_TkNHitFr);
   fChain->SetBranchAddress("TkMuI", TkMuI, &b_TkMuI);
   fChain->SetBranchAddress("ClN", &ClN, &b_ClN);
   fChain->SetBranchAddress("ClRow", ClRow, &b_ClRow);
   fChain->SetBranchAddress("ClCol", ClCol, &b_ClCol);
   fChain->SetBranchAddress("ClLx", ClLx, &b_ClLx);
   fChain->SetBranchAddress("ClLxE", ClLxE, &b_ClLxE);
   fChain->SetBranchAddress("ClLy", ClLy, &b_ClLy);
   fChain->SetBranchAddress("ClLyE", ClLyE, &b_ClLyE);
   fChain->SetBranchAddress("ClGx", ClGx, &b_ClGx);
   fChain->SetBranchAddress("ClGy", ClGy, &b_ClGy);
   fChain->SetBranchAddress("ClGz", ClGz, &b_ClGz);
   fChain->SetBranchAddress("ClSize", ClSize, &b_ClSize);
   fChain->SetBranchAddress("ClSizeX", ClSizeX, &b_ClSizeX);
   fChain->SetBranchAddress("ClSizeY", ClSizeY, &b_ClSizeY);
   fChain->SetBranchAddress("ClFlipped", ClFlipped, &b_ClFlipped);
   fChain->SetBranchAddress("ClLayer", ClLayer, &b_ClLayer);
   fChain->SetBranchAddress("ClLadder", ClLadder, &b_ClLadder);
   fChain->SetBranchAddress("ClModule", ClModule, &b_ClModule);
   fChain->SetBranchAddress("ClDisk", ClDisk, &b_ClDisk);
   fChain->SetBranchAddress("ClBlade", ClBlade, &b_ClBlade);
   fChain->SetBranchAddress("ClPanel", ClPanel, &b_ClPanel);
   fChain->SetBranchAddress("ClPlaquette", ClPlaquette, &b_ClPlaquette);
   fChain->SetBranchAddress("ClRing", ClRing, &b_ClRing);
   fChain->SetBranchAddress("ClDetId", ClDetId, &b_ClDetId);
   fChain->SetBranchAddress("ClFedId", ClFedId, &b_ClFedId);
   fChain->SetBranchAddress("ClROC", ClROC, &b_ClROC);
   fChain->SetBranchAddress("ClChannel", ClChannel, &b_ClChannel);
   fChain->SetBranchAddress("ClCharge", ClCharge, &b_ClCharge);
   fChain->SetBranchAddress("ClChargeCorr", ClChargeCorr, &b_ClChargeCorr);
   fChain->SetBranchAddress("ClType", ClType, &b_ClType);
   fChain->SetBranchAddress("ClRhLx", ClRhLx, &b_ClRhLx);
   fChain->SetBranchAddress("ClRhLxE", ClRhLxE, &b_ClRhLxE);
   fChain->SetBranchAddress("ClRhLy", ClRhLy, &b_ClRhLy);
   fChain->SetBranchAddress("ClRhLyE", ClRhLyE, &b_ClRhLyE);
   fChain->SetBranchAddress("ClRhGx", ClRhGx, &b_ClRhGx);
   fChain->SetBranchAddress("ClRhGy", ClRhGy, &b_ClRhGy);
   fChain->SetBranchAddress("ClRhGz", ClRhGz, &b_ClRhGz);
   fChain->SetBranchAddress("ClRhProb", ClRhProb, &b_ClRhProb);
   fChain->SetBranchAddress("ClRhProbX", ClRhProbX, &b_ClRhProbX);
   fChain->SetBranchAddress("ClRhProbY", ClRhProbY, &b_ClRhProbY);
   fChain->SetBranchAddress("ClRhQualWord", ClRhQualWord, &b_ClRhQualWord);
   fChain->SetBranchAddress("ClRhqBin", ClRhqBin, &b_ClRhqBin);
   fChain->SetBranchAddress("ClRhSpansTwoROCs", ClRhSpansTwoROCs, &b_ClRhSpansTwoROCs);
   fChain->SetBranchAddress("ClRhIsOnEdge", ClRhIsOnEdge, &b_ClRhIsOnEdge);
   fChain->SetBranchAddress("ClRhHasBadPixels", ClRhHasBadPixels, &b_ClRhHasBadPixels);
   fChain->SetBranchAddress("ClTkN", ClTkN, &b_ClTkN);
   fChain->SetBranchAddress("ClTkI", ClTkI, &b_ClTkI);
   fChain->SetBranchAddress("ClDgN", ClDgN, &b_ClDgN);
   fChain->SetBranchAddress("ClDgI", ClDgI, &b_ClDgI);
   fChain->SetBranchAddress("ClSimHitN", ClSimHitN, &b_ClSimHitN);
   fChain->SetBranchAddress("ClSimHitPID", ClSimHitPID, &b_ClSimHitPID);
   fChain->SetBranchAddress("ClSimHitPRC", ClSimHitPRC, &b_ClSimHitPRC);
   fChain->SetBranchAddress("ClSimHitTrkID", ClSimHitTrkID, &b_ClSimHitTrkID);
   fChain->SetBranchAddress("ClSimHitLx", ClSimHitLx, &b_ClSimHitLx);
   fChain->SetBranchAddress("ClSimHitLy", ClSimHitLy, &b_ClSimHitLy);
   fChain->SetBranchAddress("ClSimHitThe", ClSimHitThe, &b_ClSimHitThe);
   fChain->SetBranchAddress("ClSimHitPhi", ClSimHitPhi, &b_ClSimHitPhi);
   fChain->SetBranchAddress("ClSimHitMom", ClSimHitMom, &b_ClSimHitMom);
   fChain->SetBranchAddress("ClSimTrN", ClSimTrN, &b_ClSimTrN);
   fChain->SetBranchAddress("ClSimTrID", ClSimTrID, &b_ClSimTrID);
   fChain->SetBranchAddress("ClSimTrFr", ClSimTrFr, &b_ClSimTrFr);
   fChain->SetBranchAddress("ClSimTrID2", ClSimTrID2, &b_ClSimTrID2);
   fChain->SetBranchAddress("ClSimTrType", ClSimTrType, &b_ClSimTrType);
   fChain->SetBranchAddress("ClSimTrQ", ClSimTrQ, &b_ClSimTrQ);
   fChain->SetBranchAddress("ClSimTrPx", ClSimTrPx, &b_ClSimTrPx);
   fChain->SetBranchAddress("ClSimTrPy", ClSimTrPy, &b_ClSimTrPy);
   fChain->SetBranchAddress("ClSimTrPz", ClSimTrPz, &b_ClSimTrPz);
   fChain->SetBranchAddress("ClSimTrEn", ClSimTrEn, &b_ClSimTrEn);
   fChain->SetBranchAddress("ClSimTrEta", ClSimTrEta, &b_ClSimTrEta);
   fChain->SetBranchAddress("ClSimTrPhi", ClSimTrPhi, &b_ClSimTrPhi);
   fChain->SetBranchAddress("ClSimTrPt", ClSimTrPt, &b_ClSimTrPt);
   fChain->SetBranchAddress("ClSimTrVx", ClSimTrVx, &b_ClSimTrVx);
   fChain->SetBranchAddress("ClSimTrVy", ClSimTrVy, &b_ClSimTrVy);
   fChain->SetBranchAddress("ClSimTrVz", ClSimTrVz, &b_ClSimTrVz);
   fChain->SetBranchAddress("DgN", &DgN, &b_DgN);
   fChain->SetBranchAddress("DgRow", DgRow, &b_DgRow);
   fChain->SetBranchAddress("DgCol", DgCol, &b_DgCol);
   fChain->SetBranchAddress("DgDetid", DgDetid, &b_DgDetid);
   fChain->SetBranchAddress("DgRoc", DgRoc, &b_DgRoc);
   fChain->SetBranchAddress("DgRocR", DgRocR, &b_DgRocR);
   fChain->SetBranchAddress("DgRocC", DgRocC, &b_DgRocC);
   fChain->SetBranchAddress("DgLx", DgLx, &b_DgLx);
   fChain->SetBranchAddress("DgLy", DgLy, &b_DgLy);
   fChain->SetBranchAddress("DgGx", DgGx, &b_DgGx);
   fChain->SetBranchAddress("DgGy", DgGy, &b_DgGy);
   fChain->SetBranchAddress("DgGz", DgGz, &b_DgGz);
   fChain->SetBranchAddress("DgAdc", DgAdc, &b_DgAdc);
   fChain->SetBranchAddress("DgCharge", DgCharge, &b_DgCharge);
   fChain->SetBranchAddress("DgClI", DgClI, &b_DgClI);
   Notify();
}

Bool_t pixelTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void pixelTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t pixelTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef pixelTree_cxx
