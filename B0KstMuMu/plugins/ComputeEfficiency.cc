// #########################################################################
// # Program to compute the efficiency for the B0 --> K*0 mu+ mu- analysis #
// # in bins of dimuon q^2, cos(theta_K), cos(theta_l), and phi            #
// #########################################################################
// # Author: Mauro Dinardo                                                 #
// #########################################################################

#include <TROOT.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TF12.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TText.h>
#include <TFitResult.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "Utils.h"
#include "B0KstMuMuSingleCandTreeContent.h"

using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::vector;
using std::ios_base;
using std::make_pair;


// ###########################################
// # How to create the analytical efficiency #
// ###########################################
// (a) Make binned efficiency with the command "Make"

// (b) Fit the theta_l variable with the command "Fit1DEff" for every q2 bin
//     - Set command option to "thetaL"

// (c) Fit the theta_K variable with the command "Fit1DEff"
//     - Set command option to "thetaK"
//     - Set "INPUT_THETAL" to the output.txt file of point (b)

// (d) Fit the phi variable with the command "Fit1DEff" for every q2 bin
//     - Set command option to "phi"

// (e) Fit the theta_l-theta_K variables with the command "FitEff2D"
//     - Set "INPUT_THETAL_THETAK" to the output.txt file of point (c)

// (f) Fit the theta_l-theta_K-phi variables with the command "FitEff3D"
//     - Set "INPUT_THETAL_THETAK" to the output.txt file of point (c)
//     - Set "INPUT_PHI"           to the output.txt file of point (d)

// (g) If you want to look at the final result run the command "Test2DEff"
//     - Copy the output.txt file of point (e) into the parameter file

// (h) If you want to look at the final result run the command "Test3DEff"
//     - Copy the output.txt file of point (f) into the parameter file


// ####################
// # Global constants #
// ####################
#define INPUT_THETAL        "ThetaL_B0ToKstMuMu.txt"
#define INPUT_PHI           "Phi_B0ToKstMuMu.txt"
#define INPUT_THETAL_THETAK "ThetaK_B0ToKstMuMu.txt"

#define RIGHTtag        true
#define SAVEPLOT        false
#define CHECKnegEFF     false
#define EFFis2Dnot3D    true
#define NFILES          100
#define GENEFF          "/efficiency/EffRndGenAnalyFilesSign_JPsi_Psi2S/Efficiency_RndGen.txt"
// "/efficiency/EffRndGenAnalyFilesSign_JPsi_Psi2S/Efficiency_RndGen.txt"
// # OR #
// "/efficiency/EffRndGenBinFilesSign_JPsi_Psi2S/Efficiency_RndGen.txt"
#define SETBATCH        true
#define PARAMETERFILEIN "/python/ParameterFile.txt"
#define ordinateRange   1e-2


// ####################
// # Global variables #
// ####################
Utils* Utility;

TTree* theTreeGenCandidatesNoFilter;
TTree* theTreeRecoCandidates;
TTree* theTreeSingleCand;
B0KstMuMuSingleCandTreeContent* NTupleGenCandidatesNoFilter;
B0KstMuMuSingleCandTreeContent* NTupleRecoCandidates;
B0KstMuMuSingleCandTreeContent* NTupleSingleCand;

vector<TF2*> effFuncs;

vector<double> q2Bins;
vector<double> cosThetaKBins;
vector<double> cosThetaLBins;
vector<double> phiBins;

Utils::effStruct myEff;


// #######################
// # Function Definition #
// #######################
void SetStyle              ();
void ComputeEfficiency     (TTree* theTree, B0KstMuMuSingleCandTreeContent* NTuple, double* Vector, double* VectorErr2Pois, double* VectorErr2Weig, unsigned int type,
			    vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, int SignalType);
void MakeHistogramsAllBins (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, Utils::effStruct myEff, int specBin);
void ReadEfficiencies      (bool isSingleEff, vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			    string fileNameInput, bool isAnalyEff, Utils::effStruct* myEff, bool CheckNegEff, bool savePlot, int specBin, bool EffIs2Dnot3D, int SignalType = 1);
void GenerateEfficiencies  (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, string fileNameInput);
void Fit1DEfficiencies     (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			    Utils::effStruct myEff, string who, unsigned int q2BinIndx, string fileNameOut);
void Fit2DEfficiencies     (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			    int SignalType, Utils::effStruct myEff, unsigned int q2BinIndx, string fileNameOut);
void Fit3DEfficiencies     (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			    int SignalType, Utils::effStruct myEff, unsigned int q2BinIndx, string fileNameOut);
void Test2DEfficiency      (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, Utils::effStruct myEff, unsigned int q2BinIndx, int SignalType, string analyORbin, bool savePlot);
void Test3DEfficiency      (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, Utils::effStruct myEff, unsigned int q2BinIndx, bool savePlot);


// ###########################
// # Function Implementation #
// ###########################
void SetStyle ()
{
  gROOT->SetStyle("Plain");
  gROOT->ForceStyle();
  gStyle->SetTextFont(42);

  gStyle->SetOptFit(1112);
  gStyle->SetOptStat(1110);
  gStyle->SetOptTitle(0);

  gStyle->SetPadRightMargin(0.02);
  gStyle->SetPadTopMargin(0.11);
  gStyle->SetPadBottomMargin(0.12);

  gStyle->SetTitleFont(42,"x");
  gStyle->SetTitleFont(42,"y");
  gStyle->SetTitleOffset(1.05,"x");
  gStyle->SetTitleOffset(0.95,"y");
  gStyle->SetTitleSize(0.05,"x");
  gStyle->SetTitleSize(0.05,"y");

  gStyle->SetLabelFont(42,"x");
  gStyle->SetLabelFont(42,"y");
  gStyle->SetLabelSize(0.05,"x");
  gStyle->SetLabelSize(0.05,"y");

  TGaxis::SetMaxDigits(3);
  gStyle->SetStatY(0.9);
}


void ComputeEfficiency (TTree* theTree, B0KstMuMuSingleCandTreeContent* NTuple, double* Vector, double* VectorErr2Pois, double* VectorErr2Weig, unsigned int type,
			vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, int SignalType)
// ##########################################################
// # Efficiency type = 1 --> total gen events before filter #
// # Efficiency type = 2 --> total gen events after filter  #
// # Efficiency type = 3 --> total reco events              #
// # Efficiency type = 4 --> total single candidate events  #
// ##########################################################
{
  // ###################
  // # Local variables #
  // ###################
  int nEntries;
  Utils::effStruct _Counter;
  double* Counter;
  double mumuq2;
  double cosThetaK;
  double cosThetaMu;
  double phiKstMuMuPlane;
  int mumuq2BinIndx;
  int cosThetaKBinIndx;
  int cosThetaMuBinIndx;
  int phiKstMuMuPlaneBinIndx;
  // ###################


  // ###################################
  // # Initialize efficiency structure #
  // ###################################
  Utility->GenEfficiency(&_Counter,q2Bins,cosThetaKBins,cosThetaLBins,phiBins);
  Utils::effValue myEffVal;
  Utility->ResetEffValue(&myEffVal,0.0);
  Utility->InitEfficiency(myEffVal,_Counter,q2Bins,cosThetaKBins,cosThetaLBins,phiBins);
  Counter = _Counter.Num1;


  NTuple->ClearNTuple();
  NTuple->SetBranchAddresses(theTree);
  nEntries = theTree->GetEntries();
  //nEntries/=100;
  cout << "\n[ComputeEfficiency::ComputeEfficiency]\t@@@ Computing efficiency type " << type;
  if      (type == 1) cout << " (before filter) @@@" << endl;
  else if (type == 2) cout << " (after filter) @@@" << endl;
  else if (type == 3) cout << " (reco events) @@@" << endl;
  else if (type == 4) cout << " (single candidate events) @@@" << endl;
  cout << "[ComputeEfficiency::ComputeEfficiency]\t@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;


  int evPassing =0;
  for (int entry = 0; entry < nEntries; entry++) {
    theTree->GetEntry(entry);

    if ((NTuple->B0pT > Utility->GetSeleCut("B0pT")) && (fabs(NTuple->B0Eta) < Utility->GetSeleCut("B0Eta")) &&

        ((NTuple->genSignal == SignalType || NTuple->genSignal == SignalType+1)) &&

        ((type == 1 || type == 3) ||

         ((type == 2) &&
          (sqrt(NTuple->genMumPx*NTuple->genMumPx + NTuple->genMumPy*NTuple->genMumPy)   > Utility->GetPreCut("MinMupT")) &&
          (sqrt(NTuple->genMupPx*NTuple->genMupPx + NTuple->genMupPy*NTuple->genMupPy)   > Utility->GetPreCut("MinMupT")) &&
          (fabs(Utility->computeEta(NTuple->genMumPx,NTuple->genMumPy,NTuple->genMumPz)) < Utility->GetPreCut("MuEta"))   &&
          (fabs(Utility->computeEta(NTuple->genMupPx,NTuple->genMupPy,NTuple->genMupPz)) < Utility->GetPreCut("MuEta")))  ||

         ((type == 4) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == Utility->RIGHTflavorTAG) &&
          (NTuple->B0MassArb > Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str()))            &&
          (NTuple->B0MassArb < Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str()))           &&

          ((((SignalType == Utility->B0ToKstMuMu)  || (SignalType == Utility->B0ToKstMuMu+1))  && (Utility->PsiRejection(NTuple->B0MassArb,NTuple->mumuMass->at(0),NTuple->mumuMassE->at(0),"rejectPsi",true) == true)) ||
           (((SignalType == Utility->B0ToJPsiKst)  || (SignalType == Utility->B0ToJPsiKst+1))  && (Utility->PsiRejection(NTuple->B0MassArb,NTuple->mumuMass->at(0),NTuple->mumuMassE->at(0),"keepJpsi")       == true)) ||
           (((SignalType == Utility->B0ToPsi2SKst) || (SignalType == Utility->B0ToPsi2SKst+1)) && (Utility->PsiRejection(NTuple->B0MassArb,NTuple->mumuMass->at(0),NTuple->mumuMassE->at(0),"keepPsiP")       == true))))))
    {
      evPassing++;
      mumuq2 = NTuple->mumuMass->at(0)*NTuple->mumuMass->at(0);
      if ((type == 4) && (NTuple->rightFlavorTag == false))
      {
        cosThetaK  = - NTuple->CosThetaKArb;
        cosThetaMu = - NTuple->CosThetaMuArb;
      }
      else
      {
        cosThetaK  = NTuple->CosThetaKArb;
        cosThetaMu = NTuple->CosThetaMuArb;
      }
      phiKstMuMuPlane = NTuple->PhiKstMuMuPlaneArb;

      mumuq2BinIndx          = Utility->SearchBin(mumuq2,q2Bins);
      cosThetaKBinIndx       = Utility->SearchBin(cosThetaK,cosThetaKBins);
      cosThetaMuBinIndx      = Utility->SearchBin(cosThetaMu,cosThetaLBins);
      phiKstMuMuPlaneBinIndx = Utility->SearchBin(phiKstMuMuPlane,phiBins);

      if ((mumuq2BinIndx != -1) && (cosThetaKBinIndx != -1) && (cosThetaMuBinIndx != -1) && (phiKstMuMuPlaneBinIndx != -1))
      {
        if ((type == 1) || (type == 2))
          // ##################################################
          // # No event-weight for muon acceptance efficiency #
          // ##################################################
          Vector[phiKstMuMuPlaneBinIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
            cosThetaMuBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
            cosThetaKBinIndx*(q2Bins->size()-1) +
            mumuq2BinIndx]++;
        else if ((type == 3) || (type == 4))
        {
          // ########################################
          // # Add event-weight for reco efficiency #
          // ########################################
          Vector[phiKstMuMuPlaneBinIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
            cosThetaMuBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
            cosThetaKBinIndx*(q2Bins->size()-1) +
            mumuq2BinIndx] +=
            NTuple->evWeight;

          Counter[phiKstMuMuPlaneBinIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
            cosThetaMuBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
            cosThetaKBinIndx*(q2Bins->size()-1) +
            mumuq2BinIndx]++;

          VectorErr2Weig[phiKstMuMuPlaneBinIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
            cosThetaMuBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
            cosThetaKBinIndx*(q2Bins->size()-1) +
            mumuq2BinIndx] +=
            NTuple->evWeightE2;
        }
      }
    }
  }
  cout << "[ComputeEfficiency::ComputeEfficiency]\t@@@ Number of events passing selection: " << evPassing << " @@@" << endl;


  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	for (unsigned int l = 0; l < phiBins->size()-1; l++)
	  {
	    if ((type == 1) || (type == 2))
	      {
		VectorErr2Pois[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] =
		  Vector[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      }

	    else if (Counter[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] == 0.0)
	      {
		VectorErr2Pois[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = 0.0;
		VectorErr2Weig[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = 0.0;
	      }

	    else if ((type == 3) || (type == 4))
	      VectorErr2Pois[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] =
		Counter[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] *
		pow(Vector[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] /
		    Counter[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0);
	  }


  Utility->DeleteEfficiency(_Counter);
}


void MakeHistogramsAllBins (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, Utils::effStruct myEff, int specBin)
{
  // ###################
  // # Local variables #
  // ###################
  stringstream myString;
  double Eff, EffErr;
  double totalEffAll;
  double totalEffSignal;
  double totalEffJPsi;
  double totalEffPsiP;
  double* q2Bins_        = Utility->MakeBinning(q2Bins);
  double* cosThetaKBins_ = Utility->MakeBinning(cosThetaKBins);
  double* cosThetaLBins_ = Utility->MakeBinning(cosThetaLBins);
  double* phiBins_       = Utility->MakeBinning(phiBins);
  // ###################
  string NumOrDen2Plot = "N2"; // It can be: "N1", "N2", "D1", "D2"
  double Yaxes = 2e5;
  double Zaxes = 1e3;
  // ###################


  // ##########################
  // # Efficiency projections #
  // ##########################
  TCanvas* cEff = new TCanvas("cEff", "cEff", 10, 10, 1200, 800);
  cEff->Divide(2,2);

  TH1D* Hq2 = new TH1D("Hq2", "Hq2", q2Bins->size()-1, q2Bins_);
  Hq2->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
  Hq2->SetYTitle("Efficiency");

  vector<TH1D*> vecHcosThetaK;
  vector<TH1D*> vecHcosThetaL;
  vector<TH1D*> vecHphi;

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    {
      myString.clear(); myString.str("");
      myString << "vecHcosThetaK_" << i;
      vecHcosThetaK.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), cosThetaKBins->size()-1, cosThetaKBins_));
      vecHcosThetaK.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
      vecHcosThetaK.back()->SetYTitle("Efficiency");

      myString.clear(); myString.str("");
      myString << "vecHcosThetaL_" << i;
      vecHcosThetaL.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), cosThetaLBins->size()-1, cosThetaLBins_));
      vecHcosThetaL.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
      vecHcosThetaL.back()->SetYTitle("Efficiency");

      myString.clear(); myString.str("");
      myString << "vecHphi_" << i;
      vecHphi.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), phiBins->size()-1, phiBins_));
      vecHphi.back()->SetXTitle("#phi");
      vecHphi.back()->SetYTitle("Efficiency");
    }


  for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
    {
      // #############################################
      // # Fill histogram : efficiency vs dimuon q^2 #
      // #############################################
      Utility->IntegrateEffButq2(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,i,myEff,&Eff,&EffErr);
      Hq2->SetBinContent(i+1,Eff);
      Hq2->SetBinError(i+1,EffErr);


      // ########################################################################
      // # Fill histogram : efficiency vs cos(theta_K) for different dimuon q^2 #
      // ########################################################################
      for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	{
	  Utility->IntegrateEffPhiCosThetaL(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,i,j,myEff,&Eff,&EffErr);
	  vecHcosThetaK[i]->SetBinContent(j+1,Eff);
	  vecHcosThetaK[i]->SetBinError(j+1,EffErr);
	}


      // ########################################################################
      // # Fill histogram : efficiency vs cos(theta_l) for different dimuon q^2 #
      // ########################################################################
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	{
	  Utility->IntegrateEffPhiCosThetaK(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,i,k,myEff,&Eff,&EffErr);
	  vecHcosThetaL[i]->SetBinContent(k+1,Eff);
	  vecHcosThetaL[i]->SetBinError(k+1,EffErr);
	}


      // ###############################################################
      // # Fill histogram : efficiency vs phi for different dimuon q^2 #
      // ###############################################################
      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  Utility->IntegrateEffCosThetaKCosThetaL(q2Bins,cosThetaKBins,cosThetaLBins,i,l,myEff,&Eff,&EffErr);
	  vecHphi[i]->SetBinContent(l+1,Eff);
	  vecHphi[i]->SetBinError(l+1,EffErr);
	}
    }


  cEff->cd(1);
  Hq2->SetMarkerStyle(20);
  Hq2->GetYaxis()->SetRangeUser(0.0,ordinateRange);
  Hq2->Draw("e1");

  cEff->cd(2);
  TLegend* legThetaK;
  if (specBin == -1)
    {
      vecHcosThetaK[0]->SetMarkerStyle(20);
      vecHcosThetaK[0]->SetMarkerColor(1);
      vecHcosThetaK[0]->SetLineColor(1);
      vecHcosThetaK[0]->SetLineWidth(2);
      vecHcosThetaK[0]->Draw("e1");
      vecHcosThetaK[0]->GetYaxis()->SetRangeUser(0.0,ordinateRange);
      legThetaK = new TLegend(0.88, 0.65, 0.97, 0.89, "");
      myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin 0";
      legThetaK->AddEntry(vecHcosThetaK[0],myString.str().c_str());
      for (unsigned int i = 1; i < q2Bins->size()-1; i++)
	{
	  vecHcosThetaK[i]->SetMarkerStyle(20+i);
	  vecHcosThetaK[i]->SetMarkerColor(1+i);
	  vecHcosThetaK[i]->SetLineColor(1+i);
	  vecHcosThetaK[i]->SetLineWidth(2);
	  vecHcosThetaK[i]->Draw("sames e1");
	  vecHcosThetaK[i]->GetYaxis()->SetRangeUser(0.0,ordinateRange);
	  myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin " << i;
	  legThetaK->AddEntry(vecHcosThetaK[i],myString.str().c_str());
	}
    }
  else
    {
      vecHcosThetaK[specBin]->SetMarkerStyle(20);
      vecHcosThetaK[specBin]->SetMarkerColor(1);
      vecHcosThetaK[specBin]->SetLineColor(1);
      vecHcosThetaK[specBin]->SetLineWidth(2);
      vecHcosThetaK[specBin]->Draw("e1");
      vecHcosThetaK[specBin]->GetYaxis()->SetRangeUser(0.0,ordinateRange);
      legThetaK = new TLegend(0.88, 0.65, 0.97, 0.89, "");
      myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin 0";
      legThetaK->AddEntry(vecHcosThetaK[specBin],myString.str().c_str());
    }
  legThetaK->SetFillColor(0);
  legThetaK->SetBorderSize(0);
  legThetaK->Draw();

  cEff->cd(3);
  TLegend* legThetaL;
  if (specBin == -1)
    {
      vecHcosThetaL[0]->SetMarkerStyle(20);
      vecHcosThetaL[0]->SetMarkerColor(1);
      vecHcosThetaL[0]->SetLineColor(1);
      vecHcosThetaL[0]->SetLineWidth(2);
      vecHcosThetaL[0]->Draw("e1");
      vecHcosThetaL[0]->GetYaxis()->SetRangeUser(0.0,ordinateRange);
      legThetaL = new TLegend(0.88, 0.65, 0.97, 0.89, "");
      myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin 0";
      legThetaL->AddEntry(vecHcosThetaL[0],myString.str().c_str());
      for (unsigned int i = 1; i < q2Bins->size()-1; i++)
	{
	  vecHcosThetaL[i]->SetMarkerStyle(20+i);
	  vecHcosThetaL[i]->SetMarkerColor(1+i);
	  vecHcosThetaL[i]->SetLineColor(1+i);
	  vecHcosThetaL[i]->SetLineWidth(2);
	  vecHcosThetaL[i]->Draw("sames e1");
	  vecHcosThetaL[i]->GetYaxis()->SetRangeUser(0.0,ordinateRange);
	  myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin " << i;
	  legThetaL->AddEntry(vecHcosThetaL[i],myString.str().c_str());
	}
    }
  else
    {
      vecHcosThetaL[specBin]->SetMarkerStyle(20);
      vecHcosThetaL[specBin]->SetMarkerColor(1);
      vecHcosThetaL[specBin]->SetLineColor(1);
      vecHcosThetaL[specBin]->SetLineWidth(2);
      vecHcosThetaL[specBin]->Draw("e1");
      vecHcosThetaL[specBin]->GetYaxis()->SetRangeUser(0.0,ordinateRange);
      legThetaL = new TLegend(0.88, 0.65, 0.97, 0.89, "");
      myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin 0";
      legThetaL->AddEntry(vecHcosThetaL[specBin],myString.str().c_str());
    }
  legThetaL->SetFillColor(0);
  legThetaL->SetBorderSize(0);
  legThetaL->Draw();

  cEff->cd(4);
  TLegend* legPhi;
  if (specBin == -1)
    {
      vecHphi[0]->SetMarkerStyle(20);
      vecHphi[0]->SetMarkerColor(1);
      vecHphi[0]->SetLineColor(1);
      vecHphi[0]->SetLineWidth(2);
      vecHphi[0]->Draw("e1");
      vecHphi[0]->GetYaxis()->SetRangeUser(0.0,ordinateRange);
      legPhi = new TLegend(0.88, 0.65, 0.97, 0.89, "");
      myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin 0";
      legPhi->AddEntry(vecHphi[0],myString.str().c_str());
      for (unsigned int i = 1; i < q2Bins->size()-1; i++)
	{
	  vecHphi[i]->SetMarkerStyle(20+i);
	  vecHphi[i]->SetMarkerColor(1+i);
	  vecHphi[i]->SetLineColor(1+i);
	  vecHphi[i]->SetLineWidth(2);
	  vecHphi[i]->Draw("sames e1");
	  vecHphi[i]->GetYaxis()->SetRangeUser(0.0,ordinateRange);
	  myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin " << i;
	  legPhi->AddEntry(vecHphi[i],myString.str().c_str());
	}
    }
  else
    {
      vecHphi[specBin]->SetMarkerStyle(20);
      vecHphi[specBin]->SetMarkerColor(1);
      vecHphi[specBin]->SetLineColor(1);
      vecHphi[specBin]->SetLineWidth(2);
      vecHphi[specBin]->Draw("e1");
      vecHphi[specBin]->GetYaxis()->SetRangeUser(0.0,ordinateRange);
      legPhi = new TLegend(0.88, 0.65, 0.97, 0.89, "");
      myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin 0";
      legPhi->AddEntry(vecHphi[specBin],myString.str().c_str());
    }
  legPhi->SetFillColor(0);
  legPhi->SetBorderSize(0);
  legPhi->Draw();

  cEff->Modified();
  cEff->Update();


  // #############################
  // # Eff. Terms 2D projections #
  // #############################
  TCanvas* cNumCosThetaK = new TCanvas("cNumCosThetaK", "cNumCosThetaK", 10, 10, 1200, 800);
  cNumCosThetaK->Divide(2,3);

  TCanvas* cNumCosThetaL = new TCanvas("cNumCosThetaL", "cNumCosThetaL", 10, 10, 1200, 800);
  cNumCosThetaL->Divide(2,3);

  TCanvas* cNumPhi = new TCanvas("cNumPhi", "cNumPhi", 10, 10, 1200, 800);
  cNumPhi->Divide(2,3);

  vector<TH2D*> vecHq2ANDcosThetaK;
  vector<TH2D*> vecHq2ANDcosThetaL;
  vector<TH2D*> vecHq2ANDphi;

  vector<TLegend*> legHq2ANDcosThetaK;
  vector<TLegend*> legHq2ANDcosThetaL;
  vector<TLegend*> legHq2ANDphi;

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    {
      for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	{
	  myString.clear(); myString.str("");
	  myString << "vecHq2ANDcosThetaK_" << i << "_" << j;
	  vecHq2ANDcosThetaK.push_back(new TH2D(myString.str().c_str(), myString.str().c_str(), cosThetaLBins->size()-1, cosThetaLBins_, phiBins->size()-1, phiBins_));
	  vecHq2ANDcosThetaK.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
	  vecHq2ANDcosThetaK.back()->SetYTitle("#phi");
	  vecHq2ANDcosThetaK.back()->SetZTitle("Eff. Terms [#]");

	  legHq2ANDcosThetaK.push_back(new TLegend(0.88, 0.6, 0.97, 0.89, ""));
	}

      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	{
	  myString.clear(); myString.str("");
	  myString << "vecHq2ANDcosThetaL_" << i << "_" << k;
	  vecHq2ANDcosThetaL.push_back(new TH2D(myString.str().c_str(), myString.str().c_str(), cosThetaKBins->size()-1, cosThetaKBins_, phiBins->size()-1, phiBins_));
	  vecHq2ANDcosThetaL.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
	  vecHq2ANDcosThetaL.back()->SetYTitle("#phi");
	  vecHq2ANDcosThetaL.back()->SetZTitle("Eff. Terms [#]");

	  legHq2ANDcosThetaL.push_back(new TLegend(0.88, 0.6, 0.97, 0.89, ""));
	}

      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  myString.clear(); myString.str("");
	  myString << "vecHq2ANDphi_" << i << "_" << l;
	  vecHq2ANDphi.push_back(new TH2D(myString.str().c_str(), myString.str().c_str(), cosThetaLBins->size()-1, cosThetaLBins_, cosThetaKBins->size()-1, cosThetaKBins_));
	  vecHq2ANDphi.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
	  vecHq2ANDphi.back()->SetYTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
	  vecHq2ANDphi.back()->SetZTitle("Eff. Terms [#]");

	  legHq2ANDphi.push_back(new TLegend(0.88, 0.6, 0.97, 0.89, ""));
	}
    }


  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
    {
      for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
	{
	  // ###############################################################################################################
	  // # Fill histogram : numerator of the efficiency vs cos(theta_l) in bins of dimuon q^2 and cos(theta_K) and phi #
	  // ###############################################################################################################
	  for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	    for (unsigned int l = 0; l < phiBins->size()-1; l++)
	      {
		if      (NumOrDen2Plot == "N1") Eff = myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "N2") Eff = myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "D1") Eff = myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else                            Eff = myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

		vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->SetBinContent(k+1,l+1,Eff);
	      }
	  
	  cNumCosThetaK->cd(j+1)->SetLogz();
	  vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->SetMarkerStyle(1);
	  vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->SetMarkerColor(1+i);
	  vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->SetLineColor(1+i);
	  vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->SetLineWidth(2);
	  if (i == 0) vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->Draw("surf bb fb");
	  else        vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->Draw("sames surf bb fb");
	  vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i]->GetZaxis()->SetRangeUser(0.1,Zaxes);

	  myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin " << i;
	  legHq2ANDcosThetaK[j]->AddEntry(vecHq2ANDcosThetaK[j*(q2Bins->size()-1)+i],myString.str().c_str());
	}

      legHq2ANDcosThetaK[j]->SetFillColor(0);
      legHq2ANDcosThetaK[j]->SetBorderSize(0);
      legHq2ANDcosThetaK[j]->Draw();
    }
  

  for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
    {
      for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
	{
	  // #######################################################################################################
	  // # Fill histogram : numerator of the efficiency vs cos(theta_K) in bins of dimuon q^2 and cos(theta_l) #
	  // #######################################################################################################
	  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	    for (unsigned int l = 0; l < phiBins->size()-1; l++)
	      {
		if      (NumOrDen2Plot == "N1") Eff = myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "N2") Eff = myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "D1") Eff = myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else                            Eff = myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

		vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->SetBinContent(j+1,l+1,Eff);
	      }

	  cNumCosThetaL->cd(k+1)->SetLogz();
	  vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->SetMarkerStyle(1);
	  vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->SetMarkerColor(1+i);
	  vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->SetLineColor(1+i);
	  vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->SetLineWidth(2);
	  if (i == 0) vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->Draw("surf bb fb");
	  else vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->Draw("sames surf bb fb");
	  vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i]->GetZaxis()->SetRangeUser(0.1,Zaxes);

	  myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin " << i;
	  legHq2ANDcosThetaL[k]->AddEntry(vecHq2ANDcosThetaL[k*(q2Bins->size()-1)+i],myString.str().c_str());
	}

      legHq2ANDcosThetaL[k]->SetFillColor(0);
      legHq2ANDcosThetaL[k]->SetBorderSize(0);
      legHq2ANDcosThetaL[k]->Draw();
    }


  for (unsigned int l = 0; l < phiBins->size()-1; l++)
    {
      for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
	{
	  // ###############################################################################################################
	  // # Fill histogram : numerator of the efficiency vs phi in bins of dimuon q^2 and cos(theta_l) and cos(theta_k) #
	  // ###############################################################################################################
	  for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	      {
		if      (NumOrDen2Plot == "N1") Eff = myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "N2") Eff = myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "D1") Eff = myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else                            Eff = myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

		vecHq2ANDphi[l*(q2Bins->size()-1)+i]->SetBinContent(k+1,j+1,Eff);
	      }
	  
	  cNumPhi->cd(l+1)->SetLogz();
	  vecHq2ANDphi[l*(q2Bins->size()-1)+i]->SetMarkerStyle(1);
	  vecHq2ANDphi[l*(q2Bins->size()-1)+i]->SetMarkerColor(1+i);
	  vecHq2ANDphi[l*(q2Bins->size()-1)+i]->SetLineColor(1+i);
	  vecHq2ANDphi[l*(q2Bins->size()-1)+i]->SetLineWidth(2);
	  if (i == 0) vecHq2ANDphi[l*(q2Bins->size()-1)+i]->Draw("surf bb fb");
	  else        vecHq2ANDphi[l*(q2Bins->size()-1)+i]->Draw("sames surf bb fb");
	  vecHq2ANDphi[l*(q2Bins->size()-1)+i]->GetZaxis()->SetRangeUser(0.1,Zaxes);

	  myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin " << i;
	  legHq2ANDphi[l]->AddEntry(vecHq2ANDphi[l*(q2Bins->size()-1)+i],myString.str().c_str());
	}

      legHq2ANDphi[l]->SetFillColor(0);
      legHq2ANDphi[l]->SetBorderSize(0);
      legHq2ANDphi[l]->Draw();
    }


  cNumCosThetaL->Modified();
  cNumCosThetaL->Update();
  cNumCosThetaK->Modified();
  cNumCosThetaK->Update();
  cNumPhi->Modified();
  cNumPhi->Update();


  // #############################
  // # Eff. Terms 1D projections #
  // #############################
  TCanvas* cNum = new TCanvas("cNum", "cNum", 10, 10, 1200, 800);
  cNum->Divide(2,2);

  vector<TH1D*> vecHq2ANDNumcosThetaK;
  vector<TH1D*> vecHq2ANDNumcosThetaL;
  vector<TH1D*> vecHq2ANDNumPhi;

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    {
      myString.clear(); myString.str("");
      myString << "vecHq2ANDNumcosThetaK_" << i;
      vecHq2ANDNumcosThetaK.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), cosThetaKBins->size()-1, cosThetaKBins_));
      vecHq2ANDNumcosThetaK.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
      vecHq2ANDNumcosThetaK.back()->SetYTitle("Eff. Terms [#]");

      myString.clear(); myString.str("");
      myString << "vecHq2ANDNumcosThetaL_" << i;
      vecHq2ANDNumcosThetaL.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), cosThetaLBins->size()-1, cosThetaLBins_));
      vecHq2ANDNumcosThetaL.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
      vecHq2ANDNumcosThetaL.back()->SetYTitle("Eff. Terms [#]");

      myString.clear(); myString.str("");
      myString << "vecHq2ANDNumPhi_" << i;
      vecHq2ANDNumPhi.push_back(new TH1D(myString.str().c_str(), myString.str().c_str(), phiBins->size()-1, phiBins_));
      vecHq2ANDNumPhi.back()->SetXTitle("#phi");
      vecHq2ANDNumPhi.back()->SetYTitle("Eff. Terms [#]");
    }


  for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
    {
      for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	{
	  Eff = 0.0;
	  
	  // ########################################################################
	  // # Fill histogram : efficiency vs cos(theta_K) for different dimuon q^2 #
	  // ########################################################################
	  for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	    for (unsigned int l = 0; l < phiBins->size()-1; l++)
	      {
		if      (NumOrDen2Plot == "N1") Eff += myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "N2") Eff += myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "D1") Eff += myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else                            Eff += myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      }

	  Eff = Eff / (cosThetaKBins->operator[](j+1) - cosThetaKBins->operator[](j));
	  vecHq2ANDNumcosThetaK[i]->SetBinContent(j+1,Eff);
	}


      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	{
	  Eff = 0.0;

	  // ########################################################################
	  // # Fill histogram : efficiency vs cos(theta_l) for different dimuon q^2 #
	  // ########################################################################
	  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	    for (unsigned int l = 0; l < phiBins->size()-1; l++)
	      {
		if      (NumOrDen2Plot == "N1") Eff += myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "N2") Eff += myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "D1") Eff += myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else                            Eff += myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      }

	  Eff = Eff / (cosThetaLBins->operator[](k+1) - cosThetaLBins->operator[](k));
	  vecHq2ANDNumcosThetaL[i]->SetBinContent(k+1,Eff);
	}


      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  Eff = 0.0;

	  // ###############################################################
	  // # Fill histogram : efficiency vs phi for different dimuon q^2 #
	  // ###############################################################
	  for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	      {
		if      (NumOrDen2Plot == "N1") Eff += myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "N2") Eff += myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else if (NumOrDen2Plot == "D1") Eff += myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		else                            Eff += myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      }

	  Eff = Eff / (phiBins->operator[](l+1) - phiBins->operator[](l));
	  vecHq2ANDNumPhi[i]->SetBinContent(l+1,Eff);
	}
    }


  cNum->cd(2);
  TLegend* legNumThetaK;
  if (specBin == -1)
    {
      vecHq2ANDNumcosThetaK[0]->SetMarkerStyle(20);
      vecHq2ANDNumcosThetaK[0]->SetMarkerColor(1);
      vecHq2ANDNumcosThetaK[0]->SetLineColor(1);
      vecHq2ANDNumcosThetaK[0]->SetLineWidth(2);
      vecHq2ANDNumcosThetaK[0]->Draw("p");
      vecHq2ANDNumcosThetaK[0]->GetYaxis()->SetRangeUser(0.1,Yaxes);
      legNumThetaK = new TLegend(0.88, 0.65, 0.97, 0.89, "");
      myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin 0";
      legNumThetaK->AddEntry(vecHq2ANDNumcosThetaK[0],myString.str().c_str());
      for (unsigned int i = 1; i < q2Bins->size()-1; i++)
	{
	  vecHq2ANDNumcosThetaK[i]->SetMarkerStyle(20+i);
	  vecHq2ANDNumcosThetaK[i]->SetMarkerColor(1+i);
	  vecHq2ANDNumcosThetaK[i]->SetLineColor(1+i);
	  vecHq2ANDNumcosThetaK[i]->SetLineWidth(2);
	  vecHq2ANDNumcosThetaK[i]->Draw("sames p");
	  vecHq2ANDNumcosThetaK[i]->GetYaxis()->SetRangeUser(0.1,Yaxes);
	  myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin " << i;
	  legNumThetaK->AddEntry(vecHq2ANDNumcosThetaK[i],myString.str().c_str());
	}
    }
  else
    {
      vecHq2ANDNumcosThetaK[specBin]->SetMarkerStyle(20);
      vecHq2ANDNumcosThetaK[specBin]->SetMarkerColor(1);
      vecHq2ANDNumcosThetaK[specBin]->SetLineColor(1);
      vecHq2ANDNumcosThetaK[specBin]->SetLineWidth(2);
      vecHq2ANDNumcosThetaK[specBin]->Draw("p");
      vecHq2ANDNumcosThetaK[specBin]->GetYaxis()->SetRangeUser(0.1,Yaxes);
      legNumThetaK = new TLegend(0.88, 0.65, 0.97, 0.89, "");
      myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin 0";
      legNumThetaK->AddEntry(vecHq2ANDNumcosThetaK[specBin],myString.str().c_str());
    }
  legNumThetaK->SetFillColor(0);
  legNumThetaK->SetBorderSize(0);
  legNumThetaK->Draw();

  cNum->cd(3);
  TLegend* legNumThetaL;
  if (specBin == -1)
    {
      vecHq2ANDNumcosThetaL[0]->SetMarkerStyle(20);
      vecHq2ANDNumcosThetaL[0]->SetMarkerColor(1);
      vecHq2ANDNumcosThetaL[0]->SetLineColor(1);
      vecHq2ANDNumcosThetaL[0]->SetLineWidth(2);
      vecHq2ANDNumcosThetaL[0]->Draw("p");
      vecHq2ANDNumcosThetaL[0]->GetYaxis()->SetRangeUser(0.1,Yaxes);
      legNumThetaL = new TLegend(0.88, 0.65, 0.97, 0.89, "");
      myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin 0";
      legNumThetaL->AddEntry(vecHq2ANDNumcosThetaL[0],myString.str().c_str());
      for (unsigned int i = 1; i < q2Bins->size()-1; i++)
	{
	  vecHq2ANDNumcosThetaL[i]->SetMarkerStyle(20+i);
	  vecHq2ANDNumcosThetaL[i]->SetMarkerColor(1+i);
	  vecHq2ANDNumcosThetaL[i]->SetLineColor(1+i);
	  vecHq2ANDNumcosThetaL[i]->SetLineWidth(2);
	  vecHq2ANDNumcosThetaL[i]->Draw("sames p");
	  vecHq2ANDNumcosThetaL[i]->GetYaxis()->SetRangeUser(0.1,Yaxes);
	  myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin " << i;
	  legNumThetaL->AddEntry(vecHq2ANDNumcosThetaL[i],myString.str().c_str());
	}
    }
  else
    {
      vecHq2ANDNumcosThetaL[specBin]->SetMarkerStyle(20);
      vecHq2ANDNumcosThetaL[specBin]->SetMarkerColor(1);
      vecHq2ANDNumcosThetaL[specBin]->SetLineColor(1);
      vecHq2ANDNumcosThetaL[specBin]->SetLineWidth(2);
      vecHq2ANDNumcosThetaL[specBin]->Draw("p");
      vecHq2ANDNumcosThetaL[specBin]->GetYaxis()->SetRangeUser(0.1,Yaxes);
      legNumThetaL = new TLegend(0.88, 0.65, 0.97, 0.89, "");
      myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin 0";
      legNumThetaL->AddEntry(vecHq2ANDNumcosThetaL[specBin],myString.str().c_str());
    }
  legNumThetaL->SetFillColor(0);
  legNumThetaL->SetBorderSize(0);
  legNumThetaL->Draw();

  cNum->cd(4);
  TLegend* legNumPhi;
  if (specBin == -1)
    {
      vecHq2ANDNumPhi[0]->SetMarkerStyle(20);
      vecHq2ANDNumPhi[0]->SetMarkerColor(1);
      vecHq2ANDNumPhi[0]->SetLineColor(1);
      vecHq2ANDNumPhi[0]->SetLineWidth(2);
      vecHq2ANDNumPhi[0]->Draw("p");
      vecHq2ANDNumPhi[0]->GetYaxis()->SetRangeUser(0.1,Yaxes);
      legNumPhi = new TLegend(0.88, 0.65, 0.97, 0.89, "");
      myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin 0";
      legNumPhi->AddEntry(vecHq2ANDNumPhi[0],myString.str().c_str());
      for (unsigned int i = 1; i < q2Bins->size()-1; i++)
	{
	  vecHq2ANDNumPhi[i]->SetMarkerStyle(20+i);
	  vecHq2ANDNumPhi[i]->SetMarkerColor(1+i);
	  vecHq2ANDNumPhi[i]->SetLineColor(1+i);
	  vecHq2ANDNumPhi[i]->SetLineWidth(2);
	  vecHq2ANDNumPhi[i]->Draw("sames p");
	  vecHq2ANDNumPhi[i]->GetYaxis()->SetRangeUser(0.1,Yaxes);
	  myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin " << i;
	  legNumPhi->AddEntry(vecHq2ANDNumPhi[i],myString.str().c_str());
	}
    }
  else
    {
      vecHq2ANDNumPhi[specBin]->SetMarkerStyle(20);
      vecHq2ANDNumPhi[specBin]->SetMarkerColor(1);
      vecHq2ANDNumPhi[specBin]->SetLineColor(1);
      vecHq2ANDNumPhi[specBin]->SetLineWidth(2);
      vecHq2ANDNumPhi[specBin]->Draw("p");
      vecHq2ANDNumPhi[specBin]->GetYaxis()->SetRangeUser(0.1,Yaxes);
      legNumPhi = new TLegend(0.88, 0.65, 0.97, 0.89, "");
      myString.clear(); myString.str(""); myString << "q#lower[0.4]{^{2}} bin 0";
      legNumPhi->AddEntry(vecHq2ANDNumPhi[specBin],myString.str().c_str());
    }
  legNumPhi->SetFillColor(0);
  legNumPhi->SetBorderSize(0);
  legNumPhi->Draw();

  cNum->Modified();
  cNum->Update();


  Utility->IntegrateEffAll(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,myEff,&totalEffAll,&EffErr);
  Utility->IntegrateEffButPsi(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,myEff,&totalEffSignal,&EffErr);
  Utility->IntegrateEffInJPsi(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,myEff,&totalEffJPsi,&EffErr);
  Utility->IntegrateEffInPsiP(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,myEff,&totalEffPsiP,&EffErr);
  cout << "\nTotal efficiency: " << totalEffAll <<  "\tTotal signal efficiency: " << totalEffSignal;
  cout << "\tTotal J/psi efficiency: " << totalEffJPsi << "\tTotal psi(2S) efficiency: " << totalEffPsiP << endl;
}


void ReadEfficiencies (bool isSingleEff, vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
		       string fileNameInput, bool isAnalyEff, Utils::effStruct* myEff, bool CheckNegEff, bool savePlot, int specBin, bool EffIs2Dnot3D, int SignalType)
{
  // ###################
  // # Local variables #
  // ###################
  const unsigned int MAXVAL = (isSingleEff == true ? 1 : NFILES+1);
  double totalEffAll;
  double totalEffSignal;
  double totalEffJPsi;
  double totalEffPsiP;
  double Eff, EffErr;
  double* q2Bins_        = Utility->MakeBinning(q2Bins);
  double* cosThetaKBins_ = Utility->MakeBinning(cosThetaKBins);
  double* cosThetaLBins_ = Utility->MakeBinning(cosThetaLBins);
  double* phiBins_       = Utility->MakeBinning(phiBins);
  vector<TF2*> effFuncs2D;
  vector<TF3*> effFuncs3D;
  string tmpString;
  stringstream myString1;
  stringstream myString2;
  // ###################
  unsigned int nBins = 100;
  double Xaxes       = 0.1e-3;
  // ###################


  // # ###############################################################
  // # Canvas for efficiencies integrated over all variables but one #
  // # ###############################################################
  TCanvas* cEff0 = new TCanvas("cEff0", "Efficiency integrated over all variables but one", 10, 10, 700, 500);
  TCanvas* cEff1 = new TCanvas("cEff1", "Efficiency integrated over all variables but one", 20, 20, 700, 500);
  TCanvas* cEff2 = new TCanvas("cEff2", "Efficiency integrated over all variables but one", 30, 30, 700, 500);
  TCanvas* cEff3 = new TCanvas("cEff3", "Efficiency integrated over all variables but one", 40, 40, 700, 500);

  tmpString = Utility->MakeAnalysisPATH(GENEFF);
  tmpString.erase(tmpString.find(".txt"),4);

  vector<TH1D*> Hq2;
  for (unsigned int i = 0; i < MAXVAL; i++)
    {
      myString1.clear(); myString1.str("");
      myString1 << "Hq2_" << i;
      Hq2.push_back(new TH1D(myString1.str().c_str(), myString1.str().c_str(), q2Bins->size()-1, q2Bins_));
      Hq2.back()->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
      Hq2.back()->SetYTitle("Efficiency");
    }

  vector<TH1D*> HcosThetaK;
  for (unsigned int i = 0; i < MAXVAL; i++)
    {
      myString1.clear(); myString1.str("");
      myString1 << "HcosThetaK_" << i;
      HcosThetaK.push_back(new TH1D(myString1.str().c_str(), myString1.str().c_str(), cosThetaKBins->size()-1, cosThetaKBins_));
      HcosThetaK.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
      HcosThetaK.back()->SetYTitle("Efficiency");
    }

  vector<TH1D*> HcosThetaL;
  for (unsigned int i = 0; i < MAXVAL; i++)
    {
      myString1.clear(); myString1.str("");
      myString1 << "HcosThetaL_" << i;
      HcosThetaL.push_back(new TH1D(myString1.str().c_str(), myString1.str().c_str(), cosThetaLBins->size()-1, cosThetaLBins_));
      HcosThetaL.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
      HcosThetaL.back()->SetYTitle("Efficiency");
    }

  vector<TH1D*> Hphi;
  for (unsigned int i = 0; i < MAXVAL; i++)
    {
      myString1.clear(); myString1.str("");
      myString1 << "Hphi_" << i;
      Hphi.push_back(new TH1D(myString1.str().c_str(), myString1.str().c_str(), phiBins->size()-1, phiBins_));
      Hphi.back()->SetXTitle("#phi");
      Hphi.back()->SetYTitle("Efficiency");
    }


  // ###################################
  // # Initialize efficiency structure #
  // ###################################
  Utility->GenEfficiency(myEff,q2Bins,cosThetaKBins,cosThetaLBins,phiBins);
  Utils::effValue myEffVal;
  Utility->ResetEffValue(&myEffVal,0.0);
  Utility->InitEfficiency(myEffVal,*myEff,q2Bins,cosThetaKBins,cosThetaLBins,phiBins);


  for (unsigned int itF = 0; itF < MAXVAL; itF++)
    {
      if ((isSingleEff == true) || (itF == MAXVAL-1))
	{
	  if (isAnalyEff == false)
	    {
	      Utility->ReadEfficiency(fileNameInput.c_str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,myEff);

	      myString1.clear(); myString1.str("");
	      myString1 << fileNameInput;
	    }
	  else
	    {
	      if (EffIs2Dnot3D == true) Utility->ReadAnalyticalEff(fileNameInput.c_str(),q2Bins,cosThetaKBins,cosThetaLBins,&effFuncs2D,"effFuncs2D",0);
	      else                      Utility->ReadAnalyticalEff(fileNameInput.c_str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,&effFuncs3D,"effFuncs3D",0);
	    }
	}
      else
	{
	  myString1.clear(); myString1.str("");
	  myString1 << tmpString << "_" << itF+1 << ".txt";

 	  if (isAnalyEff == false) Utility->ReadEfficiency(myString1.str().c_str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,myEff);
	  else if (EffIs2Dnot3D == true)
	    {	  
	      myString2.clear(); myString2.str("");
	      myString2 << "effFuncs2D_" << itF;
	      Utility->ReadAnalyticalEff(myString1.str().c_str(),q2Bins,cosThetaKBins,cosThetaLBins,&effFuncs2D,myString2.str().c_str(),1);
	    }
	  else 
	    {
	      myString2.clear(); myString2.str("");
	      myString2 << "effFuncs3D_" << itF;
	      Utility->ReadAnalyticalEff(myString1.str().c_str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,&effFuncs3D,myString2.str().c_str(),1);
	    }

	  cout << "[ComputeEfficiency::ReadEfficiencies]\tRead randomly generated efficiency n." << itF+1 << endl;
	}


      // ##########################
      // # Save binned efficiency #
      // ##########################
      if ((isAnalyEff == false) && (savePlot == true))
	{
	  myString2.clear(); myString2.str("");
	  myString2 << tmpString;

	  tmpString = myString1.str();
	  tmpString.erase(tmpString.find(".txt"),4);

	  myString1.clear(); myString1.str("");
	  myString1 << tmpString;

	  tmpString = myString2.str();


	  for (unsigned int q2BinIndx = (specBin == -1 ? 0 : specBin); q2BinIndx < (specBin == -1 ? q2Bins->size()-1 : specBin+1); q2BinIndx++)
	    {
	      myString2.clear(); myString2.str("");
	      myString2 << Utility->GetHisto2DEffName(SignalType) << "_" << itF << "_" << q2BinIndx;
	      TH2D* hisFunc2D = Utility->Get2DEffHistoButPhi(myString2.str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2BinIndx,*myEff);
	      
	      myString2.clear(); myString2.str("");
	      myString2 << myString1.str() << "_H2Deff_q2Bin_" << q2BinIndx << ".txt";
	      Utility->Put2DEffHistoq2Bin(myString2.str().c_str(),hisFunc2D);
	      
	      
	      myString2.clear(); myString2.str("");
	      myString2 << Utility->GetHisto3DEffName(SignalType) << "_" << itF << "_" << q2BinIndx;
	      TH3D* hisFunc3D = Utility->Get3DEffHistoq2Bin(myString2.str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2BinIndx,*myEff);
	      
	      myString2.clear(); myString2.str("");
	      myString2 << myString1.str() << "_H3Deff_q2Bin_" << q2BinIndx << ".txt";
	      Utility->Put3DEffHistoq2Bin(myString2.str().c_str(),hisFunc3D);
	    }
	}


      // ################################################
      // # Check if analytical efficiency goes negative #
      // ################################################
      if ((isAnalyEff == true) && (CheckNegEff == true))
      	for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
	  if (((EffIs2Dnot3D == true)  && (Utility->EffMinValue2D(cosThetaKBins,cosThetaLBins,effFuncs2D[i])         < 0.0)) ||
	      ((EffIs2Dnot3D == false) && (Utility->EffMinValue3D(cosThetaKBins,cosThetaLBins,phiBins,effFuncs3D[i]) < 0.0)))
	    cout << "@@@ Negative efficiency function #" << itF+1 << " for q2 bin #" << i << " ! @@@" << endl;

      for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
	{
	  // #############################################
	  // # Fill histogram : efficiency vs dimuon q^2 #
	  // #############################################

	  Eff    = 0.0;
	  EffErr = 0.0;

	  if (isAnalyEff == false) Utility->IntegrateEffButq2(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,i,*myEff,&Eff,&EffErr);
	  else
	    {	      
	      if (EffIs2Dnot3D == true)
		Eff = Eff + effFuncs2D[i]->Integral(cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),
						    cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1)) /
		  ((cosThetaKBins->operator[](cosThetaKBins->size()-1)-cosThetaKBins->operator[](0)) * (cosThetaLBins->operator[](cosThetaLBins->size()-1)-cosThetaLBins->operator[](0)));
	      else
		Eff = Eff + effFuncs3D[i]->Integral(cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),
						    cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1),
						    phiBins->operator[](0),phiBins->operator[](phiBins->size()-1)) /
		  ((cosThetaKBins->operator[](cosThetaKBins->size()-1)-cosThetaKBins->operator[](0)) *
		   (cosThetaLBins->operator[](cosThetaLBins->size()-1)-cosThetaLBins->operator[](0)) *
		   (phiBins->operator[](phiBins->size()-1)-phiBins->operator[](0)));	      
	    }

	  Hq2[itF]->SetBinContent(i+1,Eff);
	  Hq2[itF]->SetBinError(i+1,EffErr);
	}
      
      for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      	{
      	  // ###############################################
      	  // # Fill histogram : efficiency vs cos(theta_K) #
      	  // ###############################################

	  Eff    = 0.0;
	  EffErr = 0.0;
	  
      	  if (isAnalyEff == false) Utility->IntegrateEffButCosThetaK(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,j,*myEff,&Eff,&EffErr);
      	  else
      	    {
	      for (unsigned int indx = (specBin == -1 ? 0 : specBin); indx < (specBin == -1 ? q2Bins->size()-1 : specBin+1); indx++)
		if (EffIs2Dnot3D == true)
		  Eff = Eff + effFuncs2D[indx]->Integral(cosThetaKBins->operator[](j),cosThetaKBins->operator[](j+1),
							 cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1)) /
		    ((cosThetaKBins->operator[](j+1)-cosThetaKBins->operator[](j)) * (cosThetaLBins->operator[](cosThetaLBins->size()-1)-cosThetaLBins->operator[](0)));
		else
		  Eff = Eff + effFuncs3D[indx]->Integral(cosThetaKBins->operator[](j),cosThetaKBins->operator[](j+1),
							 cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1),
							 phiBins->operator[](0),phiBins->operator[](phiBins->size()-1)) /
		    ((cosThetaKBins->operator[](j+1)-cosThetaKBins->operator[](j)) *
		     (cosThetaLBins->operator[](cosThetaLBins->size()-1)-cosThetaLBins->operator[](0)) *
		     (phiBins->operator[](phiBins->size()-1)-phiBins->operator[](0)));
	      Eff = Eff / (q2Bins->operator[](specBin == -1 ? q2Bins->size()-1 : specBin+1)-q2Bins->operator[](specBin == -1 ? 0 : specBin));
	      EffErr = EffErr / (q2Bins->operator[](specBin == -1 ? q2Bins->size()-1 : specBin+1)-q2Bins->operator[](specBin == -1 ? 0 : specBin));
      	    }

      	  HcosThetaK[itF]->SetBinContent(j+1,Eff);
      	  HcosThetaK[itF]->SetBinError(j+1,EffErr);
      	}

      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
      	{
      	  // ###############################################
      	  // # Fill histogram : efficiency vs cos(theta_l) #
      	  // ###############################################

	  Eff    = 0.0;
	  EffErr = 0.0;
	  
      	  if (isAnalyEff == false) Utility->IntegrateEffButCosThetaL(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,k,*myEff,&Eff,&EffErr);
      	  else
      	    {
	      for (unsigned int indx = (specBin == -1 ? 0 : specBin); indx < (specBin == -1 ? q2Bins->size()-1 : specBin+1); indx++)
		if (EffIs2Dnot3D == true)
		  Eff = Eff + effFuncs2D[indx]->Integral(cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),
							 cosThetaLBins->operator[](k),cosThetaLBins->operator[](k+1)) /
		    ((cosThetaKBins->operator[](cosThetaKBins->size()-1)-cosThetaKBins->operator[](0)) * (cosThetaLBins->operator[](k+1)-cosThetaLBins->operator[](k)));
		else
		  Eff = Eff + effFuncs3D[indx]->Integral(cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),
							 cosThetaLBins->operator[](k),cosThetaLBins->operator[](k+1),
							 phiBins->operator[](0),phiBins->operator[](phiBins->size()-1)) /
		    ((cosThetaKBins->operator[](cosThetaKBins->size()-1)-cosThetaKBins->operator[](0)) *
		     (cosThetaLBins->operator[](k+1)-cosThetaLBins->operator[](k)) *
		     (phiBins->operator[](phiBins->size()-1)-phiBins->operator[](0)));	      
	      Eff = Eff / (q2Bins->operator[](specBin == -1 ? q2Bins->size()-1 : specBin+1)-q2Bins->operator[](specBin == -1 ? 0 : specBin));
	      EffErr = EffErr / (q2Bins->operator[](specBin == -1 ? q2Bins->size()-1 : specBin+1)-q2Bins->operator[](specBin == -1 ? 0 : specBin));
      	    }

      	  HcosThetaL[itF]->SetBinContent(k+1,Eff);
      	  HcosThetaL[itF]->SetBinError(k+1,EffErr);
      	}

      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  // ######################################
	  // # Fill histogram : efficiency vs phi #
	  // ######################################

	  Eff    = 0.0;
	  EffErr = 0.0;

	  if (isAnalyEff == false) Utility->IntegrateEffButPhi(q2Bins,cosThetaKBins,cosThetaLBins,l,*myEff,&Eff,&EffErr);
	  else if (EffIs2Dnot3D == false)
	    {
	      for (unsigned int indx = (specBin == -1 ? 0 : specBin); indx < (specBin == -1 ? q2Bins->size()-1 : specBin+1); indx++)
		Eff = Eff + effFuncs3D[indx]->Integral(cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),
						       cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1),
						       phiBins->operator[](l),phiBins->operator[](l+1)) /
		  ((cosThetaKBins->operator[](cosThetaKBins->size()-1)-cosThetaKBins->operator[](0)) *
		   (cosThetaLBins->operator[](cosThetaLBins->size()-1)-cosThetaLBins->operator[](0)) *
		   (phiBins->operator[](l+1)-phiBins->operator[](l)));	      
	      Eff = Eff / (q2Bins->operator[](specBin == -1 ? q2Bins->size()-1 : specBin+1)-q2Bins->operator[](specBin == -1 ? 0 : specBin));
	      EffErr = EffErr / (q2Bins->operator[](specBin == -1 ? q2Bins->size()-1 : specBin+1)-q2Bins->operator[](specBin == -1 ? 0 : specBin));
	    }

	  Hphi[itF]->SetBinContent(l+1,Eff);
	  Hphi[itF]->SetBinError(l+1,EffErr);
	}
      
      
      if ((isSingleEff == true) || (itF == MAXVAL-1))
	{
	  cEff0->cd();
	  Hq2[itF]->SetMarkerStyle(20);
	  Hq2[itF]->GetYaxis()->SetRangeUser(0,ordinateRange);
	  if (isSingleEff == true) Hq2[itF]->Draw("pe1");
	  else                     Hq2[itF]->Draw("same pe1");
	  
	  cEff1->cd();
	  HcosThetaK[itF]->SetMarkerStyle(20);
	  HcosThetaK[itF]->GetYaxis()->SetRangeUser(0,ordinateRange);
	  if (isSingleEff == true) HcosThetaK[itF]->Draw("pe1");
	  else                     HcosThetaK[itF]->Draw("same pe1");

	  cEff2->cd();
	  HcosThetaL[itF]->SetMarkerStyle(20);
	  HcosThetaL[itF]->GetYaxis()->SetRangeUser(0,ordinateRange);
	  if (isSingleEff == true) HcosThetaL[itF]->Draw("pe1");
	  else                     HcosThetaL[itF]->Draw("same pe1");

	  cEff3->cd();
	  Hphi[itF]->SetMarkerStyle(20);
	  Hphi[itF]->GetYaxis()->SetRangeUser(0,ordinateRange);
	  if (isSingleEff == true) Hphi[itF]->Draw("pe1");
	  else                     Hphi[itF]->Draw("same pe1");
	}
      else
	{
	  cEff0->cd();
	  Hq2[itF]->SetLineColor(kBlue);
	  Hq2[itF]->GetYaxis()->SetRangeUser(0,ordinateRange);
	  if (itF == 0) Hq2[itF]->Draw();
	  else          Hq2[itF]->Draw("same");

	  cEff1->cd();
	  HcosThetaK[itF]->SetLineColor(kBlue);
	  HcosThetaK[itF]->GetYaxis()->SetRangeUser(0,ordinateRange);
	  if (itF == 0) HcosThetaK[itF]->Draw();
	  else          HcosThetaK[itF]->Draw("same");

	  cEff2->cd();
	  HcosThetaL[itF]->SetLineColor(kBlue);
	  HcosThetaL[itF]->GetYaxis()->SetRangeUser(0,ordinateRange);
	  if (itF == 0) HcosThetaL[itF]->Draw();
	  else          HcosThetaL[itF]->Draw("same");

	  cEff3->cd();
	  Hphi[itF]->SetLineColor(kBlue);
	  Hphi[itF]->GetYaxis()->SetRangeUser(0,ordinateRange);
	  if (itF == 0) Hphi[itF]->Draw();
	  else          Hphi[itF]->Draw("same");
	}
      
      effFuncs2D.clear();
      effFuncs3D.clear();
    }
 
  cEff0->Modified();
  cEff0->Update();
  cEff1->Modified();
  cEff1->Update();
  cEff2->Modified();
  cEff2->Update();
  cEff3->Modified();
  cEff3->Update();

  if (savePlot == true)
    {
      if (isAnalyEff == false) cEff0->Print("Histo1.pdf");
      cEff1->Print("Histo2.pdf");
      cEff2->Print("Histo3.pdf");
      if (isAnalyEff == false) cEff3->Print("Histo4.pdf");
    }


  TCanvas* cStatK = new TCanvas("cStatK", "cStatK", 10, 10, 1200, 800);
  cStatK->Divide(2,3);

  vector<TH1D*> HstatK;
  for (unsigned int i = 0; i < cosThetaKBins->size()-1; i++)
    {
      myString1.clear(); myString1.str("");
      myString1 << "HstatK_" << i;
      HstatK.push_back(new TH1D(myString1.str().c_str(), myString1.str().c_str(), nBins, -Xaxes, Xaxes));
      HstatK.back()->SetXTitle("Diff. to reference");
      HstatK.back()->SetYTitle("Entries [#]");

      cStatK->cd(i+1);
      HstatK.back()->Draw();
    }


  TCanvas* cStatL = new TCanvas("cStatL", "cStatL", 10, 10, 1200, 800);
  cStatL->Divide(2,3);
  
  vector<TH1D*> HstatL;
  for (unsigned int i = 0; i < cosThetaLBins->size()-1; i++)
    {
      myString1.clear(); myString1.str("");
      myString1 << "HstatL_" << i;
      HstatL.push_back(new TH1D(myString1.str().c_str(), myString1.str().c_str(), nBins, -Xaxes, Xaxes));
      HstatL.back()->SetXTitle("Diff. to reference");
      HstatL.back()->SetYTitle("Entries [#]");

      cStatL->cd(i+1);
      HstatL.back()->Draw();
    }


  TCanvas* cStatPhi = new TCanvas("cStatPhi", "cStatPhi", 10, 10, 1200, 800);
  cStatPhi->Divide(2,3);
  
  vector<TH1D*> HstatPhi;
  for (unsigned int i = 0; i < phiBins->size()-1; i++)
    {
      myString1.clear(); myString1.str("");
      myString1 << "HstatPhi_" << i;
      HstatPhi.push_back(new TH1D(myString1.str().c_str(), myString1.str().c_str(), nBins, -Xaxes, Xaxes));
      HstatPhi.back()->SetXTitle("Diff. to reference");
      HstatPhi.back()->SetYTitle("Entries [#]");

      cStatPhi->cd(i+1);
      HstatPhi.back()->Draw();
    }


  if (isSingleEff == false)
    {
      for (unsigned int i = 0; i < cosThetaKBins->size()-1; i++)
	{
	  for (unsigned int itF = 0; itF < MAXVAL-1; itF++)
	    HstatK[i]->Fill(HcosThetaK[itF]->GetBinContent(i+1) - HcosThetaK[MAXVAL-1]->GetBinContent(i+1));
	  
  	  cStatK->cd(i+1);
  	  HstatK[i]->SetFillColor(kAzure+6);
  	  HstatK[i]->Draw();
  	}
      cStatK->Modified();
      cStatK->Update();


      for (unsigned int i = 0; i < cosThetaLBins->size()-1; i++)
	{
	  for (unsigned int itF = 0; itF < MAXVAL-1; itF++)
	    HstatL[i]->Fill(HcosThetaL[itF]->GetBinContent(i+1) - HcosThetaL[MAXVAL-1]->GetBinContent(i+1));

  	  cStatL->cd(i+1);
  	  HstatL[i]->SetFillColor(kAzure+6);
  	  HstatL[i]->Draw();
  	}
      cStatL->Modified();
      cStatL->Update();


      for (unsigned int i = 0; i < phiBins->size()-1; i++)
	{
	  for (unsigned int itF = 0; itF < MAXVAL-1; itF++)
	    HstatPhi[i]->Fill(Hphi[itF]->GetBinContent(i+1) - Hphi[MAXVAL-1]->GetBinContent(i+1));
	  
  	  cStatPhi->cd(i+1);
  	  HstatPhi[i]->SetFillColor(kAzure+6);
  	  HstatPhi[i]->Draw();
  	}
      cStatPhi->Modified();
      cStatPhi->Update();
    }


  Utility->IntegrateEffAll(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,*myEff,&totalEffAll,&EffErr);
  Utility->IntegrateEffButPsi(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,*myEff,&totalEffSignal,&EffErr);
  Utility->IntegrateEffInJPsi(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,*myEff,&totalEffJPsi,&EffErr);
  Utility->IntegrateEffInPsiP(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,*myEff,&totalEffPsiP,&EffErr);

  cout << "\nTotal efficiency: " << totalEffAll <<  "\tTotal signal efficiency: " << totalEffSignal;
  cout << "\tTotal J/psi efficiency: " << totalEffJPsi << "\tTotal psi(2S) efficiency: " << totalEffPsiP << endl;


  if ((isSingleEff == true) && (isAnalyEff == false)) MakeHistogramsAllBins(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,*myEff,specBin);
}


void GenerateEfficiencies (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, string fileNameInput)
{
  // ###################
  // # Local variables #
  // ###################
  ofstream fileOutput;

  string myString;
  stringstream fileNameOutput;

  double toGenN = 0.0;
  double toGenD = 0.0;
  double weight = 0.0;

  Utils::effStruct orgEff;
  Utils::effStruct newEff;
  Utils::effValue orgEffVal;
  Utils::effValue newEffVal;

  TRandom* rnd = new TRandom(1);
  // ###################


  Utility->ReadEfficiency(fileNameInput.c_str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,&orgEff);
  Utility->GenEfficiency(&newEff,q2Bins,cosThetaKBins,cosThetaLBins,phiBins);

  myString = GENEFF;
  myString.erase(myString.find(".txt"),4);


  for (unsigned int itF = 0; itF < NFILES; itF++)
    {
      Utility->ResetEffValue(&newEffVal,0.0);
      Utility->InitEfficiency(newEffVal,newEff,q2Bins,cosThetaKBins,cosThetaLBins,phiBins);

      for (unsigned int i = 0; i < q2Bins->size()-1; i++)
	for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	  for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	    for (unsigned int l = 0; l < phiBins->size()-1; l++)
	      {
		orgEffVal.Num1 = orgEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		orgEffVal.Num2 = orgEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		orgEffVal.Den1 = orgEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		orgEffVal.Den2 = orgEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

		orgEffVal.Err2PoisNum1 = orgEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		orgEffVal.Err2PoisNum2 = orgEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		orgEffVal.Err2PoisDen1 = orgEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		orgEffVal.Err2PoisDen2 = orgEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
		

		// ##############################
		// # Generate binned efficiency #
		// ##############################		
		if (orgEffVal.Err2PoisDen1 != 0.)
		  {
		    toGenD                 = orgEffVal.Den1 * orgEffVal.Den1 / orgEffVal.Err2PoisDen1;
		    weight                 = orgEffVal.Err2PoisDen1 / orgEffVal.Den1;
		    newEffVal.Den1         = rnd->PoissonD(toGenD) * weight;
		    newEffVal.Err2PoisDen1 = newEffVal.Den1 * weight;
		  }
		else
		  {
		    newEffVal.Den1         = 0.0;
		    newEffVal.Err2PoisDen1 = 0.0;
		  }
		
		if (orgEffVal.Err2PoisNum1 != 0.)
		  {
		    toGenN                 = orgEffVal.Num1 * orgEffVal.Num1 / orgEffVal.Err2PoisNum1;
		    weight                 = orgEffVal.Err2PoisNum1 / orgEffVal.Num1;
		    newEffVal.Num1         = fabs( rnd->Gaus(newEffVal.Den1/weight * toGenN / toGenD, sqrt (newEffVal.Den1/weight * toGenN / toGenD * (1. - toGenN / toGenD))) );
		    newEffVal.Err2PoisNum1 = newEffVal.Num1 * weight;
		  }
		else
		  {
		    newEffVal.Num1         = 0.0;
		    newEffVal.Err2PoisNum1 = 0.0;
		  }

		if (orgEffVal.Err2PoisDen2 != 0.)
		  {
		    toGenD                 = orgEffVal.Den2 * orgEffVal.Den2 / orgEffVal.Err2PoisDen2;
		    weight                 = orgEffVal.Err2PoisDen2 / orgEffVal.Den2;
		    newEffVal.Den2         = rnd->PoissonD(toGenD) * weight;
		    newEffVal.Err2PoisDen2 = newEffVal.Den2 * weight;
		  }
		else
		  {
		    newEffVal.Den2         = 0.0;
		    newEffVal.Err2PoisDen2 = 0.0;
		  }

		if (orgEffVal.Err2PoisNum2 != 0.)
		  {
		    toGenN                 = orgEffVal.Num2 * orgEffVal.Num2 / orgEffVal.Err2PoisNum2;
		    weight                 = orgEffVal.Err2PoisNum2 / orgEffVal.Num2;
		    newEffVal.Num2         = fabs( rnd->Gaus(newEffVal.Den2/weight * toGenN / toGenD, sqrt (newEffVal.Den2/weight * toGenN / toGenD * (1. - toGenN / toGenD))) );
		    newEffVal.Err2PoisNum2 = newEffVal.Num2 * weight;
		  }
		else
		  {
		    newEffVal.Num2         = 0.0;
		    newEffVal.Err2PoisNum2 = 0.0;
		  }


		newEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = newEffVal.Num1;
		newEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = newEffVal.Num2;
		newEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = newEffVal.Den1;
		newEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = newEffVal.Den2;

		newEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = newEffVal.Err2PoisNum1;
		newEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = newEffVal.Err2PoisNum2;
		newEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = newEffVal.Err2PoisDen1;
		newEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = newEffVal.Err2PoisDen2;
	      }


      // ##########################
      // # Save binned efficiency #
      // ##########################
      fileNameOutput.clear(); fileNameOutput.str("");
      fileNameOutput << myString.c_str() << "_" << itF+1 << ".txt";
      Utility->SaveEfficiency(fileNameOutput.str().c_str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,newEff);
      cout << "[ComputeEfficiency::GenerateEfficiencies]\tSaved randomly generated binned effiiency n." << itF+1 << endl;
    }


  delete rnd;
}


void Fit1DEfficiencies (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			Utils::effStruct myEff, string who, unsigned int q2BinIndx, string fileNameOut)
{
  // ###################
  // # Local variables #
  // ###################
  TFitResultPtr fitResults;
  double Eff, EffErr;
  double* cosThetaKBins_ = Utility->MakeBinning(cosThetaKBins);
  double* cosThetaLBins_ = Utility->MakeBinning(cosThetaLBins);
  double* phiBins_       = Utility->MakeBinning(phiBins);
  vector<TH1D*> histFit;
  TF1* effFunc1D;
  stringstream myString;
  string tmpString;
  double coeff;
  ofstream fileOutput;
  ifstream fileInput;
  bool nullParameter;
  // ###################


  TCanvas* cEff = new TCanvas("cEff", "cEff", 10, 10, 1200, 800);
  cEff->Divide(3,2);

  if (who == "thetaL")
    {
      fileOutput.open(fileNameOut.replace(fileNameOut.find(".txt"),4,"L.txt").c_str(),ofstream::app);
      if (fileOutput.good() == false)
	{
	  cout << "[ComputeEfficiency::Fit1DEfficiencies]\tError opening file : " << fileNameOut.c_str() << endl;
	  exit (EXIT_FAILURE);
	}


      // ##############################
      // # Read analytical efficiency #
      // ##############################
      effFunc1D = new TF1("effFunc1D",Utility->TellMeEffFuncThetaL().c_str(),-1.0,1.0);


      for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	{
	  Utility->InitEffFuncThetaL(effFunc1D);

	  myString.clear(); myString.str("");
	  myString << "histFit_" << j;
	  histFit.push_back(new TH1D(myString.str().c_str(),myString.str().c_str(),cosThetaLBins->size()-1,cosThetaLBins_));

	  for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	    {
	      Utility->IntegrateEffPhi(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2Bins->operator[](q2BinIndx),cosThetaKBins->operator[](j),cosThetaLBins->operator[](k),myEff,&Eff,&EffErr);
	      
	      histFit.back()->SetBinContent(k+1,Eff);
	      histFit.back()->SetBinError(k+1,EffErr);
	    }


	  // ######################################################################
	  // # Add constraint where it is necessary to bound the function at zero #
	  // ######################################################################
	  Utility->AddConstraintThetaL(&histFit.back(),q2BinIndx,j,q2BinIndx);


	  // #########################################################################
	  // # Perform the fit of the analytical efficiency to the binned efficiency #
	  // #########################################################################
	  histFit.back()->SetMarkerStyle(20);
	  histFit.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
	  histFit.back()->SetYTitle("Efficiency");
	  histFit.back()->GetYaxis()->SetRangeUser(-ordinateRange,ordinateRange);

	  cEff->cd(j+1);
	  fitResults = histFit.back()->Fit("effFunc1D","V");
	  histFit.back()->Draw("pe1");
	  if (fitResults != 0) exit (EXIT_FAILURE);


	  // ################
	  // # Save results #
	  // ################
	  fileOutput << (q2Bins->operator[](q2BinIndx)+q2Bins->operator[](q2BinIndx+1))/2.;
	  for (int i = 0; i < effFunc1D->GetNpar(); i++) fileOutput << "   " << effFunc1D->GetParameter(i) << "   " << effFunc1D->GetParError(i);
	  fileOutput << endl;

	  cout << "\n@@@ Fit for for bin n." << j << " @@@" << endl;
	  cout << "@@@ chi2/DoF = " << effFunc1D->GetChisquare() / effFunc1D->GetNDF() << " (" << effFunc1D->GetChisquare() << "/" << effFunc1D->GetNDF() << ")";
	  cout << "\tCL : " << TMath::Prob(effFunc1D->GetChisquare(),effFunc1D->GetNDF()) << " @@@" << endl;
	  cout << "@@@ Value at " << cosThetaLBins->operator[](0) << " : " << effFunc1D->Eval(cosThetaLBins->operator[](0)) << " @@@" << endl;
	  cout << "@@@ Value at " << cosThetaLBins->operator[](cosThetaLBins->size()-1) << " : " << effFunc1D->Eval(cosThetaLBins->operator[](cosThetaLBins->size()-1)) << " @@@\n" << endl;

	  if (Utility->EffMinValue1D(cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1),effFunc1D) < 0.0)
	    {
	      cout << "NEGATIVE EFFICIENCY !" << endl;
	      if (CHECKnegEFF == true) exit (EXIT_FAILURE);
	    }
	}


      fileOutput.close();
    }
  else if (who == "thetaK")
    {
      fileOutput.open(fileNameOut.replace(fileNameOut.find(".txt"),4,"K.txt").c_str(),ofstream::app);
      if (fileOutput.good() == false)
	{
	  cout << "[ComputeEfficiency::Fit1DEfficiencies]\tError opening file : " << fileNameOut.c_str() << endl;
	  exit (EXIT_FAILURE);
	}


      // ##############################
      // # Read analytical efficiency #
      // ##############################
      effFunc1D = new TF1("effFunc1D",Utility->TellMeEffFuncThetaK().c_str(),-1.0,1.0);


      for (unsigned int k = 0; k < Utility->NcoeffThetaL; k++)
	{
	  Utility->InitEffFuncThetaK(effFunc1D);

	  myString.clear(); myString.str("");
	  myString << "histFit_" << k;
	  histFit.push_back(new TH1D(myString.str().c_str(),myString.str().c_str(),cosThetaKBins->size()-1,cosThetaKBins_));


	  // ###############################
	  // # Read coefficients from file #
	  // ###############################
	  fileInput.open(INPUT_THETAL,ofstream::in);
	  if (fileInput.good() == false)
	    {
	      cout << "[ComputeEfficiency::Fit1DEfficiencies]\tError opening file : " << INPUT_THETAL << endl;
	      exit (EXIT_FAILURE);
	    }
	  for (unsigned int j = 0; j < (cosThetaKBins->size()-1)*q2BinIndx; j++) getline(fileInput,tmpString);
	  nullParameter = true;
	  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
	    {
	      getline(fileInput,tmpString);
	      stringstream rawStringK(tmpString);
	      rawStringK >> coeff; // Discard q2 bin value
	      for (unsigned int i = 0; i < k*2; i++) rawStringK >> coeff;
	      
	      rawStringK >> coeff;
	      histFit.back()->SetBinContent(j+1,coeff);
	      rawStringK >> coeff;
	      histFit.back()->SetBinError(j+1,coeff);

	      cout << "Cos(theta_K) bin " << j << " --> reading coef. " << k << " for var. theta_l: " << histFit.back()->GetBinContent(j+1) << "+/-" <<  histFit.back()->GetBinError(j+1) << endl;

	      if (coeff == 0.0) nullParameter = nullParameter*true; // All parameter errors must be zero in order to ignore the parameter
	      else              nullParameter = false;
	    }


	  // #########################################################################
	  // # Perform the fit of the analytical efficiency to the binned efficiency #
	  // #########################################################################
	  histFit.back()->SetMarkerStyle(20);
	  histFit.back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
	  histFit.back()->SetYTitle("Efficiency");
	  histFit.back()->GetYaxis()->SetRangeUser(-ordinateRange,ordinateRange);

	  cEff->cd(k+1);
	  fitResults = histFit.back()->Fit("effFunc1D","V");
	  histFit.back()->Draw("pe1");
	  if (fitResults != 0) exit (EXIT_FAILURE);


	  // ################
	  // # Save results #
	  // ################
	  fileOutput << (q2Bins->operator[](q2BinIndx)+q2Bins->operator[](q2BinIndx+1))/2.;
	  for (int i = 0; i < effFunc1D->GetNpar(); i++) fileOutput << "   " << effFunc1D->GetParameter(i) << "   " << (nullParameter == false ? effFunc1D->GetParError(i) : 0.0);
	  fileOutput << endl;

	  cout << "\n@@@ Fit for for bin n." << k << " @@@" << endl;
	  cout << "@@@ chi2/DoF = " << effFunc1D->GetChisquare() / effFunc1D->GetNDF() << " (" << effFunc1D->GetChisquare() << "/" << effFunc1D->GetNDF() << ")";
	  cout << "\tCL : " << TMath::Prob(effFunc1D->GetChisquare(),effFunc1D->GetNDF()) << " @@@" << endl;
	  cout << "@@@ Value at " << cosThetaKBins->operator[](0) << " : " << effFunc1D->Eval(cosThetaKBins->operator[](0)) << " @@@" << endl;
	  cout << "@@@ Value at " << cosThetaKBins->operator[](cosThetaKBins->size()-1) << " : " << effFunc1D->Eval(cosThetaKBins->operator[](cosThetaKBins->size()-1)) << " @@@\n" << endl;

	  fileInput.close();
	}


      fileOutput.close();
    }
  else if (who == "phi")
    {
      fileOutput.open(fileNameOut.replace(fileNameOut.find("Theta.txt"),9,"Phi.txt").c_str(),ofstream::app);
      if (fileOutput.good() == false)
	{
	  cout << "[ComputeEfficiency::Fit1DEfficiencies]\tError opening file : " << fileNameOut.c_str() << endl;
	  exit (EXIT_FAILURE);
	}

      
      // ##############################
      // # Read analytical efficiency #
      // ##############################
      effFunc1D = new TF1("effFunc1D",Utility->TellMeEffFuncPhi().c_str(),-Utility->PI,Utility->PI);
      Utility->InitEffFuncPhi(effFunc1D);


      myString.clear(); myString.str("");
      myString << "histFit";
      histFit.push_back(new TH1D(myString.str().c_str(),myString.str().c_str(),phiBins->size()-1,phiBins_));

      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  Utility->IntegrateEffCosThetaKCosThetaL(q2Bins,cosThetaKBins,cosThetaLBins,q2BinIndx,l,myEff,&Eff,&EffErr);

	  histFit.back()->SetBinContent(l+1,Eff);
	  histFit.back()->SetBinError(l+1,EffErr);
	}


      // ##############################
      // # Read analytical efficiency #
      // ##############################
      histFit.back()->SetMarkerStyle(20);
      histFit.back()->SetXTitle("#phi");
      histFit.back()->SetYTitle("Efficiency");
      histFit.back()->GetYaxis()->SetRangeUser(-ordinateRange,ordinateRange);

      cEff->cd(1);
      fitResults = histFit.back()->Fit("effFunc1D","V");
      histFit.back()->Draw("pe1");
      if (fitResults != 0) exit (EXIT_FAILURE);
	  

      // ################
      // # Save results #
      // ################
      fileOutput << (q2Bins->operator[](q2BinIndx)+q2Bins->operator[](q2BinIndx+1))/2.;
      for (int i = 0; i < effFunc1D->GetNpar(); i++) fileOutput << "   " << effFunc1D->GetParameter(i) << "   " << effFunc1D->GetParError(i);
      fileOutput << endl;

      cout << "@@@ Value at " << phiBins->operator[](0) << " : " << effFunc1D->Eval(phiBins->operator[](0)) << " @@@" << endl;
      cout << "@@@ Value at " << phiBins->operator[](phiBins->size()-1) << " : " << effFunc1D->Eval(phiBins->operator[](phiBins->size()-1)) << " @@@\n" << endl;

      if (Utility->EffMinValue1D(phiBins->operator[](0),phiBins->operator[](phiBins->size()-1),effFunc1D) < 0.0)
	{
	  cout << "NEGATIVE EFFICIENCY !" << endl;
	  if (CHECKnegEFF == true) exit (EXIT_FAILURE);
	}
      

      fileOutput.close();
    }


  cEff->Modified();
  cEff->Update();
  histFit.clear();
}


void Fit2DEfficiencies (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			int SignalType, Utils::effStruct myEff, unsigned int q2BinIndx, string fileNameOut, bool savePlot)
{
  // ###################
  // # Local variables #
  // ###################
  TFitResultPtr fitResults;
  stringstream myString;
  vector<TF2*> effFuncs2D;
  // ###################


  TCanvas* cTestGlobalFit = new TCanvas("cTestGlobalFit", "cTestGlobalFit", 10, 10, 700, 500);


  // ##########################
  // # Read binned efficiency #
  // ##########################
  myString.clear(); myString.str("");
  myString << "Histo_q2Bin_" << q2BinIndx;
  TH2D* hisFunc2D = Utility->Get2DEffHistoButPhi(myString.str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2BinIndx,myEff);


  // ##############################
  // # Read analytical efficiency #
  // ##############################
  Utility->ReadAnalyticalEff(INPUT_THETAL_THETAK,q2Bins,cosThetaKBins,cosThetaLBins,&effFuncs2D,"effFuncs2D",0);
  effFuncs2D[q2BinIndx]->SetRange(cosThetaKBins->operator[](0),
				  cosThetaLBins->operator[](0),
				  cosThetaKBins->operator[](cosThetaKBins->size()-1),
				  cosThetaLBins->operator[](cosThetaLBins->size()-1));

  for (int i = 0; i < effFuncs2D[q2BinIndx]->GetNpar(); i++)
    if (effFuncs2D[q2BinIndx]->GetParError(i) == 0.0) effFuncs2D[q2BinIndx]->FixParameter(i,effFuncs2D[q2BinIndx]->GetParameter(i));


  // ############################################################################################
  // # Add constraint along Y (= cosThetaL) where it is necessary to bound the function at zero #
  // ############################################################################################
  Utility->AddConstraintThetaKThetaL(&hisFunc2D,cosThetaKBins,q2BinIndx,SignalType,q2BinIndx);


  // #########################################################################
  // # Perform the fit of the analytical efficiency to the binned efficiency #
  // #########################################################################
  fitResults = hisFunc2D->Fit(effFuncs2D[q2BinIndx]->GetName(),"VMS");
  cout << "[ComputeEfficiency::Fit2DEfficiencies]\tFit status: " << fitResults << endl;

  if (fitResults != 0)
    {
      fitResults = hisFunc2D->Fit(effFuncs2D[q2BinIndx]->GetName(),"VS");
      cout << "[ComputeEfficiency::Fit2DEfficiencies]\tFit status: " << fitResults << endl;
    }

  TMatrixTSym<double> covMatrix(fitResults->GetCovarianceMatrix());

  if (fitResults != 0) exit (EXIT_FAILURE);
  else if (covMatrix.Determinant() <= 0.0)
    {
      cout << "[ComputeEfficiency::Fit2DEfficiencies]\tFit status: covariance matrix with negative determinant" << endl;
      exit (EXIT_FAILURE);
    }
  else if (covMatrix.IsSymmetric() == false)
    {
      cout << "[ComputeEfficiency::Fit2DEfficiencies]\tFit status: covariance matrix not symmetric" << endl;
      exit (EXIT_FAILURE);
    }

  cTestGlobalFit->cd();
  effFuncs2D[q2BinIndx]->Draw("surf1 fb");
  hisFunc2D->Draw("lego2 fb");
  cTestGlobalFit->Modified();
  cTestGlobalFit->Update();

  cout << "@@@ chi2/DoF = " << effFuncs2D[q2BinIndx]->GetChisquare() / effFuncs2D[q2BinIndx]->GetNDF() << " (" << effFuncs2D[q2BinIndx]->GetChisquare() << "/" << effFuncs2D[q2BinIndx]->GetNDF() << ")";
  cout << "\tCL : " << TMath::Prob(effFuncs2D[q2BinIndx]->GetChisquare(),effFuncs2D[q2BinIndx]->GetNDF()) << " @@@" << endl;


  // ################################################
  // # Check if analytical efficiency goes negative #
  // ################################################
  if (Utility->EffMinValue2D(cosThetaKBins,cosThetaLBins,effFuncs2D[q2BinIndx]) < 0.0)
    {
      cout << "NEGATIVE EFFICIENCY !" << endl;
      if (CHECKnegEFF == true) exit (EXIT_FAILURE);
    }
  Utility->SaveAnalyticalEff(fileNameOut.c_str(),effFuncs2D[q2BinIndx],(q2Bins->operator[](q2BinIndx) + q2Bins->operator[](q2BinIndx+1)) / 2.);
  fileNameOut.replace(fileNameOut.find(".txt"),4,"FullCovariance.txt");
  Utility->SaveAnalyticalEffFullCovariance(fileNameOut.c_str(),&covMatrix,(q2Bins->operator[](q2BinIndx) + q2Bins->operator[](q2BinIndx+1)) / 2.);


  // ########################################
  // # Check integrity of covariance matrix #
  // ########################################
  vector<TMatrixTSym<double>*>* covMatrices = new vector<TMatrixTSym<double>*>;
  Utility->ReadAnalyticalEffFullCovariance(fileNameOut.c_str(),covMatrices,"2D",0);


  // #############
  // # Save plot #
  // #############
  if (savePlot == true)
    {
      myString.clear(); myString.str("");
      myString << "Eff2D_q2Bin_" << q2BinIndx << ".pdf";
      cTestGlobalFit->Print(myString.str().c_str());
    }


  effFuncs2D.clear();
  covMatrix.Clear();
  covMatrices->clear();
  delete covMatrices;
}


void Fit3DEfficiencies (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins,
			int SignalType, Utils::effStruct myEff, unsigned int q2BinIndx, string fileNameOut, bool savePlot)
// @TMP@
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  gStyle->SetPadTopMargin(0.02);
  gStyle->SetPadRightMargin(0.04);


  // ###################
  // # Local variables #
  // ###################
  ifstream fileInput;
  TFitResultPtr fitResults;
  stringstream myString;
  string tmpString;
  double coeff;
  TH1* tmpHist1D;
  TH2* tmpHist2D;
  TF1* effFunc1D;
  vector<TF2*> effFuncs2D;
  TF3* effFunc3D;
  TH3D* effHis3D;
  double Zmin, Zmax;
  // #######################
  double parErrorPhi = 1e-4;
  double Yaxes       = 2e-1;
  double Zscale      = 1.1;
  // #######################


  TCanvas* cShow2DAnaEff = new TCanvas("cShow2DAnaEff", "cShow2DAnaEff", 10, 10, 900, 500);
  cShow2DAnaEff->Divide(3,1);

  TCanvas* cTestGlobalFit2D = new TCanvas("cTestGlobalFit2D", "cTestGlobalFit2D", 10, 10, 900, 500);
  cTestGlobalFit2D->Divide(3,1);

  TCanvas* cTestGlobalFit1D = new TCanvas("cTestGlobalFit1D", "cTestGlobalFit11D", 10, 10, 900, 500);
  cTestGlobalFit1D->Divide(3,1);

  
  // ##########################
  // # Read binned efficiency #
  // ##########################
  myString.clear(); myString.str("");
  myString << "Histo_q2Bin_" << q2BinIndx;
  TH3D* hisFunc3D = Utility->Get3DEffHistoq2Bin(myString.str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2BinIndx,myEff);
  Zmin = hisFunc3D->GetMinimum();
 
  cShow2DAnaEff->cd(1);
  hisFunc3D->Draw("ISO");
  cTestGlobalFit2D->cd(1);
  hisFunc3D->Draw("ISO");
  cTestGlobalFit1D->cd(1);
  hisFunc3D->Draw("ISO");


  // #################################
  // # Read analytical efficiency    #
  // #################################
  // # cos(theta_k) and cos(theta_l) #
  // #################################
  Utility->ReadAnalyticalEff(INPUT_THETAL_THETAK,q2Bins,cosThetaKBins,cosThetaLBins,&effFuncs2D,"effFuncs2D",0);

  // #######
  // # phi #
  // #######
  effFunc1D = new TF1("effFunc1D",Utility->TellMeEffFuncPhi().c_str(),Utility->PI,Utility->PI);
  Utility->InitEffFuncPhi(effFunc1D);

  cout << "\n[ComputeEfficiency::Fit3DEfficiencies]\t@@@ Reading coefficients for analytical efficiency for phi from file " << INPUT_PHI << " @@@" << endl;
  fileInput.open(INPUT_PHI,ofstream::in);
  if (fileInput.good() == false)
    {
      cout << "[ComputeEfficiency::Fit3DEfficiencies]\tError opening file : " << INPUT_PHI << endl;
      exit (EXIT_FAILURE);
    }
  for (unsigned int j = 0; j < q2BinIndx; j++) getline(fileInput,tmpString);
  getline(fileInput,tmpString);
  stringstream rawStringK(tmpString);
  rawStringK >> coeff; // Discard q2 bin value
  for (unsigned int l = 0; l < Utility->NcoeffPhi; l++)
    {
      rawStringK >> coeff;
      effFunc1D->SetParameter(l,coeff);
      rawStringK >> coeff;
      effFunc1D->SetParError(l,coeff);

      cout << "Reading coef. " << l << " for var. phi: " << effFunc1D->GetParameter(l) << "+/-" << effFunc1D->GetParError(l) << endl;
    }


  // ########################################################
  // # Putting together cos(theta_k), cos(theta_l), and phi #
  // ########################################################
  effFunc3D = new TF3("effFunc3D",Utility->TellMeEffFuncThetaKThetaLPhi().c_str(),
		      cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),
		      cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1),
		      phiBins->operator[](0)      ,phiBins->operator[](phiBins->size()-1)           );
  effFunc3D->SetMarkerStyle(20);
  effFunc3D->SetMarkerColor(kBlack);

  for (int i = 0; i < effFuncs2D[q2BinIndx]->GetNpar(); i++)
    {
      effFunc3D->SetParameter(i,effFuncs2D[q2BinIndx]->GetParameter(i));
      effFunc3D->SetParError(i,effFuncs2D[q2BinIndx]->GetParError(i));
    }

  for (int i = 0; i < effFunc1D->GetNpar(); i++)
    {
      effFunc3D->SetParameter(effFuncs2D[q2BinIndx]->GetNpar()+i,effFunc1D->GetParameter(i));
      effFunc3D->SetParError(effFuncs2D[q2BinIndx]->GetNpar()+i,effFunc1D->GetParError(i));
    }


  // ##################################################################
  // # Add constraint along X or Y or Z to bound the function at zero #
  // ##################################################################
  Utility->AddConstraintThetaKThetaLPhi(SignalType);
  hisFunc3D->SetMarkerStyle(20);
  hisFunc3D->SetMarkerColor(kBlack);
  hisFunc3D->GetXaxis()->SetTitleOffset(1.35);
  hisFunc3D->GetYaxis()->SetTitleOffset(1.35);
  hisFunc3D->GetZaxis()->SetTitleOffset(1.35);


  // #########################################################################
  // # Perform the fit of the analytical efficiency to the binned efficiency #
  // #########################################################################
  fitResults = hisFunc3D->Fit(effFunc3D->GetName(),"VMS");
  cout << "[ComputeEfficiency::Fit3DEfficiencies]\tFit status: " << fitResults << endl;

  if (fitResults != 0)
    {
      fitResults = hisFunc3D->Fit(effFunc3D->GetName(),"VS");
      cout << "[ComputeEfficiency::Fit3DEfficiencies]\tFit status: " << fitResults << endl;
    }

  if (fitResults != 0)
    {
      for (int i = 1; i < effFunc1D->GetNpar(); i++)
	{
	  effFunc3D->SetParameter(effFuncs2D[q2BinIndx]->GetNpar()+i,0.0);
	  effFunc3D->SetParError(effFuncs2D[q2BinIndx]->GetNpar()+i,parErrorPhi);
	}

      fitResults = hisFunc3D->Fit(effFunc3D->GetName(),"VS");
      cout << "[ComputeEfficiency::Fit3DEfficiencies]\tFit status: " << fitResults << endl;
    }

  TMatrixTSym<double> covMatrix(fitResults->GetCovarianceMatrix());

  if (fitResults != 0) exit (EXIT_FAILURE);
  else if (covMatrix.Determinant() <= 0.0)
    {
      cout << "[ComputeEfficiency::Fit3DEfficiencies]\tFit status: covariance matrix with negative determinant" << endl;
      exit (EXIT_FAILURE);
    }
  else if (covMatrix.IsSymmetric() == false)
    {
      cout << "[ComputeEfficiency::Fit3DEfficiencies]\tFit status: covariance matrix not symmetric" << endl;
      exit (EXIT_FAILURE);
    }

  effHis3D = Utility->Get3DEffHistoq2Bin("effHis3D",q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2BinIndx,myEff);
  effHis3D->SetMarkerStyle(22);
  effHis3D->SetMarkerColor(kRed);
  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
    for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  coeff = effFunc3D->Eval((cosThetaKBins->operator[](j)+cosThetaKBins->operator[](j+1))/2.,
				  (cosThetaLBins->operator[](k)+cosThetaLBins->operator[](k+1))/2.,
				  (phiBins->operator[](l)+phiBins->operator[](l+1))/2.);
	  effHis3D->SetBinContent(j+1,k+1,l+1,coeff);
	  effHis3D->SetBinError(j+1,k+1,l+1,0.0);
	}

  // #######################
  // # Make 2D projections #
  // #######################
  cTestGlobalFit2D->cd(1);
  tmpHist2D = static_cast<TH2D*>(hisFunc3D->Project3D("xy"));
  Zmax = tmpHist2D->GetMaximum()*Zscale;
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("lego2 fb");

  tmpHist2D = static_cast<TH2D*>(effHis3D->Project3D("xy"));
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("same surf fb");

  cShow2DAnaEff->cd(1);
  tmpHist2D = static_cast<TH2D*>(effHis3D->Project3D("xy"));
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("surf1 fb");

  cTestGlobalFit2D->cd(2);
  tmpHist2D = static_cast<TH2D*>(hisFunc3D->Project3D("xz"));
  Zmax = tmpHist2D->GetMaximum()*Zscale;
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("lego2 fb");

  tmpHist2D = static_cast<TH2D*>(effHis3D->Project3D("xz"));
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("same surf fb");

  cShow2DAnaEff->cd(2);
  tmpHist2D = static_cast<TH2D*>(effHis3D->Project3D("xz"));
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("surf1 fb");

  cTestGlobalFit2D->cd(3);
  tmpHist2D = static_cast<TH2D*>(hisFunc3D->Project3D("yz"));
  Zmax = tmpHist2D->GetMaximum()*Zscale;
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("lego2 fb");

  tmpHist2D = static_cast<TH2D*>(effHis3D->Project3D("yz"));
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("same surf fb");

  cShow2DAnaEff->cd(3);
  tmpHist2D = static_cast<TH2D*>(effHis3D->Project3D("yz"));
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("surf1 fb");

  cTestGlobalFit2D->Modified();
  cTestGlobalFit2D->Update();
  cShow2DAnaEff->Modified();
  cShow2DAnaEff->Update();

  // #######################
  // # Make 1D projections #
  // #######################
  cTestGlobalFit1D->cd(1);
  tmpHist1D = hisFunc3D->Project3D("x");
  tmpHist1D->GetYaxis()->SetRangeUser(0.0,Yaxes);
  tmpHist1D->Draw("pe1");
  effHis3D->Project3D("x")->Draw("same pe1");

  cTestGlobalFit1D->cd(2);
  tmpHist1D = hisFunc3D->Project3D("y");
  tmpHist1D->GetYaxis()->SetRangeUser(0.0,Yaxes);
  tmpHist1D->Draw("pe1");
  effHis3D->Project3D("y")->Draw("same pe1");

  cTestGlobalFit1D->cd(3);
  tmpHist1D = hisFunc3D->Project3D("z");
  tmpHist1D->GetYaxis()->SetRangeUser(0.0,Yaxes);
  tmpHist1D->Draw("pe1");
  effHis3D->Project3D("z")->Draw("same pe1");

  TLegend* leg = new TLegend(0.2, 0.8, 0.85, 0.95, "");
  leg->AddEntry(hisFunc3D,"Binned efficiency");
  leg->AddEntry(effHis3D,"Analytical efficiency");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  cTestGlobalFit1D->Modified();
  cTestGlobalFit1D->Update();

  cout << "@@@ chi2/DoF = " << fitResults->Chi2() / fitResults->Ndf() << " (" << fitResults->Chi2() << "/" << fitResults->Ndf() << ")";
  cout << "\tCL : " << TMath::Prob(fitResults->Chi2(),fitResults->Ndf()) << " @@@" << endl;


  // ################################################
  // # Check if analytical efficiency goes negative #
  // ################################################
  if (Utility->EffMinValue3D(cosThetaKBins,cosThetaLBins,phiBins,effFunc3D) < 0.0)
    {
      cout << "NEGATIVE EFFICIENCY !" << endl;
      if (CHECKnegEFF == true) exit (EXIT_FAILURE);
    }
  Utility->SaveAnalyticalEff(fileNameOut.c_str(),effFunc3D,(q2Bins->operator[](q2BinIndx) + q2Bins->operator[](q2BinIndx+1)) / 2.);
  fileNameOut.replace(fileNameOut.find(".txt"),4,"FullCovariance.txt");
  Utility->SaveAnalyticalEffFullCovariance(fileNameOut.c_str(),&covMatrix,(q2Bins->operator[](q2BinIndx) + q2Bins->operator[](q2BinIndx+1)) / 2.);
  

  // ########################################
  // # Check integrity of covariance matrix #
  // ########################################
  vector<TMatrixTSym<double>*>* covMatrices = new vector<TMatrixTSym<double>*>;
  Utility->ReadAnalyticalEffFullCovariance(fileNameOut.c_str(),covMatrices,"3D",0);


  // #############
  // # Save plot #
  // #############
  if (savePlot == true)
    {
      myString.clear(); myString.str("");
      myString << "Eff1D_q2Bin_" << q2BinIndx << ".pdf";
      cTestGlobalFit1D->Print(myString.str().c_str());

      myString.clear(); myString.str("");
      myString << "Eff2D_q2Bin_" << q2BinIndx << ".pdf";
      cTestGlobalFit2D->Print(myString.str().c_str());

      myString.clear(); myString.str("");
      myString << "Eff2Danaly_q2Bin_" << q2BinIndx << ".pdf";
      cShow2DAnaEff->Print(myString.str().c_str());
    }


  covMatrix.Clear();
  covMatrices->clear();
  delete covMatrices;
  delete effHis3D;
}


void Test2DEfficiency (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, Utils::effStruct myEff, unsigned int q2BinIndx, int SignalType, string analyORbin, bool savePlot)
{

  // ###################
  // # Local variables #
  // ###################
  stringstream myString;
  vector<TF2*> effFuncs2D[3];
  vector<TF12*> effFuncSlice[3];
  vector<TH1D*> histoSliceOrg[3];
  vector<TH1D*> histoSliceNew[3];
  TF2* effFunc2D[3]     = {NULL,NULL,NULL};
  TH2D* hisFunc2Dnew[3] = {NULL,NULL,NULL};
  TH2D* hisFunc2Dorg[3];
  //TH3D* hisFunc3D;
  // ###################

  /*myString.clear(); myString.str("");
  myString << "H3Deff_q2Bin_" << q2BinIndx;
  hisFunc3D = Utility->Get3DEffHistoq2Bin (myString.str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2BinIndx,myEff);
  TFile* effHistoFile1=TFile::Open("effHisto3DFile.root","UPDATE");
  hisFunc3D->Write();
  effHistoFile1->Close();*/

  TCanvas* cEff[3];
  int cv_cnt[3] = {3,3,3};

  TCanvas* cEffComp = new TCanvas("cEffComp", "cEffComp");
  cEffComp->Divide(2,1);

  // ################################
  // # Loop over couple of variable #
  // ################################
  for (int coup=0; coup<3; coup++)
    {
      cEff[coup] = new TCanvas(Form("cEff%i",coup), Form("cEff%i",coup), 10, 10, 1500, 900);
      cEff[coup]->Divide(5,3);


      // ##########################
      // # Read binned efficiency #
      // ##########################
      myString.clear(); myString.str("");
      myString << Form("H2Deff_%i_q2Bin_Org_",coup) << q2BinIndx;
      if (coup==0) hisFunc2Dorg[coup] = Utility->Get2DEffHistoButPhi(myString.str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2BinIndx,myEff);
      else if (coup==1) hisFunc2Dorg[coup] = Utility->Get2DEffHistoButCosThetaK(myString.str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2BinIndx,myEff);
      else hisFunc2Dorg[coup] = Utility->Get2DEffHistoButCosThetaL(myString.str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2BinIndx,myEff);

      cEff[coup]->cd(1);
      hisFunc2Dorg[coup]->GetXaxis()->SetLabelSize(0.05);
      hisFunc2Dorg[coup]->GetYaxis()->SetLabelSize(0.05);
      hisFunc2Dorg[coup]->GetZaxis()->SetLabelSize(0.05);
      hisFunc2Dorg[coup]->GetXaxis()->SetTitleSize(0.06);
      hisFunc2Dorg[coup]->GetYaxis()->SetTitleSize(0.06);
      hisFunc2Dorg[coup]->GetXaxis()->SetTitleOffset(1.25);
      hisFunc2Dorg[coup]->GetYaxis()->SetTitleOffset(1.25);
      hisFunc2Dorg[coup]->Draw("lego2 fb");

      // These are the 2D histogram of the efficiency. Try to save them somewhere
      TFile* effHistoFile=TFile::Open("effHisto2DFile.root","UPDATE");
      hisFunc2Dorg[coup]->Write();
      effHistoFile->Close();


      // #######################################################
      // # Read analytical OR binned (interpolated) efficiency #
      // #######################################################
      if (coup == 0 && strcmp(analyORbin.c_str(),"ANALY") == 0)
	{
	  Utility->ReadAnalyticalEff(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),q2Bins,cosThetaKBins,cosThetaLBins,&effFuncs2D[coup],"effFuncs2D",Utility->ParFileBlockN("analyEffokTag"));
	  effFunc2D[coup] = effFuncs2D[coup][q2BinIndx];
	  cEff[coup]->cd(2);
	  effFunc2D[coup]->GetXaxis()->SetLabelSize(0.05);
	  effFunc2D[coup]->GetYaxis()->SetLabelSize(0.05);
	  effFunc2D[coup]->GetZaxis()->SetLabelSize(0.05);
	  effFunc2D[coup]->GetXaxis()->SetTitleSize(0.06);
	  effFunc2D[coup]->GetYaxis()->SetTitleSize(0.06);
	  effFunc2D[coup]->GetXaxis()->SetTitleOffset(1.25);
	  effFunc2D[coup]->GetYaxis()->SetTitleOffset(1.25);
	  effFunc2D[coup]->GetZaxis()->SetTitleOffset(1.25);
	  effFunc2D[coup]->Draw("surf2 fb");
	  Utility->EffMinValue2D(cosThetaKBins,cosThetaLBins,effFunc2D[coup]);
	}
      else if (coup == 0 && strcmp(analyORbin.c_str(),"BIN") == 0)
	{
	  myString.clear(); myString.str("");
	  myString << "H2Deff_q2Bin_new_" << q2BinIndx;
	  hisFunc2Dnew[coup] = Utility->Get2DEffHistoq2Bin(cosThetaKBins,cosThetaLBins,q2BinIndx,SignalType,true,make_pair(-1.0,1.0),make_pair(-1.0,1.0));
	  cEff[coup]->cd(2);
	  hisFunc2Dnew[coup]->GetXaxis()->SetLabelSize(0.05);
	  hisFunc2Dnew[coup]->GetYaxis()->SetLabelSize(0.05);
	  hisFunc2Dnew[coup]->GetZaxis()->SetLabelSize(0.05);
	  hisFunc2Dnew[coup]->GetXaxis()->SetTitleSize(0.06);
	  hisFunc2Dnew[coup]->GetYaxis()->SetTitleSize(0.06);
	  hisFunc2Dnew[coup]->GetXaxis()->SetTitleOffset(1.25);
	  hisFunc2Dnew[coup]->GetYaxis()->SetTitleOffset(1.25);
	  hisFunc2Dnew[coup]->GetZaxis()->SetTitleOffset(1.25);
	  hisFunc2Dnew[coup]->Draw("surf2 fb");
	}
      else if (coup == 0)
	{
	  cout << "[ComputeEfficiency::Test2DEfficiency]\tWrong parameter option : "  << analyORbin.c_str() << endl;
	  exit (EXIT_FAILURE);
	}
    }

  // ##############################
  // # Show slices of cos(thetaK) #
  // ##############################
  for (int coup=0; coup<3; coup++)
    {
      if (coup == 1) continue;

      for (unsigned int binIndx = 0; binIndx < cosThetaKBins->size()-1; binIndx++)
	{
	  cEff[coup]->cd(cv_cnt[coup]);
	  cv_cnt[coup]++;

	  /*if (strcmp(analyORbin.c_str(),"ANALY") == 0)
	    {
	      myString.clear(); myString.str("");
	      myString << "effFuncSlice_binK_" << binIndx;
	      effFuncSlice[coup].push_back(new TF12(myString.str().c_str(),effFunc2D[coup],(cosThetaKBins->operator[](binIndx)+cosThetaKBins->operator[](binIndx+1))/2.,"y"));
	      effFuncSlice[coup].back()->GetXaxis()->SetTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
	      effFuncSlice[coup].back()->GetYaxis()->SetTitle("Efficiency");
	      effFuncSlice[coup].back()->GetXaxis()->SetLabelSize(0.05);
	      effFuncSlice[coup].back()->GetYaxis()->SetLabelSize(0.05);
	      effFuncSlice[coup].back()->GetXaxis()->SetTitleSize(0.06);
	      effFuncSlice[coup].back()->SetLineWidth(2);
	      effFuncSlice[coup].back()->SetLineColor(kRed);
	      effFuncSlice[coup].back()->Draw();
	    }
	  else
	    {
	      int binx, biny, binz;
	      myString.clear(); myString.str("");
	      myString << "histoSliceNew_binK_" << binIndx;
	      hisFunc2Dnew[coup]->GetBinXYZ(hisFunc2Dnew[coup]->FindBin((cosThetaKBins->operator[](binIndx)+cosThetaKBins->operator[](binIndx+1))/2.,0.),binx,biny,binz);
	      histoSliceNew[coup].push_back(hisFunc2Dnew[coup]->ProjectionY(myString.str().c_str(),binx,binx));
	      histoSliceNew[coup].back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
	      histoSliceNew[coup].back()->SetYTitle("Efficiency");
	      histoSliceNew[coup].back()->GetXaxis()->SetLabelSize(0.05);
	      histoSliceNew[coup].back()->GetYaxis()->SetLabelSize(0.05);
	      histoSliceNew[coup].back()->GetXaxis()->SetTitleSize(0.06);
	      histoSliceNew[coup].back()->GetXaxis()->SetTitleOffset(0.7);
	      histoSliceNew[coup].back()->SetLineWidth(2);
	      histoSliceNew[coup].back()->SetLineColor(kRed);
	      histoSliceNew[coup].back()->Draw("hist");
	      }*/

	  myString.clear(); myString.str("");
	  myString << "histoSliceOrg_binK_" << coup << "_" << binIndx;
	  histoSliceOrg[coup].push_back(hisFunc2Dorg[coup]->ProjectionY(myString.str().c_str(),binIndx+1,binIndx+1));
	  if (coup == 0) histoSliceOrg[coup].back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
	  else histoSliceOrg[coup].back()->SetXTitle("#phi");
	  histoSliceOrg[coup].back()->GetXaxis()->SetLabelSize(0.05);
	  histoSliceOrg[coup].back()->GetYaxis()->SetLabelSize(0.05);
	  histoSliceOrg[coup].back()->GetXaxis()->SetTitleSize(0.06);
	  histoSliceOrg[coup].back()->GetXaxis()->SetTitleOffset(0.7);
	  histoSliceOrg[coup].back()->SetYTitle("Efficiency");
	  histoSliceOrg[coup].back()->SetMarkerStyle(20);
	  histoSliceOrg[coup].back()->Draw("pe1");
	  if (coup == 2) {
	    TString opt = "lpe1";
	    if (binIndx > 0) opt = "same lpe1";
	    cEffComp->cd(1);
	    histoSliceOrg[coup].back()->SetLineColor(binIndx+1);
	    histoSliceOrg[coup].back()->SetMarkerColor(binIndx+1);
	    histoSliceOrg[coup].back()->DrawNormalized(opt);
	    histoSliceOrg[coup].back()->SetLineColor(1);
	    histoSliceOrg[coup].back()->SetMarkerColor(1);
	  }
	}
    }

  // ##############################
  // # Show slices of cos(thetaL) #
  // ##############################
  for (int coup=0; coup<2; coup++)
    {
      for (unsigned int binIndx = 0; binIndx < cosThetaLBins->size()-1; binIndx++)
	{
	  cEff[coup]->cd(cv_cnt[coup]);
	  cv_cnt[coup]++;

	  /*if (strcmp(analyORbin.c_str(),"ANALY") == 0)
	    {
	      myString.clear(); myString.str("");
	      myString << "effFuncSlice_binL_" << binIndx;
	      effFuncSlice[coup].push_back(new TF12(myString.str().c_str(),effFunc2D[coup],(cosThetaLBins->operator[](binIndx)+cosThetaLBins->operator[](binIndx+1))/2.,"x"));
	      effFuncSlice[coup].back()->GetXaxis()->SetTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
	      effFuncSlice[coup].back()->GetYaxis()->SetTitle("Efficiency");
	      effFuncSlice[coup].back()->GetXaxis()->SetLabelSize(0.05);
	      effFuncSlice[coup].back()->GetYaxis()->SetLabelSize(0.05);
	      effFuncSlice[coup].back()->GetXaxis()->SetTitleSize(0.06);
	      effFuncSlice[coup].back()->SetLineWidth(2);
	      effFuncSlice[coup].back()->SetLineColor(kRed);
	      effFuncSlice[coup].back()->Draw();
	    }
	  else
	    {
	      int binx, biny, binz;
	      myString.clear(); myString.str("");
	      myString << "histoSliceNew_binL_" << binIndx;
	      hisFunc2Dnew[coup]->GetBinXYZ(hisFunc2Dnew[coup]->FindBin(0.,(cosThetaLBins->operator[](binIndx)+cosThetaLBins->operator[](binIndx+1))/2.),binx,biny,binz);
	      histoSliceNew[coup].push_back(hisFunc2Dnew[coup]->ProjectionX(myString.str().c_str(),biny,biny));
	      histoSliceNew[coup].back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
	      histoSliceNew[coup].back()->SetYTitle("Efficiency");
	      histoSliceNew[coup].back()->GetXaxis()->SetLabelSize(0.05);
	      histoSliceNew[coup].back()->GetYaxis()->SetLabelSize(0.05);
	      histoSliceNew[coup].back()->GetXaxis()->SetTitleSize(0.06);
	      histoSliceNew[coup].back()->GetXaxis()->SetTitleOffset(0.7);
	      histoSliceNew[coup].back()->SetLineWidth(2);
	      histoSliceNew[coup].back()->SetLineColor(kRed);
	      histoSliceNew[coup].back()->Draw("hist");
	      }*/

	  myString.clear(); myString.str("");
	  myString << "histoSliceOrg_binL_" << coup << "_" << binIndx;
	  if (coup == 0) histoSliceOrg[coup].push_back(hisFunc2Dorg[coup]->ProjectionX(myString.str().c_str(),binIndx+1,binIndx+1));
	  else histoSliceOrg[coup].push_back(hisFunc2Dorg[coup]->ProjectionY(myString.str().c_str(),binIndx+1,binIndx+1));
	  if (coup == 0) histoSliceOrg[coup].back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
	  else histoSliceOrg[coup].back()->SetXTitle("#phi");
	  histoSliceOrg[coup].back()->GetXaxis()->SetLabelSize(0.05);
	  histoSliceOrg[coup].back()->GetYaxis()->SetLabelSize(0.05);
	  histoSliceOrg[coup].back()->GetXaxis()->SetTitleSize(0.06);
	  histoSliceOrg[coup].back()->GetXaxis()->SetTitleOffset(0.7);
	  histoSliceOrg[coup].back()->SetMarkerStyle(20);
	  histoSliceOrg[coup].back()->SetYTitle("Efficiency");
	  histoSliceOrg[coup].back()->Draw("pe1");
	  if (coup == 1) {
	    TString opt = "lpe1";
	    if (binIndx > 0) opt = "same lpe1";
	    cEffComp->cd(2);
	    histoSliceOrg[coup].back()->SetLineColor(binIndx+1);
	    histoSliceOrg[coup].back()->SetMarkerColor(binIndx+1);
	    histoSliceOrg[coup].back()->DrawNormalized(opt);
	    histoSliceOrg[coup].back()->SetLineColor(1);
	    histoSliceOrg[coup].back()->SetMarkerColor(1);
	  }
	}
    }
      
  // ######################
  // # Show slices of phi #
  // ######################
  for (int coup=1; coup<3; coup++)
    {
      for (unsigned int binIndx = 0; binIndx < phiBins->size()-1; binIndx++)
	{
	  cEff[coup]->cd(cv_cnt[coup]);
	  cv_cnt[coup]++;

	  /*if (strcmp(analyORbin.c_str(),"ANALY") == 0)
	    {
	      myString.clear(); myString.str("");
	      myString << "effFuncSlice_binL_" << binIndx;
	      effFuncSlice[coup].push_back(new TF12(myString.str().c_str(),effFunc2D[coup],(cosThetaLBins->operator[](binIndx)+cosThetaLBins->operator[](binIndx+1))/2.,"x"));
	      effFuncSlice[coup].back()->GetXaxis()->SetTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
	      effFuncSlice[coup].back()->GetYaxis()->SetTitle("Efficiency");
	      effFuncSlice[coup].back()->GetXaxis()->SetLabelSize(0.05);
	      effFuncSlice[coup].back()->GetYaxis()->SetLabelSize(0.05);
	      effFuncSlice[coup].back()->GetXaxis()->SetTitleSize(0.06);
	      effFuncSlice[coup].back()->SetLineWidth(2);
	      effFuncSlice[coup].back()->SetLineColor(kRed);
	      effFuncSlice[coup].back()->Draw();
	    }
	  else
	    {
	      int binx, biny, binz;
	      myString.clear(); myString.str("");
	      myString << "histoSliceNew_binL_" << binIndx;
	      hisFunc2Dnew[coup]->GetBinXYZ(hisFunc2Dnew[coup]->FindBin(0.,(cosThetaLBins->operator[](binIndx)+cosThetaLBins->operator[](binIndx+1))/2.),binx,biny,binz);
	      histoSliceNew[coup].push_back(hisFunc2Dnew[coup]->ProjectionX(myString.str().c_str(),biny,biny));
	      histoSliceNew[coup].back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
	      histoSliceNew[coup].back()->SetYTitle("Efficiency");
	      histoSliceNew[coup].back()->GetXaxis()->SetLabelSize(0.05);
	      histoSliceNew[coup].back()->GetYaxis()->SetLabelSize(0.05);
	      histoSliceNew[coup].back()->GetXaxis()->SetTitleSize(0.06);
	      histoSliceNew[coup].back()->GetXaxis()->SetTitleOffset(0.7);
	      histoSliceNew[coup].back()->SetLineWidth(2);
	      histoSliceNew[coup].back()->SetLineColor(kRed);
	      histoSliceNew[coup].back()->Draw("hist");
	      }*/

	  myString.clear(); myString.str("");
	  myString << "histoSliceOrg_binP_" << coup << "_" << binIndx;
	  histoSliceOrg[coup].push_back(hisFunc2Dorg[coup]->ProjectionX(myString.str().c_str(),binIndx+1,binIndx+1));
	  if (coup == 2) histoSliceOrg[coup].back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
	  else histoSliceOrg[coup].back()->SetXTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
	  histoSliceOrg[coup].back()->GetXaxis()->SetLabelSize(0.05);
	  histoSliceOrg[coup].back()->GetYaxis()->SetLabelSize(0.05);
	  histoSliceOrg[coup].back()->GetXaxis()->SetTitleSize(0.06);
	  histoSliceOrg[coup].back()->GetXaxis()->SetTitleOffset(0.7);
	  histoSliceOrg[coup].back()->SetMarkerStyle(20);
	  histoSliceOrg[coup].back()->SetYTitle("Efficiency");
	  histoSliceOrg[coup].back()->Draw("pe1");
	}
    }
	  
  for (int coup=0; coup<3; coup++)
    {
      cEff[coup]->cd(1);
      hisFunc2Dorg[coup]->Draw("lego2 fb");
      cEff[coup]->Modified();
      cEff[coup]->Update();
      


      // #############
      // # Save plot #
      // #############
      if (savePlot == true || 1)
	{
	  myString.clear(); myString.str("");
	  myString << "Test2DEff" << "_" << q2BinIndx << "_" << coup << ".pdf";
	  cEff[coup]->Print(myString.str().c_str());
	  myString.clear(); myString.str("");
	  myString << "Test2DEffComp" << "_" << q2BinIndx << "_" << coup << ".pdf";
	  if (coup==0) cEffComp->Print(myString.str().c_str());
	}
    }
}


void Test3DEfficiency (vector<double>* q2Bins, vector<double>* cosThetaKBins, vector<double>* cosThetaLBins, vector<double>* phiBins, Utils::effStruct myEff, unsigned int q2BinIndx, bool savePlot)
// @TMP@
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  gStyle->SetPadRightMargin(0.04);
  gStyle->SetPadBottomMargin(0.12);


  // ###################
  // # Local variables #
  // ###################
  stringstream myString;
  double coeff;
  TH1* tmpHist1D;
  TH2* tmpHist2D;
  vector<TF3*> effFuncs3D;
  TF3 *EffFunc3D;
  TH3D* effHis3D;
  double Zmin, Zmax;
  // ##################
  double Yaxes  = 2e-1;
  double Zscale =  1.1;
  // ##################


  TCanvas* cEff = new TCanvas("cEff", "cEff", 10, 10, 1000, 1200);
  cEff->Divide(3,3);


  // ##########################
  // # Read binned efficiency #
  // ##########################
  myString.clear(); myString.str("");
  myString << "hisFunc3D_" << q2BinIndx;
  TH3D* hisFunc3D = Utility->Get3DEffHistoq2Bin(myString.str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2BinIndx,myEff);
  Zmin = hisFunc3D->GetMinimum();
  hisFunc3D->SetMarkerStyle(20);
  hisFunc3D->SetMarkerColor(kBlack);
  hisFunc3D->GetXaxis()->SetTitleOffset(1.4);
  hisFunc3D->GetYaxis()->SetTitleOffset(1.4);
  hisFunc3D->GetZaxis()->SetTitleOffset(1.4);

 
  // ##############################
  // # Read analytical efficiency #
  // ##############################
  // Utility->ReadAnalyticalEff(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),q2Bins,cosThetaKBins,cosThetaLBins,phiBins,&effFuncs3D,"effFuncs3D",Utility->ParFileBlockN("analyEffokTag"));
  // EffFunc3D = effFuncs3D[q2BinIndx];
  // Utility->EffMinValue3D(cosThetaKBins,cosThetaLBins,phiBins,EffFunc3D);

  effHis3D = Utility->Get3DEffHistoq2Bin("effHis3D",q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2BinIndx,myEff);
  effHis3D->SetMarkerStyle(22);
  effHis3D->SetMarkerColor(kRed);
  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
    for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
        // coeff = EffFunc3D->Eval((cosThetaKBins->operator[](j)+cosThetaKBins->operator[](j+1))/2.,
        //                   (cosThetaLBins->operator[](k)+cosThetaLBins->operator[](k+1))/2.,
        //                   (phiBins->operator[](l)+phiBins->operator[](l+1))/2.);
	  
	  effHis3D->SetBinContent(j+1,k+1,l+1,coeff);
	  effHis3D->SetBinError(j+1,k+1,l+1,0.0);
	}


  // #######################
  // # Make 2D projections #
  // #######################
  cEff->cd(1);
  tmpHist2D = static_cast<TH2D*>(hisFunc3D->Project3D("xy"));
  Zmax = tmpHist2D->GetMaximum()*Zscale;
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("lego2 fb");
  tmpHist2D = static_cast<TH2D*>(effHis3D->Project3D("xy"));
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("same surf fb");

  cEff->cd(4);
  tmpHist2D->Draw("surf2 fb");

  cEff->cd(2);
  tmpHist2D = static_cast<TH2D*>(hisFunc3D->Project3D("xz"));
  Zmax = tmpHist2D->GetMaximum()*Zscale;
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("lego2 fb");
  tmpHist2D = static_cast<TH2D*>(effHis3D->Project3D("xz"));
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("same surf fb");

  cEff->cd(5);
  tmpHist2D->Draw("surf2 fb");

  cEff->cd(3);
  tmpHist2D = static_cast<TH2D*>(hisFunc3D->Project3D("yz"));
  Zmax = tmpHist2D->GetMaximum()*Zscale;
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("lego2 fb");
  tmpHist2D = static_cast<TH2D*>(effHis3D->Project3D("yz"));
  tmpHist2D->GetZaxis()->SetRangeUser(Zmin,Zmax);
  tmpHist2D->Draw("same surf fb");

  cEff->cd(6);
  tmpHist2D->Draw("surf2 fb");


  // #######################
  // # Make 1D projections #
  // #######################
  cEff->cd(7);
  tmpHist1D = hisFunc3D->Project3D("x");
  tmpHist1D->GetYaxis()->SetRangeUser(0.0,Yaxes);
  tmpHist1D->Draw("pe1");
  effHis3D->Project3D("x")->Draw("same pe1");

  cEff->cd(8);
  tmpHist1D = hisFunc3D->Project3D("y");
  tmpHist1D->GetYaxis()->SetRangeUser(0.0,Yaxes);
  tmpHist1D->Draw("pe1");
  effHis3D->Project3D("y")->Draw("same pe1");

  cEff->cd(9);
  tmpHist1D = hisFunc3D->Project3D("z");
  tmpHist1D->GetYaxis()->SetRangeUser(0.0,Yaxes);
  tmpHist1D->Draw("pe1");
  effHis3D->Project3D("z")->Draw("same pe1");


  cEff->Modified();
  cEff->Update();


  // #############
  // # Save plot #
  // #############
  if (savePlot == true)
    {
      myString.clear(); myString.str("");
      myString << "Test2DEff" << "_" << q2BinIndx << ".pdf";
      cEff->Print(myString.str().c_str());
    }

  // Save histos
  TFile* outfile=TFile::Open("effHisto3DFile.root","UPDATE");
  hisFunc3D->Write(Form("hisFunc3D_%d",q2BinIndx));
  effHis3D->Write(Form("effHis3D_%d",q2BinIndx));
  outfile->Close();

  delete effHis3D;
}


int main (int argc, char** argv)
{
  if (argc >= 2)
    {
      string option = argv[1];


      cout << "\n[ComputeEfficiency::main]\t@@@ Settings @@@" << endl;
      cout << "INPUT_THETAL: "           << INPUT_THETAL << endl;
      cout << "INPUT_PHI: "              << INPUT_PHI << endl;
      cout << "INPUT_THETAL_THETAK: "    << INPUT_THETAL_THETAK << endl;

      cout << "\nRIGHTtag: "      << RIGHTtag << endl;
      cout << "SAVEPLOT: "        << SAVEPLOT << endl;
      cout << "CHECKnegEFF: "     << CHECKnegEFF << endl;
      cout << "EFFis2Dnot3D: "    << EFFis2Dnot3D << endl;
      cout << "NFILES: "          << NFILES << endl;
      cout << "GENEFF: "          << GENEFF << endl;
      cout << "SETBATCH: "        << SETBATCH << endl;
      cout << "PARAMETERFILEIN: " << PARAMETERFILEIN << endl;
      cout << "ordinateRange: "   << ordinateRange << endl;


      // ##########################
      // # Set histo layout style #
      // ##########################
      SetStyle();
      gStyle->SetOptStat(0);


      if ((option == "Make") && (argc == 7))
	{
	  string fileNameGenCandidatesNoFilter = argv[2];
	  string fileNameRecoCandidates        = argv[3];
	  string fileNameSingleCand            = argv[4];
	  string fileNameOutput                = argv[5];
	  string SignalType                    = argv[6];

	  
	  if (SETBATCH == true)
	    {
	      cout << "\n[ComputeEfficiency::main]\t@@@ Setting batch mode @@@" << endl;
	      gROOT->SetBatch(true);
	    }
	  TApplication theApp ("Applications", &argc, argv);

      
	  TFile* NtplFileGenCandidatesNoFilter = TFile::Open(fileNameGenCandidatesNoFilter.c_str(), "READ");
	  theTreeGenCandidatesNoFilter         = (TTree*) NtplFileGenCandidatesNoFilter->Get("B0KstMuMu/B0KstMuMuNTuple");
	  NTupleGenCandidatesNoFilter          = new B0KstMuMuSingleCandTreeContent();
	  NTupleGenCandidatesNoFilter->Init();

	  TFile* NtplFileRecoCandidates = TFile::Open(fileNameRecoCandidates.c_str(), "READ");
	  theTreeRecoCandidates         = (TTree*) NtplFileRecoCandidates->Get("B0KstMuMu/B0KstMuMuNTuple");
	  NTupleRecoCandidates          = new B0KstMuMuSingleCandTreeContent();
	  NTupleRecoCandidates->Init();

	  TFile* NtplFileSingleCand = TFile::Open(fileNameSingleCand.c_str(), "READ");
	  theTreeSingleCand         = (TTree*) NtplFileSingleCand->Get("B0KstMuMu/B0KstMuMuNTuple");
	  NTupleSingleCand          = new B0KstMuMuSingleCandTreeContent();
	  NTupleSingleCand->Init();


	  Utility = new Utils(RIGHTtag);
	  if (Utility->RIGHTflavorTAG == true) Utility->ReadAllBins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"goodTag");
	  else                                 Utility->ReadAllBins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"misTag");
	  Utility->ReadPreselectionCut(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str());
	  Utility->ReadSelectionCuts(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str());
	  Utility->ReadGenericParam(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str());


	  // ###################################
	  // # Initialize efficiency structure #
	  // ###################################
	  Utility->GenEfficiency(&myEff,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins);
	  Utils::effValue myEffVal;
	  Utility->ResetEffValue(&myEffVal,0.0);
	  Utility->InitEfficiency(myEffVal,myEff,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins);


	  ComputeEfficiency(theTreeGenCandidatesNoFilter,NTupleGenCandidatesNoFilter,myEff.Den1,myEff.Err2PoisDen1,myEff.Err2WeigDen1, 1 ,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,atoi(SignalType.c_str()));
	  ComputeEfficiency(theTreeGenCandidatesNoFilter,NTupleGenCandidatesNoFilter,myEff.Num1,myEff.Err2PoisNum1,myEff.Err2WeigNum1, 2 ,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,atoi(SignalType.c_str()));
	  ComputeEfficiency(theTreeRecoCandidates,       NTupleRecoCandidates,       myEff.Den2,myEff.Err2PoisDen2,myEff.Err2WeigDen2, 3 ,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,atoi(SignalType.c_str()));
	  ComputeEfficiency(theTreeSingleCand,           NTupleSingleCand,           myEff.Num2,myEff.Err2PoisNum2,myEff.Err2WeigNum2, 4 ,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,atoi(SignalType.c_str()));


	  Utility->SaveEfficiency(fileNameOutput.c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,myEff);
	  if (SETBATCH == false) MakeHistogramsAllBins(&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,myEff,-1);
	  cout << "\n[ComputeEfficiency::main]\t@@@ Efficiency computation is done @@@" << endl;


	  delete Utility;
	  if (SETBATCH == false) theApp.Run (); // Eventloop on air
	  return EXIT_SUCCESS;
	}
      else if (((option == "ReadBin") || (option == "ReadGenBin")) && (argc >= 4))
	{
	  string fileNameInput = argv[2];
	  string SignalType    = argv[3];
	  int specBin = -1;
	  if (argc == 5) specBin = atoi(argv[4]);

	  
	  if (SETBATCH == true)
	    {
	      cout << "\n[ComputeEfficiency::main]\t@@@ Setting batch mode @@@" << endl;
	      gROOT->SetBatch(true);
	    }
	  TApplication theApp ("Applications", &argc, argv);

    
	  Utility = new Utils(RIGHTtag);
	  if (Utility->RIGHTflavorTAG == true) Utility->ReadAllBins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"goodTag");
	  else                                 Utility->ReadAllBins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"misTag");
	  if (option == "ReadBin") ReadEfficiencies(true, &q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,fileNameInput,false,&myEff,CHECKnegEFF,SAVEPLOT,specBin,EFFis2Dnot3D,atoi(SignalType.c_str()));
	  else                     ReadEfficiencies(false,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,fileNameInput,false,&myEff,CHECKnegEFF,SAVEPLOT,specBin,EFFis2Dnot3D,atoi(SignalType.c_str()));


	  delete Utility;
	  if (SETBATCH == false) theApp.Run (); // Eventloop on air
	  return EXIT_SUCCESS;
	}
      else if (((option == "ReadAnaly") || (option == "ReadGenAnaly")) && (argc >= 3))
	{
	  string fileNameInput = argv[2];
	  int specBin = -1;
	  if (argc == 4) specBin = atoi(argv[3]);

	  
	  if (SETBATCH == true)
	    {
	      cout << "\n[ComputeEfficiency::main]\t@@@ Setting batch mode @@@" << endl;
	      gROOT->SetBatch(true);
	    }
	  TApplication theApp ("Applications", &argc, argv);


	  Utility = new Utils(RIGHTtag);
	  if (Utility->RIGHTflavorTAG == true) Utility->ReadAllBins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"goodTag");
	  else                                 Utility->ReadAllBins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"misTag");
	  if (option == "ReadAnaly") ReadEfficiencies(true, &q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,fileNameInput,true,&myEff,CHECKnegEFF,SAVEPLOT,specBin,EFFis2Dnot3D);
	  else                       ReadEfficiencies(false,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,fileNameInput,true,&myEff,CHECKnegEFF,SAVEPLOT,specBin,EFFis2Dnot3D);


	  delete Utility;
	  if (SETBATCH == false) theApp.Run (); // Eventloop on air
	  return EXIT_SUCCESS;
	}
      else if ((option == "GenBin") && (argc = 3))
	{
	  string fileNameInput = argv[2];


	  Utility = new Utils(RIGHTtag);
	  if (Utility->RIGHTflavorTAG == true) Utility->ReadAllBins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"goodTag");
	  else                                 Utility->ReadAllBins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"misTag");
	  GenerateEfficiencies(&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,fileNameInput);


	  delete Utility;
	  return EXIT_SUCCESS;
	}
      else if (((option == "Fit1DEff") && (argc == 6)) || ((option == "Fit2DEff") && (argc == 5)) || ((option == "Fit3DEff") && (argc == 5)))
	{
	  string SignalType      = argv[2];
	  string fileNameInput   = argv[3];
	  unsigned int q2BinIndx = atoi(argv[4]);
	  string whichVar2Fit    = "";
	  if (option == "Fit1DEff") whichVar2Fit = argv[5];

	  
	  if (SETBATCH == true)
	    {
	      cout << "\n[ComputeEfficiency::main]\t@@@ Setting batch mode @@@" << endl;
	      gROOT->SetBatch(true);
	    }
	  TApplication theApp ("Applications", &argc, argv);


	  Utility = new Utils(RIGHTtag);
	  if (Utility->RIGHTflavorTAG == true) Utility->ReadAllBins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"goodTag");
	  else                                 Utility->ReadAllBins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"misTag");
	  Utility->ReadEfficiency(fileNameInput.c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,&myEff);

	  if      (option == "Fit1DEff") Fit1DEfficiencies(&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,myEff,whichVar2Fit,q2BinIndx,"Theta.txt");
	  else if (option == "Fit2DEff") Fit2DEfficiencies(&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,atoi(SignalType.c_str()),myEff,q2BinIndx,"ThetaKThetaL.txt",SAVEPLOT);
	  else if (option == "Fit3DEff") Fit3DEfficiencies(&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,atoi(SignalType.c_str()),myEff,q2BinIndx,"ThetaKThetaLPhi.txt",SAVEPLOT);


	  delete Utility;
	  if (SETBATCH == false) theApp.Run (); // Eventloop on air
	  return EXIT_SUCCESS;
	}
      else if (((option == "Test2DEff") || (option == "Test3DEff")) && (argc == 6))
	{
	  string SignalType      = argv[2];
	  string fileNameInput   = argv[3];
	  string analyORbin      = argv[4];
	  unsigned int q2BinIndx = atoi(argv[5]);

	  
	  if (SETBATCH == true)
	    {
	      cout << "\n[ComputeEfficiency::main]\t@@@ Setting batch mode @@@" << endl;
	      gROOT->SetBatch(true);
	    }
	  TApplication theApp ("Applications", &argc, argv);


	  Utility = new Utils(RIGHTtag);
	  if (Utility->RIGHTflavorTAG == true) Utility->ReadAllBins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"goodTag");
	  else                                 Utility->ReadAllBins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"misTag");
	  Utility->ReadEfficiency(fileNameInput.c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,&myEff);
	  
	  if      (option == "Test2DEff") Test2DEfficiency(&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,myEff,q2BinIndx,atoi(SignalType.c_str()),analyORbin,SAVEPLOT);
	  else if (option == "Test3DEff") Test3DEfficiency(&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,myEff,q2BinIndx,SAVEPLOT);


	  delete Utility;
	  if (SETBATCH == false) theApp.Run (); // Eventloop on air
	  return EXIT_SUCCESS;
	}
      else
	{
	  cout << "Efficiency = [#ev. after filter (GEN-level) / #ev. generated (GEN-level before filter)] * [#ev. in single cand. (reco-level truth-matched) / #ev. after filter (GEN-level in multi cand.)]" << endl;
	  cout << "1. process inputFileGenCandidatesNoFilter.root with AddVars2Candidates nvGen" << endl;
	  cout << "2. process inputRecoCandidates.root with AddVars2Candidates nvGen" << endl;
	  cout << "3. make inputFileSingleCand.root from inputRecoCandidates.root" << endl; 

	  cout << "\nParameter missing: " << endl;
	  cout << "./ComputeEfficiency [Make ReadBin ReadGenBin ReadAnaly ReadGenAnaly GenBin Fit1DEff Fit2DEff Fit3DEff Test2DEff Test3DEff] " << endl;
	  cout << "[inputFileGenCandidatesNoFilter.root inputRecoCandidates.root inputFileSingleCand.root] " << endl;
	  cout << "[out/in]putFile.txt [SignalType] [q2 bin indx.]" << endl;
	  cout << "SignalType : if B0 --> K*0 mumu : 1; if B0 --> J/psi K*0 : 3; if B0 --> psi(2S) K*0 : 5" << endl;

	  cout << "\nMake           --> [root files for efficiency computation] AND [outputFile.txt] AND [SignalType]" << endl;

	  cout << "ReadBin        --> [file with binned efficiency] AND [SignalType] AND [q2 bin indx.(optional)]" << endl;
	  cout << "ReadGenBin     --> [file with binned efficiency] AND [SignalType] AND [q2 bin indx.(optional)]" << endl;

	  cout << "ReadAnaly      --> [file with analytical efficiency] AND [q2 bin indx.(optional)]" << endl;	  
	  cout << "ReadGenAnaly   --> [file with analytical efficiency] AND [q2 bin indx.(optional)]" << endl;

	  cout << "GenBin         --> [file with binned efficiency]" << endl;

	  cout << "Fit1DEff       --> [SignalType] AND [file with binned efficiency] AND [q2 bin indx.] [thetaL thetaK phi]" << endl;
	  cout << "Fit2DEff       --> [SignalType] AND [file with binned efficiency] AND [q2 bin indx.]" << endl;
	  cout << "Fit3DEff       --> [SignalType] AND [file with binned efficiency] AND [q2 bin indx.]" << endl;

	  cout << "Test2DEff      --> [SignalType] AND [file with binned efficiency] AND [ANALY or BIN] AND [q2 bin indx.]" << endl;
	  cout << "Test3DEff      --> [SignalType] AND [file with binned efficiency] AND [ANALY or BIN] AND [q2 bin indx.]" << endl;

	  return EXIT_FAILURE;
	}
    }
  else
    {
      cout << "Efficiency = [#ev. after filter (GEN-level) / #ev. generated (GEN-level before filter)] * [#ev. in single cand. (reco-level truth-matched) / #ev. after filter (GEN-level in multi cand.)]" << endl;
      cout << "1. process inputFileGenCandidatesNoFilter.root with AddVars2Candidates nvGen" << endl;
      cout << "2. process inputRecoCandidates.root with AddVars2Candidates nvGen" << endl;
      cout << "3. make inputFileSingleCand.root from inputRecoCandidates.root" << endl; 

      cout << "\nParameter missing: " << endl;
      cout << "./ComputeEfficiency [Make ReadBin ReadGenBin ReadAnaly ReadGenAnaly GenBin Fit1DEff Fit2DEff Fit3DEff Test2DEff Test3DEff] " << endl;
      cout << "[inputFileGenCandidatesNoFilter.root inputRecoCandidates.root inputFileSingleCand.root] " << endl;
      cout << "[out/in]putFile.txt [SignalType] [q2 bin indx.]" << endl;
      cout << "SignalType : if B0 --> K*0 mumu : 1; if B0 --> J/psi K*0 : 3; if B0 --> psi(2S) K*0 : 5" << endl;

      cout << "\nMake           --> [root files for efficiency computation] AND [outputFile.txt] AND [SignalType]" << endl;

      cout << "ReadBin        --> [file with binned efficiency] AND [SignalType] AND [q2 bin indx.(optional)]" << endl;
      cout << "ReadGenBin     --> [file with binned efficiency] AND [SignalType] AND [q2 bin indx.(optional)]" << endl;
      
      cout << "ReadAnaly      --> [file with analytical efficiency] AND [q2 bin indx.(optional)]" << endl;	  
      cout << "ReadGenAnaly   --> [file with analytical efficiency] AND [q2 bin indx.(optional)]" << endl;
      
      cout << "GenBin         --> [file with binned efficiency]" << endl;
      
      cout << "Fit1DEff       --> [SignalType] AND [file with binned efficiency] AND [q2 bin indx.] [thetaL thetaK phi]" << endl;
      cout << "Fit2DEff       --> [SignalType] AND [file with binned efficiency] AND [q2 bin indx.]" << endl;
      cout << "Fit3DEff       --> [SignalType] AND [file with binned efficiency] AND [q2 bin indx.]" << endl;
      
      cout << "Test2DEff      --> [SignalType] AND [file with binned efficiency] AND [ANALY or BIN] AND [q2 bin indx.]" << endl;
      cout << "Test3DEff      --> [SignalType] AND [file with binned efficiency] AND [ANALY or BIN] AND [q2 bin indx.]" << endl;
      
      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
