// ########################################################################
// # Program to plot final histograms for the B0 --> K*0 mu+ mu- analysis #
// ########################################################################
// # Author: Mauro Dinardo                                                #
// ########################################################################

#include <TROOT.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TCutG.h>
#include <TMath.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TPaveStats.h>
#include <TExec.h>
#include <TGraphBentErrors.h>
#include <TKey.h>
#include <TFitResult.h>

#include <RooRealVar.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include "ReadParameters.h"
#include "Utils.h"
#include "B0KstMuMuSingleCandTreeContent.h"

using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::vector;


// ####################
// # Global constants #
// ####################
#define PARAMETERFILEIN      "/python/ParameterFile.txt"
#define ParameterFILE_MCGEN  "/results/ParameterFile_Sig_Psi_MCGEN.txt"
#define ParameterFILE_MCRECO "/results/ParameterFile_Sig_MCRECO.txt"

#define YvalueOutsideLimits 20.0 // Value given to bins with zero error in order not to show them
#define FORPAPER true            // "true" = make special layout for publication in "MakePhysicsPlots" member function
#define q0SM  4.0                // Standard Model value of AFB zero crossing point
#define q0SME 0.2                // Error on q0SM
#define SIGNLUMI 5951.1
#define JPSILUMI   52.5
#define PSIPLUMI   58.5
#define DATALUMI   20.5

// ##################
// # SM predictions #
// ##################
#define SMBINFL_nL  "Data2012B0KstMuMuResults/PredictionSM/FLBinned_noLattice.data"
#define SMBINAFB_nL "Data2012B0KstMuMuResults/PredictionSM/AFBBinned_noLattice.data"
#define SMBINBF_nL  "Data2012B0KstMuMuResults/PredictionSM/dBRBinned_noLattice.data"

#define SMBINFL_yL  "Data2012B0KstMuMuResults/PredictionSM/FLBinned_yesLattice.data"
#define SMBINAFB_yL "Data2012B0KstMuMuResults/PredictionSM/AFBBinned_yesLattice.data"
#define SMBINBF_yL  "Data2012B0KstMuMuResults/PredictionSM/dBRBinned_yesLattice.data"

// ######################
// # Data/MC comparison #
// ######################
#define SingleCand_MCkstMuMu  "Data2012B0KstMuMuResults/MonteCarlo2012/SingleCand/singleCand_B0ToKstMuMu_MC_NTuple.root"
#define SingleCand_MCkstJPsi  "Data2012B0KstMuMuResults/MonteCarlo2012/SingleCand/singleCand_B0ToJPsiKst_MC_NTuple.root"
#define SingleCand_MCkstPsi2S "Data2012B0KstMuMuResults/MonteCarlo2012/SingleCand/singleCand_B0ToPsi2SKst_MC_NTuple.root"
#define SingleCand_Data       "Data2012B0KstMuMuResults/Data2012/SingleCand/singleCand_B0ToKstMuMu_Data2012ABCD_NTuples.root"


// ####################
// # Global variables #
// ####################
Utils* Utility;


// #######################
// # Function Definition #
// #######################
void SetStyle              ();
void DrawString            (double Lumi);
void MakeComparisonDataMC  (unsigned int plotType);
TCutG* DrawExclusion       (double Xlow, double Xhigh, double Ylow, double Yhigh, string cutName, unsigned int fillStyle, unsigned int color);
TGraphAsymmErrors* ReadFromASCII (string fileName, unsigned int PlotType, vector<double>* q2Bins, vector<double>* vxs, vector<double>* vys, vector<double>* vxel, vector<double>* vxeh, vector<double>* vyel, vector<double>* vyeh);
void CheckPhysicsRegion    ();
void MakePhysicsPlots      (unsigned int PlotType);
void GenNTupleFromMultyRun (string fileName, unsigned int q2BinIndx);
void PlotMuMu              (string fileName, bool bkgSub);
void PlotKst               (string fileName, bool bkgSub, bool fitParamAreFixed);
void PlotKK                (string fileName, bool bkgSub, string RECOorGEN);
void PlotMuHadMass         (string fileName);
void MakeFitResPlots       (string fileName, string plotType, int specBin, string varName, double lowBound, double highBound);
void MakePvaluePlot        (string fileName, int specBin);


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


void DrawString (double Lumi)
{
  stringstream myString;
  double scaleRespect2CMS = 0.75;


  myString.clear(); myString.str("");
  myString << "CMS";
  TLatex* LumiTex1 = new TLatex(0.11,0.9,myString.str().c_str());
  LumiTex1->SetTextFont(61);
  LumiTex1->SetTextSize(0.05);
  LumiTex1->SetTextColor(kBlack);
  LumiTex1->SetNDC(true);
  LumiTex1->DrawLatex(0.11,0.9,myString.str().c_str());


  myString.clear(); myString.str("");
  myString << "#it{Preliminary}";
  TLatex* LumiTex2 = new TLatex(0.19,0.9,myString.str().c_str());
  LumiTex2->SetTextFont(42);
  LumiTex2->SetTextSize(0.05 * scaleRespect2CMS);
  LumiTex2->SetTextColor(kBlack);
  LumiTex2->SetNDC(true);
  LumiTex2->DrawLatex(0.19,0.9,myString.str().c_str());


  myString.clear(); myString.str("");
  myString << Lumi <<  " fb#lower[0.4]{^{#font[122]{\55}1}} (8 TeV)";
  TLatex* LumiTex3 = new TLatex(0.8,0.9,myString.str().c_str());
  LumiTex3->SetTextFont(42);
  LumiTex3->SetTextSize(0.05 * scaleRespect2CMS);
  LumiTex3->SetTextColor(kBlack);
  LumiTex3->SetNDC(true);
  LumiTex3->DrawLatex(0.8,0.9,myString.str().c_str());
}


void MakeComparisonDataMC (unsigned int plotType)
// #############################################
// # plotType =  0 --> B0 pT                   #
// # plotType =  1 --> B0 eta                  #
// # plotType =  2 --> B0 vtx CL               #
// # plotType =  3 --> B0 L/sigma              #
// # plotType =  4 --> B0 cos(alpha)           #
// #############################################
// # plotType =  5 --> mu+ pT                  #
// # plotType =  6 --> mu- pT                  #
// # plotType =  7 --> mu+ eta                 #
// # plotType =  8 --> mu- eta                 #
// # plotType =  9 --> mu+ phi, eta range      #
// # plotType = 10 --> mu+ phi, eta range      #
// # plotType = 11 --> mu+ phi, eta range      #
// # plotType = 12 --> mu+ phi, eta range      #
// # plotType = 13 --> mu- phi, eta range      #
// # plotType = 14 --> mu- phi, eta range      #
// # plotType = 15 --> mu- phi, eta range      #
// # plotType = 16 --> mu- phi, eta range      #
// #############################################
// # plotType = 17 --> K*0 trk+ pT             #
// # plotType = 18 --> K*0 trk- pT             #
// # plotType = 19 --> K*0 trk+ eta            #
// # plotType = 20 --> K*0 trk+ eta            #
// # plotType = 21 --> K*0 trk+ phi, eta range #
// # plotType = 22 --> K*0 trk+ phi, eta range #
// # plotType = 23 --> K*0 trk+ phi, eta range #
// # plotType = 24 --> K*0 trk+ phi, eta range #
// # plotType = 25 --> K*0 trk- phi, eta range #
// # plotType = 26 --> K*0 trk- phi, eta range #
// # plotType = 27 --> K*0 trk- phi, eta range #
// # plotType = 28 --> K*0 trk- phi, eta range #
// # plotType = 29 --> K*0 trk+ DCA/sigma      #
// # plotType = 30 --> K*0 trk- DCA/sigma      #
// #############################################
// # plotType = 31 --> cos(theta_K)            #
// # plotType = 32 --> cos(theta_l)            #
// # plotType = 33 --> phi                     #
// #############################################
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  SetStyle();

  
  stringstream myString;
  const unsigned int NHisto = 3;
  TLegend* leg;
  vector<TFile*> Vfiles;
  vector<TTree*> TreeMC;
  vector<TH1D*> h1DVecSig;
  vector<TH1D*> h1DVecBkg;
  vector<string> queryMC;
  string weightVar = "evWeight";
  string fileName;


  // ##################
  // # Read the trees #
  // ##################
  Vfiles.push_back(TFile::Open(SingleCand_MCkstMuMu,"READ"));
  TreeMC.push_back((TTree*)Vfiles.back()->Get("B0KstMuMu/B0KstMuMuNTuple"));

  Vfiles.push_back(TFile::Open(SingleCand_MCkstJPsi,"READ"));
  TreeMC.push_back((TTree*)Vfiles.back()->Get("B0KstMuMu/B0KstMuMuNTuple"));

  Vfiles.push_back(TFile::Open(SingleCand_MCkstPsi2S,"READ"));
  TreeMC.push_back((TTree*)Vfiles.back()->Get("B0KstMuMu/B0KstMuMuNTuple"));


  Vfiles.push_back(TFile::Open(SingleCand_Data,"READ"));
  TTree* TreeData = (TTree*)Vfiles.back()->Get("B0KstMuMu/B0KstMuMuNTuple");


  double minX = 0.0;
  double maxX = 40.0;
  unsigned int nBinsX = 100;
  string Xtitle = "";

  double minY = 0.0;
  double maxY = 0.0;

  double signalSigma = sqrt( atof(Utility->GetGenericParam("FRACMASSS").c_str()) * atof(Utility->GetGenericParam("SIGMAS1").c_str()) * atof(Utility->GetGenericParam("SIGMAS1").c_str()) +
			     (1. - atof(Utility->GetGenericParam("FRACMASSS").c_str())) * atof(Utility->GetGenericParam("SIGMAS2").c_str()) * atof(Utility->GetGenericParam("SIGMAS2").c_str()) );
  cout << "\n[MakePlots::MakeComparisonDataMC]\t@@@ Signal sigma: " << signalSigma << " @@@" << endl;


  // ######
  // # B0 #
  // ######
  if (plotType == 0)
    {
      Xtitle = "B^{0} p_{T} (GeV)";
      maxX = 100.0;

      fileName = "B0pT.pdf";
    }
  else if (plotType == 1)
    {
      Xtitle = "B^{0} #eta";
      minX = -3.0;
      maxX = 3.0;

      fileName = "B0eta.pdf";
    }
  else if (plotType == 2)
    {
      Xtitle = "B^{0} vtx. CL";
      minX = 0.0;
      maxX = 1.0;

      fileName = "B0vtxCL.pdf";
   }
  else if (plotType == 3)
    {
      Xtitle = "B^{0} L/#sigma";
      minX = 0.0;
      maxX = 100.0;

      fileName = "B0vtxLS.pdf";
   }
  else if (plotType == 4)
    {
      Xtitle = "B^{0} cos(#alpha)";
      minX = 0.9994;
      maxX = 1.0;

      nBinsX = 28;

      fileName = "B0cosAlphaBS.pdf";
   }
  
  // #########
  // # Muons #
  // #########
  else if (plotType == 5)
    {
      Xtitle = "#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}} p_{T} (GeV)";
      maxX = 40.0;

      fileName = "MuppT.pdf";
    }
  else if (plotType == 6)
    {
      Xtitle = "#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{\55}}}} p_{T} (GeV)";
      maxX = 40.0;

      fileName = "MumpT.pdf";
    }
  else if (plotType == 7)
    {
      Xtitle = "#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}} #eta";
      minX = -2.4;
      maxX = 2.4;

      fileName = "Mupeta.pdf";
    }
  else if (plotType == 8)
    {
      Xtitle = "#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{\55}}}} #eta";
      minX = -2.4;
      maxX = 2.4;

      fileName = "Mumeta.pdf";
    }
  else if (plotType == 9)
    {
      Xtitle = "#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = -2.4;
      maxY = -1.2;

      nBinsX = 15;

      fileName = "Mupphi_eta0.pdf";
    }
  else if (plotType == 10)
    {
      Xtitle = "#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = -1.2;
      maxY = 0;

      nBinsX = 15;

      fileName = "Mupphi_eta1.pdf";
    }
  else if (plotType == 11)
    {
      Xtitle = "#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = 0.0;
      maxY = 1.2;

      nBinsX = 15;

      fileName = "Mupphi_eta2.pdf";
    }
  else if (plotType == 12)
    {
      Xtitle = "#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = 1.2;
      maxY = 2.4;

      nBinsX = 15;

      fileName = "Mupphi_eta3.pdf";
    }
  else if (plotType == 13)
    {
      Xtitle = "#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{\55}}}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = -2.4;
      maxY = -1.2;

      nBinsX = 15;

      fileName = "Mumphi_eta0.pdf";
    }
  else if (plotType == 14)
    {
      Xtitle = "#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{\55}}}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = -1.2;
      maxY = 0;

      nBinsX = 15;

      fileName = "Mumphi_eta1.pdf";
    }
  else if (plotType == 15)
    {
      Xtitle = "#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{\55}}}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = 0.0;
      maxY = 1.2;

      nBinsX = 15;

      fileName = "Mumphi_eta2.pdf";
    }
  else if (plotType == 16)
    {
      Xtitle = "#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{\55}}}} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = 1.2;
      maxY = 2.4;

      nBinsX = 15;

      fileName = "Mumphi_eta3.pdf";
    }

  // ###########
  // # Hadrons #
  // ###########
  else if (plotType == 17)
    {
      Xtitle = "#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}} trk#font[122]{+} p_{T} (GeV)";
      maxX = 20.0;

      // @TMP@ : special binning for hadron pT
      // nBinsX = 16;

      fileName = "KstTrkppT.pdf";
    }
  else if (plotType == 18)
    {
      Xtitle = "#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}} trk#font[122]{\55} p_{T} (GeV)";
      maxX = 20.0;

      // @TMP@ : special binning for hadron pT
      // nBinsX = 16;

      fileName = "KstTrkmpT.pdf";
    }
  else if (plotType == 19)
    {
      Xtitle = "#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}} trk#font[122]{+} #eta";
      minX = -3.0;
      maxX = 3.0;

      fileName = "KstTrkpeta.pdf";
    }
  else if (plotType == 20)
    {
      Xtitle = "#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}} trk#font[122]{\55} #eta";
      minX = -3.0;
      maxX = 3.0;

      fileName = "KstTrkmeta.pdf";
    }
  else if (plotType == 21)
    {
      Xtitle = "#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}} trk#font[122]{+} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = -3.0;
      maxY = -1.2;

      nBinsX = 15;

      fileName = "KstTrkpphi_eta0.pdf";
    }
  else if (plotType == 22)
    {
      Xtitle = "#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}} trk#font[122]{+} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = -1.2;
      maxY = 0.0;

      nBinsX = 15;

      fileName = "KstTrkpphi_eta1.pdf";
    }
  else if (plotType == 23)
    {
      Xtitle = "#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}} trk#font[122]{+} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = 0.0;
      maxY = 1.2;

      nBinsX = 15;

      fileName = "KstTrkpphi_eta2.pdf";
    }
  else if (plotType == 24)
    {
      Xtitle = "#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}} trk#font[122]{+} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = 1.2;
      maxY = 3.0;
 
      nBinsX = 15;

      fileName = "KstTrkpphi_eta3.pdf";
   }
  else if (plotType == 25)
    {
      Xtitle = "#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}} trk#font[122]{\55} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = -3.0;
      maxY = -1.2;
 
      nBinsX = 15;

      fileName = "KstTrkmphi_eta0.pdf";
   }
  else if (plotType == 26)
    {
      Xtitle = "#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}} trk#font[122]{\55} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = -1.2;
      maxY = 0.0;
 
      nBinsX = 15;

      fileName = "KstTrkmphi_eta1.pdf";
   }
  else if (plotType == 27)
    {
      Xtitle = "#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}} trk#font[122]{\55} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = 0.0;
      maxY = 1.2;

      nBinsX = 15;

      fileName = "KstTrkmphi_eta2.pdf";
    }
  else if (plotType == 28)
    {
      Xtitle = "#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}} trk#font[122]{\55} #phi";
      minX = -3.15;
      maxX = 3.15;

      // Eta range
      minY = 1.2;
      maxY = 3.0;
 
      nBinsX = 15;

      fileName = "KstTrkmphi_eta3.pdf";
   }
  else if (plotType == 29)
    {
      Xtitle = "#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}} trk#font[122]{+} DCA/#sigma";
      minX = 0.0;
      maxX = 40;

      nBinsX = 40;

      fileName = "KstTrkpDCABS.pdf";
   }
  else if (plotType == 30)
    {
      Xtitle = "#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}} trk#font[122]{\55} DCA/#sigma";
      minX = 0.0;
      maxX = 40;

      nBinsX = 40;

      fileName = "KstTrkmDCABS.pdf";
   }

  // ##########
  // # Angles #
  // ##########
  else if (plotType == 31)
    {
      Xtitle = "cos(#theta#lower[-0.4]{_{#font[122]{K}}})";
      minX = -1.0;
      maxX = 1.0;

      nBinsX = 20;

      fileName = "CosThetaK_dataMC.pdf";
   }
  else if (plotType == 32)
    {
      Xtitle = "cos(#theta#lower[-0.4]{_{#font[12]{l}}})";
      minX = -1.0;
      maxX = 1.0;

      nBinsX = 20;

      fileName = "CosThetaL_dataMC.pdf";
   }
  else if (plotType == 33)
    {
      Xtitle = "#phi";
      minX = -3.15;
      maxX = 3.15;

      nBinsX = 20;

      fileName = "phi_dataMC.pdf";
   }


  // #################
  // # 1D histograms #
  // #################
  double* binning = new double[nBinsX+1];
  // @TMP@ : special binning for hadron pT
  // if ((plotType == 17) || (plotType == 18))
  //   {
  //     unsigned int it = 0;

  //     binning[it++] =  0.0;
  //     binning[it++] =  0.5;
  //     binning[it++] =  1.0;
  //     binning[it++] =  1.5;
  //     binning[it++] =  2.0;
  //     binning[it++] =  2.5;
  //     binning[it++] =  3.0;
  //     binning[it++] =  3.5;
  //     binning[it++] =  4.0;
  //     binning[it++] =  5.0;
  //     binning[it++] =  6.0;
  //     binning[it++] =  7.0;
  //     binning[it++] =  8.0;
  //     binning[it++] = 10.0;
  //     binning[it++] = 12.0;
  //     binning[it++] = 16.0;
  //     binning[it++] = 20.0;
  //   }
  // else
    for (unsigned int it = 0; it < nBinsX+1; it++) binning[it] = minX + (maxX-minX) / static_cast<double>(nBinsX) * static_cast<double>(it);

  for (unsigned int i = 0; i < NHisto; i++)
    {
      myString.clear(); myString.str("");
      myString << "h1Dsig" << i;
      h1DVecSig.push_back(new TH1D(myString.str().c_str(),myString.str().c_str(),nBinsX,binning));
      h1DVecSig.back()->Sumw2();
      h1DVecSig.back()->SetXTitle(Xtitle.c_str());
      h1DVecSig.back()->SetYTitle("Norm. entries");
    }

  for (unsigned int i = 0; i < NHisto; i++)
    {
      myString.clear(); myString.str("");
      myString << "h1Dbkg" << i;
      h1DVecBkg.push_back(new TH1D(myString.str().c_str(),myString.str().c_str(),nBinsX,binning));
      h1DVecBkg.back()->Sumw2();
      h1DVecBkg.back()->SetXTitle(Xtitle.c_str());
      h1DVecBkg.back()->SetYTitle("Norm. entries");
    }

  TH1D* hDsig1D = new TH1D("hDsig1D","hDsig1D",nBinsX,binning);
  hDsig1D->Sumw2();
  hDsig1D->SetXTitle(Xtitle.c_str());
  hDsig1D->SetYTitle("Norm. entries");
  hDsig1D->SetMarkerStyle(20);

  TH1D* hDbkg1D = new TH1D("hDbkg1D","hDbkg1D",nBinsX,binning);
  hDbkg1D->Sumw2();
  hDbkg1D->SetXTitle(Xtitle.c_str());
  hDbkg1D->SetYTitle("Norm. entries");
  hDbkg1D->SetMarkerStyle(20);


  // #################
  // # Query NTuples #
  // #################
  string sigMassQuery = "";
  string bkgMassQuery = "";
  string query        = "";
  string tmpstring    = "";
  string selection    = "";
  string aVar         = "";
  string bVar         = "";
  
  myString.clear(); myString.str("");
  myString << "((abs(B0MassArb - " << Utility->B0Mass << ") < " << atof(Utility->GetGenericParam("NSigmaB0S").c_str())*signalSigma << ") && ";
  myString << "((abs(mumuMass - " << Utility->JPsiMass << ") < " << atof(Utility->GetGenericParam("NSigmaPsi").c_str()) << " * mumuMassE) ||";
  myString << " (abs(mumuMass - " << Utility->PsiPMass << ") < " << atof(Utility->GetGenericParam("NSigmaPsi").c_str()) << " * mumuMassE)))";
  sigMassQuery = myString.str();

  myString.clear(); myString.str("");
  myString << "(((B0MassArb > " << Utility->B0Mass + atof(Utility->GetGenericParam("NSigmaB0B").c_str())*signalSigma << " && B0MassArb < "
	   << Utility->B0Mass + (atof(Utility->GetGenericParam("NSigmaB0B").c_str()) + atof(Utility->GetGenericParam("NSigmaB0S").c_str()))*signalSigma << ") || ";
  myString << "(B0MassArb < " << Utility->B0Mass - atof(Utility->GetGenericParam("NSigmaB0B").c_str())*signalSigma << " && B0MassArb > "
	   << Utility->B0Mass - (atof(Utility->GetGenericParam("NSigmaB0B").c_str()) + atof(Utility->GetGenericParam("NSigmaB0S").c_str()))*signalSigma << ")) && ";
  myString << "((abs(mumuMass - " << Utility->JPsiMass << ") < " << atof(Utility->GetGenericParam("NSigmaPsi").c_str()) << " * mumuMassE) ||";
  myString << " (abs(mumuMass - " << Utility->PsiPMass << ") < " << atof(Utility->GetGenericParam("NSigmaPsi").c_str()) << " * mumuMassE)))";
  bkgMassQuery = myString.str();


  // ######
  // # B0 #
  // ######
  if      (plotType == 0) query = "B0pT";        // B0 pT
  else if (plotType == 1) query = "B0Eta";       // B0 eta
  else if (plotType == 2) query = "bVtxCL";      // B0 vtx CL
  else if (plotType == 3) query = "bLBS/bLBSE";  // B0 L/sigma
  else if (plotType == 4) query = "bCosAlphaBS"; // B0 cos(alpha)

  // #########
  // # Muons #
  // #########
  else if (plotType == 5) query = "sqrt(mupPx*mupPx+mupPy*mupPy)"; // mu+ pT
  else if (plotType == 6) query = "sqrt(mumPx*mumPx+mumPy*mumPy)"; // mu- pT
  else if (plotType == 7) query = "0.5*log((sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) + mupPz) / (sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) - mupPz))"; // mu+ eta
  else if (plotType == 8) query = "0.5*log((sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) + mumPz) / (sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) - mumPz))"; // mu- eta
  else if (plotType == 9)
    {
      query = "atan(mupPy / mupPx)"; // mu+ phi in eta range
      tmpstring = "0.5*log((sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) + mupPz) / (sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) - mupPz))";
      myString.clear(); myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // mu+ eta range
      aVar = "mupPx";
      bVar = "mupPy";
    }
  else if (plotType == 10)
    {
      query = "atan(mupPy / mupPx)"; // mu+ phi in eta range
      tmpstring = "0.5*log((sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) + mupPz) / (sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) - mupPz))";
      myString.clear(); myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // mu+ eta range
      aVar = "mupPx";
      bVar = "mupPy";
    }
  else if (plotType == 11)
    {
      query = "atan(mupPy / mupPx)"; // mu+ phi in eta range
      tmpstring = "0.5*log((sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) + mupPz) / (sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) - mupPz))";
      myString.clear(); myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // mu+ eta range
      aVar = "mupPx";
      bVar = "mupPy";
    }
  else if (plotType == 12)
    {
      query = "atan(mupPy / mupPx)"; // mu+ phi in eta range
      tmpstring = "0.5*log((sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) + mupPz) / (sqrt(mupPx*mupPx + mupPy*mupPy + mupPz*mupPz) - mupPz))";
      myString.clear(); myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // mu+ eta range
      aVar = "mupPx";
      bVar = "mupPy";
    }
  else if (plotType == 13)
    {
      query = "atan(mumPy / mumPx)"; // mu- phi in eta range
      tmpstring = "0.5*log((sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) + mumPz) / (sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) - mumPz))";
      myString.clear(); myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // mu- eta range
      aVar = "mumPx";
      bVar = "mumPy";
    }
  else if (plotType == 14)
    {
      query = "atan(mumPy / mumPx)"; // mu- phi in eta range
      tmpstring = "0.5*log((sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) + mumPz) / (sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) - mumPz))";
      myString.clear(); myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // mu- eta range
      aVar = "mumPx";
      bVar = "mumPy";
    }
  else if (plotType == 15)
    {
      query = "atan(mumPy / mumPx)"; // mu- phi in eta range
      tmpstring = "0.5*log((sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) + mumPz) / (sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) - mumPz))";
      myString.clear(); myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // mu- eta range
      aVar = "mumPx";
      bVar = "mumPy";
    }
  else if (plotType == 16)
    {
      query = "atan(mumPy / mumPx)"; // mu- phi in eta range
      tmpstring = "0.5*log((sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) + mumPz) / (sqrt(mumPx*mumPx + mumPy*mumPy + mumPz*mumPz) - mumPz))";
      myString.clear(); myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // mu- eta range
      aVar = "mumPx";
      bVar = "mumPy";
    }

  // ###########
  // # Hadrons #
  // ########### 
  else if (plotType == 17) query = "sqrt(kstTrkpPx*kstTrkpPx+kstTrkpPy*kstTrkpPy)"; // K*0 trk+ pT
  else if (plotType == 18) query = "sqrt(kstTrkmPx*kstTrkmPx+kstTrkmPy*kstTrkmPy)"; // K*0 trk- pT
  else if (plotType == 19) query = "0.5*log((sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) + kstTrkpPz) / (sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) - kstTrkpPz))"; // K*0 trk+ eta
  else if (plotType == 20) query = "0.5*log((sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) + kstTrkmPz) / (sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) - kstTrkmPz))"; // K*0 trk- eta
  else if (plotType == 21)
    {
      query = "atan(kstTrkpPy / kstTrkpPx)"; // K*0 trk+ phi in eta range
      tmpstring = "0.5*log((sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) + kstTrkpPz) / (sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) - kstTrkpPz))";
      myString.clear(); myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // K*0 trk+ eta range
      aVar = "kstTrkpPx";
      bVar = "kstTrkpPy";
    }
  else if (plotType == 22)
    {
      query = "atan(kstTrkpPy / kstTrkpPx)"; // K*0 trk+ phi in eta range
      tmpstring = "0.5*log((sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) + kstTrkpPz) / (sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) - kstTrkpPz))";
      myString.clear(); myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // K*0 trk+ eta range
      aVar = "kstTrkpPx";
      bVar = "kstTrkpPy";
    }
  else if (plotType == 23)
    {
      query = "atan(kstTrkpPy / kstTrkpPx)"; // K*0 trk+ phi in eta range
      tmpstring = "0.5*log((sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) + kstTrkpPz) / (sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) - kstTrkpPz))";
      myString.clear(); myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // K*0 trk+ eta range
      aVar = "kstTrkpPx";
      bVar = "kstTrkpPy";
    }
  else if (plotType == 24)
    {
      query = "atan(kstTrkpPy / kstTrkpPx)"; // K*0 trk+ phi in eta range
      tmpstring = "0.5*log((sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) + kstTrkpPz) / (sqrt(kstTrkpPx*kstTrkpPx + kstTrkpPy*kstTrkpPy + kstTrkpPz*kstTrkpPz) - kstTrkpPz))";
      myString.clear(); myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // K*0 trk+ eta range
      aVar = "kstTrkpPx";
      bVar = "kstTrkpPy";
    }
  else if (plotType == 25)
    {
      query = "atan(kstTrkmPy / kstTrkmPx)"; // K*0 trk- phi in eta range
      tmpstring = "0.5*log((sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) + kstTrkmPz) / (sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) - kstTrkmPz))";
      myString.clear(); myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // K*0 trk- eta range
      aVar = "kstTrkmPx";
      bVar = "kstTrkmPy";
    }
  else if (plotType == 26)
    {
      query = "atan(kstTrkmPy / kstTrkmPx)"; // K*0 trk- phi in eta range
      tmpstring = "0.5*log((sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) + kstTrkmPz) / (sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) - kstTrkmPz))";
      myString.clear(); myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // K*0 trk- eta range
      aVar = "kstTrkmPx";
      bVar = "kstTrkmPy";
    }
  else if (plotType == 27)
    {
      query = "atan(kstTrkmPy / kstTrkmPx)"; // K*0 trk- phi in eta range
      tmpstring = "0.5*log((sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) + kstTrkmPz) / (sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) - kstTrkmPz))";
      myString.clear(); myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // K*0 trk- eta range
      aVar = "kstTrkmPx";
      bVar = "kstTrkmPy";
    }
  else if (plotType == 28)
    {
      query = "atan(kstTrkmPy / kstTrkmPx)"; // K*0 trk- phi in eta range
      tmpstring = "0.5*log((sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) + kstTrkmPz) / (sqrt(kstTrkmPx*kstTrkmPx + kstTrkmPy*kstTrkmPy + kstTrkmPz*kstTrkmPz) - kstTrkmPz))";
      myString.clear(); myString.str("");
      myString << "(" << tmpstring << " > " << minY << " && " << tmpstring << " < " << maxY << ")";
      selection = myString.str(); // K*0 trk- eta range
      aVar = "kstTrkmPx";
      bVar = "kstTrkmPy";
    }
  else if (plotType == 29) query = "kstTrkpDCABS/kstTrkpDCABSE"; // K*0 trk+ DCA/sigma
  else if (plotType == 30) query = "kstTrkmDCABS/kstTrkmDCABSE"; // K*0 trk- DCA/sigma
  
  // ##########
  // # Angles #
  // ##########
  else if (plotType == 31) query = "CosThetaKArb";       // cos(theta_K)
  else if (plotType == 32) query = "CosThetaMuArb";      // cos(theta_l)
  else if (plotType == 33) query = "PhiKstMuMuPlaneArb"; // phi



  TCanvas* c0 = new TCanvas("c0","c0",10,10,700,500);


  // #######################
  // # Query the MC NTuple #
  // #######################
  TH1D* hM1D;
  cout << "\n\n[MakePlots::MakeComparisonDataMC]\t@@@ Query the MC @@@" << endl;
  for (unsigned int i = 0; i < NHisto; i++)
    {
      if (((plotType >= 9) && (plotType <= 16)) || ((plotType >= 21) && (plotType <= 28)))
	{
	  TH1D* hTmp = new TH1D("hTmp","hTmp",nBinsX,minX,maxX);


	  tmpstring = (weightVar != "" ? weightVar + "*" : "") + "(" + selection + " && " + aVar + " > 0 && " + sigMassQuery + ")";
	  myString.clear(); myString.str("");
	  myString << query + ">>hTmp";
	  cout << "\nPlot: " << myString.str().c_str() << endl;
	  cout << "Selection: " << tmpstring.c_str() << endl;
	  TreeMC[i]->Draw(myString.str().c_str(),tmpstring.c_str());
	  h1DVecSig[i]->Add(hTmp);

	  tmpstring = (weightVar != "" ? weightVar + "*" : "") + "(" + selection + " && " + aVar + " < 0 && " + bVar + " < 0 && " + sigMassQuery + ")";
	  myString.clear(); myString.str("");
	  myString << query + "-TMath::Pi()>>hTmp";
	  cout << "\nPlot: " << myString.str().c_str() << endl;
	  cout << "Selection: " << tmpstring.c_str() << endl;
	  TreeMC[i]->Draw(myString.str().c_str(),tmpstring.c_str());
	  h1DVecSig[i]->Add(hTmp);

	  tmpstring = (weightVar != "" ? weightVar + "*" : "") + "(" + selection + " && " + aVar + " < 0 && " + bVar + " > 0 && " + sigMassQuery + ")";
	  myString.clear(); myString.str("");
	  myString << query + "+TMath::Pi()>>hTmp";
	  cout << "\nPlot: " << myString.str().c_str() << endl;
	  cout << "Selection: " << tmpstring.c_str() << endl;
	  TreeMC[i]->Draw(myString.str().c_str(),tmpstring.c_str());
	  h1DVecSig[i]->Add(hTmp);


	  tmpstring = (weightVar != "" ? weightVar + "*" : "") + "(" + selection + " && " + aVar + " > 0 && " + bkgMassQuery + ")";
	  myString.clear(); myString.str("");
	  myString << query + ">>hTmp";
	  cout << "\nPlot: " << myString.str().c_str() << endl;
	  cout << "Selection: " << tmpstring.c_str() << endl;
	  TreeMC[i]->Draw(myString.str().c_str(),tmpstring.c_str());
	  h1DVecBkg[i]->Add(hTmp);

	  tmpstring = (weightVar != "" ? weightVar + "*" : "") + "(" + selection + " && " + aVar + " < 0 && " + bVar + " < 0 && " + bkgMassQuery + ")";
	  myString.clear(); myString.str("");
	  myString << query + "-TMath::Pi()>>hTmp";
	  cout << "\nPlot: " << myString.str().c_str() << endl;
	  cout << "Selection: " << tmpstring.c_str() << endl;
	  TreeMC[i]->Draw(myString.str().c_str(),tmpstring.c_str());
	  h1DVecBkg[i]->Add(hTmp);

	  tmpstring = (weightVar != "" ? weightVar + "*" : "") + "(" + selection + " && " + aVar + " < 0 && " + bVar + " > 0 && " + bkgMassQuery + ")";
	  myString.clear(); myString.str("");
	  myString << query + "+TMath::Pi()>>hTmp";
	  cout << "\nPlot: " << myString.str().c_str() << endl;
	  cout << "Selection: " << tmpstring.c_str() << endl;
	  TreeMC[i]->Draw(myString.str().c_str(),tmpstring.c_str());
	  h1DVecBkg[i]->Add(hTmp);


	  h1DVecSig[i]->Add(h1DVecBkg[i], -1.0);
	  delete hTmp;
	}
      else
	{
	  tmpstring = (weightVar != "" ? weightVar + "*" : "") + "(" + sigMassQuery + ")";
	  myString.clear(); myString.str("");
	  myString << "h1Dsig" << i;
	  queryMC.push_back(query + ">>" + myString.str().c_str());
	  cout << "\nPlot: " << queryMC.back().c_str() << endl;
	  cout << "Selection: " << tmpstring.c_str() << endl;
	  TreeMC[i]->Draw(queryMC.back().c_str(),tmpstring.c_str());


	  tmpstring = (weightVar != "" ? weightVar + "*" : "") + "(" + bkgMassQuery + ")";
	  myString.clear(); myString.str("");
	  myString << "h1Dbkg" << i;
	  queryMC.push_back(query + ">>" + myString.str().c_str());
	  cout << "\nPlot: " << queryMC.back().c_str() << endl;
	  cout << "Selection: " << tmpstring.c_str() << endl;
	  TreeMC[i]->Draw(queryMC.back().c_str(),tmpstring.c_str());


	  h1DVecSig[i]->Add(h1DVecBkg[i], -1.0);
	}

      cout << "\nSelected signal events: " << h1DVecSig[i]->Integral() << "\tSelected background events: " << h1DVecBkg[i]->Integral() << endl;
    }
  
  
  // #########################
  // # Query the Data NTuple #
  // #########################
  string queryData;
  cout << "\n\n[MakePlots::MakeComparisonDataMC]\t@@@ Query the data @@@" << endl;
  if (((plotType >= 9) && (plotType <= 16)) || ((plotType >= 21) && (plotType <= 28)))
    {
      TH1D* hTmp = new TH1D("hTmp","hTmp",nBinsX,minX,maxX);
      

      tmpstring = selection + " && " + aVar + " > 0 && " + sigMassQuery;
      myString.clear(); myString.str("");
      myString << query + ">>hTmp";
      cout << "\nPlot: " << myString.str().c_str() << endl;
      cout << "Selection: " << tmpstring.c_str() << endl;
      TreeData->Draw(myString.str().c_str(),tmpstring.c_str());
      hDsig1D->Add(hTmp);
      
      tmpstring = selection + " && " + aVar + " < 0 && " + bVar + " < 0 && " + sigMassQuery;
      myString.clear(); myString.str("");
      myString << query + "-TMath::Pi()>>hTmp";
      cout << "\nPlot: " << myString.str().c_str() << endl;
      cout << "Selection: " << tmpstring.c_str() << endl;
      TreeData->Draw(myString.str().c_str(),tmpstring.c_str());
      hDsig1D->Add(hTmp);
      
      tmpstring = selection + " && " + aVar + " < 0 && " + bVar + " > 0 && " + sigMassQuery;
      myString.clear(); myString.str("");
      myString << query + "+TMath::Pi()>>hTmp";
      cout << "\nPlot: " << myString.str().c_str() << endl;
      cout << "Selection: " << tmpstring.c_str() << endl;
      TreeData->Draw(myString.str().c_str(),tmpstring.c_str());
      hDsig1D->Add(hTmp);


      tmpstring = selection + " && " + aVar + " > 0 && " + bkgMassQuery;
      myString.clear(); myString.str("");
      myString << query + ">>hTmp";
      cout << "\nPlot: " << myString.str().c_str() << endl;
      cout << "Selection: " << tmpstring.c_str() << endl;
      TreeData->Draw(myString.str().c_str(),tmpstring.c_str());
      hDbkg1D->Add(hTmp);
      
      tmpstring = selection + " && " + aVar + " < 0 && " + bVar + " < 0 && " + bkgMassQuery;
      myString.clear(); myString.str("");
      myString << query + "-TMath::Pi()>>hTmp";
      cout << "\nPlot: " << myString.str().c_str() << endl;
      cout << "Selection: " << tmpstring.c_str() << endl;
      TreeData->Draw(myString.str().c_str(),tmpstring.c_str());
      hDbkg1D->Add(hTmp);
      
      tmpstring = selection + " && " + aVar + " < 0 && " + bVar + " > 0 && " + bkgMassQuery;
      myString.clear(); myString.str("");
      myString << query + "+TMath::Pi()>>hTmp";
      cout << "\nPlot: " << myString.str().c_str() << endl;
      cout << "Selection: " << tmpstring.c_str() << endl;
      TreeData->Draw(myString.str().c_str(),tmpstring.c_str());
      hDbkg1D->Add(hTmp);


      hDsig1D->Add(hDbkg1D, -1.0);
      delete hTmp;
    }
  else
    {
      queryData = query + ">>hDsig1D";
      cout << "\nPlot: " << queryData.c_str() << endl;
      cout << "Selection: " << sigMassQuery.c_str() << endl;
      TreeData->Draw(queryData.c_str(),sigMassQuery.c_str());


      queryData = query + ">>hDbkg1D";
      cout << "\nPlot: " << queryData.c_str() << endl;
      cout << "Selection: " << bkgMassQuery.c_str() << endl;
      TreeData->Draw(queryData.c_str(),bkgMassQuery.c_str());


      hDsig1D->Add(hDbkg1D, -1.0);
    }
  
  cout << "\nSelected signal events: " << hDsig1D->Integral() << "\tSelected background events: " << hDbkg1D->Integral() << endl;
  
  
  // #########################################
  // # Rescale MCs by the Branching Fraction #
  // #########################################

  // # B0 --> K* mu mu #
  h1DVecSig[0]->Scale(1. / SIGNLUMI * DATALUMI);

  // # B0 --> K* J/psi #
  h1DVecSig[1]->Scale(1. / JPSILUMI * DATALUMI);

  // # B0 --> K* psi(2S) #
  h1DVecSig[2]->Scale(1. / PSIPLUMI * DATALUMI);

  hM1D = (TH1D*)h1DVecSig[0]->Clone("hM1D");
  hM1D->SetXTitle(Xtitle.c_str());
  hM1D->SetYTitle("Norm. entries");
  hM1D->SetLineColor(kBlack);
  hM1D->SetFillColor(kAzure+6);

  for (unsigned int i = 1; i < NHisto; i++) hM1D->Add(h1DVecSig[i]);
      
  hM1D->Scale(1./hM1D->Integral());
  hDsig1D->Scale(1./hDsig1D->Integral());
  for (unsigned int it = 0; it < nBinsX; it++)
    {
      hM1D->SetBinContent(it+1,hM1D->GetBinContent(it+1) / hM1D->GetBinWidth(it+1));
      hM1D->SetBinError(it+1,hM1D->GetBinError(it+1) / hM1D->GetBinWidth(it+1));
      
      hDsig1D->SetBinContent(it+1,hDsig1D->GetBinContent(it+1) / hDsig1D->GetBinWidth(it+1));
      hDsig1D->SetBinError(it+1,hDsig1D->GetBinError(it+1) / hDsig1D->GetBinWidth(it+1));
    }

  c0->cd();
  hM1D->Draw("e3");
  hDsig1D->Draw("e1p sames");
  hM1D->GetYaxis()->SetRangeUser(0.0,(hM1D->GetBinContent(hM1D->GetMaximumBin()) > hDsig1D->GetBinContent(hDsig1D->GetMaximumBin()) ?
  				      hM1D->GetBinContent(hM1D->GetMaximumBin()) : hDsig1D->GetBinContent(hDsig1D->GetMaximumBin()))*1.1);
  hDsig1D->GetYaxis()->SetRangeUser(0.0,(hM1D->GetBinContent(hM1D->GetMaximumBin()) > hDsig1D->GetBinContent(hDsig1D->GetMaximumBin()) ?
  					 hM1D->GetBinContent(hM1D->GetMaximumBin()) : hDsig1D->GetBinContent(hDsig1D->GetMaximumBin()))*1.1);

  if ((plotType == 3)  ||
      (plotType == 4)  ||
      ((plotType >= 9) && (plotType <= 16))  ||
      ((plotType >= 21) && (plotType <= 28)) ||
      (plotType == 29) ||
      (plotType == 30) ||
      (plotType == 31) ||
      (plotType == 32) ||
      (plotType == 33))
    leg = new TLegend(0.15, 0.15, 0.34, 0.25, "");
  else if ((plotType == 1)  ||
	   (plotType == 7)  ||
	   (plotType == 8)  ||
	   (plotType == 19) ||
	   (plotType == 20)) leg = new TLegend(0.15, 0.75, 0.34, 0.85, "");
  else leg = new TLegend(0.8, 0.15, 0.97, 0.25, "");
  leg->AddEntry(hM1D,"MC");
  leg->AddEntry(hDsig1D,"Data");
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();

  c0->Modified();
  c0->Update();

  TPaveStats* stM = (TPaveStats*)hM1D->FindObject("stats");
  stM->SetFillColor(kAzure+6);
  TPaveStats* stD = (TPaveStats*)hDsig1D->FindObject("stats");
  if ((plotType == 0)  ||
      (plotType == 3)  ||
      (plotType == 5)  ||
      (plotType == 6)  ||
      (plotType == 17) ||
      (plotType == 18) ||
      (plotType == 29) ||
      (plotType == 30) ||
      (plotType == 31))
    {
      stD->SetX1NDC(0.55);
      stD->SetX2NDC(0.75);
      stD->SetY2NDC(0.89);

      stM->SetY2NDC(0.89);
    }
  else if ((plotType == 1)  ||
	   (plotType == 7)  ||
	   (plotType == 8)  ||
	   (plotType == 19) ||
	   (plotType == 20))
    {
      stD->SetX1NDC(0.83);
      stD->SetY1NDC(0.3);
      stD->SetY2NDC(0.6);

      stM->SetX1NDC(0.83);
      stM->SetY2NDC(0.89);
    }
  else if ((plotType == 2) || (plotType == 4))
    {
      stD->SetX1NDC(0.15);
      stD->SetX2NDC(0.35);
      stD->SetY2NDC(0.89);

      stM->SetX1NDC(0.4);
      stM->SetX2NDC(0.6);
      stM->SetY2NDC(0.89);
    }
  else
    {
      stD->SetX1NDC(0.4);
      stD->SetX2NDC(0.6);
      stD->SetY1NDC(0.15);
      stD->SetY2NDC(0.45);

      stM->SetX1NDC(0.65);
      stM->SetX2NDC(0.85);
      stM->SetY1NDC(0.15);
      stM->SetY2NDC(0.45);
    }


  // ######################################
  // # Estimate the percentage difference #
  // ######################################
  TH1D* hdiff = (TH1D*)hM1D->Clone("hdiff");
  hdiff->Add(hDsig1D,-1);
  for (int i = 0; i < hdiff->GetNbinsX(); i++) hdiff->SetBinContent(i+1,fabs(hdiff->GetBinContent(i+1)));
  cout << "\n\n[MakePlots::MakeComparisonDataMC]\t@@@ Percentage difference: " << hdiff->Integral() / 2. << " @@@" << endl;


  // ###################
  // # data / MC ratio #
  // ###################
  TCanvas* c1 = new TCanvas("c1","c1",20,20,700,500);

  c1->cd();
  TH1D* hratio = (TH1D*)hDsig1D->Clone("hratio");
  hratio->GetYaxis()->SetRangeUser(0.6,1.4);
  hratio->SetStats(false);
  hratio->Divide(hM1D);
  hratio->Draw("pe1");
  c1->Modified();
  c1->Update();


  c0->Print(fileName.c_str());
  c1->Print(fileName.replace(fileName.find(".pdf"),4,"_ratio.pdf").c_str());


  TreeMC.clear();
  h1DVecSig.clear();
  h1DVecBkg.clear();
  queryMC.clear();
}


TCutG* DrawExclusion (double Xlow, double Xhigh, double Ylow, double Yhigh, string cutName, unsigned int fillStyle, unsigned int color)
{
  TCutG* ExclusionZone = new TCutG(cutName.c_str(),5);
  ExclusionZone->SetVarX("");
  ExclusionZone->SetVarY("");
  ExclusionZone->SetPoint(0,Xlow,Ylow);
  ExclusionZone->SetPoint(1,Xhigh,Ylow);
  ExclusionZone->SetPoint(2,Xhigh,Yhigh);
  ExclusionZone->SetPoint(3,Xlow,Yhigh);
  ExclusionZone->SetPoint(4,Xlow,Ylow);
  ExclusionZone->SetFillColor(color);
  ExclusionZone->SetFillStyle(fillStyle);
  ExclusionZone->Draw("F");

  return ExclusionZone;
}


TGraphAsymmErrors* ReadFromASCII (string fileName, unsigned int PlotType, vector<double>* q2Bins, vector<double>* vxs, vector<double>* vys, vector<double>* vxel, vector<double>* vxeh, vector<double>* vyel, vector<double>* vyeh)
{
  ifstream inputFile;

  // ############################################
  // # Variables to read values from ASCII file #
  // ############################################
  double xs  = 0.0;
  double ys  = 0.0;

  double xel = 0.0;
  double xeh = 0.0;

  double yel = 0.0;
  double yeh = 0.0;

  vxs->clear();
  vys->clear();
  vxel->clear();
  vxeh->clear();
  vyel->clear();
  vyeh->clear();

  
  inputFile.open(fileName.c_str(), ifstream::in);
  if (inputFile.good() == false)
    {
      cout << "[MakePlots::ReadFromASCII]\tError opening file : " << fileName.c_str() << endl;
      exit (EXIT_FAILURE);
    }
  cout << "\n[MakePlots::ReadFromASCII]\tReding SM values from file : " << fileName.c_str() << endl; 
  inputFile >> xs >> ys >> xel >> xeh >> yel >> yeh;
  while (inputFile.eof() == false)
    {
      if ((Utility->ValIsInPsi(q2Bins,xs) == false) &&
	  (xs >= q2Bins->operator[](0)) && (xs <= q2Bins->operator[](q2Bins->size()-1)) &&
	  (yel != 0.0) && (yeh != 0.0))
	{
	  if ((PlotType == 0) || (PlotType == 10)) // Fl
	    {
	      vxs->push_back(xs);
	      vys->push_back(ys);
	      vxel->push_back(xel);
	      vxeh->push_back(xeh);
	      vyel->push_back(yel);
	      vyeh->push_back(yeh);
	    }
	  else if ((PlotType == 1) || (PlotType == 11)) // Afb
	    {
	      vxs->push_back(xs);
	      vys->push_back(ys);
	      vxel->push_back(xel);
	      vxeh->push_back(xeh);
	      vyel->push_back(yeh);
	      vyeh->push_back(yel);
	    }
	  else if ((PlotType == 2) || (PlotType == 12)) // Branching fraction
	    {
	      vxs->push_back(xs);
	      vys->push_back(ys);
	      vxel->push_back(xel);
	      vxeh->push_back(xeh);
	      vyel->push_back(yel);
	      vyeh->push_back(yeh);
	    }
	  cout << "xs: " << xs << "\tys: " << ys << "\txel: " << xel << "\txeh: " << xeh << "\tyel: " << yel << "\tyeh: " << yeh << endl;
	}
      else
	{
	  vxs->push_back(xs);
	  vys->push_back(YvalueOutsideLimits);
	  vxel->push_back(xel);
	  vxeh->push_back(xeh);
	  vyel->push_back(yel);
	  vyeh->push_back(yeh);
	}
      inputFile >> xs >> ys >> xel >> xeh >> yel >> yeh;
    }
  inputFile.close();

  TGraphAsymmErrors* ge = new TGraphAsymmErrors(vxs->size(), &(*vxs)[0], &(*vys)[0], &(*vxel)[0], &(*vxeh)[0], &(*vyel)[0], &(*vyeh)[0]);

  return ge;
}


void CheckPhysicsRegion ()
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  SetStyle();

  
  stringstream myString;
  double offset = 0.02;
  TCanvas* canv0 = new TCanvas("canv0","canv0",10,10,1200,600);
  TGraphAsymmErrors* ge;
  TGraphAsymmErrors* geTMP;
  TH1D* histo = new TH1D("histo","histo",100,-1.0,1.0);
  histo->GetYaxis()->SetRangeUser(0.0,1.0);
  histo->SetXTitle("A_{FB}");
  histo->SetYTitle("F_{L}");
  TLine* line1;
  TLine* line2;
  TLatex* binIndx;

  vector<double> q2Bins;
  Utility->Readq2Bins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins);
  double LUMI = Utility->ReadLumi(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str());
  Utility->MakeGraphVar(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&ge,"Fl");
  ge->SetMarkerColor(kBlack);
  ge->SetMarkerStyle(20);
  Utility->MakeGraphVar(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&geTMP,"Afb");

  cout << "\n[MakePlots::CheckPhysicsRegion]\t@@@ I've found " << geTMP->GetN() << " data points @@@" << endl;

  canv0->cd();
  histo->Draw();

  for (int i = 0; i < geTMP->GetN(); i++)
    {
      ge->SetPoint(i,geTMP->GetY()[i],ge->GetY()[i]);
      ge->SetPointEXhigh(i,geTMP->GetErrorYhigh(i));
      ge->SetPointEXlow(i,geTMP->GetErrorYlow(i));

      myString.clear(); myString.str("");
      myString << i;
      binIndx = new TLatex(0.5,0.5,myString.str().c_str());
      binIndx->SetTextFont(61);
      binIndx->SetTextSize(0.04);
      binIndx->SetTextColor(kBlue);
      binIndx->SetNDC(false);
      binIndx->DrawLatex(geTMP->GetY()[i] + offset,ge->GetY()[i] + offset,myString.str().c_str());
    }


  ge->Draw("pe1");

  line1 = new TLine(-3.0/4.0,0.0,0.0,1.0);
  line1->SetLineColor(kRed);
  line1->SetLineWidth(2);
  line1->SetLineStyle(2);
  line1->Draw("same");

  line2 = new TLine(+3.0/4.0,0.0,0.0,1.0);
  line2->SetLineColor(kRed);
  line2->SetLineWidth(2);
  line2->SetLineStyle(2);
  line2->Draw("same");

  DrawString(LUMI);
  canv0->Modified();
  canv0->Update();
}


void MakePhysicsPlots (unsigned int PlotType)
// ####################
// # PlotType:        #
// ####################
// # 0 = GEN-RECO-FL  #
// # 1 = GEN-RECO-AFB #
// # 2 = GEN-RECO-BF  #
// ####################
// # 10 = DATA-FL     #
// # 11 = DATA-AFB    #
// # 12 = DATA-BF     #
// # 13 = DATA-FS     #
// # 14 = DATA-AS     #
// ####################
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  SetStyle();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(8);


  unsigned int DoF = 0;
  double myGlobalChi2 = 0.0;


  // ##################################
  // # Read q^2 bins from config file #
  // ##################################
  vector<double> q2Bins;
  Utility->Readq2Bins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins);
  unsigned int JPsibin = Utility->GetJPsiBin(&q2Bins);
  unsigned int PsiPbin = Utility->GetPsiPBin(&q2Bins);
  double LUMI          = Utility->ReadLumi(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str());
  double* q2Bins_      = Utility->MakeBinning(&q2Bins);


  TCanvas* canv0;
  if (FORPAPER == false) canv0 = new TCanvas("canv0","canv0",10,10,700,900);
  else                   canv0 = new TCanvas("canv0","canv0",10,10,700,500);
  TPad *pad1, *pad2, *pad3;
  TPaveText* paveText = NULL;
  TGraphAsymmErrors* ge0  = NULL;
  TGraphAsymmErrors* ge1  = NULL;
  TGraphAsymmErrors* ge00 = NULL;
  TGraphAsymmErrors* ge11 = NULL;
  TGraphBentErrors*  geb  = NULL;
  TLine* line;


  // ############################################
  // # Variables to read values from ASCII file #
  // ############################################
  vector<double> vxs, vys;
  vector<double> vxel, vxeh;
  vector<double> vyel, vyeh;

 
  if (PlotType == 0) // Fl
    {
      Utility->MakeGraphVar(Utility->MakeAnalysisPATH(ParameterFILE_MCGEN).c_str(),&ge0,"Fl");
      ge0->SetMarkerColor(kBlue);
      ge0->SetMarkerStyle(22);
      ge0->SetMarkerSize(1.2);
      ge0->SetFillColor(kWhite);
      ge0->SetFillStyle(0);
      ge0->SetLineColor(kBlue);
      ge0->SetLineWidth(2);
      ge0->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge0->GetYaxis()->SetRangeUser(-0.02,1.0);
      ge0->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});F_{L}");

      Utility->MakeGraphVar(Utility->MakeAnalysisPATH(ParameterFILE_MCRECO).c_str(),&ge1,"Fl");
      ge1->SetMarkerColor(kBlack);
      ge1->SetMarkerStyle(20);
      ge1->SetMarkerSize(1.2);
      ge1->SetFillColor(kGreen-7);
      ge1->SetFillStyle(3001);
      ge1->SetLineColor(kBlack);
      ge1->SetLineWidth(2);
      ge1->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge1->GetYaxis()->SetRangeUser(-0.02,1.0);
      ge1->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});F_{L}");
     }
  else if (PlotType == 1) // Afb
    {
      Utility->MakeGraphVar(Utility->MakeAnalysisPATH(ParameterFILE_MCGEN).c_str(),&ge0,"Afb");
      ge0->SetMarkerColor(kBlue);
      ge0->SetMarkerStyle(22);
      ge0->SetMarkerSize(1.2);
      ge0->SetFillColor(kWhite);
      ge0->SetFillStyle(0);
      ge0->SetLineColor(kBlue);
      ge0->SetLineWidth(2);
      ge0->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge0->GetYaxis()->SetRangeUser(-1.04,1.0);
      ge0->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});A_{FB}");
 
      Utility->MakeGraphVar(Utility->MakeAnalysisPATH(ParameterFILE_MCRECO).c_str(),&ge1,"Afb");
      ge1->SetMarkerColor(kBlack);
      ge1->SetMarkerStyle(20);
      ge1->SetMarkerSize(1.2);
      ge1->SetFillColor(kGreen-7);
      ge1->SetFillStyle(3001);
      ge1->SetLineColor(kBlack);
      ge1->SetLineWidth(2);
      ge1->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge1->GetYaxis()->SetRangeUser(-1.04,1.0);
      ge1->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});A_{FB}");
    }
  else if (PlotType == 2) // Branching fraction
    {
      Utility->MakeGraphVar(Utility->MakeAnalysisPATH(ParameterFILE_MCGEN).c_str(),&ge0,"BF");
      ge0->SetMarkerColor(kBlue);
      ge0->SetMarkerStyle(22);
      ge0->SetMarkerSize(1.2);
      ge0->SetFillColor(kWhite);
      ge0->SetFillStyle(0);
      ge0->SetLineColor(kBlue);
      ge0->SetLineWidth(2);
      ge0->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge0->GetYaxis()->SetRangeUser(-0.3,12.);
      ge0->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}8}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");
 
      Utility->MakeGraphVar(Utility->MakeAnalysisPATH(ParameterFILE_MCRECO).c_str(),&ge1,"BF");
      ge1->SetMarkerColor(kBlack);
      ge1->SetMarkerStyle(20);
      ge1->SetMarkerSize(1.2);
      ge1->SetFillColor(kGreen-7);
      ge1->SetFillStyle(3001);
      ge1->SetLineColor(kBlack);
      ge1->SetLineWidth(2);
      ge1->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge1->GetYaxis()->SetRangeUser(-0.3,12.);
      ge1->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}8}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");
    }
  else if (PlotType == 10) // Fl
    {
      Utility->MakeGraphVar(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&ge0,"Fl");
      ge0->SetMarkerColor(kBlack);
      ge0->SetMarkerStyle(20);
      ge0->SetMarkerSize(1.2);
      ge0->SetLineColor(kBlack);
      ge0->SetLineWidth(2);
      ge0->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge0->GetYaxis()->SetRangeUser(0.0,1.0);
      ge0->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^2});F_{L}");
    }
  else if (PlotType == 11) // Afb
    {
      Utility->MakeGraphVar(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&ge0,"Afb");
      ge0->SetMarkerColor(kBlack);
      ge0->SetMarkerStyle(20);
      ge0->SetMarkerSize(1.2);
      ge0->SetLineColor(kBlack);
      ge0->SetLineWidth(2);
      ge0->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge0->GetYaxis()->SetRangeUser(-1.0,1.0);
      ge0->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});A_{FB}");
    }
  else if (PlotType == 12) // Branching fraction
    {
      Utility->MakeGraphVar(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&ge0,"BF");
      ge0->SetMarkerColor(kBlack);
      ge0->SetMarkerStyle(20);
      ge0->SetMarkerSize(1.2);
      ge0->SetLineColor(kBlack);
      ge0->SetLineWidth(2);
      ge0->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge0->GetYaxis()->SetRangeUser(0.0,12.);
      ge0->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}8}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");
    }
  else if (PlotType == 13) // Fs
    {
      Utility->MakeGraphVar(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&ge0,"FS");
      ge0->SetMarkerColor(kBlack);
      ge0->SetMarkerStyle(20);
      ge0->SetMarkerSize(1.2);
      ge0->SetLineColor(kBlack);
      ge0->SetLineWidth(2);
      ge0->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge0->GetYaxis()->SetRangeUser(0.0,1.0);
      ge0->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});F_{S}");
    }
  else if (PlotType == 14) // As
    {
      Utility->MakeGraphVar(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&ge0,"AS");
      ge0->SetMarkerColor(kBlack);
      ge0->SetMarkerStyle(20);
      ge0->SetMarkerSize(1.2);
      ge0->SetLineColor(kBlack);
      ge0->SetLineWidth(2);
      ge0->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);
      ge0->GetYaxis()->SetRangeUser(-1.0,1.0);
      ge0->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});A_{S}");
    }
  else
    {
      cout << "[MakePlots::MakePhysicsPlots]\tWrong option number" << endl;
      exit (EXIT_FAILURE);
    }


  // ################################################
  // # Average theory over q^2 bins from ASCII file #
  // ################################################
  TGraphAsymmErrors* geStepTh_nL = NULL;
  if      ((PlotType == 0)  || (PlotType == 10)) geStepTh_nL = ReadFromASCII(SMBINFL_nL,PlotType,&q2Bins,&vxs,&vys,&vxel,&vxeh,&vyel,&vyeh);  // Fl
  else if ((PlotType == 1)  || (PlotType == 11)) geStepTh_nL = ReadFromASCII(SMBINAFB_nL,PlotType,&q2Bins,&vxs,&vys,&vxel,&vxeh,&vyel,&vyeh); // Afb
  else if ((PlotType == 2)  || (PlotType == 12)) geStepTh_nL = ReadFromASCII(SMBINBF_nL,PlotType,&q2Bins,&vxs,&vys,&vxel,&vxeh,&vyel,&vyeh);  // Branching fraction
  else if ((PlotType == 13) || (PlotType == 14)) geStepTh_nL = ReadFromASCII(SMBINBF_nL,12,&q2Bins,&vxs,&vys,&vxel,&vxeh,&vyel,&vyeh);        // Fs OR As
  geStepTh_nL->SetMarkerColor(kBlack);
  geStepTh_nL->SetMarkerStyle(1);
  geStepTh_nL->SetFillColor(kBlue);
  geStepTh_nL->SetFillStyle(3002);
  geStepTh_nL->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);


  TGraphAsymmErrors* geStepTh_yL = NULL;
  if      ((PlotType == 0)  || (PlotType == 10)) geStepTh_yL = ReadFromASCII(SMBINFL_yL,PlotType,&q2Bins,&vxs,&vys,&vxel,&vxeh,&vyel,&vyeh);  // Fl
  else if ((PlotType == 1)  || (PlotType == 11)) geStepTh_yL = ReadFromASCII(SMBINAFB_yL,PlotType,&q2Bins,&vxs,&vys,&vxel,&vxeh,&vyel,&vyeh); // Afb
  else if ((PlotType == 2)  || (PlotType == 12)) geStepTh_yL = ReadFromASCII(SMBINBF_yL,PlotType,&q2Bins,&vxs,&vys,&vxel,&vxeh,&vyel,&vyeh);  // Branching fraction
  else if ((PlotType == 13) || (PlotType == 14)) geStepTh_yL = ReadFromASCII(SMBINBF_yL,12,&q2Bins,&vxs,&vys,&vxel,&vxeh,&vyel,&vyeh);        // Fs OR As
  geStepTh_yL->SetMarkerColor(kBlack);
  geStepTh_yL->SetMarkerStyle(1);
  geStepTh_yL->SetFillColor(kBlack);
  geStepTh_yL->SetFillStyle(0);
  geStepTh_yL->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);


  // ############################
  // # Adding systematic errors #
  // ############################
  vector<vector<double>*> vecObs; // Vector containing the pointers to the vectors containing the fit-observable systematic errors
  Utility->ReadFitSystematics(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&vecObs);
  if (PlotType == 0) // Fl
    {
      // ############################################
      // # Graph containing only statistical errors #
      // ############################################
      ge11 = new TGraphAsymmErrors(*ge1);
      ge11->SetMarkerColor(kBlack);
      ge11->SetMarkerStyle(20);
      ge11->SetMarkerSize(1.2);
      ge11->SetFillColor(kGreen-7);
      ge11->SetFillStyle(3001);
      ge11->SetLineColor(kBlack);
      ge11->SetLineWidth(2);
      ge11->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);


      // ######################################################
      // # Graph containing statistical and systematic errors #
      // ######################################################
      for (int i = 0; i < ge1->GetN(); i++)
      	if (Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false)
      	  {
      	    ge1->SetPointEYhigh(i,sqrt(ge1->GetErrorYhigh(i)*ge1->GetErrorYhigh(i) + vecObs[0]->operator[](i)*vecObs[0]->operator[](i)));
      	    ge1->SetPointEYlow(i, sqrt(ge1->GetErrorYlow(i)*ge1->GetErrorYlow(i)   + vecObs[1]->operator[](i)*vecObs[1]->operator[](i)));
      	  }
    }
  else if (PlotType == 1) // Afb
    {
      // ############################################
      // # Graph containing only statistical errors #
      // ############################################
      ge11 = new TGraphAsymmErrors(*ge1);
      ge11->SetMarkerColor(kBlack);
      ge11->SetMarkerStyle(20);
      ge11->SetMarkerSize(1.2);
      ge11->SetFillColor(kGreen-7);
      ge11->SetFillStyle(3001);
      ge11->SetLineColor(kBlack);
      ge11->SetLineWidth(2);
      ge11->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);


      // ######################################################
      // # Graph containing statistical and systematic errors #
      // ######################################################
      for (int i = 0; i < ge1->GetN(); i++)
	if (Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false)
	  {
	    ge1->SetPointEYhigh(i,sqrt(ge1->GetErrorYhigh(i)*ge1->GetErrorYhigh(i) + vecObs[2]->operator[](i)*vecObs[2]->operator[](i)));
	    ge1->SetPointEYlow(i, sqrt(ge1->GetErrorYlow(i)*ge1->GetErrorYlow(i)   + vecObs[3]->operator[](i)*vecObs[3]->operator[](i)));
	  }
    }
  else if (PlotType == 2) // Branching fraction
    {
      // ############################################
      // # Graph containing only statistical errors #
      // ############################################
      ge11 = new TGraphAsymmErrors(*ge1);
      ge11->SetMarkerColor(kBlack);
      ge11->SetMarkerStyle(20);
      ge11->SetMarkerSize(1.2);
      ge11->SetFillColor(kGreen-7);
      ge11->SetFillStyle(3001);
      ge11->SetLineColor(kBlack);
      ge11->SetLineWidth(2);
      ge11->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);


      // ######################################################
      // # Graph containing statistical and systematic errors #
      // ######################################################
      for (int i = 0; i < ge1->GetN(); i++)
	if (Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false)
	  {
	    ge1->SetPointEYhigh(i,sqrt(ge1->GetErrorYhigh(i)*ge1->GetErrorYhigh(i) + vecObs[8]->operator[](i)*vecObs[8]->operator[](i)));
	    ge1->SetPointEYlow(i, sqrt(ge1->GetErrorYlow(i)*ge1->GetErrorYlow(i)   + vecObs[9]->operator[](i)*vecObs[9]->operator[](i)));
	  }
    }
  else if (PlotType == 10) // Fl
    {
      // ############################################
      // # Graph containing only statistical errors #
      // ############################################
      ge00 = new TGraphAsymmErrors(*ge0);
      ge00->SetMarkerColor(kBlack);
      ge00->SetMarkerStyle(20);
      ge00->SetMarkerSize(1.2);
      ge00->SetLineColor(kBlack);
      ge00->SetLineWidth(2);
      ge00->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);


      // ######################################################
      // # Graph containing statistical and systematic errors #
      // ######################################################
      for (int i = 0; i < ge0->GetN(); i++)
	if (Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false)
	  {
	    ge0->SetPointEXlow(i,0.0);
	    ge0->SetPointEXhigh(i,0.0);

	    ge0->SetPointEYhigh(i,sqrt(ge0->GetErrorYhigh(i)*ge0->GetErrorYhigh(i) + vecObs[0]->operator[](i)*vecObs[0]->operator[](i)));
	    ge0->SetPointEYlow(i, sqrt(ge0->GetErrorYlow(i)*ge0->GetErrorYlow(i)   + vecObs[1]->operator[](i)*vecObs[1]->operator[](i)));
	  }


      geStepTh_nL->GetYaxis()->SetRangeUser(0.0,1.0);
      geStepTh_nL->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});F_{L}");

      geStepTh_yL->GetYaxis()->SetRangeUser(0.0,1.0);
      geStepTh_yL->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});F_{L}");
    }
  else if (PlotType == 11) // Afb
    {
      // ############################################
      // # Graph containing only statistical errors #
      // ############################################
      ge00 = new TGraphAsymmErrors(*ge0);
      ge00->SetMarkerColor(kBlack);
      ge00->SetMarkerStyle(20);
      ge00->SetMarkerSize(1.2);
      ge00->SetLineColor(kBlack);
      ge00->SetLineWidth(2);
      ge00->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);


      // ######################################################
      // # Graph containing statistical and systematic errors #
      // ######################################################
      for (int i = 0; i < ge0->GetN(); i++)
	if (Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false)
	  {
	    ge0->SetPointEXlow(i,0.0);
	    ge0->SetPointEXhigh(i,0.0);

	    ge0->SetPointEYhigh(i,sqrt(ge0->GetErrorYhigh(i)*ge0->GetErrorYhigh(i) + vecObs[2]->operator[](i)*vecObs[2]->operator[](i)));
	    ge0->SetPointEYlow(i, sqrt(ge0->GetErrorYlow(i)*ge0->GetErrorYlow(i)   + vecObs[3]->operator[](i)*vecObs[3]->operator[](i)));
	  }

      
      geStepTh_nL->GetYaxis()->SetRangeUser(-1.0,1.0);
      geStepTh_nL->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});A_{FB}");

      geStepTh_yL->GetYaxis()->SetRangeUser(-1.0,1.0);
      geStepTh_yL->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});A_{FB}");
    }
  else if (PlotType == 12) // Branching fraction
    {
      // ############################################
      // # Graph containing only statistical errors #
      // ############################################
      ge00 = new TGraphAsymmErrors(*ge0);
      ge00->SetMarkerColor(kBlack);
      ge00->SetMarkerStyle(20);
      ge00->SetMarkerSize(1.2);
      ge00->SetLineColor(kBlack);
      ge00->SetLineWidth(2);
      ge00->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);


      // ######################################################
      // # Graph containing statistical and systematic errors #
      // ######################################################
      for (int i = 0; i < ge0->GetN(); i++)
	if (Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false)
	  {
	    ge0->SetPointEXlow(i,0.0);
	    ge0->SetPointEXhigh(i,0.0);

	    ge0->SetPointEYhigh(i,sqrt(ge0->GetErrorYhigh(i)*ge0->GetErrorYhigh(i) + vecObs[8]->operator[](i)*vecObs[8]->operator[](i)));
	    ge0->SetPointEYlow(i, sqrt(ge0->GetErrorYlow(i)*ge0->GetErrorYlow(i)   + vecObs[9]->operator[](i)*vecObs[9]->operator[](i)));
	  }


      // #################################################################
      // # Divide the theoretical branching fraction by the q2 bin width #
      // #################################################################
      for (unsigned int i = 0; i < q2Bins.size()-1; i++)
        if ((Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false) && (geStepTh_nL->GetErrorYlow(i) != 0.0) && (geStepTh_nL->GetErrorYhigh(i) != 0.0))
      	  {
      	    geStepTh_nL->SetPoint(i,geStepTh_nL->GetX()[i],geStepTh_nL->GetY()[i] /  (q2Bins[i+1] - q2Bins[i]));
      	    geStepTh_nL->SetPointEYlow(i,geStepTh_nL->GetErrorYlow(i) / (q2Bins[i+1] - q2Bins[i]));
      	    geStepTh_nL->SetPointEYhigh(i,geStepTh_nL->GetErrorYhigh(i) / (q2Bins[i+1] - q2Bins[i]));

	  }
      
      geStepTh_nL->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);
      geStepTh_nL->GetYaxis()->SetRangeUser(0.0,12.0);
      geStepTh_nL->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}8}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");


      // #################################################################
      // # Divide the theoretical branching fraction by the q2 bin width #
      // #################################################################
      for (unsigned int i = 0; i < q2Bins.size()-1; i++)
        if ((Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false) && (geStepTh_yL->GetErrorYlow(i) != 0.0) && (geStepTh_yL->GetErrorYhigh(i) != 0.0))
      	  {
	    geStepTh_yL->SetPoint(i,geStepTh_yL->GetX()[i],geStepTh_yL->GetY()[i] /  (q2Bins[i+1] - q2Bins[i]));
	    geStepTh_yL->SetPointEYlow(i,geStepTh_yL->GetErrorYlow(i) / (q2Bins[i+1] - q2Bins[i]));
	    geStepTh_yL->SetPointEYhigh(i,geStepTh_yL->GetErrorYhigh(i) / (q2Bins[i+1] - q2Bins[i]));
	  }


      geStepTh_yL->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);
      geStepTh_yL->GetYaxis()->SetRangeUser(0.0,12.0);
      geStepTh_yL->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});dBF/dq#lower[0.4]{^{2}} (10#lower[0.4]{^{#font[122]{\55}8}} #times GeV#lower[0.4]{^{#font[122]{\55}2}})");
    }
  else if (PlotType == 13) // Fs
    {
      // ############################################
      // # Graph containing only statistical errors #
      // ############################################
      ge00 = new TGraphAsymmErrors(*ge0);
      ge00->SetMarkerColor(kBlack);
      ge00->SetMarkerStyle(20);
      ge00->SetMarkerSize(1.2);
      ge00->SetLineColor(kBlack);
      ge00->SetLineWidth(2);
      ge00->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);
      
      
      for (unsigned int i = 0; i < q2Bins.size()-1; i++) geStepTh_nL->SetPoint(i,geStepTh_nL->GetX()[i],YvalueOutsideLimits);
      geStepTh_nL->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);
      geStepTh_nL->GetYaxis()->SetRangeUser(0.0,1.0);
      geStepTh_nL->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});F_{S}");
    }
  else if (PlotType == 14) // As
    {
      // ############################################
      // # Graph containing only statistical errors #
      // ############################################
      ge00 = new TGraphAsymmErrors(*ge0);
      ge00->SetMarkerColor(kBlack);
      ge00->SetMarkerStyle(20);
      ge00->SetMarkerSize(1.2);
      ge00->SetLineColor(kBlack);
      ge00->SetLineWidth(2);
      ge00->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);

      
      for (unsigned int i = 0; i < q2Bins.size()-1; i++) geStepTh_nL->SetPoint(i,geStepTh_nL->GetX()[i],YvalueOutsideLimits);
      geStepTh_nL->GetXaxis()->SetLimits(q2Bins[0],q2Bins[q2Bins.size()-1]);
      geStepTh_nL->GetYaxis()->SetRangeUser(-1.0,1.0);
      geStepTh_nL->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});A_{S}");
    }
  

  // ################################
  // # Pad for the actual histogram #
  // ################################
  canv0->cd();
  if (FORPAPER == false)
    {
      pad1 = new TPad("pad1","pad1",0,0.5,1,1);
      pad1->SetBottomMargin(0);
    }
  else pad1 = new TPad("pad1","pad1",0,0,1,1);
  pad1->Draw();
  pad1->cd();


  TLegend* leg = new TLegend(0.8, 0.7, 0.97, 0.89, "");
  if ((PlotType == 0) || (PlotType == 1) || (PlotType == 2)) // Fl OR Afb OR Branching fraction
    {
      ge0->Draw("ape1");
      ge11->Draw("same e2");
      ge1->Draw("same pe1");

      leg->AddEntry(ge0,"GEN-MC","PL");
      leg->AddEntry(ge11,"RECO-MC","EPFL");
    }
  else if ((PlotType == 10) || (PlotType == 11) || (PlotType == 12)) // Fl OR Afb OR Branching fraction
    {
      geStepTh_nL->Draw("ae2");
      geStepTh_yL->Draw("same e2");
      ge00->Draw("same pe1");
      // ###################################
      // # Code to make slanted error bars #
      // ###################################
      vector<double> exld;
      vector<double> exhd;
      vector<double> eyld;
      vector<double> eyhd;
      for (int i = 0; i < ge0->GetN(); i++)
	{
	  exld.push_back(0.0);
	  exhd.push_back(0.0);
	  eyld.push_back(0.4);
	  eyhd.push_back(-0.4);
	}
      geb = new TGraphBentErrors(ge0->GetN(), ge0->GetX(), ge0->GetY(), ge0->GetEXlow(), ge0->GetEXhigh(), ge0->GetEYlow(), ge0->GetEYhigh(), &exld[0], &exhd[0], &eyld[0], &eyhd[0]);
      geb->SetMarkerColor(kBlack);
      geb->SetMarkerStyle(20);
      geb->SetMarkerSize(1.2);
      geb->SetLineColor(kBlack);
      geb->SetLineWidth(2);
      geb->SetLineStyle(kDashed);
      exld.clear();
      exhd.clear();
      eyld.clear();
      eyhd.clear();
      // #####################
      // # Bended error bars #
      // #####################
      // @TMP@
      // geb->Draw("same pe1");
      ge0->Draw("same pez");

      leg->AddEntry(ge00,"Data","EPL");
      leg->AddEntry(geStepTh_nL,"#LT SM #GT (LCSR)","F");
      leg->AddEntry(geStepTh_yL,"#LT SM #GT (Lattice)","F");
    }
  else if ((PlotType == 13) || (PlotType == 14)) // Fs OR As
    {
      geStepTh_nL->Draw("ae2");
      ge00->Draw("same pe1");

      leg->AddEntry(ge00,"Data","EPL");
    }
  leg->SetFillStyle(0);
  leg->SetFillColor(0);
  leg->SetBorderSize(0);
  leg->Draw();


  // ###################################
  // # Extract AFB zero crossing point #
  // ###################################
  if ((PlotType == 1) || (PlotType == 11))
    {      
      TF1* ZeroCrox = new TF1("ZeroCrox","[0]*x + [1]",q2Bins[0],q2Bins[JPsibin]);
      ZeroCrox->SetParName(0,"a");
      ZeroCrox->SetParName(1,"b");
      ZeroCrox->SetLineColor(kBlack);
      ZeroCrox->SetLineWidth(2);

      ZeroCrox->SetParameter(0,4.0);
      ZeroCrox->SetParameter(1,-0.3);

      TFitResultPtr fitResults = ge0->Fit("ZeroCrox","S V E MR0");
      TMatrixTSym<double> covMatrix(fitResults->GetCovarianceMatrix());

      double q0  = -ZeroCrox->GetParameter(1) / ZeroCrox->GetParameter(0);
      double q0E = q0 * sqrt(pow(ZeroCrox->GetParError(0) / ZeroCrox->GetParameter(0),2.) + pow(ZeroCrox->GetParError(1) / ZeroCrox->GetParameter(1),2.) - 2./(ZeroCrox->GetParameter(0) * ZeroCrox->GetParameter(1)) * covMatrix(0,1));
      cout << "\n@@@ Zero crossing point: " << q0 << " +/- " << q0E << " @@@" << endl;
      cout << "@@@ p-value (for SM compatibility): " << TMath::Erfc(fabs(q0 - q0SM) / sqrt(q0E*q0E + q0SME*q0SME)) << " @@@" << endl;


      // #################################
      // # Print AFB zero crossing point #
      // #################################
      if (FORPAPER == false)
	{
	  ZeroCrox->Draw("same");

	  paveText = new TPaveText(0.12,0.05,0.27,0.18,"NDC");
	  paveText->SetTextAlign(11);
	  paveText->SetBorderSize(0.0);
	  paveText->SetFillColor(kWhite);
	  paveText->AddText(Form("Zero crossing: %.1f#pm%.1f",q0,q0E));
	  paveText->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",ZeroCrox->GetChisquare()/static_cast<double>(ZeroCrox->GetNDF())));
	  paveText->AddText(Form("%s%.2f","p-value = ",TMath::Prob(ZeroCrox->GetChisquare(),ZeroCrox->GetNDF())));
	  paveText->SetFillStyle(0);
	  paveText->SetTextSize(0.035);
	  paveText->Draw();
	}
    }


  // ########################
  // # Draw exclusion zones #
  // ########################
  DrawExclusion(q2Bins[JPsibin],q2Bins[JPsibin+1],-1.2,30.,"RejectJPsi1",3001,kGray);
  DrawExclusion(q2Bins[PsiPbin],q2Bins[PsiPbin+1],-1.2,30.,"RejectPsiP1",3001,kGray);


  // ####################
  // # Pad for the chi2 #
  // ####################
  canv0->cd();
  pad2 = new TPad("pad2","pad2",0,0.25,1,0.5);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0);
  pad2->SetLogy();
  if (FORPAPER == false) pad2->Draw();
  pad2->cd();


  TH1D* chi2Histo = new TH1D("chi2Histo","chi2Histo",q2Bins.size()-1,q2Bins_);
  TLegend* ratioLeg = new TLegend(0.8, 0.85, 0.97, 0.95, "");
  if ((PlotType == 0) || (PlotType == 1) || (PlotType == 2)) // Fl OR Afb OR Branching fraction
    {
      double tmpVar;
      for (int i = 0; i < chi2Histo->GetNbinsX(); i++)
	{
	  if (Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false)
	    {
	      tmpVar = pow(ge0->GetY()[i] - ge1->GetY()[i],2.) /
		(ge0->GetY()[i] > ge1->GetY()[i] ? pow(ge0->GetErrorYlow(i),2.) + pow(ge1->GetErrorYhigh(i),2.) : pow(ge0->GetErrorYhigh(i),2.) + pow(ge1->GetErrorYlow(i),2.));

	      myGlobalChi2 = myGlobalChi2 + tmpVar;
	      chi2Histo->SetBinContent(i+1,tmpVar);

	      DoF++;
	    }
	  else chi2Histo->SetBinContent(i+1,0.0);
	}

      ratioLeg->AddEntry(chi2Histo,"#chi#lower[0.4]{^{2}}(RECO, GEN)");
    }
  else if ((PlotType == 10) || (PlotType == 11) || (PlotType == 12)) // Fl OR Afb OR Branching fraction
    {      
      double tmpVar;
      for (int i = 0; i < chi2Histo->GetNbinsX(); i++)
	{
	  if ((Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false) && (geStepTh_nL->GetErrorYlow(i) != 0.0) && (geStepTh_nL->GetErrorYhigh(i) != 0.0))
	    {
	      tmpVar = pow(geStepTh_nL->GetY()[i] - ge0->GetY()[i],2.) /
		(geStepTh_nL->GetY()[i] > ge0->GetY()[i] ? pow(geStepTh_nL->GetErrorYlow(i),2.) + pow(ge0->GetErrorYhigh(i),2.) : pow(geStepTh_nL->GetErrorYhigh(i),2.) + pow(ge0->GetErrorYlow(i),2.));
	      
	      myGlobalChi2 = myGlobalChi2 + tmpVar;
	      chi2Histo->SetBinContent(i+1,tmpVar);

	      DoF++;
	    }
	  else chi2Histo->SetBinContent(i+1,0.0);
	}
      
      ratioLeg->AddEntry(chi2Histo,"#chi#lower[0.4]{^{2}}(#LT SM #GT, Data)");
    }
  cout << "\n@@@ Global chi2 = " << myGlobalChi2 / static_cast<double>(DoF) << " (" << myGlobalChi2 << "/" << static_cast<double>(DoF) << ") @@@" << endl;
  myGlobalChi2 = myGlobalChi2 / static_cast<double>(DoF);

  chi2Histo->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});#chi#lower[0.4]{^{2}}");
  chi2Histo->GetXaxis()->SetLabelSize(0.06);
  chi2Histo->GetXaxis()->SetTitleOffset(0.8);
  chi2Histo->GetXaxis()->SetTitleSize(0.07);
  chi2Histo->GetYaxis()->SetLabelSize(0.06);
  chi2Histo->GetYaxis()->SetTitleOffset(0.7);
  chi2Histo->GetYaxis()->SetTitleSize(0.07);
  chi2Histo->SetFillColor(kAzure+6);
  chi2Histo->SetFillStyle(1001);
  chi2Histo->GetYaxis()->SetRangeUser(0.008,30.0);
  chi2Histo->Draw();

  ratioLeg->SetFillColor(0);
  ratioLeg->SetBorderSize(1);
  ratioLeg->Draw();


  // ########################
  // # Draw exclusion zones #
  // ########################
  DrawExclusion(q2Bins[JPsibin],q2Bins[JPsibin+1],0.008,30.0,"RejectJPsi2",3001,kGray);
  DrawExclusion(q2Bins[PsiPbin],q2Bins[PsiPbin+1],0.008,30.0,"RejectPsiP2",3001,kGray);


  // ###########################
  // # Pad for the probability #
  // ###########################
  canv0->cd();
  pad3 = new TPad("pad3","pad3",0,0.0,1,0.25);
  pad3->SetTopMargin(0);
  pad3->SetBottomMargin(0.15);
  pad3->SetLogy();
  if (FORPAPER == false) pad3->Draw();
  pad3->cd();

  TH1D* probHisto = new TH1D("probHisto","probHisto",q2Bins.size()-1,q2Bins_);
  TLegend* probLeg = new TLegend(0.8, 0.85, 0.97, 0.95, "");
  for (int i = 0; i < chi2Histo->GetNbinsX(); i++)
    if (Utility->ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false) probHisto->SetBinContent(i+1,TMath::Prob(chi2Histo->GetBinContent(i+1),1));
  probLeg->AddEntry(probHisto,"p-value");

  probHisto->SetTitle(";q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}});p-value");
  probHisto->GetXaxis()->SetLabelSize(0.06);
  probHisto->GetXaxis()->SetTitleOffset(0.8);
  probHisto->GetXaxis()->SetTitleSize(0.07);
  probHisto->GetYaxis()->SetLabelSize(0.06);
  probHisto->GetYaxis()->SetTitleOffset(0.7);
  probHisto->GetYaxis()->SetTitleSize(0.07);
  probHisto->SetFillColor(kGreen-7);
  probHisto->SetFillStyle(1001);
  probHisto->GetYaxis()->SetRangeUser(1e-3,1.2);
  probHisto->Draw();

  probLeg->SetFillColor(0);
  probLeg->SetBorderSize(1);
  probLeg->Draw();


  // ########################
  // # Draw exclusion zones #
  // ########################
  DrawExclusion(q2Bins[JPsibin],q2Bins[JPsibin+1],1e-3,1.13,"RejectJPsi3",3001,kGray);
  DrawExclusion(q2Bins[PsiPbin],q2Bins[PsiPbin+1],1e-3,1.13,"RejectPsiP3",3001,kGray);
  DrawExclusion(q2Bins[0],q2Bins[q2Bins.size()-1],1e-3,0.05,"CL95",3002,kRed-9);


  // #################################
  // # Write global chi2 and p-value #
  // #################################
  pad1->cd();
  if (FORPAPER == false)
    {
      paveText = new TPaveText(0.75,0.05,0.95,0.18,"NDC");
      paveText->SetTextAlign(11);
      paveText->SetBorderSize(0.0);
      paveText->SetFillColor(kWhite);
      paveText->AddText("Global compatibility");
      paveText->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myGlobalChi2));
      paveText->AddText(Form("%s%.2f","p-value = ",TMath::Prob(myGlobalChi2*static_cast<double>(DoF),DoF)));
      paveText->SetFillStyle(0);
      paveText->SetTextSize(0.035);
      paveText->Draw();
    }
  DrawString(LUMI);


  // #################################
  // # Draw horizontal line at y = 0 #
  // #################################
  if ((PlotType == 1) || (PlotType == 11) || (PlotType == 14)) // Afb OR As
    {
      line = new TLine(q2Bins[0],0.0,q2Bins[q2Bins.size()-1],0.0);
      line->SetLineStyle(kDashed);
      line->Draw();
    }
  

  canv0->cd();
  canv0->Modified();
  canv0->Update();
}


void GenNTupleFromMultyRun (string fileName, unsigned int q2BinIndx)
{
  ifstream inputFile;

  unsigned int nPar = 18;

  double ID;
  double fit_Fl,  errorHi_Fl,  errorLo_Fl,  pdf_Fl;
  double fit_Afb, errorHi_Afb, errorLo_Afb, pdf_Afb;
  double fit_BF,  error_BF, pdf_BF;
  double effMuMuGoodTag, effMuMuMisTag;
  double fitOrg_Fl,  errorHiOrg_Fl,  errorLoOrg_Fl;
  double fitOrg_Afb, errorHiOrg_Afb, errorLoOrg_Afb;
  double nll;

  vector<double>                vecVar;
  vector<string>                ParVector;
  vector<vector<string>*>       fitParam;
  vector<vector<unsigned int>*> configParam;

  ReadParameters* ParameterFile = new ReadParameters(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str());


  // ###################
  // # Read parameters #
  // ###################
  Utility->ReadFitStartingValues(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&fitParam,&configParam,Utility->ParFileBlockN("fitValBins"));
  ParameterFile->ReadFromFile(Utility->ParFileBlockN("BF"),&ParVector);
  pdf_Fl  = atof(fitParam[Utility->GetFitParamIndx("FlS")]->operator[](q2BinIndx).c_str());
  pdf_Afb = atof(fitParam[Utility->GetFitParamIndx("AfbS")]->operator[](q2BinIndx).c_str());
  pdf_BF  = atof(ParVector[q2BinIndx].c_str());


  // #####################
  // # NTuple definition #
  // #####################
  TTree* FitResults = new TTree("FitResults","Systematic-error fit results");

  FitResults->Branch("ID",&ID,"ID/D");

  FitResults->Branch("fit_Fl",    &fit_Fl,    "fit_Fl/D");
  FitResults->Branch("errorHi_Fl",&errorHi_Fl,"errorHi_Fl/D");
  FitResults->Branch("errorLo_Fl",&errorLo_Fl,"errorLo_Fl/D");
  FitResults->Branch("pdf_Fl",    &pdf_Fl,    "pdf_Fl/D");

  FitResults->Branch("fit_Afb",    &fit_Afb,    "fit_Afb/D");
  FitResults->Branch("errorHi_Afb",&errorHi_Afb,"errorHi_Afb/D");
  FitResults->Branch("errorLo_Afb",&errorLo_Afb,"errorLo_Afb/D");
  FitResults->Branch("pdf_Afb",    &pdf_Afb,    "pdf_Afb/D");

  FitResults->Branch("fit_BF",  &fit_BF,  "fit_BF/D");
  FitResults->Branch("error_BF",&error_BF,"error_BF/D");
  FitResults->Branch("pdf_BF",  &pdf_BF,  "pdf_BF/D");
  
  FitResults->Branch("effMuMuGoodTag",&effMuMuGoodTag,"effMuMuGoodTag/D");
  FitResults->Branch("effMuMuMisTag", &effMuMuMisTag, "effMuMuMisTag/D");

  FitResults->Branch("fitOrg_Fl",    &fitOrg_Fl,    "fitOrg_Fl/D");
  FitResults->Branch("errorHiOrg_Fl",&errorHiOrg_Fl,"errorHiOrg_Fl/D");
  FitResults->Branch("errorLoOrg_Fl",&errorLoOrg_Fl,"errorLoOrg_Fl/D");

  FitResults->Branch("fitOrg_Afb",    &fitOrg_Afb,    "fitOrg_Afb/D");
  FitResults->Branch("errorHiOrg_Afb",&errorHiOrg_Afb,"errorHiOrg_Afb/D");
  FitResults->Branch("errorLoOrg_Afb",&errorLoOrg_Afb,"errorLoOrg_Afb/D");

  FitResults->Branch("nll",&nll,"nll/D");


  // #########################
  // # Read values from file #
  // #########################
  cout << "\n[MakePlots::GenNTupleFromMultyRun]\t@@@ Opening systematics file: " << fileName.c_str() << " @@@" << endl;
  inputFile.open(fileName.c_str(), ifstream::in);
  if (inputFile.good() == false)
    {
      cout << "[MakePlots::GenNTupleFromMultyRun]\tError opening file : " << fileName.c_str() << endl;
      exit (EXIT_FAILURE);
    }


  for (unsigned int i = 0; i < nPar; i++)
    {
      vecVar.push_back(0.0);
      inputFile >> vecVar.back();
      cout << "var" << i << ": " << vecVar.back() << "\t";
    }
  cout << endl;

  while (inputFile.eof() == false)
    {
      if (vecVar[1] != -2.0)
	{
	  ID          = vecVar[0];

	  fit_Fl      = vecVar[1];
	  errorHi_Fl  = vecVar[2];
	  errorLo_Fl  = vecVar[3];

	  fit_Afb     = vecVar[4];
	  errorHi_Afb = vecVar[5];
	  errorLo_Afb = vecVar[6];

	  fit_BF      = vecVar[7];
	  error_BF    = vecVar[8];

	  effMuMuGoodTag = vecVar[9];
	  effMuMuMisTag  = vecVar[10];

	  fitOrg_Fl      = vecVar[11];
	  errorHiOrg_Fl  = vecVar[12];
	  errorLoOrg_Fl  = vecVar[13];

	  fitOrg_Afb     = vecVar[14];
	  errorHiOrg_Afb = vecVar[15];
	  errorLoOrg_Afb = vecVar[16];

	  nll         = vecVar[17];

	  FitResults->Fill();
	}

      for (unsigned int i = 0; i < nPar; i++)
	{
	  vecVar[i];
	  inputFile >> vecVar[i];
	  cout << "var" << i << ": " << vecVar[i] << "\t";
	}
      cout << endl;
    }


  // ###############
  // # Save ntuple #
  // ###############
  fileName.replace(fileName.find(".txt"),4,".root");
  TFile* fileID = new TFile(fileName.c_str(),"RECREATE");
  fileID->cd();
  FitResults->Write();
  fileID->Close();
  delete fileID;

  inputFile.close();
  ParVector.clear();
  delete ParameterFile;
}


void PlotMuMu (string fileName, bool bkgSub)
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  SetStyle();

  
  int nEntries;
  double minX = 0.8;
  double maxX = 5.0;
  unsigned int nBins = 300;
  stringstream myString;
  string sigMassQuery = "";
  string bkgMassQuery = "";

  double signalSigma = sqrt( atof(Utility->GetGenericParam("FRACMASSS").c_str()) * atof(Utility->GetGenericParam("SIGMAS1").c_str()) * atof(Utility->GetGenericParam("SIGMAS1").c_str()) +
			     (1. - atof(Utility->GetGenericParam("FRACMASSS").c_str())) * atof(Utility->GetGenericParam("SIGMAS2").c_str()) * atof(Utility->GetGenericParam("SIGMAS2").c_str()) );
  cout << "\n[MakePlots::PlotMuMu]\t@@@ Signal sigma: " << signalSigma << " @@@" << endl;


  if (bkgSub == true)
    {
      myString.clear(); myString.str("");
      myString << "(abs(B0MassArb - " << Utility->B0Mass << ") < " << atof(Utility->GetGenericParam("NSigmaB0S").c_str())*signalSigma << ")";
      sigMassQuery = myString.str();
      
      myString.clear(); myString.str("");
      myString << "((B0MassArb > " << Utility->B0Mass+atof(Utility->GetGenericParam("NSigmaB0B").c_str())*signalSigma << " && B0MassArb < "
	       << Utility->B0Mass+(atof(Utility->GetGenericParam("NSigmaB0B").c_str())+atof(Utility->GetGenericParam("NSigmaB0S").c_str()))*signalSigma << ") || ";
      myString << "(B0MassArb < " << Utility->B0Mass-atof(Utility->GetGenericParam("NSigmaB0B").c_str())*signalSigma << " && B0MassArb > "
	       << Utility->B0Mass-(atof(Utility->GetGenericParam("NSigmaB0B").c_str())+atof(Utility->GetGenericParam("NSigmaB0S").c_str()))*signalSigma << "))";
      bkgMassQuery = myString.str(); 
    }
  else
    {
      myString.clear(); myString.str("");
      myString << "(B0MassArb > " << Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str()) << " && B0MassArb < " << Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str()) << ")";
      sigMassQuery = myString.str();

      myString.clear(); myString.str("");
      myString << "(B0MassArb > " << Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str()) << " && B0MassArb < " << Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str()) << ")";
      bkgMassQuery = myString.str(); 
    }


  TFile* file0 = TFile::Open(fileName.c_str(),"READ");
  TTree* theTree = (TTree*)file0->Get("B0KstMuMu/B0KstMuMuNTuple");

  nEntries = theTree->GetEntries();
  cout << "\n[MakePlots::PlotMuMu]\t@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;


  TCanvas* c0 = new TCanvas("c0","c0",10,10,700,500);
  c0->cd();

  TH1D* hDsig = new TH1D("hDsig","hDsig",nBins,minX,maxX);
  hDsig->SetXTitle("m(#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{\55}}}}) (GeV)");
  hDsig->SetYTitle("Entries / (0.014 GeV)");
  hDsig->SetMarkerStyle(20);

  TH1D* hDbkg = new TH1D("hDbkg","hDbkg",nBins,minX,maxX);
  hDbkg->SetXTitle("m(#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{\55}}}}) (GeV)");
  hDbkg->SetYTitle("Entries / (0.014 GeV)");
  hDbkg->SetMarkerStyle(20);

  theTree->Draw("mumuMass>>hDsig",sigMassQuery.c_str(),"goff");
  theTree->Draw("mumuMass>>hDbkg",bkgMassQuery.c_str(),"goff");

  if (bkgSub == true) hDsig->Add(hDbkg, -1.0);
 
  hDsig->Draw("e1p");
  c0->Modified();
  c0->Update();
}


void PlotKst (string fileName, bool bkgSub, bool fitParamAreFixed)
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  SetStyle();

  
  int nEntries;
  double minX  = 0.81;
  double maxX  = 0.98;
  double extra = 0.015;
  unsigned int nBins = 50;
  stringstream myString;
  string sigMassQuery = "";
  string bkgMassQuery = "";
  string tmpstring    = "";

  double signalSigma = sqrt( atof(Utility->GetGenericParam("FRACMASSS").c_str()) * atof(Utility->GetGenericParam("SIGMAS1").c_str()) * atof(Utility->GetGenericParam("SIGMAS1").c_str()) +
			     (1. - atof(Utility->GetGenericParam("FRACMASSS").c_str())) * atof(Utility->GetGenericParam("SIGMAS2").c_str()) * atof(Utility->GetGenericParam("SIGMAS2").c_str()) );
  cout << "\n[MakePlots::PlotKst]\t@@@ Signal sigma: " << signalSigma << " @@@" << endl;


  if (bkgSub == true)
    {
      myString.clear(); myString.str("");
      myString << "(abs(B0MassArb - " << Utility->B0Mass << ") < " << atof(Utility->GetGenericParam("NSigmaB0S").c_str())*signalSigma << ")";
      sigMassQuery = myString.str();
      
      myString.clear(); myString.str("");
      myString << "((B0MassArb > " << Utility->B0Mass+atof(Utility->GetGenericParam("NSigmaB0B").c_str())*signalSigma << " && B0MassArb < "
	       << Utility->B0Mass+(atof(Utility->GetGenericParam("NSigmaB0B").c_str())+atof(Utility->GetGenericParam("NSigmaB0S").c_str()))*signalSigma << ") || ";
      myString << "(B0MassArb < " << Utility->B0Mass-atof(Utility->GetGenericParam("NSigmaB0B").c_str())*signalSigma << " && B0MassArb > "
	       << Utility->B0Mass-(atof(Utility->GetGenericParam("NSigmaB0B").c_str())+atof(Utility->GetGenericParam("NSigmaB0S").c_str()))*signalSigma << "))";
      bkgMassQuery = myString.str(); 
    }
  else
    {
      myString.clear(); myString.str("");
      myString << "(B0MassArb > " << Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str()) << " && B0MassArb < " << Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str()) << ")";
      sigMassQuery = myString.str();

      myString.clear(); myString.str("");
      myString << "(B0MassArb > " << Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str()) << " && B0MassArb < " << Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str()) << ")";
      bkgMassQuery = myString.str(); 
    }


  TCanvas* c0 = new TCanvas("c0","c0",10,10,700,500);
  TCanvas* c1 = new TCanvas("c1","c1",20,20,700,500);

  TFile* file0 = TFile::Open(fileName.c_str(),"READ");
  TTree* theTree = (TTree*)file0->Get("B0KstMuMu/B0KstMuMuNTuple");

  nEntries = theTree->GetEntries();
  cout << "\n[MakePlots::PlotKst]\t@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;

  TH1D* h1Dsig = new TH1D("h1Dsig","h1Dsig",nBins,minX - extra,maxX + extra);
  h1Dsig->SetXTitle("m(#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}}) (GeV)");
  h1Dsig->SetYTitle("Entries / (0.004 GeV)");
  h1Dsig->SetMarkerStyle(20);

  TH1D* h2Dsig = new TH1D("h2Dsig","h2Dsig",nBins,minX - extra,maxX + extra);
  h2Dsig->SetXTitle("m(#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}}) (GeV)");
  h2Dsig->SetYTitle("Entries / (0.004 GeV)");
  h2Dsig->SetMarkerStyle(21);
  h2Dsig->SetMarkerColor(kRed);


  TH1D* h1Dbkg = new TH1D("h1Dbkg","h1Dbkg",nBins,minX - extra,maxX + extra);
  h1Dbkg->SetXTitle("m(#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}}) (GeV)");
  h1Dbkg->SetYTitle("Entries / (0.004 GeV)");
  h1Dbkg->SetMarkerStyle(20);

  TH1D* h2Dbkg = new TH1D("h2Dbkg","h2Dbkg",nBins,minX - extra,maxX + extra);
  h2Dbkg->SetXTitle("m(#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}}) (GeV)");
  h2Dbkg->SetYTitle("Entries / (0.004 GeV)");
  h2Dbkg->SetMarkerStyle(21);
  h2Dbkg->SetMarkerColor(kRed);

  TF1* fsig = new TF1("fsig","[6]*([2]*TMath::Gaus(x,[0],[1]) + [5]*TMath::Gaus(x,[3],[4]))",minX,maxX);
  TF1* fbkg = new TF1("fbkg","[0] + [1]*x",minX,maxX);
  TF1* f0 = new TF1("f0","fsig + fbkg",minX,maxX);

  fsig->SetLineColor(kBlue);
  fsig->SetFillColor(kAzure+6);
  fsig->SetFillStyle(3345);

  fbkg->SetLineColor(kRed);
  fbkg->SetFillColor(kRed);
  fbkg->SetFillStyle(3354);

  f0->SetParName(0,"#mu-1");
  f0->SetParName(1,"#sigma-1");
  f0->SetParName(2,"Ampli-S1");
  f0->SetParName(3,"#mu-2");
  f0->SetParName(4,"#sigma-2");
  f0->SetParName(5,"Ampli-S2");
  f0->SetParName(6,"Ampli-S");

  f0->SetParName(7,"Ampli-B");
  f0->SetParName(8,"Slope-B");


  if (fitParamAreFixed == true)
    {
      // #######################################
      // # From fit to RECO B0 --> J/psi K* MC #
      // #######################################
      f0->FixParameter(0, 8.94440e-01);
      f0->FixParameter(1, 1.27969e-02);
      f0->FixParameter(2, 8.36339e+03);
      f0->FixParameter(3, 8.94996e-01);
      f0->FixParameter(4, 2.99200e-02);
      f0->FixParameter(5, 1.01805e+04);
    }
  else
    {
      // ####################
      // # Starting valiues #
      // ####################
      f0->SetParameter(0, 0.892);
      f0->SetParameter(1, 0.01);
      f0->SetParameter(2, 6000);
      f0->SetParameter(3, 0.892);
      f0->SetParameter(4, 0.03);
      f0->SetParameter(5, 5000);
    }
  f0->SetParameter(6, 1.0);
  f0->SetParameter(7, -4000);
  f0->SetParameter(8,  6000);


  tmpstring = "(B0notB0bar == 1)";
  tmpstring += " && " + sigMassQuery;
  theTree->Draw("kstMass>>h1Dsig",tmpstring.c_str());

  tmpstring = "(B0notB0bar == 0)";
  tmpstring += " && " + sigMassQuery;
  theTree->Draw("kstBarMass>>h2Dsig",tmpstring.c_str());

  TH1D* h3Dsig = (TH1D*)h1Dsig->Clone("h3Dsig");
  h3Dsig->SetXTitle("m(#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}}) (GeV)");
  h3Dsig->SetYTitle("Entries / (0.004 GeV)");
  h3Dsig->SetMarkerStyle(20);

  h3Dsig->Add(h2Dsig, 1.0);


  tmpstring = "(B0notB0bar == 1)";
  tmpstring += " && " + bkgMassQuery;
  theTree->Draw("kstMass>>h1Dbkg",tmpstring.c_str());

  tmpstring = "(B0notB0bar == 0)";
  tmpstring += " && " + bkgMassQuery;
  theTree->Draw("kstBarMass>>h2Dbkg",tmpstring.c_str());

  TH1D* h3Dbkg = (TH1D*)h1Dbkg->Clone("h3Dbkg");
  h3Dbkg->SetXTitle("m(#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}}) (GeV)");
  h3Dbkg->SetYTitle("Entries / (0.004 GeV)");
  h3Dbkg->SetMarkerStyle(20);

  h3Dbkg->Add(h2Dbkg, 1.0);

  if (bkgSub == true)
    {
      h1Dsig->Add(h1Dbkg, -1.0);
      h2Dsig->Add(h2Dbkg, -1.0);
      h3Dsig->Add(h3Dbkg, -1.0);
    }


  c0->cd();
  h3Dsig->Fit("f0","R0");
  h3Dsig->Draw("e1");
  h3Dsig->GetFunction("f0")->Draw("same");

  fsig->SetParameter(0,f0->GetParameter(0));
  fsig->SetParameter(1,f0->GetParameter(1));
  fsig->SetParameter(2,f0->GetParameter(2));
  fsig->SetParameter(3,f0->GetParameter(3));
  fsig->SetParameter(4,f0->GetParameter(4));
  fsig->SetParameter(5,f0->GetParameter(5));
  fsig->SetParameter(6,f0->GetParameter(6));
  fsig->Draw("same");

  fbkg->SetParameter(0,f0->GetParameter(7));
  fbkg->SetParameter(1,f0->GetParameter(8));
  fbkg->Draw("same");

  TLegend* leg0 = new TLegend(0.15, 0.6, 0.35, 0.85, "");
  leg0->AddEntry(h3Dsig,"Data");
  leg0->AddEntry(f0,"Total fit");
  leg0->AddEntry(fsig,"Signal");
  leg0->AddEntry(fbkg,"Background");
  leg0->SetFillColor(0);
  leg0->SetBorderSize(0);
  leg0->Draw();

  c0->Modified();
  c0->Update();


  c1->cd();
  h1Dsig->Draw("e1p");
  h2Dsig->Draw("e1p sames");

  TLegend* leg1 = new TLegend(0.15, 0.6, 0.25, 0.85, "");
  leg1->AddEntry(h1Dsig,"#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}}");
  leg1->AddEntry(h2Dsig,"#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{*0}}}}#lower[-0.4]{_{bar}}");
  leg1->SetFillColor(0);
  leg1->SetBorderSize(0);
  leg1->Draw();

  c1->Modified();
  c1->Update();


  cout << "\nSignal integral [" << minX << ";" << maxX << "] = "  << fsig->Integral(minX,maxX) / (static_cast<double>(maxX + extra - (minX - extra)) / static_cast<double>(nBins)) << endl;
  cout << "Backgrund integral [" << minX << ";" << maxX << "] = " << fbkg->Integral(minX,maxX) / (static_cast<double>(maxX + extra - (minX - extra)) / static_cast<double>(nBins)) << endl;
}


void PlotKK (string fileName, bool bkgSub, string RECOorGEN)
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  SetStyle();
  gStyle->SetPadRightMargin(0.14);


  int nEntries;

  double KKminX      = 0.95;
  double KKmaxX      = 1.45;
  double KstminX     = 0.800;
  double KstmaxX     = 0.990;
  double Kstextra    = 0.015;
  double DalitzminX  = 0.65;
  double DalitzmaxX  = 0.95;
  double DalitzminY  = 14.5;
  double DalitzmaxY  = 27.0;
  unsigned int nBins = 100;
  double KKmass      = 0.0;
  double massPsiK    = 0.0;
  double massKpi     = 0.0;

  double signalSigma = sqrt( atof(Utility->GetGenericParam("FRACMASSS").c_str()) * atof(Utility->GetGenericParam("SIGMAS1").c_str()) * atof(Utility->GetGenericParam("SIGMAS1").c_str()) +
			     (1. - atof(Utility->GetGenericParam("FRACMASSS").c_str())) * atof(Utility->GetGenericParam("SIGMAS2").c_str()) * atof(Utility->GetGenericParam("SIGMAS2").c_str()) );
  cout << "\n[MakePlots::PlotKK]\t@@@ Signal sigma: " << signalSigma << " @@@" << endl;


  TFile* file0 = TFile::Open(fileName.c_str(),"READ");
  TTree* theTree = (TTree*)file0->Get("B0KstMuMu/B0KstMuMuNTuple");

  B0KstMuMuSingleCandTreeContent* NTuple = new B0KstMuMuSingleCandTreeContent();
  NTuple->Init();
  NTuple->ClearNTuple();
  NTuple->SetBranchAddresses(theTree);
  
  nEntries = theTree->GetEntries();
  cout << "\n[MakePlots::PlotKK]\t@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;


  TCanvas* c0 = new TCanvas("c0","c0",10,10,700,500);
  TCanvas* c1 = new TCanvas("c1","c1",20,20,700,500);
  TCanvas* c2 = new TCanvas("c2","c2",30,30,700,500);

  TH1D* hKKSig = new TH1D("hKKSig","hKKSig",nBins,KKminX,KKmaxX);
  hKKSig->SetXTitle("m(#font[122]{K}#lower[0.6]{^{#font[122]{+}}}#font[122]{K}#lower[0.6]{^{#font[122]{\55}}}) (GeV)");
  hKKSig->SetYTitle("Entries / (0.005 GeV)");
  hKKSig->SetFillColor(kAzure+6);

  TH1D* hKKBkg = new TH1D("hKKBkg","hKKBkg",nBins,KKminX,KKmaxX);
  hKKBkg->SetXTitle("m(#font[122]{K}#lower[0.6]{^{#font[122]{+}}}#font[122]{K}#lower[0.6]{^{#font[122]{\55}}}) (GeV)");
  hKKBkg->SetYTitle("Entries / (0.005 GeV)");
  hKKBkg->SetFillColor(kAzure+6);

  TH1D* hKstSig = new TH1D("hKstSig","hKstSig",nBins,KstminX - Kstextra,KstmaxX + Kstextra);
  hKstSig->SetXTitle("m(K #pi) (GeV)");
  hKstSig->SetYTitle("Entries / (0.022 GeV)");
  hKstSig->SetMarkerStyle(20);

  TH1D* hKstBkg = new TH1D("hKstBkg","hKstBkg",nBins,KstminX - Kstextra,KstmaxX + Kstextra);
  hKstBkg->SetXTitle("m(K #pi) (GeV)");
  hKstBkg->SetYTitle("Entries / (0.022 GeV)");
  hKstSig->SetMarkerStyle(20);

  TH2D* hDalitzSig = new TH2D("hDalitzSig","hDalitzSig",nBins,DalitzminX,DalitzmaxX,nBins,DalitzminY,DalitzmaxY);
  hDalitzSig->SetXTitle("m#lower[0.4]{^{2}}(K #pi) (GeV#lower[0.4]{^{2}})");
  hDalitzSig->SetYTitle("m#lower[0.4]{^{2}}(J/#psi K) (GeV#lower[0.4]{^{2}})");
  hDalitzSig->SetZTitle("Entries / (0.03x0.125 (GeV#lower[0.4]{^{4}}))");

  TH2D* hDalitzBkg = new TH2D("hDalitzBkg","hDalitzBkg",nBins,DalitzminX,DalitzmaxX,nBins,DalitzminY,DalitzmaxY);
  hDalitzBkg->SetXTitle("m#lower[0.4]{^{2}}(K #pi) (GeV#lower[0.4]{^{2}})");
  hDalitzBkg->SetYTitle("m#lower[0.4]{^{2}}(J/#psi K) (GeV#lower[0.4]{^{2}})");
  hDalitzBkg->SetZTitle("Entries / (0.03x0.125 (GeV#lower[0.4]{^{4}}))");


  if      (RECOorGEN == "RECO") cout << "\n@@@ Using RECO information @@@" << endl;
  else if (RECOorGEN == "GEN")  cout << "\n@@@ Using GEN information @@@" << endl;
  else
    {
      cout << "[MakePlots::PlotKK]\tWrong parameter: " << RECOorGEN << endl;
      exit (EXIT_FAILURE);
    }


  for (int entry = 0; entry < nEntries; entry++)
    {
      theTree->GetEntry(entry);


      if (RECOorGEN == "RECO")
	{
	  // #################
	  // # Use RECO info #
	  // #################
	  KKmass = Utility->computeInvMass(NTuple->kstTrkmPx->at(0),NTuple->kstTrkmPy->at(0),NTuple->kstTrkmPz->at(0),Utility->kaonMass,
					   NTuple->kstTrkpPx->at(0),NTuple->kstTrkpPy->at(0),NTuple->kstTrkpPz->at(0),Utility->kaonMass);
	  
	  if (NTuple->B0notB0bar == true)
	    {
	      massPsiK = Utility->computeInvMass(NTuple->mumPx->at(0) + NTuple->mupPx->at(0),NTuple->mumPy->at(0) + NTuple->mupPy->at(0),NTuple->mumPz->at(0) + NTuple->mupPz->at(0),NTuple->mumuMass->at(0),
						 NTuple->kstTrkpPx->at(0),NTuple->kstTrkpPy->at(0),NTuple->kstTrkpPz->at(0),Utility->kaonMass);
	      massKpi = Utility->computeInvMass(NTuple->kstTrkmPx->at(0),NTuple->kstTrkmPy->at(0),NTuple->kstTrkmPz->at(0),Utility->pionMass,
						NTuple->kstTrkpPx->at(0),NTuple->kstTrkpPy->at(0),NTuple->kstTrkpPz->at(0),Utility->kaonMass);
	    }
	  else
	    {
	      massPsiK = Utility->computeInvMass(NTuple->mumPx->at(0) + NTuple->mupPx->at(0),NTuple->mumPy->at(0) + NTuple->mupPy->at(0),NTuple->mumPz->at(0) + NTuple->mupPz->at(0),NTuple->mumuMass->at(0),
						 NTuple->kstTrkmPx->at(0),NTuple->kstTrkmPy->at(0),NTuple->kstTrkmPz->at(0),Utility->kaonMass);
	      massKpi = Utility->computeInvMass(NTuple->kstTrkmPx->at(0),NTuple->kstTrkmPy->at(0),NTuple->kstTrkmPz->at(0),Utility->kaonMass,
						NTuple->kstTrkpPx->at(0),NTuple->kstTrkpPy->at(0),NTuple->kstTrkpPz->at(0),Utility->pionMass);
	    }
	}
      else if (RECOorGEN == "GEN")
	{
	  // ################
	  // # Use GEN info #
	  // ################
	  KKmass = Utility->computeInvMass(NTuple->genKstTrkmPx,NTuple->genKstTrkmPy,NTuple->genKstTrkmPz,Utility->kaonMass,
					   NTuple->genKstTrkpPx,NTuple->genKstTrkpPy,NTuple->genKstTrkpPz,Utility->kaonMass);
	      
	  if (NTuple->B0notB0bar == true)
	    {
	      massPsiK = Utility->computeInvMass(NTuple->genMumPx + NTuple->genMupPx,NTuple->genMumPy + NTuple->genMupPy,NTuple->genMumPz + NTuple->genMupPz,NTuple->mumuMass->at(0),
						 NTuple->genKstTrkpPx,NTuple->genKstTrkpPy,NTuple->genKstTrkpPz,Utility->kaonMass);
	      massKpi = Utility->computeInvMass(NTuple->genKstTrkmPx,NTuple->genKstTrkmPy,NTuple->genKstTrkmPz,Utility->pionMass,
						NTuple->genKstTrkpPx,NTuple->genKstTrkpPy,NTuple->genKstTrkpPz,Utility->kaonMass);
	    }
	  else
	    {
	      massPsiK = Utility->computeInvMass(NTuple->genMumPx + NTuple->genMupPx,NTuple->genMumPy + NTuple->genMupPy,NTuple->genMumPz + NTuple->genMupPz,NTuple->mumuMass->at(0),
						 NTuple->genKstTrkmPx,NTuple->genKstTrkmPy,NTuple->genKstTrkmPz,Utility->kaonMass);
	      massKpi = Utility->computeInvMass(NTuple->genKstTrkmPx,NTuple->genKstTrkmPy,NTuple->genKstTrkmPz,Utility->kaonMass,
						NTuple->genKstTrkpPx,NTuple->genKstTrkpPy,NTuple->genKstTrkpPz,Utility->pionMass);
	    }
	}


      if ((bkgSub == false) || (fabs(NTuple->B0MassArb - Utility->B0Mass) < atof(Utility->GetGenericParam("NSigmaB0S").c_str())*signalSigma))
	{
	  // ####################
	  // # Make signal plot #
	  // ####################
	  hKKSig->Fill(KKmass);
	  hKstSig->Fill(massKpi);
	  hDalitzSig->Fill(massKpi*massKpi,massPsiK*massPsiK);
	}
      else if (((NTuple->B0MassArb > Utility->B0Mass + atof(Utility->GetGenericParam("NSigmaB0B").c_str())*signalSigma) &&
		(NTuple->B0MassArb < Utility->B0Mass + (atof(Utility->GetGenericParam("NSigmaB0B").c_str())+atof(Utility->GetGenericParam("NSigmaB0S").c_str()))*signalSigma)) ||
	       ((NTuple->B0MassArb < Utility->B0Mass - atof(Utility->GetGenericParam("NSigmaB0B").c_str())*signalSigma) &&
		(NTuple->B0MassArb > Utility->B0Mass - (atof(Utility->GetGenericParam("NSigmaB0B").c_str())+atof(Utility->GetGenericParam("NSigmaB0S").c_str()))*signalSigma)))
	{
	  // ########################
	  // # Make background plot #
	  // ########################
	  hKKBkg->Fill(KKmass);
	  hKstBkg->Fill(massKpi);
	  hDalitzBkg->Fill(massKpi*massKpi,massPsiK*massPsiK);
	}
    }


  if (bkgSub == true) hKKSig->Add(hKKBkg, -1.0);

  c0->cd();
  hKKSig->Draw();
  c0->Modified();
  c0->Update();


  if (bkgSub == true) hKstSig->Add(hKstBkg, -1.0);

  c1->cd();
  hKstSig->Draw("e1p");
  c1->Modified();
  c1->Update();


  if (bkgSub == true) hDalitzSig->Add(hDalitzBkg, -1.0);

  c2->cd();
  hDalitzSig->Draw("cont4z");
  c2->Modified();
  c2->Update();
}


void PlotMuHadMass (string fileName)
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  SetStyle();

  
  double minX = 0.0;
  double maxX = 5.0;
  double MuHadMass;
  int nEntries;
  unsigned int nBins = 100;


  TCanvas* c0 = new TCanvas("c0","c0",10,10,700,500);
  TCanvas* c1 = new TCanvas("c1","c1",20,20,700,500);


  TFile* file0 = TFile::Open(fileName.c_str(),"READ");
  TTree* theTree = (TTree*)file0->Get("B0KstMuMu/B0KstMuMuNTuple");

  B0KstMuMuSingleCandTreeContent* NTuple = new B0KstMuMuSingleCandTreeContent();
  NTuple->Init();
  NTuple->ClearNTuple();
  NTuple->SetBranchAddresses(theTree);
 
  nEntries = theTree->GetEntries();
  cout << "\n[MakePlots::PlotMuHadMass]\t@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;

  TH1D* histoK = new TH1D("histoK","histoK",nBins,minX,maxX);
  histoK->SetXTitle("m(K #mu) (GeV)");
  histoK->SetYTitle("Entries / (0.056 GeV)");
  histoK->SetFillColor(kAzure+6);

  TH1D* histoPi = new TH1D("histoPi","histoPi",nBins,minX,maxX);
  histoPi->SetXTitle("m(#pi #mu) (GeV)");
  histoPi->SetYTitle("Entries / (0.056 GeV)");
  histoPi->SetFillColor(kAzure+6);

  for (int entry = 0; entry < nEntries; entry++)
    {
     theTree->GetEntry(entry);

     if ((NTuple->B0MassArb > Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str())) && (NTuple->B0MassArb < Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str())))
       {
	 if (NTuple->B0notB0bar == true)
	   {
	     MuHadMass = Utility->computeInvMass(NTuple->mumPx->at(0),NTuple->mumPy->at(0),NTuple->mumPz->at(0),Utility->muonMass,
						 NTuple->kstTrkpPx->at(0),NTuple->kstTrkpPy->at(0),NTuple->kstTrkpPz->at(0),Utility->kaonMass);
	 
	     histoK->Fill(MuHadMass);
	 
	     MuHadMass = Utility->computeInvMass(NTuple->mupPx->at(0),NTuple->mupPy->at(0),NTuple->mupPz->at(0),Utility->muonMass,
						 NTuple->kstTrkmPx->at(0),NTuple->kstTrkmPy->at(0),NTuple->kstTrkmPz->at(0),Utility->pionMass);
	 
	     histoPi->Fill(MuHadMass);
	   }
	 else
	   {
	     MuHadMass = Utility->computeInvMass(NTuple->mumPx->at(0),NTuple->mumPy->at(0),NTuple->mumPz->at(0),Utility->muonMass,
						 NTuple->kstTrkpPx->at(0),NTuple->kstTrkpPy->at(0),NTuple->kstTrkpPz->at(0),Utility->pionMass);
	 
	     histoPi->Fill(MuHadMass);
	 
	     MuHadMass = Utility->computeInvMass(NTuple->mupPx->at(0),NTuple->mupPy->at(0),NTuple->mupPz->at(0),Utility->muonMass,
						 NTuple->kstTrkmPx->at(0),NTuple->kstTrkmPy->at(0),NTuple->kstTrkmPz->at(0),Utility->kaonMass);
	 
	     histoK->Fill(MuHadMass);
	   }
       }
    }


  c0->cd();
  histoK->Draw();
  c0->Modified();
  c0->Update();

  c1->cd();
  histoPi->Draw();
  c1->Modified();
  c1->Update();
}


void MakeFitResPlots (string fileName, string plotType, int specBin, string varName, double lowBound, double highBound)
// ##################################
// # varName is one of variables in #
// # the ExtractYield ntuple output #
// ##################################
// # plotType = "Fl"                #
// # plotType = "Afb"               #
// # plotType = "P1"                #
// # plotType = "P2"                #
// # plotType = "BF"                #
// ##################################
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  SetStyle();
  gStyle->SetPadRightMargin(0.08);


  stringstream myString;
  double val;

  vector<vector<double>*> vecNLL;
  Utility->ReadNLLval(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&vecNLL);
  val = Utility->GetNLLval(&vecNLL,plotType,specBin);

  TCanvas* c0 = new TCanvas("c0","c0",10,10,700,500);

  TFile* fileID = TFile::Open(fileName.c_str(),"READ");
  TTree* theTree = (TTree*)fileID->Get("FitResults");
  TH1D* histo = new TH1D("histo","histo",100,lowBound,highBound);
  histo->SetXTitle(varName.c_str());
  histo->SetYTitle("Enstries [#]");
  histo->SetFillColor(kGreen-7);


  myString.clear(); myString.str("");
  myString << varName << " >>histo";
  theTree->Draw(myString.str().c_str(),"","goff");

  
  c0->cd();
  histo->Draw();
  DrawExclusion(val,histo->GetBinLowEdge(histo->GetNbinsX())+histo->GetBinWidth(1),0.0,histo->GetMaximum()*1.1,"NLL",3001,kGray);

  c0->Modified();
  c0->Update();

  TPaveStats* st = (TPaveStats*)histo->FindObject("stats");
  st->SetX1NDC(0.15);
  st->SetX2NDC(0.35);
  st->SetY1NDC(0.85);
  st->SetY2NDC(0.55);

  c0->Modified();
  c0->Update();
}


void ScatterPlotB0MuMu (string fileName, unsigned int cutType)
// ##############################################################  
// # cutType = 0 --> no psi rejection                           #
// # cutType = 1 --> psi rejection                              #
// # cutType = 2 --> "1" + B0&psi cut for J/psi                 #
// # cutType = 3 --> "1" + B0&psi cut for psi(2S)               #
// # cutType = 4 --> "1" + B0&psi cut between J/psi and psi(2S) #
// # cutType = 5 --> "1" + "2" + "3" + "4"                      #
// ##############################################################  
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  SetStyle();
  gStyle->SetOptStat(0);

  
  stringstream myString;
  double lowB0Mass  = Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str());
  double highB0Mass = Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str());
  vector<double> q2Bins;
  Utility->Readq2Bins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins);

  TFile* fileID = TFile::Open(fileName.c_str(),"READ");
  TTree* B0KstMuMuNTuple = (TTree*)fileID->Get("B0KstMuMu/B0KstMuMuNTuple");

  
  TCanvas* c0 = new TCanvas("c0","c0",1200,800);
  c0->cd();
  TH2D* hs = new TH2D("hs","hs",200,sqrt(q2Bins[0]),sqrt(q2Bins[q2Bins.size()-1]),200,lowB0Mass,highB0Mass);
  hs->SetXTitle("m(#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#kern[-0.1]{#mu}#kern[-1.3]{#lower[0.6]{^{#font[122]{\55}}}}) (GeV)");
  hs->SetYTitle("m(K#pi#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#kern[-0.1]{#mu}#kern[-1.3]{#lower[0.6]{^{#font[122]{\55}}}}) (GeV)");
  hs->SetZTitle("Entries [#]");

  
  myString.clear(); myString.str("");
  if (cutType == 0) myString << "(B0MassArb > " << lowB0Mass << ") && (B0MassArb < " << highB0Mass << ")";
  else if (cutType == 1)
    {
      myString << "(B0MassArb > " << lowB0Mass << ") && (B0MassArb < " << highB0Mass << ")";
      myString << " && ";
      myString << "(abs(mumuMass-" << Utility->JPsiMass << ") > " << atof(Utility->GetGenericParam("NSigmaPsi").c_str()) << "*mumuMassE) && (abs(mumuMass-" << Utility->PsiPMass << ") > " << atof(Utility->GetGenericParam("NSigmaPsi").c_str()) << "*mumuMassE)";
    }
  else if (cutType == 2)
    {
      myString << "(B0MassArb > " << lowB0Mass << ") && (B0MassArb < " << highB0Mass << ")";
      myString << " && ";
      myString << "(abs(mumuMass-" << Utility->JPsiMass << ") > " << atof(Utility->GetGenericParam("NSigmaPsi").c_str()) << "*mumuMassE) && (abs(mumuMass-" << Utility->PsiPMass << ") > " << atof(Utility->GetGenericParam("NSigmaPsi").c_str()) << "*mumuMassE)";
      myString << " && ";
      myString << "((mumuMass < " << Utility->JPsiMass << ") && !((abs((B0MassArb - " << Utility->B0Mass << ") - (mumuMass - " << Utility->JPsiMass << ")) < " << atof(Utility->GetGenericParam("B&psiMassJpsiLo").c_str()) << ") || (abs((B0MassArb - " << Utility->B0Mass << ") - (mumuMass - " << Utility->PsiPMass << ")) < " << atof(Utility->GetGenericParam("B&psiMassPsiPLo").c_str()) << ")))";
    }
  else if (cutType == 3)
    {
      myString << "(B0MassArb > " << lowB0Mass << ") && (B0MassArb < " << highB0Mass << ")";
      myString << " && ";
      myString << "(abs(mumuMass-" << Utility->JPsiMass << ") > " << atof(Utility->GetGenericParam("NSigmaPsi").c_str()) << "*mumuMassE) && (abs(mumuMass-" << Utility->PsiPMass << ") > " << atof(Utility->GetGenericParam("NSigmaPsi").c_str()) << "*mumuMassE)";
      myString << " && ";
      myString << "((mumuMass > " << Utility->PsiPMass << ") && !((abs((B0MassArb - " << Utility->B0Mass << ") - (mumuMass - " << Utility->PsiPMass << ")) < " << atof(Utility->GetGenericParam("B&psiMassPsiPHi").c_str()) << ") || (abs((B0MassArb - " << Utility->B0Mass << ") - (mumuMass - " << Utility->JPsiMass << ")) < " << atof(Utility->GetGenericParam("B&psiMassJpsiHi").c_str()) << ")))";
    }
  else if (cutType == 4)
    {
      myString << "(B0MassArb > " << lowB0Mass << ") && (B0MassArb < " << highB0Mass << ")";
      myString << " && ";
      myString << "(abs(mumuMass-" << Utility->JPsiMass << ") > " << atof(Utility->GetGenericParam("NSigmaPsi").c_str()) << "*mumuMassE) && (abs(mumuMass-" << Utility->PsiPMass << ") > " << atof(Utility->GetGenericParam("NSigmaPsi").c_str()) << "*mumuMassE)";
      myString << " && ";
      myString << "((mumuMass > " << Utility->JPsiMass << ") && (mumuMass < " << Utility->PsiPMass << ") && !((abs((B0MassArb - " << Utility->B0Mass << ") - (mumuMass - " << Utility->JPsiMass << ")) < " << atof(Utility->GetGenericParam("B&psiMassJpsiHi").c_str()) << ") || (abs((B0MassArb - " << Utility->B0Mass << ") - (mumuMass - " << Utility->PsiPMass << ")) < " << atof(Utility->GetGenericParam("B&psiMassPsiPLo").c_str()) << ")))";
    }
  else if (cutType == 5)
    {
      myString << "(B0MassArb > " << lowB0Mass << ") && (B0MassArb < " << highB0Mass << ")";
      myString << " && ";
      myString << "(abs(mumuMass-" << Utility->JPsiMass << ") > " << atof(Utility->GetGenericParam("NSigmaPsi").c_str()) << "*mumuMassE) && (abs(mumuMass-" << Utility->PsiPMass << ") > " << atof(Utility->GetGenericParam("NSigmaPsi").c_str()) << "*mumuMassE)";
      myString << " && (";
      myString << "((mumuMass < " << Utility->JPsiMass << ") && !((abs((B0MassArb - " << Utility->B0Mass << ") - (mumuMass - " << Utility->JPsiMass << ")) < " << atof(Utility->GetGenericParam("B&psiMassJpsiLo").c_str()) << ") || (abs((B0MassArb - " << Utility->B0Mass << ") - (mumuMass - " << Utility->PsiPMass << ")) < " << atof(Utility->GetGenericParam("B&psiMassPsiPLo").c_str()) << ")))";
      myString << " || ";
      myString << "((mumuMass > " << Utility->PsiPMass << ") && !((abs((B0MassArb - " << Utility->B0Mass << ") - (mumuMass - " << Utility->PsiPMass << ")) < " << atof(Utility->GetGenericParam("B&psiMassPsiPHi").c_str()) << ") || (abs((B0MassArb - " << Utility->B0Mass << ") - (mumuMass - " << Utility->JPsiMass << ")) < " << atof(Utility->GetGenericParam("B&psiMassJpsiHi").c_str()) << ")))";
      myString << " || ";
      myString << "((mumuMass > " << Utility->JPsiMass << ") && (mumuMass < " << Utility->PsiPMass << ") && !((abs((B0MassArb - " << Utility->B0Mass << ") - (mumuMass - " << Utility->JPsiMass << ")) < " << atof(Utility->GetGenericParam("B&psiMassJpsiHi").c_str()) << ") || (abs((B0MassArb - " << Utility->B0Mass << ") - (mumuMass - " << Utility->PsiPMass << ")) < " << atof(Utility->GetGenericParam("B&psiMassPsiPLo").c_str()) << ")))";
      myString << ")";
    }
  cout << "Query string : " << myString.str().c_str() << endl;
  B0KstMuMuNTuple->Draw("B0MassArb:mumuMass>>hs",myString.str().c_str(),"goff");

  hs->Draw();

  
  // ######################
  // # Draw signal region #
  // ######################
  TCutG* rejectJPsi = new TCutG("rejectJPsi",5);
  rejectJPsi->SetVarX("");
  rejectJPsi->SetVarY("");
  rejectJPsi->SetPoint(0,sqrt(q2Bins[Utility->GetJPsiBin(&q2Bins)]),lowB0Mass);
  rejectJPsi->SetPoint(1,sqrt(q2Bins[Utility->GetJPsiBin(&q2Bins)+1]),lowB0Mass);
  rejectJPsi->SetPoint(2,sqrt(q2Bins[Utility->GetJPsiBin(&q2Bins)+1]),highB0Mass);
  rejectJPsi->SetPoint(3,sqrt(q2Bins[Utility->GetJPsiBin(&q2Bins)]),highB0Mass);
  rejectJPsi->SetPoint(4,sqrt(q2Bins[Utility->GetJPsiBin(&q2Bins)]),lowB0Mass);
  rejectJPsi->SetFillColor(kRed);
  rejectJPsi->SetFillStyle(3002);
  rejectJPsi->Draw("F");

  TCutG* rejectPsiP = new TCutG("rejectPsiP",5);
  rejectPsiP->SetVarX("");
  rejectPsiP->SetVarY("");
  rejectPsiP->SetPoint(0,sqrt(q2Bins[Utility->GetPsiPBin(&q2Bins)]),lowB0Mass);
  rejectPsiP->SetPoint(1,sqrt(q2Bins[Utility->GetPsiPBin(&q2Bins)+1]),lowB0Mass);
  rejectPsiP->SetPoint(2,sqrt(q2Bins[Utility->GetPsiPBin(&q2Bins)+1]),highB0Mass);
  rejectPsiP->SetPoint(3,sqrt(q2Bins[Utility->GetPsiPBin(&q2Bins)]),highB0Mass);
  rejectPsiP->SetPoint(4,sqrt(q2Bins[Utility->GetPsiPBin(&q2Bins)]),lowB0Mass);
  rejectPsiP->SetFillColor(kRed);
  rejectPsiP->SetFillStyle(3002);
  rejectPsiP->Draw("F");


  c0->Modified();
  c0->Update();
}


void MakePvaluePlot (string fileName, int specBin)
// #################################################
// # fileName must of be of the form: *_q2bin.root #
// #################################################
// # If specBin == -1 then loop over all q^2 bins  #
// #################################################
{
  // ##########################
  // # Set histo layout style #
  // ##########################
  SetStyle();
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  

  stringstream myString;
  double val;
  TFile* fileID;

  vector<double> q2Bins;
  Utility->Readq2Bins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins);
  double* q2Bins_ = Utility->MakeBinning(&q2Bins);

  vector<vector<double>*> vecNLL;
  Utility->ReadNLLval(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&vecNLL);

  TCanvas* c0 = new TCanvas("c0","c0",10,10,700,500);
  c0->cd();
  TH1D* pval = new TH1D("pval","pval",q2Bins.size()-1,q2Bins_);
  pval->SetMarkerStyle(20);
  pval->SetXTitle("q#lower[0.4]{^{2}} (GeV#lower[0.4]{^{2}})");
  pval->SetYTitle("p-value");
  pval->SetMinimum(0.0);
  pval->SetMaximum(1.05);


  cout << "\n[MakePlots::MakePvaluePlot]\t@@@ p-value computaion of NLL distribution from paramter file: " << PARAMETERFILEIN << " @@@" << endl;

  for (int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? static_cast<int>(q2Bins.size()-1) : specBin+1); i++)
    {
      if ((i == Utility->GetJPsiBin(&q2Bins)) || (i == Utility->GetPsiPBin(&q2Bins))) continue;

      val = Utility->GetNLLval(&vecNLL,"Afb",i);

      myString.clear(); myString.str("");
      myString << i << ".root";
      fileName.replace(fileName.find(".root")-1,6,myString.str());
      cout << "\nReading NLL distribution from file: " << fileName.c_str() << endl;
      fileID = new TFile(fileName.c_str(),"READ");

      if (fileID->IsZombie() == false)
	{
	  TTree* theTree = (TTree*)fileID->Get("FitResults");
	  int nEvents = theTree->GetEntries();
	  myString.clear(); myString.str("");
	  myString << "nll > " << val;
	  int integral = theTree->Draw("nll",myString.str().c_str(),"goff");

	  pval->SetBinContent(i+1,static_cast<double>(integral) / static_cast<double>(nEvents));
	  pval->SetBinError(i+1,pval->GetBinContent(i+1)*1e-3);
	  cout << "p-value for q^2 bin #" << i << " --> " << pval->GetBinContent(i+1) << endl;

	  fileID->Close();
	  delete fileID;
	}
    }
  

  c0->cd();
  pval->Draw("e1p");
  DrawExclusion(pval->GetBinLowEdge(1),pval->GetBinLowEdge(pval->GetNbinsX())+pval->GetBinWidth(pval->GetNbinsX()),0.0,0.05,"p-val",3001,kRed-9);

  c0->Modified();
  c0->Update();
}


int main (int argc, char** argv)
{
  if (argc >= 2)
    {
      string option = argv[1];
      unsigned int intVal = 0;
      double realVal1 = 0;
      double realVal2 = 0;
      string fileName;
      string tmpStr1;
      string tmpStr2;
      

      if ((option == "GenMultyRun") && (argc == 4))
	{
	  fileName = argv[2];
	  intVal   = atoi(argv[3]);
	}
      else if (((option == "Phy") || (option == "DataMC")) && (argc == 3)) intVal = atoi(argv[2]);
      else if ((option == "Pval") && (argc == 4))
	{
	  fileName = argv[2];
	  intVal   = atoi(argv[3]);
	}
      else if ((option == "FitRes") && (argc == 8))
	{
	  fileName = argv[2];
	  tmpStr1  = argv[3];
	  intVal   = atoi(argv[4]);
	  tmpStr2  = argv[5];
	  realVal1 = atof(argv[6]);
	  realVal2 = atof(argv[7]);
	}
      else if (((option == "MuMuMass") || (option == "KstMass")) && (argc == 4))
	{
	  fileName = argv[2];
	  intVal   = atoi(argv[3]);
	}
      else if ((option == "KKMass") && (argc == 5))
	{
	  fileName = argv[2];
	  intVal   = atoi(argv[3]);
	  tmpStr1  = argv[4];
	}
      else if ((option == "MuHadMass") && (argc == 3)) fileName = argv[2];
      else if ((option == "ScatB0MuMu") && (argc == 4))
	{
	  fileName = argv[2];
	  intVal   = atoi(argv[3]);
	}
      else if (option != "PhyRegion")
	{
	  cout << "./MakePlots [Phy GenMultyRun DataMC PhyRegion Pval FitRes MuMuMass KKMass KstMass MuHadMass ScatB0MuMu]" << endl;
	  cout << "            [Phy: 0-2||10-12]" << endl;
	  cout << "            [GenMultyRun: fileName q^2_bin_index]" << endl;
	  cout << "            [DataMC: 0-27]" << endl;
	  cout << "            [Pval: toyFileName q^2_bin_index]" << endl;
	  cout << "            [FitRes: toyFileName plotType q^2_bin_index varName lowBound highBound]" << endl;
	  cout << "            [MuMuMass OR KstMass: dataFileName bkgSub]" << endl;
	  cout << "            [KKMass: dataFileName bkgSub RECOorGEN]" << endl;
	  cout << "            [MuHadMass: dataFileName]" << endl;
	  cout << "            [ScatB0MuMu: dataFileName option]" << endl;

	  return EXIT_FAILURE;
	}


      cout << "\n[MakePlots::main]\t@@@ Settings @@@" << endl;
      cout << "PARAMETERFILEIN : "     << PARAMETERFILEIN << endl;
      cout << "ParameterFILE_MCGEN: "  << ParameterFILE_MCGEN << endl;
      cout << "ParameterFILE_MCRECO: " << ParameterFILE_MCRECO << endl;

      cout << "\nSMBINFL_nL: " << SMBINFL_nL << endl;
      cout << "SMBINAFB_nL: "  << SMBINAFB_nL << endl;
      cout << "SMBINBF_nL: "   << SMBINBF_nL << endl;

      cout << "\nSMBINFL_yL: " << SMBINFL_yL << endl;
      cout << "SMBINAFB_yL: "  << SMBINAFB_yL << endl;
      cout << "SMBINBF_yL: "   << SMBINBF_yL << endl;

      cout << "\nSingleCand_MCkstJPsi: " << SingleCand_MCkstJPsi << endl;
      cout << "SingleCand_MCkstPsi2S: "  << SingleCand_MCkstPsi2S << endl;
      cout << "SingleCand_Data: "        << SingleCand_Data << endl;

      cout << "\nYvalueOutsideLimits: " << YvalueOutsideLimits << endl;
      cout << "FORPAPER: "              << FORPAPER << endl;

      cout << "\noption: " << option << endl;
      cout << "intVal: "   << intVal << endl;
      cout << "realVal1: " << realVal1 << endl;
      cout << "realVal2: " << realVal2 << endl;
      cout << "fileName: " << fileName << endl;
      cout << "tmpStr1: "  << tmpStr1 << endl; 
      cout << "tmpStr2: "  << tmpStr2 << endl; 

      if (option == "GenMultyRun") gROOT->SetBatch(true);
      TApplication theApp ("Applications", &argc, argv);

      Utility = new Utils();
      Utility->ReadGenericParam(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str());

      if      (option == "Phy")         MakePhysicsPlots(intVal);
      else if (option == "GenMultyRun") GenNTupleFromMultyRun(fileName,intVal);
      else if (option == "DataMC")      MakeComparisonDataMC(intVal);
      else if (option == "PhyRegion")   CheckPhysicsRegion();
      else if (option == "Pval")        MakePvaluePlot(fileName,intVal);
      else if (option == "FitRes")      MakeFitResPlots(fileName,tmpStr1,intVal,tmpStr2,realVal1,realVal2);
      else if (option == "MuMuMass")    PlotMuMu(fileName,intVal);
      else if (option == "KKMass")      PlotKK(fileName,intVal,tmpStr1);
      else if (option == "KstMass")     PlotKst(fileName,intVal,true);
      else if (option == "MuHadMass")   PlotMuHadMass(fileName);
      else if (option == "ScatB0MuMu")  ScatterPlotB0MuMu(fileName,intVal);
  
      delete Utility;
      if (option != "GenMultyRun") theApp.Run (); // Eventloop on air
      return EXIT_SUCCESS;
    }
  else
    {
      cout << "./MakePlots [Phy GenMultyRun DataMC PhyRegion Pval FitRes MuMuMass KKMass KstMass MuHadMass ScatB0MuMu]" << endl;
      cout << "            [Phy: 0-2||10-12]" << endl;
      cout << "            [GenMultyRun: fileName q^2_bin_index]" << endl;
      cout << "            [DataMC: 0-27]" << endl;
      cout << "            [Pval: toyFileName q^2_bin_index]" << endl;
      cout << "            [FitRes: toyFileName plotType q^2_bin_index varName lowBound highBound]" << endl;
      cout << "            [MuMuMass OR KstMass: dataFileName bkgSub]" << endl;
      cout << "            [KKMass: dataFileName bkgSub RECOorGEN]" << endl;
      cout << "            [MuHadMass: dataFileName]" << endl;
      cout << "            [ScatB0MuMu: dataFileName cutType]" << endl;

      cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "For [MuMuMass KKMass KstMass]:" << endl;
      cout << "bkgSub = 0 (do not subtract background)" << endl;
      cout << "bkgSub = 1 (subtract background)" << endl;

      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "For [KKMass]:" << endl;
      cout << "RECOorGEN = RECO (use RECO information)" << endl;
      cout << "RECOorGEN = GEN (use GEN information)" << endl;

      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "For [Phy]:" << endl;
      cout << "0 = Fl GEN vs RECO-MC" << endl;
      cout << "1 = Afb GEN vs RECO-MC" << endl;
      cout << "2 = BF GEN vs RECO-MC" << endl;
      cout << "10 = Fl Data vs Theory" << endl;
      cout << "11 = Afb Data vs Theory" << endl;
      cout << "12 = BF Data vs Theory" << endl;
      cout << "13 = Fs Data" << endl;
      cout << "14 = As Data" << endl;

      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "For [FitRes]:" << endl;
      cout << "Fl" << endl;
      cout << "Afb" << endl;
      cout << "P1" << endl;
      cout << "P2" << endl;
      cout << "BF" << endl;

      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "For [ScatB0MuMu]:" << endl;
      cout << "cutType = 0 --> no psi rejection" << endl;
      cout << "cutType = 1 --> psi rejection" << endl;
      cout << "cutType = 2 --> 1 + B0&psi cut for J/psi" << endl;
      cout << "cutType = 3 --> 1 + B0&psi cut for psi(2S)";
      cout << "cutType = 4 --> 1 + B0&psi cut between J/psi and psi(2S)" << endl;
      cout << "cutType = 5 --> 1 + 2 + 3 + 4" << endl;
      
      cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
      cout << "For [DataMC]:" << endl;
      cout << "0 = B0 pT" << endl;
      cout << "1 = B0 eta" << endl;
      cout << "2 = B0 vtx CL" << endl;
      cout << "3 = B0 L/sigma" << endl;
      cout << "4 = B0 cos(alpha)" << endl;

      cout << "5 = mu+ pT" << endl;
      cout << "6 = mu- pT" << endl;
      cout << "7 = mu+ eta" << endl;
      cout << "8 = mu- eta" << endl;
      cout << "9 = mu+ phi, eta range" << endl;
      cout << "10 = mu+ phi, eta range" << endl;
      cout << "11 = mu+ phi, eta range" << endl;
      cout << "12 = mu+ phi, eta range" << endl;
      cout << "13 = mu- phi, eta range" << endl;
      cout << "14 = mu- phi, eta range" << endl;
      cout << "15 = mu- phi, eta range" << endl;
      cout << "16 = mu- phi, eta range" << endl;

      cout << "17 = K*0 trk+ pT" << endl;
      cout << "18 = K*0 trk- pT" << endl;
      cout << "19 = K*0 trk+ eta" << endl;
      cout << "20 = K*0 trk- eta" << endl;
      cout << "21 = K*0 trk+ phi, eta range" << endl;
      cout << "22 = K*0 trk+ phi, eta range" << endl;
      cout << "23 = K*0 trk+ phi, eta range" << endl;
      cout << "24 = K*0 trk+ phi, eta range" << endl;
      cout << "25 = K*0 trk- phi, eta range" << endl;
      cout << "26 = K*0 trk- phi, eta range" << endl;
      cout << "27 = K*0 trk- phi, eta range" << endl;
      cout << "28 = K*0 trk- phi, eta range" << endl;

      cout << "29 = K*0 trk+ DCA/sigma" << endl;
      cout << "30 = K*0 trk- DCA/sigma" << endl;

      cout << "31 = cos(theta_K)" << endl;
      cout << "32 = cos(theta_l)" << endl;
      cout << "33 = phi" << endl;

      return EXIT_FAILURE;
    }
  
  return EXIT_SUCCESS;
}
