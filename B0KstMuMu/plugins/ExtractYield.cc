// ################################################################################
// # Program to perform the full angular analysis of the decay B0 --> K*0 mu+ mu- #
// ################################################################################

#include <TROOT.h>
#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH3.h>
#include <TF2.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TLatex.h>
#include <TCutG.h>
#include <Math/Functor.h>

#include <TMinuit.h>
#include <RooRealVar.h>
#include <RooGaussian.h>
#include <RooAbsPdf.h>
#include <RooHistPdf.h>
#include <RooAddPdf.h>
#include <RooDataSet.h>
#include <RooGenericPdf.h>
#include <RooPlot.h>
#include <RooArgSet.h>
#include <RooFitResult.h>
#include <RooPolynomial.h>
#include <RooMCStudy.h>
#include <RooMinuit.h>
#include <RooWorkspace.h>
#include <RooConstVar.h>
#include <RooRandom.h>
#include <RooDataHist.h>
#include <RooFunctorBinding.h>
#include <RooStats/RooStatsUtils.h>
#include <RooMinimizer.h>

#include <ctime>
#include <iostream>
#include <utility>
#include <sstream>

#include "ReadParameters.h"
#include "Utils.h"
#include "B0KstMuMuSingleCandTreeContent.h"

using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::vector;
using std::ios_base;
using std::pair;
using std::make_pair;
using namespace RooFit;


// ####################
// # Global constants #
// ####################
#define NBINS        20
#define MULTYIELD     1. // Multiplication factor to the number of entry in toy-MC
#define NCOEFFPOLYBKG 5  // Maximum number of coefficients (= degree) of the polynomial describing the combinatorial bkg

#define nJPSIS 230000.0
#define nJPSIB   2500.0
#define nPSIPS  15000.0
#define nPSIPB   1500.0

#define _USE_MATH_DEFINES

// ##########################################
// # Internal flags to control the workflow #
// ##########################################
#define MAKEmumuPLOTS false
#define SETBATCH      true
#define PROFILENLL    true //[true= plot the profile likelihood]
#define PLOT          false //[true= plot the results]
#define SAVEPOLY      false // ["true" = save bkg polynomial coefficients in new parameter file; "false" = save original values]
#define SAVEPLOT      false   //2015-08-20
#define RESETsigANG   false // Reset signal angular parameters before starting the fit
#define RESETcomANG   false // Reset combinatorial bkg angular parameters before starting the fit
#define FULLTOYS      false // Run generation-and-fit toys
#define FUNCERRBAND   false // Show the p.d.f. error band
#define MINIMIZER     "Minuit2" // Minimizer type for 3D MODEL actual fit ["Minuit"; "Minuit2"]
#define GENPARAMS     "All" // Option to generate parameters for parameter file: "All" "misTagFrac" "FlP5pFsAs" "combBkgAng"

// ##################
// # External files #
// ##################
#define PARAMETERFILEIN  "/python/ParameterFile.txt"
#define PARAMETERFILEOUT "ParameterFileOut.txt"


// ############################################
// # Global variables from configuration file #
// ############################################
double PsiYieldGoodTag, PsiYieldGoodTagErr;
double PsiYieldMisTag,  PsiYieldMisTagErr;
double LUMI;

string CTRLfitWRKflow;
string ParameterFILE;

vector<double> q2Bins;
vector<double> cosThetaKBins;
vector<double> cosThetaLBins;
vector<double> phiBins;

vector<vector<string>*>             fitParam;    // Vector containing the pointers to the vectors containing the starting values for the fit
vector<vector<unsigned int>*>       configParam; // Vector containing the pointers to the vectors containing the configuration parameters for the fit

// ####################
// # Global variables #
// ####################
TTree* theTree;
Utils* Utility;
B0KstMuMuSingleCandTreeContent* NTuple;

ofstream fileFitResults;
ofstream fileFitSystematics;

double* q2BinsHisto;

bool useToyDataset;

bool scanInitVal;
int scanIndx;
double scanIndxMax = 800.;

// double maxAs5[9] = {1,1,1,,,,,,};
// double minAs5[9] = {-1,-1,,,,,,,};
double maxP1 [9] = { 0.85, 0.35, 1.00, 0.00,0,-0.25,0, 0.10,-0.17};
double minP1 [9] = {-0.65,-1.00,-0.15,-0.90,0,-0.80,0,-0.75,-0.85};
double maxP5p[9] = { 0.65, 0.50,-0.60,-0.30,0,-0.45,0,-0.35,-0.32};
double minP5p[9] = {-0.45,-1.10,-1.30,-0.95,0,-0.85,0,-0.95,-0.82};

double scanStep1[9] = {0.0159847,0.0110658,0.0140337,0.00486634,0,0.00471914,0,0.00420572,0.00854125};
double scanStep2[9] = {0.42,0.67,0.33,0.24,0,0.17,0,0.07,0.22};
double scanStep3[9] = {0.34,0.38,0.34,0.08,0,0.18,0,0.12,0.19};
double scanStep4[9] = {0.11,0.08,0.10,0.07,0,0.04,0,0.07,0.07};

// double scanStep1[9] = {0.0159847,0.0110658,0.0140337,0.00486634,0,0.00471914,0,0.00420572,0.00854125};
// double scanStep2[9] = {0.0730085,0.321776,0.0728803,0.0316799,0,0.0130877,0,0.1,0.0597653};
// // double scanStep2[9] = {0.0730085,0.121776,0.0728803,0.0316799,0,0.0130877,0,0.0309526,0.0597653};
// double scanStep3[9] = {0.0874036,0.105432,0.0627559,0.0272426,0,0.00773635,0,0.0328062,0.0461656};
// double scanStep4[9] = {0.1,0.1,0.1,0.1,0,0.1,0,0.1,0.1};

// double scanStep1[9] = {0.0163855,0.0106618,0.0148608,0.00756547,0,0.00846779,0,0.00867193,0.0131687};
// double scanStep2[9] = {0.0739349,0.120751,0.0766692,0.0320697,0,0.0136168,0,0.0521343,0.0497501};
// double scanStep3[9] = {0.089052,0.113716,0.0625183,0.0328239,0,0.00776573,0,0.0535013,0.0442178};

// ####################################
// # Useful variables from the NTuple #
// ####################################
RooDataSet* SingleCandNTuple_RejectPsi;
RooDataSet* SingleCandNTuple;
RooRealVar* B0MassArb;
RooRealVar* mumuMass;
RooRealVar* mumuMassE;
RooRealVar* ctK;                //CosThetaKArb
RooRealVar* ctL;                //CosThetaMuArb
RooRealVar* phi;                //PhiKstMuMuPlaneArb
RooRealVar* truthMatchSignal;
RooRealVar* rightFlavorTag;


// #################################
// # Variables and pdf for the fit #
// #################################

// ##################
// # Signal B0 mass #
// ##################
RooRealVar* meanS;

RooRealVar* sigmaS1;
RooAbsPdf*  MassS1;

RooRealVar* sigmaS2;
RooAbsPdf*  MassS2;

RooRealVar* fracMassS;
RooAbsPdf*  MassSignal;

// #################
// # Signal angles #
// #################
RooRealVar* FlS;
RooRealVar* P5pS;
RooRealVar* FsS;
RooRealVar* AsS;
RooRealVar* P1S;              
RooRealVar* As5S;             
RooAbsPdf*  AngleS;


// ####################
// # Total signal pdf #
// ####################
RooAbsPdf* Signal;
RooAbsPdf* SignalT;

// ####################################
// # Combinatorial background B0 mass #
// ####################################
RooRealVar* var1;
RooRealVar* var2;
RooAbsPdf*  BkgMassExp1;
RooAbsPdf*  BkgMassExp2;

RooRealVar* fracMassBExp;
RooAbsPdf*  BkgMassComb;

// #########################
// # Mistag signal B0 mass #
// #########################
RooRealVar* sigmaMisTag1;
RooAbsPdf*  MassMisTag1;
RooRealVar* sigmaMisTag2;
RooAbsPdf*  MassMisTag2;
RooRealVar* fracMisTag;
RooAbsPdf*  MassMisTag;
RooAbsPdf*  AngleMisTag;

// ##############################
// # Peaking background B0 mass #
// ##############################
RooRealVar* meanR1;
RooRealVar* sigmaR1;
RooRealVar* meanR2;
RooRealVar* sigmaR2;
RooAbsPdf*  BkgMassRPeak1;
RooAbsPdf*  BkgMassRPeak2;

RooRealVar* fracMassBRPeak;
RooAbsPdf*  BkgMassRPeak;

RooRealVar* meanL1;
RooRealVar* sigmaL1;
RooRealVar* meanL2;
RooRealVar* sigmaL2;
RooAbsPdf*  BkgMassLPeak1;
RooAbsPdf*  BkgMassLPeak2;

RooRealVar* fracMassBLPeak;
RooAbsPdf*  BkgMassLPeak;

RooRealVar* fracMassBPeak;
RooAbsPdf*  BkgMassPeak;

// #####################
// # Background angles #
// #####################
RooRealVar* p1Poly[NCOEFFPOLYBKG];
RooAbsPdf*  BkgAngleP1;
RooRealVar* c1Poly[NCOEFFPOLYBKG];
RooAbsPdf*  BkgAngleC1;

RooRealVar* p2Poly[NCOEFFPOLYBKG];
RooAbsPdf*  BkgAngleP2;
RooRealVar* c2Poly[NCOEFFPOLYBKG];
RooAbsPdf*  BkgAngleC2;

RooRealVar* p3Poly[NCOEFFPOLYBKG];
RooAbsPdf*  BkgAngleP3;
RooRealVar* c3Poly[NCOEFFPOLYBKG];
RooAbsPdf*  BkgAngleC3;

RooAbsPdf* BkgAnglesC;
RooAbsPdf* BkgAnglesP;

// ########################
// # Total background pdf #
// ########################
RooAbsPdf* BkgMassAngleComb;
RooAbsPdf* MassAngleMisTag;
RooAbsPdf* BkgMassAnglePeak;

RooRealVar* nSig;
RooRealVar* nMisTagFrac;
RooRealVar* nBkgPeak;
RooRealVar* nBkgComb;

// ##################################
// # Total pdf for B0 --> K*0 mu mu #
// ##################################
RooAbsPdf* TotalPDFRejectPsi;

// ##############################################
// # Vector containing the Gaussian constraints #
// ##############################################
RooArgSet vecConstr;


// #######################
// # Function Definition #
// #######################

// ################################################
// # Structure to make the product of two p.d.f.s #
// ################################################
struct MyProdPdf
{
  public:
    MyProdPdf (RooAbsPdf& pdf1, RooAbsPdf& pdf2) : _pdf1(pdf1), _pdf2(pdf2)
  { 
    const RooArgSet* allvar1 = pdf1.getVariables();
    const RooArgSet* allvar2 = pdf2.getVariables();

    _vars.add(*allvar1);
    _vars.add(*allvar2,false);

    delete allvar1;
    delete allvar2;
  }

    int ndim ()
    {
      return _vars.getSize();
    }

    const RooArgList& vars() const
    {
      return _vars;
    }

    double operator() (const double* v)
    {
      for (int i = 0; i < ndim(); ++i) ((RooRealVar&)_vars[i]).setVal(v[i]);

      return _pdf1.getVal() * _pdf2.getVal();
    }


  private:
    RooAbsPdf& _pdf1; 
    RooAbsPdf& _pdf2; 
    RooArgList _vars;

};

bool CheckGoodFit               (RooFitResult* fitResult, TPaveText* paveText = NULL);
RooRealVar* GetVar              (RooAbsPdf* pdf, string varName);
bool GetValueAndErrors          (RooAbsPdf* pdf, string varName, stringstream* myString, double& val, double& errLo, double& errHi);
void SetValueAndErrors          (RooAbsPdf* pdf, string varName, double multi, stringstream* myString, double* val, double* errLo, double* errHi);
void PrintVariables             (RooArgSet* setVar, string type);
void ClearVars                  (RooArgSet* vecVars);
void CloseAllAndQuit            (TApplication* theApp, TFile* NtplFile);

void SetStyle                   ();
string MakeName                 (const RooDataSet* data, unsigned int ID);
void DrawString                 (double Lumi, RooPlot* myFrame = NULL);
void AddGaussConstraint         (RooArgSet* vecConstr, RooAbsPdf* pdf, string varName);
void BuildMassConstraints       (RooArgSet* vecConstr, RooAbsPdf* pdf, string varName);
void BuildAngularConstraints    (RooArgSet* vecConstr, RooAbsPdf* pdf, string varName,vector<vector<string>*>* fitParam = NULL, unsigned int q2BinIndx = 0);
RooAbsPdf* MakeAngWithEffPDF (unsigned int q2BinIndx, RooRealVar* y, RooRealVar* z,RooRealVar* p, unsigned int FitType, bool useEffPDF, RooArgSet* VarsAng);
unsigned int CopyFitResults     (RooAbsPdf* pdf, unsigned int q2BinIndx, vector<vector<string>*>* fitParam,unsigned int countMisTag = 0,unsigned int countGoodTag =0);

string GeneratePolynomial       (RooRealVar* var, unsigned int nCoef, string sCoef);

void LoadDatasets               (TFile* NtplFile, int specBin, int FitType, int toyIndx);
void MakeDatasets               (B0KstMuMuSingleCandTreeContent* NTuple, unsigned int FitType);

//###############
// 3-D models   #
//###############

void InstantiateGen3AnglesFit     (RooAbsPdf** TotalPDF,
                                   bool useEffPDF,
                                   RooRealVar* y, RooRealVar* z,RooRealVar* p,
                                   string fitName, unsigned int FitType,
                                   vector<vector<unsigned int>*>* configParam,
                                   vector<vector<string>*>* fitParam,
                                   unsigned int q2BinIndx);
void InstantiateReco3AnglesFit (RooAbsPdf** TotalPDF,
                                bool useEffPDF,
                                RooRealVar* y, RooRealVar* z,RooRealVar* p,
                                string fitName, unsigned int FitType,
                                vector<vector<unsigned int>*>* configParam,
                                vector<vector<string>*>* fitParam,
                                unsigned int q2BinIndx);
RooFitResult* Make3AnglesFit (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* y, RooRealVar* z,RooRealVar* p,  unsigned int FitType, RooArgSet* vecConstr, TCanvas* Canv, unsigned int ID);

void Iterative3AnglesFitq2Bins (RooDataSet* dataSet,
                                bool useEffPDF,
                                RooRealVar* y, RooRealVar* z,RooRealVar* p,
                                int specBin,
                                unsigned int FitType,
                                vector<double>* q2Bins,
                                vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam,
                                RooArgSet* vecConstr,
                                unsigned int ID = 0);



//##############################
//    #  4D model ##############
//#############################
void InstantiateMass3AnglesFit (RooAbsPdf** TotalPDF,
                                bool useEffPDF,
                                RooRealVar* x, RooRealVar* y, RooRealVar* z,RooRealVar* p,
                                string fitName, unsigned int FitType,
                                vector<vector<unsigned int>*>* configParam,
                                vector<vector<string>*>* fitParam,
                                unsigned int q2BinIndx);
RooFitResult* MakeMass3AnglesFit (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* x, RooRealVar* y, RooRealVar* z,RooRealVar* p,  unsigned int FitType, RooArgSet* vecConstr, TPaveText* extText, unsigned int ID);

void IterativeMass3AnglesFitq2Bins (RooDataSet* dataSet,
                                    bool useEffPDF,
                                    RooRealVar* x, RooRealVar* y, RooRealVar* z,RooRealVar* p,
                                    int specBin,
                                    unsigned int FitType,
                                    vector<double>* q2Bins,
                                    vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam,
                                    RooArgSet* vecConstr,
                                    unsigned int ID = 0);



// ###########################
// # Function Implementation #
// ###########################
bool CheckGoodFit (RooFitResult* fitResult, TPaveText* paveText)
  // ####################################################
  // # Covariance matrix quality:                       #
  // # -1 : "Unknown, matrix was externally provided"   #
  // #  0 : "Not calculated at all"                     #
  // #  1 : "Approximation only, not accurate"          #
  // #  2 : "Full matrix, but forced positive-definite" #
  // #  3 : "Full, accurate covariance matrix"          #
  // ####################################################
{
  if (fitResult != NULL)
  {
    if ((fitResult->covQual() == 3) && (fitResult->status() == 0))
    {
      if (paveText != NULL) paveText->AddText("Fit status: GOOD");
      return true;
    }
    else
    {
      if (paveText != NULL) paveText->AddText("Fit status: BAD");
      return false;
    }
  }


  return false;
}


RooRealVar* GetVar (RooAbsPdf* pdf, string varName)
{
  return (RooRealVar*)(pdf->getVariables()->find(varName.c_str()));
}


bool GetValueAndErrors (RooAbsPdf* pdf, string varName, stringstream* myString, double& val, double& errLo, double& errHi)
{
  if (GetVar(pdf,varName.c_str()) != NULL)
  {
    val   = GetVar(pdf,varName.c_str())->getVal();
    errLo = GetVar(pdf,varName.c_str())->getErrorLo();
    errHi = GetVar(pdf,varName.c_str())->getErrorHi();

    (*myString) << val << "   " << errHi << "   " << errLo << "   ";
    return true;
  }
  else (*myString) << "0.0   0.0   0.0";


  return false;
}


void SetValueAndErrors (RooAbsPdf* pdf, string varName, double multi, stringstream* myString, double* val, double* errLo, double* errHi)
  // #############################################################################
  // # If the error is an empty string --> setAsymError = -1/+1 and setError = 1 #
  // # If both errLo and errHi are 0.0 --> setAsymError = -1/+1 and setError = 1 #
  // #############################################################################
{
  string tmpStr;


  if (myString->str().empty() == false)
  {
    tmpStr.clear();
    (*myString) >> tmpStr;
    *val = atof(tmpStr.c_str()) * multi;

    tmpStr.clear();
    (*myString) >> tmpStr;
    if (tmpStr.empty() == true) *errLo = -1.0;
    else *errLo = atof(tmpStr.c_str()) * multi;

    tmpStr.clear();
    (*myString) >> tmpStr;
    if (tmpStr.empty() == true) *errHi = 1.0;
    else *errHi = atof(tmpStr.c_str()) * multi;
  }

  if ((pdf != NULL) && (GetVar(pdf,varName) != NULL))
  {
    pdf->getVariables()->setRealValue(varName.c_str(),*val);
    if ((*errLo == 0.0) && (*errHi == 0.0))
    {
      GetVar(pdf,varName)->setAsymError(0.0,0.0);
      GetVar(pdf,varName)->setError(0.0);
    }
    else
    {
      GetVar(pdf,varName)->setAsymError(*errLo,*errHi);
      GetVar(pdf,varName)->setError((*errHi - *errLo) / 2.);
    }
  }
}


void PrintVariables (RooArgSet* setVar, string type)
{
  RooRealVar* tmpVar;
  int nEleSet = setVar->getSize();


  if (type == "vars")
  {
    cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
    cout <<   "@@@ Printing variables @@@" << endl;
    cout <<   "@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

    if (setVar != NULL)
    {
      TIterator* it = setVar->createIterator();
      for (int i = 0; i < nEleSet; i++)
      {
        tmpVar = (RooRealVar*)it->Next();
        cout << "Variable: " << i;
        cout << "\tname: "   << tmpVar->GetName();
        cout << "\tvalue: "  << tmpVar->getVal();
        cout << "\terr: "    << tmpVar->getError();
        cout << "\tErrLo: "  << tmpVar->getErrorLo();
        cout << "\tErrHi: "  << tmpVar->getErrorHi() << endl;
      }
    }
  }
  else if (type == "cons")
  {
    cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
    cout <<   "@@@@@@@@@ Printing constraints @@@@@@@@@" << endl;
    cout <<   "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

    if (setVar != NULL)
    {
      TIterator* it = setVar->createIterator();
      for (int i = 0; i < nEleSet; i++) PrintVariables(((RooAbsPdf*)it->Next())->getVariables(),"vars");
    }
  }
  else
  {
    cout << "[ExtractYield::PrintVariables]\tWrong parameter: " << type << endl;
    exit (EXIT_FAILURE);
  }
}


void ClearVars (RooArgSet* vecVars)
{
  if (vecVars != NULL)
  {
    int nEle = vecVars->getSize();

    TIterator* it = vecVars->createIterator();
    for (int i = 0; i < nEle; i++) delete it->Next();

    vecVars->removeAll();
  }
}

void CloseAllAndQuit (TApplication* theApp, TFile* NtplFile)
{
  fileFitResults.close();

  if (NtplFile                   != NULL) NtplFile->Close("R");

  gROOT->CloseFiles();
  theApp->Terminate(0);
}


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


string MakeName (const RooDataSet* data, unsigned int ID)
{
  stringstream myString;


  myString.clear(); myString.str("");
  myString << data->GetName() << "_" << ID;


  return myString.str();
}


void DrawString (double Lumi, RooPlot* myFrame)
{
  stringstream myString;
  double scaleRespect2CMS = 0.75;


  myString.clear(); myString.str("");
  myString << "CMS";
  TLatex* LumiTex1 = new TLatex(0.18,0.9,myString.str().c_str());
  LumiTex1->SetTextFont(61);
  LumiTex1->SetTextSize(0.05);
  LumiTex1->SetTextColor(kBlack);
  LumiTex1->SetNDC(true);
  if (myFrame == NULL) LumiTex1->DrawLatex(0.18,0.9,myString.str().c_str());
  else
  {
    LumiTex1->Paint();
    myFrame->addObject(LumiTex1);
  }


  myString.clear(); myString.str("");
  myString << "#it{Preliminary}";
  TLatex* LumiTex2 = new TLatex(0.265,0.9,myString.str().c_str());
  LumiTex2->SetTextFont(42);
  LumiTex2->SetTextSize(0.05 * scaleRespect2CMS);
  LumiTex2->SetTextColor(kBlack);
  LumiTex2->SetNDC(true);
  if (myFrame == NULL) LumiTex2->DrawLatex(0.265,0.9,myString.str().c_str());
  else
  {
    LumiTex2->Paint();
    myFrame->addObject(LumiTex2);
  }


  myString.clear(); myString.str("");
  myString << Lumi <<  " fb#lower[0.4]{^{#font[122]{\55}1}} (8 TeV)";
  TLatex* LumiTex3 = new TLatex(0.8,0.9,myString.str().c_str());
  LumiTex3->SetTextFont(42);
  LumiTex3->SetTextSize(0.05 * scaleRespect2CMS);
  LumiTex3->SetTextColor(kBlack);
  LumiTex3->SetNDC(true);
  if (myFrame == NULL) LumiTex3->DrawLatex(0.8,0.9,myString.str().c_str());
  else
  {
    LumiTex3->Paint();
    myFrame->addObject(LumiTex3);
  }
}

void AddGaussConstraint (RooArgSet* vecConstr, RooAbsPdf* pdf, string varName)
{
  stringstream myString;


  myString.clear(); myString.str("");
  myString << varName;

  if (GetVar(pdf,myString.str().c_str())->isConstant() == false)
  {
    RooRealVar* varConstr = GetVar(pdf,myString.str().c_str());
    double mean  = GetVar(pdf,myString.str().c_str())->getVal();
    double sigma = (varConstr->getErrorHi() - varConstr->getErrorLo()) / 2.;

    myString << "_constr";

    RooGaussian* newConstr = new RooGaussian(myString.str().c_str(), myString.str().c_str(), *varConstr, RooConst(mean), RooConst(sigma));
    vecConstr->add(*newConstr);
  }
}

void BuildMassConstraints (RooArgSet* vecConstr, RooAbsPdf* pdf, string varName)
  // ##############################################################
  // // # All:    apply constraint to all physical angular variables #
  // // # sign:   apply constraint to signal                         #
  // // # mistag: apply constraint to mistag                         #
  // // # peak:   apply constraint to peaking bkg                    #
  // // ##############################################################
{
  if (atoi(Utility->GetGenericParam("ApplyConstr").c_str()) == true)
  {
    if ((GetVar(pdf,"sigmaS1")        != NULL) && ((varName == "All") || (varName == "sign"))) AddGaussConstraint(vecConstr, pdf, "sigmaS1");
    if ((GetVar(pdf,"sigmaS2")        != NULL) && ((varName == "All") || (varName == "sign"))) AddGaussConstraint(vecConstr, pdf, "sigmaS2");
    if ((GetVar(pdf,"fracMassS")      != NULL) && ((varName == "All") || (varName == "sign"))) AddGaussConstraint(vecConstr, pdf, "fracMassS");

    if ((GetVar(pdf,"sigmaMisTag1")   != NULL) && ((varName == "All") || (varName == "mistag"))) AddGaussConstraint(vecConstr, pdf, "sigmaMisTag1");
    if ((GetVar(pdf,"sigmaMisTag2")   != NULL) && ((varName == "All") || (varName == "mistag"))) AddGaussConstraint(vecConstr, pdf, "sigmaMisTag2");
    if ((GetVar(pdf,"fracMisTag")     != NULL) && ((varName == "All") || (varName == "mistag"))) AddGaussConstraint(vecConstr, pdf, "fracMisTag");

    if ((GetVar(pdf,"meanR1")         != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, pdf, "meanR1");
    if ((GetVar(pdf,"sigmaR1")        != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, pdf, "sigmaR1");
    if ((GetVar(pdf,"meanR2")         != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, pdf, "meanR2");
    if ((GetVar(pdf,"sigmaR2")        != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, pdf, "sigmaR2");
    if ((GetVar(pdf,"fracMassBRPeak") != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, pdf, "fracMassBRPeak");

    if ((GetVar(pdf,"meanL1")         != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, pdf, "meanL1");
    if ((GetVar(pdf,"sigmaL1")        != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, pdf, "sigmaL1");
    if ((GetVar(pdf,"meanL2")         != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, pdf, "meanL2");
    if ((GetVar(pdf,"sigmaL2")        != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, pdf, "sigmaL2");
    if ((GetVar(pdf,"fracMassBLPeak") != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, pdf, "fracMassBLPeak");

    if ((GetVar(pdf,"fracMassBPeak")  != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, pdf, "fracMassBPeak");

    if (((strcmp(CTRLfitWRKflow.c_str(),"trueAll") == 0) || (strcmp(CTRLfitWRKflow.c_str(),"allEvts") == 0)) && (atoi(Utility->GetGenericParam("CtrlMisTagWrkFlow").c_str()) == 1))
      if ((GetVar(pdf,"nMisTagFrac") != NULL) && ((varName == "All") || (varName == "mistag"))) AddGaussConstraint(vecConstr, pdf, "nMisTagFrac");

    if ((GetVar(pdf,"nBkgPeak") != NULL) && ((varName == "All") || (varName == "peak"))) AddGaussConstraint(vecConstr, pdf, "nBkgPeak");
  }
}


void BuildAngularConstraints (RooArgSet* vecConstr, RooAbsPdf* pdf, string varName, vector<vector<string>*>* fitParam, unsigned int q2BinIndx)
  // #######################################################
  // // # sign: apply constraint to angular S-wave signal     #
  // // # comb: apply constraint to angular combinatorial bkg #
  // // # peak: apply constraint to angular peaking bkg       #
  // // #######################################################
{
  stringstream myString;


  if (atoi(Utility->GetGenericParam("ApplyConstr").c_str()) == true)
  {
    if (((varName == "peak") && ((GetVar(pdf,"nSig") != NULL) || (GetVar(pdf,"nMisTagFrac") != NULL))) || (varName == "comb"))
    {
      string polyType;
      if (varName == "peak") polyType = "p";
      else                   polyType = "c";


      for (unsigned int i = 0; i < NCOEFFPOLYBKG; i++)
      {
        myString.clear(); myString.str("");
        myString << polyType.c_str() << "1Poly" << i;
        if (GetVar(pdf,myString.str().c_str()) != NULL) AddGaussConstraint(vecConstr, pdf, myString.str().c_str());
        myString.clear(); myString.str("");
        myString << polyType.c_str() << "2Poly" << i;
        if (GetVar(pdf,myString.str().c_str()) != NULL) AddGaussConstraint(vecConstr, pdf, myString.str().c_str());
        myString.clear(); myString.str("");
        myString << polyType.c_str() << "3Poly" << i;
        if (GetVar(pdf,myString.str().c_str()) != NULL) AddGaussConstraint(vecConstr, pdf, myString.str().c_str());
      }
    }
    else if (varName == "sign")
    {
      double value, errLo, errHi;

      myString.clear(); myString.str("");
      myString << (atof(fitParam->operator[](Utility->GetFitParamIndx("FsS"))->operator[](q2BinIndx).c_str())) << "   " << errLo << "   " << errHi;
      SetValueAndErrors(pdf,"FsS",1.0,&myString,&value,&errLo,&errHi);

      if (GetVar(pdf,"FsS") != NULL) GetVar(pdf,"FsS")->setConstant(true);

      myString.clear(); myString.str("");
      myString << (atof(fitParam->operator[](Utility->GetFitParamIndx("AsS"))->operator[](q2BinIndx).c_str())) << "   " << errLo << "   " << errHi;
      SetValueAndErrors(pdf,"AsS",1.0,&myString,&value,&errLo,&errHi);

      if (GetVar(pdf,"AsS") != NULL) GetVar(pdf,"AsS")->setConstant(true);
    }
    else if ((varName != "sign") && (varName == "comb") && (varName == "peak"))
    {
      cout << "[ExtractYield::BuildAngularConstraints]\tWrong parameter: " << varName << endl;
      exit (EXIT_FAILURE);
    }
  }
}

RooAbsPdf* MakeAngWithEffPDF (unsigned int q2BinIndx, RooRealVar* y, RooRealVar* z, RooRealVar* p, unsigned int FitType, bool useEffPDF, RooArgSet* VarsAng)
  // ###################
  // # y: cos(theta_l) #
  // # z: cos(theta_K) #
  // # p: phi
  // ###################
{
  cout << "MakeAngWithEffPDF" << endl;
  stringstream myString;
  vector<RooRealVar*> vecParam;
  RooAbsPdf* AnglesPDF = 0;
  double a, b, c;

  if ((FitType == 1) ||(FitType == 6) || (FitType == 206) || (FitType == 106)) 
  {
    // #####################################
    // # Make 3D angular*efficiency p.d.f. #
    // # For correctly tagged events       #
    // #####################################
    FlS  = new RooRealVar("FlS","F_{L}",0.5,0.0,1.0);
    P5pS = new RooRealVar("P5pS","P_{5p}",0.5,-10.0,10.);
    P1S = new RooRealVar("P1S","P_{1}",0.0,-1.0,1.);    
    VarsAng->add(*FlS);
    VarsAng->add(*P5pS);
    VarsAng->add(*P1S);          
    VarsAng->add(*y);
    VarsAng->add(*z);
    VarsAng->add(*p);     

    myString.clear(); myString.str("");
    if (atoi(Utility->GetGenericParam("UseSPwave").c_str()) == false)
    {
      // #####################
      // # P-wave decay rate #
      // #####################
      myString << "(9/(8*3.14159265) * (2 * " << "FlS" << " * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * " << z->getPlotLabel() << "*" << z->getPlotLabel() << " + ";
      myString << "1/2 * (1-" << "FlS" << ") * (1+" << y->getPlotLabel() << "*" << y->getPlotLabel() << ")* " << "(1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") +";
      myString << "1/2*" << "P1S" <<"*(1-" << "FlS" << ") * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() <<")*(1-" << y->getPlotLabel() << "* "<<y->getPlotLabel() << ")* cos(2*" << p->getPlotLabel()<< " ) + " ;
      myString << "2*" << "P5pS" <<"* cos(" << p->getPlotLabel()<< ")*" << z->getPlotLabel() << " * sqrt("<< "FlS" << "* (1-" << "FlS" << ")*(1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")*(1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << "))))";    
    }
    else
    {
      FsS = new RooRealVar("FsS","F_{S}",0.0,0.0,1.0);
      AsS = new RooRealVar("AsS","A_{S}",0.0,-1.0,1.0);
      As5S = new RooRealVar("As5S","A_{s5}",0.0,-1.0,1.0);    
      VarsAng->add(*FsS);
      VarsAng->add(*AsS);
      VarsAng->add(*As5S);
      // ###########################
      // # S and P-wave decay rate #
      // ###########################
      myString << "(9/(8 * 3.14159265) * (2/3 *((" <<" FsS" << " + " <<" AsS" << "*" << z->getPlotLabel() << ") * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") +" << Utility->Transformer("As5S",false,a,b,c) << "*sqrt((1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ")*(1-" <<z->getPlotLabel() << "*" << z->getPlotLabel()<< "))*cos("<<p->getPlotLabel() <<")) +";
      myString << "(1-" << "FsS" << ") * ";
      myString << "(2 * " <<" FlS" << " * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * " << z->getPlotLabel() << "*" << z->getPlotLabel() << " + ";
      myString << "1/2 * (1-" <<" FlS"<< ") * (1+" << y->getPlotLabel() << "*" << y->getPlotLabel() << ")* " << "(1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") +";
      myString << "1/2*" << "P1S" <<"*(1-" << "FlS" << ") * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() <<")*(1-" << y->getPlotLabel() << "* "<<y->getPlotLabel() << ")* cos(2*" << p->getPlotLabel()
        << " ) + " ;
      myString << "2*" << "P5pS" <<"* cos(" << p->getPlotLabel()<< ")*" << z->getPlotLabel() << " * sqrt("<< "FlS" << "* (1-" << "FlS" << ")*(1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")*(1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ")))))";
    }

    cout << "\n[ExtractYield::MakeAngWithEffPDF]\t@@@ 3D angular p.d.f. @@@" << endl;
    cout << myString.str().c_str() << endl;
    if (useEffPDF == true)
    {
      RooGenericPdf* _AnglesPDF = new RooGenericPdf("_AnglesPDF",myString.str().c_str(),RooArgSet(*VarsAng));
      _AnglesPDF->Print("v");
      // #############################
      // # Make 3D efficiency p.d.f. #
      // #############################
      RooAbsPdf*  EffPdf_R = Utility->ReadRTEffPDF(q2BinIndx, z, y,p);
      MyProdPdf* myprodpdf = new MyProdPdf (*_AnglesPDF, *EffPdf_R);
      ROOT::Math::Functor* prodFunctor = new ROOT::Math::Functor(*myprodpdf,myprodpdf->ndim());
      AnglesPDF  = new RooFunctorPdfBinding("AngleS","Signal * Efficiency",*prodFunctor,myprodpdf->vars());
      cout << "\n[ExtractYield::MakeAngWithEffPDF]\t@@@ 3D angular*efficiency p.d.f. @@@" << endl;
      cout << myString.str().c_str() << endl;
    }
    else 
      AnglesPDF = new RooGenericPdf("AngleS",myString.str().c_str(),RooArgSet(*VarsAng));
  }
  else if ((FitType == 1*10) ||(FitType == 6*10) || (FitType == 206*10) || (FitType == 106*10)) 
  {
    cout << "qui" << endl;
    // ####################################
    // # Make 3D angular*fficiency p.d.f. #
    // # For incorrectly tagged events    #
    // ####################################

    myString.clear(); myString.str("");
    if (atoi(Utility->GetGenericParam("UseSPwave").c_str()) == false)
    {
      // #####################
      // # P-wave decay rate #
      // #####################
      myString << "(9/(8*3.14159265) * (2 * " << "FlS" << " * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * " << z->getPlotLabel() << "*" << z->getPlotLabel() << " + ";
      myString << "1/2 * (1-" << "FlS" << ") * (1+" << y->getPlotLabel() << "*" << y->getPlotLabel() << ")* " << "(1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") +";
      myString << "1/2*" << "P1S" <<"*(1-" <<"FlS" << ") * (1-"<< z->getPlotLabel() << "*" << z->getPlotLabel() <<")*(1-" << y->getPlotLabel() << "* "<<y->getPlotLabel() << ")* cos(2*" << p->getPlotLabel() << " ) - " ;
      myString << "2*" << "P5pS" <<"* cos(" << p->getPlotLabel()<< ")*" << z->getPlotLabel() << " * sqrt("<< "FlS" << "* (1-" << "FlS" << ")*(1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")*(1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << "))))";
    }
    else
    {
      // ###########################
      // # S and P-wave decay rate #
      // ###########################
      myString << "(9/(8 * 3.14159265) * (2/3 *((" <<" FsS" << " - " << "AsS" << "*" << z->getPlotLabel() << ") * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") + "<< Utility->Transformer("As5S",false,a,b,c) << "*sqrt((1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ")*(1-" <<z->getPlotLabel() << "*" << z->getPlotLabel()<< "))*cos("<< p->getPlotLabel() <<")) +";
      myString << "(1-" << "FsS" << ") * ";
      myString << "(2 * " << "FlS" << " * (1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ") * " << z->getPlotLabel() << "*" << z->getPlotLabel() << " + ";
      myString << "1/2 * (1-" << "FlS" << ") * (1+" << y->getPlotLabel() << "*" << y->getPlotLabel() << ")* " << "(1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ") +";
      myString << "1/2*" << "P1S" <<"*(1-" << "FlS" << ") * (1-" << z->getPlotLabel() << "*" << z->getPlotLabel() <<")*(1-" << y->getPlotLabel() << "* "<<y->getPlotLabel() << ")* cos(2*" << p->getPlotLabel() << " ) - " ;
      myString << "2*" << "P5pS" <<"* cos(" << p->getPlotLabel()<< ")*" << z->getPlotLabel() << " * sqrt("<<"FlS" << "* (1-" << "FlS" << ")*(1-" << z->getPlotLabel() << "*" << z->getPlotLabel() << ")*(1-" << y->getPlotLabel() << "*" << y->getPlotLabel() << ")))))";
    } 

    cout << "\n[ExtractYield::MakeAngWithEffPDF]\t@@@ 3D angular p.d.f. @@@" << endl;
    cout << myString.str().c_str() << endl;
    if (useEffPDF == true)
    {
      RooGenericPdf*_AnglesPDF = new RooGenericPdf("_AnglesPDF",myString.str().c_str(),RooArgSet(*VarsAng));
      // #############################
      // # Make 3D efficiency p.d.f. #
      // #############################
      RooAbsPdf* EffPdf_W = Utility->ReadWTEffPDF( q2BinIndx, z , y,p);
      MyProdPdf* myprodpdf = new MyProdPdf(*_AnglesPDF, *EffPdf_W);
      ROOT::Math::Functor* prodFunctor = new ROOT::Math::Functor(*myprodpdf,myprodpdf->ndim());
      AnglesPDF = new RooFunctorPdfBinding("AnglesM","MisTag * Efficiency",*prodFunctor,myprodpdf->vars());
      cout << "\n[ExtractYield::MakeAngWithEffPDF]\t@@@ 3D angular*efficiency p.d.f. @@@" << endl;
      cout << myString.str().c_str() << endl;
    }
    else 
      AnglesPDF = new RooGenericPdf("AngleM",myString.str().c_str(),RooArgSet(*VarsAng));
  }


  vecParam.clear();
  return AnglesPDF;
}


unsigned int CopyFitResults (RooAbsPdf* pdf, unsigned int q2BinIndx, vector<vector<string>*>* fitParam,unsigned int countMisTag, unsigned int countGoodTag)
  // ####################################################################
  // # Only for polynomial coeficients:                                 #
  // # If errLo and errHi are 0.0 --> set poly. coefficient as constant #
  // ####################################################################
{
  stringstream myString;
  stringstream myCoeff;
  double value, errLo, errHi;
  unsigned int NCoeffPolyBKGpeak1;
  unsigned int NCoeffPolyBKGcomb1;
  unsigned int NCoeffPolyBKGpeak2;
  unsigned int NCoeffPolyBKGcomb2;
  unsigned int NCoeffPolyBKGpeak3;
  unsigned int NCoeffPolyBKGcomb3;

  TRandom RG(scanIndx+1);


  if (GetVar(pdf,"meanS") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("meanS"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"meanS",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"meanS")->setConstant(true);
  }
  if (GetVar(pdf,"sigmaS1") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaS1"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"sigmaS1",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"sigmaS1")->setConstant(true);
  }
  if (GetVar(pdf,"sigmaS2") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaS2"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"sigmaS2",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"sigmaS2")->setConstant(true);
  }
  if (GetVar(pdf,"fracMassS") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("fracMassS"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"fracMassS",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"fracMassS")->setConstant(true); 
  }


  if (GetVar(pdf,"var1") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("var1"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"var1",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"var1")->setConstant(true);
  }
  if (GetVar(pdf,"var2") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("var2"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"var2",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"var2")->setConstant(true);
  }
  if (GetVar(pdf,"fracMassBExp") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("fracMassBExp"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"fracMassBExp",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"fracMassBExp")->setConstant(true);
  }

  if (GetVar(pdf,"nBkgComb") != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("nBkgPeak"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(pdf,"nBkgComb",MULTYIELD,&myString,&value,&errLo,&errHi);
      GetVar(pdf,"nBkgComb")->setConstant(false);
    }
  

  if (GetVar(pdf,"sigmaMisTag1") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag1"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"sigmaMisTag1",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"sigmaMisTag1")->setConstant(true);
  }
  if (GetVar(pdf,"sigmaMisTag2") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaMisTag2"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"sigmaMisTag2",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"sigmaMisTag2")->setConstant(true);
  }
  if (GetVar(pdf,"fracMisTag") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("fracMisTag"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"fracMisTag",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"fracMisTag")->setConstant(true);
  }


  if (GetVar(pdf,"meanR1") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("meanR1"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"meanR1",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"meanR1")->setConstant(true);
  }

  if (GetVar(pdf,"sigmaR1") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaR1"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"sigmaR1",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"sigmaR1")->setConstant(true);
  }
  if (GetVar(pdf,"meanR2") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("meanR2"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"meanR2",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"meanR2")->setConstant(true);
  }
  if (GetVar(pdf,"sigmaR2") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaR2"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"sigmaR2",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"sigmaR2")->setConstant(true);
  }
  if (GetVar(pdf,"fracMassBRPeak") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("fracMassBRPeak"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"fracMassBRPeak",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"fracMassBRPeak")->setConstant(true);
  }
  if (GetVar(pdf,"meanL1") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("meanL1"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"meanL1",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"meanL1")->setConstant(true);
  }
  if (GetVar(pdf,"sigmaL1") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaL1"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"sigmaL1",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"sigmaL1")->setConstant(true);
  }
  if (GetVar(pdf,"meanL2") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("meanL2"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"meanL2",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"meanL2")->setConstant(true);
  }
  if (GetVar(pdf,"sigmaL2") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("sigmaL2"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"sigmaL2",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"sigmaL2")->setConstant(true);
  }
  if (GetVar(pdf,"fracMassBLPeak") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("fracMassBLPeak"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"fracMassBLPeak",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"fracMassBLPeak")->setConstant(true);
  }
  if (GetVar(pdf,"fracMassBPeak") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("fracMassBPeak"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"fracMassBPeak",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"fracMassBPeak")->setConstant(true);
  }

  if (GetVar(pdf,"nMisTagFrac") != NULL)
  {
    if (atoi(Utility->GetGenericParam("CtrlMisTagWrkFlow").c_str()) != 3 ||(countMisTag==0 && countGoodTag==0))
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx("nMisTagFrac"))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(pdf,"nMisTagFrac",1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,"nMisTagFrac")->setConstant(true); 
    }
    else
    {
      myString.clear(); myString.str("");
      myString << static_cast<double>(countMisTag) / static_cast<double>(countMisTag + countGoodTag);
      SetValueAndErrors(pdf,"nMisTagFrac",1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,"nMisTagFrac")->setConstant(true);
    }                    
  }
  if (GetVar(pdf,"nBkgPeak") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("nBkgPeak"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"nBkgPeak",MULTYIELD,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"nBkgPeak")->setConstant(true);
  }
  if (GetVar(pdf,"nSig") != NULL)
  {
    // #if you just want to fit the data#
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("nSig"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"nSig",MULTYIELD,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"nSig")->setConstant(false);
    //#if you want to scan values#
    // myString.clear(); myString.str("");
    // myString << RooRandom::uniform(&RG);
    // SetValueAndErrors(pdf,"nSig",MULTYIELD,&myString,&value,&errLo,&errHi);
    // GetVar(pdf,"nSig")->setConstant(false);
    // cout<<"nsig=" << myString.str()<<endl;
    // fileFitResults << "nsig=" << myString.str() << endl;
  }


  NCoeffPolyBKGpeak1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP1"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGcomb1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC1"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGpeak2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP2"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGcomb2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC2"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGpeak3 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP3"))->operator[](q2BinIndx).c_str());
  NCoeffPolyBKGcomb3 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC3"))->operator[](q2BinIndx).c_str());
  if ((NCoeffPolyBKGpeak1 > NCOEFFPOLYBKG) || (NCoeffPolyBKGcomb1 > NCOEFFPOLYBKG) ||
      (NCoeffPolyBKGpeak2 > NCOEFFPOLYBKG) || (NCoeffPolyBKGcomb2 > NCOEFFPOLYBKG) ||
      (NCoeffPolyBKGpeak3 > NCOEFFPOLYBKG) || (NCoeffPolyBKGcomb3 > NCOEFFPOLYBKG))
  {
    cout << "[ExtractYield::CopyFitResults]\tDegree of poly bkg is not within allowed limits : ";
    cout << NCoeffPolyBKGpeak1 << "\t" << NCoeffPolyBKGcomb1 << "\t" << NCoeffPolyBKGpeak2 << "\t" << NCoeffPolyBKGcomb2 << "\t" << NCoeffPolyBKGpeak3 << "\t" << NCoeffPolyBKGcomb3 << endl;
    exit (EXIT_FAILURE);
  }


  for (unsigned int i = 0; i < NCoeffPolyBKGpeak1; i++)
  {
    myCoeff.clear(); myCoeff.str("");
    myCoeff << "p1Poly" << i;
    if (GetVar(pdf,myCoeff.str().c_str()) != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(pdf,myCoeff.str(),1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,myCoeff.str().c_str())->setConstant(true);
    }
  }
  for (unsigned int i = 0; i < NCoeffPolyBKGcomb1; i++)
  {
    myCoeff.clear(); myCoeff.str("");
    myCoeff << "c1Poly" << i;
    if (GetVar(pdf,myCoeff.str().c_str()) != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(pdf,myCoeff.str(),1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,myCoeff.str().c_str())->setConstant(true);
    }
  }
  for (unsigned int i = 0; i < NCoeffPolyBKGpeak2; i++)
  {
    myCoeff.clear(); myCoeff.str("");
    myCoeff << "p2Poly" << i;
    if (GetVar(pdf,myCoeff.str().c_str()) != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(pdf,myCoeff.str(),1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,myCoeff.str().c_str())->setConstant(true);
    }
  }
  for (unsigned int i = 0; i < NCoeffPolyBKGcomb2; i++)
  {
    myCoeff.clear(); myCoeff.str("");
    myCoeff << "c2Poly" << i;
    if (GetVar(pdf,myCoeff.str().c_str()) != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(pdf,myCoeff.str(),1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,myCoeff.str().c_str())->setConstant(true);
    }
  }
  for (unsigned int i = 0; i < NCoeffPolyBKGpeak3; i++)
  {
    myCoeff.clear(); myCoeff.str("");
    myCoeff << "p3Poly" << i;
    if (GetVar(pdf,myCoeff.str().c_str()) != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(pdf,myCoeff.str(),1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,myCoeff.str().c_str())->setConstant(true);
    }
  }
  for (unsigned int i = 0; i < NCoeffPolyBKGcomb3; i++)
  {
    myCoeff.clear(); myCoeff.str("");
    myCoeff << "c3Poly" << i;
    if (GetVar(pdf,myCoeff.str().c_str()) != NULL)
    {
      myString.clear(); myString.str("");
      myString << fitParam->operator[](Utility->GetFitParamIndx(myCoeff.str().c_str()))->operator[](q2BinIndx).c_str();
      SetValueAndErrors(pdf,myCoeff.str(),1.0,&myString,&value,&errLo,&errHi);
      GetVar(pdf,myCoeff.str().c_str())->setConstant(true);
    }
  }

  if (GetVar(pdf,"FlS") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("FlS"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"FlS",1.0,&myString,&value,&errLo,&errHi);
    // if (scanInitVal) pdf->getVariables()->setRealValue("FlS",pdf->getVariables()->getRealValue("FlS")+RooRandom::gaussian(&RG)*3*scanStep1[q2BinIndx]);
    GetVar(pdf,"FlS")->setConstant(true);
  }
  if (GetVar(pdf,"P5pS") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("P5pS"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"P5pS",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"P5pS")->setConstant(false);
    if (!PROFILENLL && scanInitVal) {
      // double inVal = 100;
      // do {
      // 	inVal = pdf->getVariables()->getRealValue("P5pS")+RooRandom::gaussian(&RG)*scanStep3[q2BinIndx];
      // } while (inVal<-10. || inVal>10.); 
      // pdf->getVariables()->setRealValue("P5pS",inVal);
      pdf->getVariables()->setRealValue("P5pS",-0.8*RooRandom::uniform(&RG));
    }
    if (PROFILENLL) {
      if (scanIndx%2==1) {
	pdf->getVariables()->setRealValue("P5pS", minP5p[q2BinIndx] + (scanIndx-1.)/(scanIndxMax-2.)*(maxP5p[q2BinIndx]-minP5p[q2BinIndx]) );
	GetVar(pdf,"P5pS")->setConstant(true);
      }
    }

    //### scan ##
    // myString.clear(); myString.str("");
    // myString << 2 * RooRandom::uniform(&RG)  - 1.;
    // SetValueAndErrors(pdf,"P5pS",1.0,&myString,&value,&errLo,&errHi);
    // GetVar(pdf,"P5pS")->setConstant(false);
    // cout<<"P5p="<<myString.str()<<endl;
    // fileFitResults << "P5p="<<myString.str() << endl;
    
  }
  if (GetVar(pdf,"P1S") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("P1S"))->operator[](q2BinIndx).c_str();
    SetValueAndErrors(pdf,"P1S",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"P1S")->setConstant(false);
    if (!PROFILENLL && scanInitVal) {
      // double inVal = 100;
      // do {
      // 	inVal = pdf->getVariables()->getRealValue("P1S")+RooRandom::gaussian(&RG)*scanStep2[q2BinIndx];
      // } while (inVal<-1. || inVal>1.); 
      // pdf->getVariables()->setRealValue("P1S",inVal);
      pdf->getVariables()->setRealValue("P1S",-0.8-pdf->getVariables()->getRealValue("P5pS"));
    }
    if (PROFILENLL) {
      if (scanIndx%2==0) {
	pdf->getVariables()->setRealValue("P1S", minP1[q2BinIndx] + scanIndx/(scanIndxMax-2.)*(maxP1[q2BinIndx]-minP1[q2BinIndx]) );
	GetVar(pdf,"P1S")->setConstant(true);
	pdf->getVariables()->setRealValue("P5pS", -0.5-0.7*pdf->getVariables()->getRealValue("P1S"));
      } else {
	double inVal = -0.8-1.3*pdf->getVariables()->getRealValue("P5pS");
	if (inVal < -0.8) inVal = -0.8;
	pdf->getVariables()->setRealValue("P1S", inVal);
      }
    }
    //#scan
    // myString.clear(); myString.str("");
    // myString << 2 * RooRandom::uniform(&RG)  - 1.;
    // SetValueAndErrors(pdf,"P1S",1.0,&myString,&value,&errLo,&errHi);
    // GetVar(pdf,"P1S")->setConstant(false);
    // cout<<"P1="<<myString.str()<<endl;
    // fileFitResults << "P1="<<myString.str() << endl;
    
  }
  if (GetVar(pdf,"FsS") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("FsS"))->operator[](q2BinIndx);
    SetValueAndErrors(pdf,"FsS",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"FsS")->setConstant(true);
  }
  if (GetVar(pdf,"AsS") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("AsS"))->operator[](q2BinIndx);
    SetValueAndErrors(pdf,"AsS",1.0,&myString,&value,&errLo,&errHi);
    GetVar(pdf,"AsS")->setConstant(true);
  }
  if (GetVar(pdf,"As5S") != NULL)
  {
    myString.clear(); myString.str("");
    myString << fitParam->operator[](Utility->GetFitParamIndx("As5S"))->operator[](q2BinIndx);
    SetValueAndErrors(pdf,"As5S",1.0,&myString,&value,&errLo,&errHi);
    if (scanInitVal) {
      // double inVal = 100;
      // do {
      // 	inVal = pdf->getVariables()->getRealValue("As5S")+RooRandom::gaussian(&RG)*scanStep4[q2BinIndx];
      // } while (inVal<-1. || inVal>1.); 
      // pdf->getVariables()->setRealValue("As5S",inVal);
      pdf->getVariables()->setRealValue("As5S",0);
    }
    GetVar(pdf,"As5S")->setConstant(false);


    // myString.clear(); myString.str("");
    // myString << RooRandom::uniform() * 2. - 1.;
    // SetValueAndErrors(pdf,"As5S",1.0,&myString,&value,&errLo,&errHi);
    // GetVar(pdf,"As5S")->setConstant(false); 
    // cout<<"As5="<<myString.str()<<endl;
    // fileFitResults << "As5="<<myString.str() << endl;
    
  }
  value = 0.0;
  if (GetVar(pdf,"nSig")         != NULL)  value = value + GetVar(pdf,"nSig")->getVal();
  if (GetVar(pdf,"nMisTagFrac") != NULL) value = value +  GetVar(pdf,"nMisTagFrac")->getVal();
  if (GetVar(pdf,"nBkgPeak")     != NULL)  value = value + GetVar(pdf,"nBkgPeak")->getVal();
  if (GetVar(pdf,"nBkgComb")     != NULL)  value = value + GetVar(pdf,"nBkgComb")->getVal();
  return value;
}

string GeneratePolynomial (RooRealVar* var, unsigned int nCoef, string sCoef)
{
  stringstream myString;


  myString.clear(); myString.str("");
  myString << "1";
  for (unsigned int i = 0; i < nCoef; i++)
  {
    myString << " + ";
    for (unsigned int j = 0; j < i+1; j++) myString << var->getPlotLabel() << "*";
    myString << sCoef << i;
  }
  cout << "[ExtractYield::GeneratePolynomial]\tI've generated the polynomial: " << myString.str().c_str() << endl;

  return myString.str();
}



void LoadDatasets (TFile* NtplFile, int specBin, int FitType, int toyIndx)
{
  cout << "LoadDatasets " << Form("toy%i_sig%i",toyIndx,specBin) << endl;
  RooWorkspace* w = (RooWorkspace*)NtplFile->Get("ws");
  cout << "RooWorkspace " << w << endl;
  SingleCandNTuple_RejectPsi = SingleCandNTuple = (RooDataSet*)w->data(Form("toy%i_sig%i",toyIndx,specBin));
  //if (FitType == 6 ||FitType == 1) SingleCandNTuple_RejectPsi->append(*(RooDataSet*)w->data(Form("toy0_bkg%i",specBin+1)));
  cout << "ntuple " << SingleCandNTuple->numEntries() << endl;

  // ###########################
  // # Define useful variables #
  // ###########################
  B0MassArb = w->var("B0MassArb");
  // B0MassArb = w->var("B0mass");
  ctK       = w->var("ctK");
  ctL       = w->var("ctL");
  phi       = w->var("phi");

  // ####################################################
  // # Setting initial values for independent variables #
  // ####################################################
  B0MassArb->setVal(Utility->B0Mass);
  ctK ->setVal(0.0);
  ctL->setVal(0.0);
  phi->setVal(0.0);
}


void MakeDatasets (B0KstMuMuSingleCandTreeContent* NTuple, unsigned int FitType)
{
  stringstream myString;

  //  RooDataSet* SingleCandNTuple;


  // ###########################
  // # Define useful variables #
  // ###########################
  B0MassArb          = new RooRealVar("B0MassArb",
                                      "#font[12]{m}(#font[122]{K}#kern[0.1]{#lower[0.4]{^{#font[122]{+}}}}#kern[-0.3]{#pi}#kern[-0.3]{#lower[0.6]{^{#font[122]{\55}}}}#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#kern[-0.1]{#mu}#kern[-1.3]{#lower[0.6]{^{#font[122]{\55}}}})",
                                      Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str()),
                                      Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str()),
                                      "GeV");
  mumuMass           = new RooRealVar("mumuMass","#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#kern[-0.1]{#mu}#kern[-1.3]{#lower[0.6]{^{#font[122]{\55}}}} inv. mass",0.0,6.0,"GeV");
  mumuMassE          = new RooRealVar("mumuMassE","#mu#kern[-0.9]{#lower[0.6]{^{#font[122]{+}}}}#kern[-0.1]{#mu}#kern[-1.3]{#lower[0.6]{^{#font[122]{\55}}}} inv. mass error",0.0,0.5,"GeV");
  ctK                = new RooRealVar("ctK","cos#theta_{K}",-1.0,1.0);
  ctL                = new RooRealVar("ctL","cos#theta_{L}",0.00,1.0);
  phi                = new RooRealVar("phi","#phi",0.0,Utility->PI);
  truthMatchSignal   = new RooRealVar("truthMatchSignal","Truth matching",0.0,1.0,"bool");
  rightFlavorTag     = new RooRealVar("rightFlavorTag","Right flavor tag",0.0,1.0,"bool");


  if ((FitType == 1) ||(FitType == 6) || (FitType == 106) || (FitType == 206))
  {
    RooArgSet Vars;
    Vars.add(*B0MassArb);
    Vars.add(*mumuMass);
    Vars.add(*mumuMassE);
    Vars.add(*ctK);
    Vars.add(*ctL);
    Vars.add(*phi);
    Vars.add(*truthMatchSignal);
    Vars.add(*rightFlavorTag);

    SingleCandNTuple           = new RooDataSet("SingleCandNTuple"          ,"SingleCandNTuple"          ,Vars);
    SingleCandNTuple_RejectPsi = new RooDataSet("SingleCandNTuple_RejectPsi","SingleCandNTuple_RejectPsi",Vars);



    // #############################
    // # Load values from the tree #
    // #############################
    NTuple->ClearNTuple();
    NTuple->SetBranchAddresses(theTree);
    int nEntries = theTree->GetEntries();
    cout << "[ExtractYield::MakeDatasets]\t@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;
    for (int entry = 0; entry < nEntries; entry++)
      {
	theTree->GetEntry(entry);
	
	if ((((FitType == 1) ||(FitType == 6)||(FitType == 106))  &&
	     (NTuple->B0MassArb > Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str())) &&
	     (NTuple->B0MassArb < Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str())))
	    ||(FitType == 206)
	    )
	  {
	    Vars.setRealValue("B0MassArb",         NTuple->B0MassArb);
	    Vars.setRealValue("mumuMass",          NTuple->mumuMass->at(0));
	    Vars.setRealValue("mumuMassE",         NTuple->mumuMassE->at(0));
	    Vars.setRealValue("ctK",      NTuple->CosThetaKArb);
	    Vars.setRealValue("ctL",   fabs(NTuple->CosThetaMuArb));
	    Vars.setRealValue("phi",fabs(NTuple->PhiKstMuMuPlaneArb));
	    Vars.setRealValue("truthMatchSignal",  NTuple->truthMatchSignal->at(0));
	    Vars.setRealValue("rightFlavorTag",    NTuple->rightFlavorTag);

	    
	    // ########################
	    // # NTuple with all data #
	    // ########################
	    SingleCandNTuple->add(Vars);
	    
	    
	    // #############################################################################
	    // # J/psi and psi(2S) rejection based on the event-by-event dimuon mass error #
	    // #############################################################################
	    if ((((FitType == 1) || (FitType == 6)) &&
		 
		 (((strcmp(CTRLfitWRKflow.c_str(),"trueGoodtag") == 0) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == true)) ||
		  ((strcmp(CTRLfitWRKflow.c_str(),"trueMistag")  == 0) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == false)) ||
		  ((strcmp(CTRLfitWRKflow.c_str(),"trueAll")     == 0) && (NTuple->truthMatchSignal->at(0) == true)) ||
		  (strcmp(CTRLfitWRKflow.c_str(),"allEvts")      == 0)) &&
		 (Utility->PsiRejection(NTuple->B0MassArb,NTuple->mumuMass->at(0),NTuple->mumuMassE->at(0),"rejectPsi",true) == true)) ||
		
		(((FitType == 106))  &&
		 (((strcmp(CTRLfitWRKflow.c_str(),"trueAll") == 0) && (NTuple->truthMatchSignal->at(0) == true))||
		  ((strcmp(CTRLfitWRKflow.c_str(),"trueGoodtag") == 0) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == true)) ||
		  ((strcmp(CTRLfitWRKflow.c_str(),"trueMistag")  == 0) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag == false))||
		  (strcmp(CTRLfitWRKflow.c_str(),"allEvts")      == 0)) &&
		 
		 (Utility->PsiRejection(NTuple->B0MassArb,NTuple->mumuMass->at(0),NTuple->mumuMassE->at(0),"rejectPsi",true) == true))
		)	

	    SingleCandNTuple_RejectPsi->add(Vars);
	  }
      }
    cout << "\n[ExtractYield::MakeDatasets]\t@@@ NTuple with all data @@@" << endl;
    SingleCandNTuple->Print("v");

    cout << "\n[ExtractYield::MakeDatasets]\t@@@ NTuple without J/psi and psi(2S) regions @@@" << endl;
    SingleCandNTuple_RejectPsi->Print("v");

  }
  // ####################################################
  // # Setting initial values for independent variables #
  // ####################################################
  B0MassArb->setVal(Utility->B0Mass);
  ctK ->setVal(0.0);
  ctL->setVal(0.0);
  phi->setVal(0.0);
}





//#############################
//# GEN level fitting  PDF #
//#############################
void InstantiateGen3AnglesFit    (RooAbsPdf** TotalPDF,
                                  bool useEffPDF,
                                  RooRealVar* y, RooRealVar* z,RooRealVar* p,
                                  unsigned int FitType,
                                  vector<vector<unsigned int>*>* configParam,
                                  unsigned int q2BinIndx)
// #########################
// # y: angle cos(theta_l) #
// # z: angle cos(theta_K) #
// # p: angle phi          #
// #########################
{
  // ################################
  // # Read configuration variables #
  // ################################
  unsigned int useSignal = configParam->operator[](Utility->GetConfigParamIndx("SigType"))->operator[](q2BinIndx);


  // ###################
  // # Local variables #
  // ###################
  stringstream myString;
  RooArgSet* VarsAng = new RooArgSet("VarsAng");

  if ((FitType == 206)) *TotalPDF = MakeAngWithEffPDF(q2BinIndx,y,z,p,FitType,useEffPDF,VarsAng);
  else
  {
    cout << "[ExtractYield::InstantiateGen3AnglesFit]\tIncorrect configuration sequence : useSignal = " << useSignal <<  endl;
    exit (EXIT_FAILURE);
  }
}

//#############################
//# Reco level fitting  PDF #
//#############################
void InstantiateReco3AnglesFit (RooAbsPdf** TotalPDF,
                                bool useEffPDF,
                                RooRealVar* y, RooRealVar* z,RooRealVar* p,
                                unsigned int FitType,
                                vector<vector<unsigned int>*>* configParam,
                                unsigned int q2BinIndx)
// #########################
// # y: angle cos(theta_l) #
// # z: angle cos(theta_K) #
// # p: angle phi          #
// #########################
{
  // ################################
  // # Read configuration variables #
  // ################################
  unsigned int useSignal = configParam->operator[](Utility->GetConfigParamIndx("SigType"))->operator[](q2BinIndx);


  // ###################
  // # Local variables #
  // ###################
  stringstream myString;
  RooArgSet* VarsAng = new RooArgSet("VarsAng");

  // ##################################################################
  // # Define angle fit variables and pdf for correctly tagged signal #
  // ##################################################################
  AngleS = MakeAngWithEffPDF(q2BinIndx,y,z,p,FitType,useEffPDF,VarsAng);

  // ############################################################
  // # Define angle fit variables and pdf for mis-tagged signal #
  // ############################################################
  AngleMisTag = MakeAngWithEffPDF(q2BinIndx,y,z,p,FitType*10,useEffPDF,VarsAng);

  // ###########################
  // # Define pdf coefficients #
  // ###########################
  nMisTagFrac = new RooRealVar("nMisTagFrac","Fraction of mistag",0.0,0.0,1.0);


  if ((FitType == 106)) *TotalPDF = new RooAddPdf("SignalRECO","AngleS + nMisTag * (AngleMisTag)",RooArgSet(*AngleMisTag,*AngleS),*nMisTagFrac);
  else
  {
    cout << "[ExtractYield::InstantiateReco3AnglesFit]\tIncorrect configuration sequence : useSignal = " << useSignal <<  endl;
    exit (EXIT_FAILURE);
  }
}



RooFitResult* Make3AnglesFit (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* y, RooRealVar* z,RooRealVar* p,  unsigned int FitType, RooArgSet* vecConstr, TCanvas* Canv, unsigned int ID)
{
  cout << "FitType " << FitType<< endl;
  // ###################
  // # Local variables #
  // ###################
  RooFitResult* fitResult = NULL;
  stringstream myString;

  const unsigned int nCanv = 4;
  TCanvas* localCanv[nCanv];
  for (unsigned int i = 0; i < nCanv; i++)
  {
    myString.clear(); myString.str("");
    myString << "localCanv" << i;
    localCanv[i] = new TCanvas(myString.str().c_str(),myString.str().c_str(),20,20,500,500);
  }

  unsigned int nElements   = 0;
  unsigned int it          = 0;
  TString legNames[6]; 
  TLegend*   legY             = NULL;
  TLegend*   legZ             = NULL;
  TLegend*   legP             = NULL;

  if ((FitType == 206) || (FitType == 106))
  {
    // ###################
    // # Make actual fit #
    // ###################
    dataSet->Print("v");

    fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),Save(true),Minimizer(MINIMIZER),NumCPU(8)); 

    // ###################################################
    // # Set p.d.f. independent variables to known point #
    // ###################################################
    if (GetVar(*TotalPDF,y->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(y->getPlotLabel(),0.0);
    if (GetVar(*TotalPDF,z->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(z->getPlotLabel(),0.0);
    if (GetVar(*TotalPDF,p->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(p->getPlotLabel(),0.0);
    if (fitResult != NULL) fitResult->Print("v");

    // ###########################
    // # Costheta-l plot results #
    // ###########################
    if (PLOT)
    {
      Canv->cd(1);
      RooPlot* myFrameY = y->frame(NBINS);

      dataSet->plotOn(myFrameY, Name(MakeName(dataSet,ID).c_str()));
      if ((FitType == 206))  legNames[nElements++] = "GEN-MC";
      else                   legNames[nElements++] = "RECO-MC";
      (*TotalPDF)->plotOn(myFrameY, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), Project(RooArgSet(*z,*p)));
      if ((FitType == 206))  legNames[nElements++] = "Gen p.d.f.";
      else                   legNames[nElements++] = "Reco p.d.f.";

      TPaveText* paveTextY = new TPaveText(0.12,0.82,0.4,0.86,"NDC");
      paveTextY->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
      CheckGoodFit(fitResult,paveTextY);
      paveTextY->SetTextAlign(11);
      paveTextY->SetBorderSize(0.0);
      paveTextY->SetFillStyle(0);
      paveTextY->SetTextSize(0.04);
      paveTextY->Paint();
      myFrameY->addObject(paveTextY);

      legY = new TLegend(0.78, 0.65, 0.97, 0.88, "");
      it = 0;
      for (unsigned int i = 0; i < 2*nElements; i++)
      {
        TString objName = myFrameY->nameOf(i);
        if ((objName == "") || (objName.Contains("paramBox") == true) || (objName.Contains("TPave") == true) || ((i > 0) && (objName == myFrameY->nameOf(i-1)))) continue;
        TObject* obj = myFrameY->findObject(objName.Data());
        legY->AddEntry(obj,legNames[it++],"PL");
        legY->SetTextFont(42);
      }     
      legY->SetFillStyle(0);
      legY->SetFillColor(0);
      legY->SetTextSize(0.04);
      legY->SetBorderSize(0);

      myFrameY->Draw();
      legY->Draw("same");

      // ###########################
      // # Costheta-k plot results #
      // ###########################
      Canv->cd(2);
      RooPlot* myFrameZ = z->frame(NBINS);

      dataSet->plotOn(myFrameZ, Name(MakeName(dataSet,ID).c_str()));
      if ((FitType == 206))  legNames[nElements++] = "GEN-MC";
      else                   legNames[nElements++] = "RECO-MC";
      (*TotalPDF)->plotOn(myFrameZ, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), Project(RooArgSet(*y,*p)));
      if ((FitType == 206))  legNames[nElements++] = "Gen p.d.f.";
      else                   legNames[nElements++] = "Reco p.d.f.";

      TPaveText* paveTextZ = new TPaveText(0.12,0.82,0.4,0.86,"NDC");
      paveTextZ->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
      paveTextZ->SetTextAlign(11);
      paveTextZ->SetBorderSize(0.0);
      paveTextZ->SetFillStyle(0);
      paveTextZ->SetTextSize(0.04);
      paveTextZ->Paint();
      myFrameZ->addObject(paveTextZ);

      legZ = new TLegend(0.78, 0.65, 0.97, 0.88, "");
      it = 0;
      for (unsigned int i = 0; i < 2*nElements; i++)
      {
        TString objName = myFrameZ->nameOf(i);
        if ((objName == "") || (objName.Contains("paramBox") == true) || (objName.Contains("TPave") == true) || ((i > 0) && (objName == myFrameZ->nameOf(i-1)))) continue;
        TObject* obj = myFrameZ->findObject(objName.Data());
        legZ->AddEntry(obj,legNames[it++],"PL");
        legZ->SetTextFont(42);
      }     
      legZ->SetFillStyle(0);
      legZ->SetFillColor(0);
      legZ->SetTextSize(0.04);
      legZ->SetBorderSize(0);

      myFrameZ->Draw();
      legZ->Draw("same");


      // ###########################
      // # Costheta-l plot results #
      // ###########################
      Canv->cd(3);
      RooPlot* myFrameP = p->frame(NBINS);

      dataSet->plotOn(myFrameP, Name(MakeName(dataSet,ID).c_str()));
      if ((FitType == 206))  legNames[nElements++] = "GEN-MC";
      else                   legNames[nElements++] = "RECO-MC";
      (*TotalPDF)->plotOn(myFrameP, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), Project(RooArgSet(*y,*z)));
      if ((FitType == 206))  legNames[nElements++] = "Gen p.d.f.";
      else                   legNames[nElements++] = "Reco p.d.f.";

      TPaveText* paveTextP = new TPaveText(0.12,0.82,0.4,0.86,"NDC");
      paveTextP->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameP->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
      paveTextP->SetTextAlign(11);
      paveTextP->SetBorderSize(0.0);
      paveTextP->SetFillStyle(0);
      paveTextP->SetTextSize(0.04);
      paveTextP->Paint();
      myFrameP->addObject(paveTextP);

      legP = new TLegend(0.78, 0.65, 0.97, 0.88, "");
      it = 0;
      for (unsigned int i = 0; i < 2*nElements; i++)
      {
        TString objName = myFrameP->nameOf(i);
        if ((objName == "") || (objName.Contains("paramBox") == true) || (objName.Contains("TPave") == true) || ((i > 0) && (objName == myFrameY->nameOf(i-1)))) continue;
        TObject* obj = myFrameP->findObject(objName.Data());
        legP->AddEntry(obj,legNames[it++],"PL");
        legP->SetTextFont(42);
      }     
      legP->SetFillStyle(0);
      legP->SetFillColor(0);
      legP->SetTextSize(0.04);
      legP->SetBorderSize(0);

      myFrameP->Draw();
      legP->Draw("same");
    }
    // ##############
    // # Save plots #
    // ##############
    Canv->Modified();
    Canv->Update();

    if (SAVEPLOT == true)
    {
      myString.clear(); myString.str("");
      myString << (*TotalPDF)->getPlotLabel() << "_Canv" << ID << ".pdf";
      Canv->Print(myString.str().c_str());

      myString.clear(); myString.str("");
      myString << (*TotalPDF)->getPlotLabel() << "_Canv" << ID << ".root";
      Canv->Print(myString.str().c_str());
    }

    // ###################
    // # Profiling scan  #
    // ###################
    if (PROFILENLL == true  &&  CheckGoodFit(fitResult) == true)
    {

      string varName;
      RooAbsReal* NLL;
      if (atoi(Utility->GetGenericParam("ApplyConstr").c_str()) == true) NLL = (*TotalPDF)->createNLL(*dataSet,Extended(true),ExternalConstraints(*vecConstr));
      else                                                               NLL = (*TotalPDF)->createNLL(*dataSet,Extended(true));
      //       RooMinuit(*NLL).migrad();
      //       RooMinuit RooMin(*NLL); 
      /*     varName = "FlS";
             GetVar(*TotalPDF,varName.c_str())->setRange(GetVar(*TotalPDF,varName.c_str())->getVal() + GetVar(*TotalPDF,varName.c_str())->getErrorLo() * atof(Utility->GetGenericParam("NSigmaB0").c_str()),
             GetVar(*TotalPDF,varName.c_str())->getVal() + GetVar(*TotalPDF,varName.c_str())->getErrorHi() * atof(Utility->GetGenericParam("NSigmaB0").c_str()));
             localCanv[0]->cd();
             RooPlot* myFrameNLLVar1 = GetVar(*TotalPDF,varName.c_str())->frame();
             NLL->plotOn(myFrameNLLVar1, ShiftToZero());
             RooAbsReal* var1Profile = NLL->createProfile(*(GetVar(*TotalPDF,varName.c_str())));
             var1Profile->plotOn(myFrameNLLVar1, LineColor(kRed));
             myFrameNLLVar1->SetMinimum(0);
             myFrameNLLVar1->SetMaximum(5);
             DrawString(LUMI,myFrameNLLVar1);
             myFrameNLLVar1->Draw();
             */
      varName = "P5pS";
      GetVar(*TotalPDF,varName.c_str())->setRange(GetVar(*TotalPDF,varName.c_str())->getVal() + GetVar(*TotalPDF,varName.c_str())->getErrorLo() * atof(Utility->GetGenericParam("NSigmaB0").c_str()),
                                                  GetVar(*TotalPDF,varName.c_str())->getVal() + GetVar(*TotalPDF,varName.c_str())->getErrorHi() * atof(Utility->GetGenericParam("NSigmaB0").c_str()));
      localCanv[1]->cd();
      RooPlot* myFrameNLLVar2 = GetVar(*TotalPDF,varName.c_str())->frame();
      NLL->plotOn(myFrameNLLVar2, ShiftToZero());
      RooAbsReal* var2Profile = NLL->createProfile(*(GetVar(*TotalPDF,varName.c_str())));
      var2Profile->plotOn(myFrameNLLVar2, LineColor(kRed));
      myFrameNLLVar2->SetMinimum(0);
      myFrameNLLVar2->SetMaximum(5);
      DrawString(LUMI,myFrameNLLVar2);
      myFrameNLLVar2->DrawClone();

      varName = "P1S";
      GetVar(*TotalPDF,varName.c_str())->setRange(GetVar(*TotalPDF,varName.c_str())->getVal() + GetVar(*TotalPDF,varName.c_str())->getErrorLo() * atof(Utility->GetGenericParam("NSigmaB0").c_str()),
                                                  GetVar(*TotalPDF,varName.c_str())->getVal() + GetVar(*TotalPDF,varName.c_str())->getErrorHi() * atof(Utility->GetGenericParam("NSigmaB0").c_str()));
      localCanv[2]->cd();
      RooPlot* myFrameNLLVar3 = GetVar(*TotalPDF,varName.c_str())->frame();
      NLL->plotOn(myFrameNLLVar3, ShiftToZero());
      RooAbsReal* var3Profile = NLL->createProfile(*(GetVar(*TotalPDF,varName.c_str())));
      var3Profile->plotOn(myFrameNLLVar3, LineColor(kRed));
      myFrameNLLVar3->SetMinimum(0);
      myFrameNLLVar3->SetMaximum(5);
      DrawString(LUMI,myFrameNLLVar3);
      myFrameNLLVar3->DrawClone();


      delete NLL;
    }
    // ##############
    // # Save plots #
    // ##############
    /*    Canv->Modified();
          Canv->Update();
          */  
    if (SAVEPLOT == true)
    {
      /*     myString.clear(); myString.str("");
             myString << (*TotalPDF)->getPlotLabel() << "_Canv" << ID << ".pdf";
             Canv->Print(myString.str().c_str());

             myString.clear(); myString.str("");
             myString << (*TotalPDF)->getPlotLabel() << "_Canv" << ID << ".root";
             Canv->Print(myString.str().c_str());
             */
      for (unsigned int i = 0; i < nCanv; i++)
      {
        myString.clear(); myString.str("");
        myString << (*TotalPDF)->getPlotLabel() << "_localCanv" << i << "_" << ID << ".pdf";
        localCanv[i]->Print(myString.str().c_str());

        myString.clear(); myString.str("");
        myString << (*TotalPDF)->getPlotLabel() << "_localCanv" << i << "_" << ID << ".root";
        localCanv[i]->Print(myString.str().c_str());
      }
    }
  } 
  return fitResult;
}


void Iterative3AnglesFitq2Bins (RooDataSet* dataSet,
                                bool useEffPDF,
                                RooRealVar* y, RooRealVar* z,RooRealVar* p,
                                int specBin,
                                unsigned int FitType,
                                vector<double>* q2Bins,
                                vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam,
                                RooArgSet* vecConstr,
                                unsigned int ID)
{
  // ###################
  // # Local variables #
  // ###################
  stringstream myString;   

  RooDataSet* dataSet_q2Bins[q2Bins->size()-1];
  RooAbsPdf*  TotalPDFq2Bins[q2Bins->size()-1];
  RooFitResult* fitResult=0;


  cout << "[ExtractYield::IterativeAnglesFitq2Bins]\tInitial number of events :" << dataSet->sumEntries() << endl;
  string outname = "Fitresult_std.root";
  if (!useToyDataset) {
    if (FitType==106) outname = "Fitresult.root";
    if (FitType==206) outname = "Fitresult2.root";
  } else {
    if (FitType==106) outname = "Fitresult1.root";
    if (FitType==206) outname = "Fitresult3.root";
  }      
  TFile f(outname.c_str(),"UPDATE") ;

  for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
  {
    myString.clear(); myString.str("");
    myString << "(mumuMass*mumuMass) > " << q2Bins->operator[](i) << " && (mumuMass*mumuMass) <= " << q2Bins->operator[](i+1);
    if (!useToyDataset) {
      cout << "\n[ExtractYield::IterativeAnglesFitq2Bins]\tCut string: " << myString.str() << endl;
      dataSet_q2Bins[i] = (RooDataSet*)dataSet->reduce(myString.str().c_str());
    } else {
      dataSet_q2Bins[i] = (RooDataSet*)dataSet;
    }
    cout << "[ExtractYield::IterativeAnglesFitq2Bins]\tNumber of events : " << dataSet_q2Bins[i]->sumEntries() << endl;

    if ((FitType == 206))
    {
      TCanvas*    Gen[q2Bins->size()-1];
      Gen[i] = new TCanvas("Gen","Gen",10,10,700,500);
      Gen[i]->Divide(2,2);

      myString.clear(); myString.str("");
      myString << "GenTotalPDFq2Bin_" << i;
      InstantiateGen3AnglesFit(&TotalPDFq2Bins[i],useEffPDF,y,z,p,FitType,configParam,i);
      // #####################
      // # Initialize p.d.f. #
      // #####################
      CopyFitResults(TotalPDFq2Bins[i],i,fitParam);

      // #####################
      // # Apply constraints #
      // #####################
      ClearVars(vecConstr);


      // ###################
      // # Perform the fit #
      // ##################
      fitResult = Make3AnglesFit(dataSet_q2Bins[i],&TotalPDFq2Bins[i],y,z,p,FitType,vecConstr,Gen[i],i);
      if (CheckGoodFit(fitResult) == true) cout << "\n[ExtractYield::Iterative3AnglesFitq2Bins]\t@@@ Fit converged ! @@@" << endl;
      else                                 cout << "\n[ExtractYield::Iterative3AnglesFitq2Bins]\t@@@ Fit didn't converge ! @@@" << endl;

    }

    if ((FitType == 106))
    {
      TCanvas*    Reco[q2Bins->size()-1];
      Reco[i] = new TCanvas("Reco","Reco",10,10,700,500);
      Reco[i]->Divide(2,2);

      unsigned int countMisTag  = 0;
      unsigned int countGoodTag = 0;
      for (int j = 0; j < static_cast<int>(dataSet_q2Bins[i]->sumEntries()); j++)
      {
        if (!useToyDataset && dataSet_q2Bins[i]->get(j)->getRealValue("truthMatchSignal") == true)
        {
          if (dataSet_q2Bins[i]->get(j)->getRealValue("rightFlavorTag") == 0.0) countMisTag++;
          else                                                                  countGoodTag++;
        }
      }
      cout << "[ExtractYield::IterativeMassFitq2Bins]\tDynamic mis-tag fraction : " << static_cast<double>(countMisTag) / static_cast<double>(countMisTag + countGoodTag) << " = (" << countMisTag << "/(" << countMisTag << "+" << countGoodTag << "))" << endl;

      myString.clear(); myString.str("");
      myString << "GenTotalPDFq2Bin_" << i;
      InstantiateReco3AnglesFit(&TotalPDFq2Bins[i],useEffPDF,y,z,p,FitType,configParam,i);

      // #####################
      // # Initialize p.d.f. #
      // #####################
      CopyFitResults(TotalPDFq2Bins[i],i,fitParam,countMisTag,countGoodTag);
      // #####################
      // # Apply constraints #
      // #####################
      ClearVars(vecConstr);
      // ###################
      // # Perform the fit #
      // ##################
      fitResult = Make3AnglesFit(dataSet_q2Bins[i],&TotalPDFq2Bins[i],y,z,p,FitType,vecConstr,Reco[i],i);
      if (CheckGoodFit(fitResult) == true) {
        cout << "\n[ExtractYield::Iterative3AnglesFitq2Bins]\t@@@ Fit converged ! @@@" << endl;
      }
      else
        cout << "\n[ExtractYield::Iterative3AnglesFitq2Bins]\t@@@ Fit didn't converge ! @@@" << endl;
    }

    // Open new ROOT file save save result 
    //if (fitResult && CheckGoodFit(fitResult) == true) fitResult->Write(Form("fitResult_Bin%d",i)) ;
    if (fitResult) {
      f.cd();
      fitResult->Write(Form("fitResult_Bin%d",i)) ;
    }
  }
  f.Close() ;
}




// ==================
// ===> 4D MODEL <===
// ==================

void InstantiateMass3AnglesFit (RooAbsPdf** TotalPDF,
                                bool useEffPDF,
                                RooRealVar* x, RooRealVar* y, RooRealVar* z,RooRealVar* p,
                                string fitName, unsigned int FitType,
                                vector<vector<unsigned int>*>* configParam,
                                vector<vector<string>*>* fitParam,
                                unsigned int q2BinIndx)
// #########################
// # x: mass               #
// # y: angle cos(theta_l) #
// # z: angle cos(theta_K) #
// # p: angle phi          #
// #########################
{
  // ################################
  // # Read configuration variables #
  // ################################
  unsigned int useSignal = configParam->operator[](Utility->GetConfigParamIndx("SigType"))->operator[](q2BinIndx);
  unsigned int usePeakB  = configParam->operator[](Utility->GetConfigParamIndx("PeakBkgType"))->operator[](q2BinIndx);
  unsigned int useCombB  = configParam->operator[](Utility->GetConfigParamIndx("CombBkgType"))->operator[](q2BinIndx);
  unsigned int useMisTag = configParam->operator[](Utility->GetConfigParamIndx("MistagType"))->operator[](q2BinIndx);


  // ###########################
  // # Read polynomial degrees #
  // ###########################
  unsigned int NCoeffPolyBKGpeak1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP1"))->operator[](q2BinIndx).c_str());
  unsigned int NCoeffPolyBKGcomb1 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC1"))->operator[](q2BinIndx).c_str());
  unsigned int NCoeffPolyBKGpeak2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyP2"))->operator[](q2BinIndx).c_str());
  unsigned int NCoeffPolyBKGcomb2 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC2"))->operator[](q2BinIndx).c_str());
  unsigned int NCoeffPolyBKGcomb3 = atoi(fitParam->operator[](Utility->GetFitParamIndx("nPolyC3"))->operator[](q2BinIndx).c_str());

  // ###################
  // # Local variables #
  // ###################
  stringstream myString;
  RooArgSet VarsC1, VarsC2, VarsC3;
  RooArgSet VarsP1, VarsP2;
  RooArgSet* VarsAng = new RooArgSet("VarsAng");


  // ################################################
  // # Define mass fit variables and pdf for signal #
  // ################################################
  meanS   = new RooRealVar("meanS","Signal mean",Utility->B0Mass,"GeV");

  sigmaS1 = new RooRealVar("sigmaS1","Signal sigma-1",0.0,"GeV");
  MassS1  = new RooGaussian("MassS1","Signal Gaussian-1",*x,*meanS,*sigmaS1);

  sigmaS2 = new RooRealVar("sigmaS2","Signal sigma-2",0.0,"GeV");
  MassS2  = new RooGaussian("MassS2","Signal Gaussian-2",*x,*meanS,*sigmaS2);

  fracMassS = new RooRealVar("fracMassS","Fraction of signal Gaussian",0.0,0.0,1.0);

  if      (useSignal == 1) MassSignal = new RooGaussian(*((RooGaussian*)MassS1),"MassSignal");
  else if (useSignal == 2) MassSignal = new RooAddPdf("MassSignal","Signal mass pdf",RooArgSet(*MassS1,*MassS2),RooArgSet(*fracMassS));
  else                     MassSignal = NULL;


  // ##################################################################
  // # Define angle fit variables and pdf for correctly tagged signal #
  // ##################################################################
  AngleS= MakeAngWithEffPDF(q2BinIndx,y,z,p,FitType,useEffPDF,VarsAng);

  Signal = new RooProdPdf("Signal","Signal Mass*Angle",RooArgSet(*MassSignal,*AngleS));


  // ###################################################################
  // # Define angle fit variables and pdf for combinatorial background #
  // ###################################################################
  for (unsigned int i = 0; i < NCoeffPolyBKGcomb1; i++)
  {
    myString.clear(); myString.str("");
    myString << "c1Poly" << i;
    c1Poly[i] = new RooRealVar(myString.str().c_str(),"Comb. bkg poly. coef.",0.0);
    VarsC1.add(*c1Poly[i]);
  }
  for (unsigned int i = 0; i < NCoeffPolyBKGcomb2; i++)
  {
    myString.clear(); myString.str("");
    myString << "c2Poly" << i;
    c2Poly[i] = new RooRealVar(myString.str().c_str(),"Comb. bkg poly. coef.",0.0);
    VarsC2.add(*c2Poly[i]);
  }
  for (unsigned int i = 0; i < NCoeffPolyBKGcomb3; i++)
  {
    myString.clear(); myString.str("");
    myString << "c3Poly" << i;
    c3Poly[i] = new RooRealVar(myString.str().c_str(),"Comb. bkg poly. coef.",0.0);
    VarsC3.add(*c3Poly[i]);
  }
  myString.clear(); myString.str("");
  myString << GeneratePolynomial(y,NCoeffPolyBKGcomb1,"c1Poly");
  BkgAngleC1 = new RooGenericPdf("BkgAngleC1",myString.str().c_str(),RooArgSet(*y,VarsC1));
  myString.clear(); myString.str("");
  myString << GeneratePolynomial(z,NCoeffPolyBKGcomb2,"c2Poly");
  BkgAngleC2 = new RooGenericPdf("BkgAngleC2",myString.str().c_str(),RooArgSet(*z,VarsC2));
  myString.clear(); myString.str("");
  myString << GeneratePolynomial(p,NCoeffPolyBKGcomb3,"c3Poly");
  BkgAngleC3 = new RooGenericPdf("BkgAngleC3",myString.str().c_str(),RooArgSet(*p,VarsC3));

  BkgAnglesC = new RooProdPdf("BkgAnglesC","Background Angle1*Angle2*Angle3",RooArgSet(*BkgAngleC1,*BkgAngleC2,*BkgAngleC3));


  // ##################################################################
  // # Define mass fit variables and pdf for combinatorial background #
  // ##################################################################
  var1 = new RooRealVar("var1","First background variables",0.0,"GeV");
  myString.clear(); myString.str("");
  myString << "exp(-(" << x->getPlotLabel() << " - meanS) / var1)";
  BkgMassExp1 = new RooGenericPdf("BkgMassExp1",myString.str().c_str(),RooArgSet(*x,*meanS,*var1));

  var2 = new RooRealVar("var2","Second background variable",0.0,"GeV");
  myString.clear(); myString.str("");
  myString << "exp(-(" << x->getPlotLabel() << " - meanS) / var2)";
  BkgMassExp2 = new RooGenericPdf("BkgMassExp2",myString.str().c_str(),RooArgSet(*x,*meanS,*var2));

  fracMassBExp = new RooRealVar("fracMassBExp","Fraction of background Exponential",0.0,0.0,1.0);

  if      (useCombB == 1) BkgMassComb = new RooGenericPdf(*((RooGenericPdf*)BkgMassExp1),"BkgMassComb");
  else if (useCombB == 2) BkgMassComb = new RooAddPdf("BkgMassComb","Background mass comb. bkg pdf",RooArgSet(*BkgMassExp1,*BkgMassExp2),RooArgSet(*fracMassBExp));
  else if (useCombB == 3)
  {
    myString.clear(); myString.str("");
    myString << "TMath::Erfc((" << x->getPlotLabel() << " - var1) / sqrt(2*var2*var2))";
    BkgMassComb = new RooGenericPdf("BkgMass",myString.str().c_str(),RooArgSet(*x,*var1,*var2));
  }
  else                    BkgMassComb = NULL;

  BkgMassAngleComb = new RooProdPdf("BkgMassAngleComb","Combinatorial bkg Mass*Angle",RooArgSet(*BkgMassComb,*BkgAnglesC));


  // #######################################################
  // # Define mass fit variables and pdf for mistag signal #
  // #######################################################
  sigmaMisTag1 = new RooRealVar("sigmaMisTag1","Mistag sigma-1",0.0,"GeV");
  MassMisTag1  = new RooGaussian("MassMisTag1","Mistag-1",*x,*meanS,*sigmaMisTag1);

  sigmaMisTag2 = new RooRealVar("sigmaMisTag2","Mistag sigma-2",0.0,"GeV");
  MassMisTag2  = new RooGaussian("MassMisTag2","Mistag-2",*x,*meanS,*sigmaMisTag2);

  fracMisTag = new RooRealVar("fracMisTag","Fraction mistag Gaussian",0.0,0.0,1.0);

  if      (useMisTag == 1) MassMisTag = new RooGaussian(*((RooGaussian*)MassMisTag1),"MassMisTag");
  else if (useMisTag == 2) MassMisTag = new RooAddPdf("MassMisTag","Mistag mass pdf",RooArgSet(*MassMisTag1,*MassMisTag2),RooArgSet(*fracMisTag));
  else                     MassMisTag = NULL;


  // ############################################################
  // # Define angle fit variables and pdf for mis-tagged signal #
  // ############################################################
  AngleMisTag = MakeAngWithEffPDF(q2BinIndx,y,z,p,FitType*10,useEffPDF,VarsAng);

  MassAngleMisTag = new RooProdPdf("MassAngleMisTag","Mistag bkg Mass*Angle",RooArgSet(*MassMisTag,*AngleMisTag));


  // #############################################################
  // # Define angle fit variables and pdf for peaking background #
  // #############################################################
  for (unsigned int i = 0; i < NCoeffPolyBKGpeak1; i++)
  {
    myString.clear(); myString.str("");
    myString << "p1Poly" << i;
    p1Poly[i] = new RooRealVar(myString.str().c_str(),"Peak. bkg poly. coef.",0.0);
    VarsP1.add(*p1Poly[i]);
  }
  for (unsigned int i = 0; i < NCoeffPolyBKGpeak2; i++)
  {
    myString.clear(); myString.str("");
    myString << "p2Poly" << i;
    p2Poly[i] = new RooRealVar(myString.str().c_str(),"Peak. bkg poly. coef.",0.0);
    VarsP2.add(*p2Poly[i]);
  }
  myString.clear(); myString.str("");
  myString << GeneratePolynomial(y,NCoeffPolyBKGpeak1,"p1Poly");
  BkgAngleP1 = new RooGenericPdf("BkgAngleP1",myString.str().c_str(),RooArgSet(*y,VarsP1));
  myString.clear(); myString.str("");
  myString << GeneratePolynomial(z,NCoeffPolyBKGpeak2,"p2Poly");
  BkgAngleP2 = new RooGenericPdf("BkgAngleP2",myString.str().c_str(),RooArgSet(*z,VarsP2));
  BkgAnglesP = new RooProdPdf("BkgAnglesP","Background Angle1*Angle2",RooArgSet(*BkgAngleP1,*BkgAngleP2));


  // ############################################################
  // # Define mass fit variables and pdf for peaking background #
  // ############################################################
  meanR1         = new RooRealVar("meanR1","Bkg right peak mean-1",0.0,"GeV");
  sigmaR1        = new RooRealVar("sigmaR1","Bkg right peak sigma-1",0.0,"GeV");
  BkgMassRPeak1  = new RooGaussian("BkgMassRPeak1","Bkg right peak-1",*x,*meanR1,*sigmaR1);

  meanR2         = new RooRealVar("meanR2","Bkg right peak mean-2",0.0,"GeV");
  sigmaR2        = new RooRealVar("sigmaR2","Bkg right peak sigma-2",0.0,"GeV");
  BkgMassRPeak2  = new RooGaussian("BkgMassRPeak2","Bkg right peak-2",*x,*meanR2,*sigmaR2);

  fracMassBRPeak = new RooRealVar("fracMassBRPeak","Fraction of background right Peak",0.0,0.0,1.0);
  BkgMassRPeak   = new RooAddPdf("BkgMassRPeak","Right peaking bkg mass pdf",RooArgSet(*BkgMassRPeak1,*BkgMassRPeak2),RooArgSet(*fracMassBRPeak));

  meanL1         = new RooRealVar("meanL1","Bkg left peak mean-1",0.0,"GeV");
  sigmaL1        = new RooRealVar("sigmaL1","Bkg left peak sigma-1",0.0,"GeV");
  BkgMassLPeak1  = new RooGaussian("BkgMassLPeak1","Bkg left peak-1",*x,*meanL1,*sigmaL1);

  meanL2         = new RooRealVar("meanL2","Bkg left peak mean-2",0.0,"GeV");
  sigmaL2        = new RooRealVar("sigmaL2","Bkg left peak sigma-2",0.0,"GeV");
  BkgMassLPeak2  = new RooGaussian("BkgMassLPeak2","Bkg left peak-2",*x,*meanL2,*sigmaL2);

  fracMassBLPeak = new RooRealVar("fracMassBLPeak","Fraction of background left Peak",0.0,0.0,1.0);
  BkgMassLPeak   = new RooAddPdf("BkgMassLPeak","Left peaking bkg mass pdf",RooArgSet(*BkgMassLPeak1,*BkgMassLPeak2),RooArgSet(*fracMassBLPeak));

  fracMassBPeak  = new RooRealVar("fracMassBPeak","Fraction of background right-left Peak",0.0,0.0,1.0);

  if      (usePeakB == 1)  BkgMassPeak = new RooGaussian(*((RooGaussian*)BkgMassRPeak1),"BkgMassPeak");
  else if (usePeakB == 2)  BkgMassPeak = new RooAddPdf(*((RooAddPdf*)BkgMassRPeak),"BkgMassPeak");
  else if (usePeakB == 11) BkgMassPeak = new RooAddPdf("BkgMassPeak","Peaking bkg mass pdf",RooArgSet(*BkgMassRPeak1,*BkgMassLPeak1),RooArgSet(*fracMassBPeak));
  else if (usePeakB == 12) BkgMassPeak = new RooAddPdf("BkgMassPeak","Peaking bkg mass pdf",RooArgSet(*BkgMassRPeak,*BkgMassLPeak),RooArgSet(*fracMassBPeak));
  else                     BkgMassPeak = NULL;

  BkgMassAnglePeak = new RooProdPdf("BkgMassAnglePeak","Peaking bkg Mass*Angle",RooArgSet(*BkgMassPeak,*BkgAnglesP));


  // ###########################
  // # Define pdf coefficients #
  // ###########################
  nSig     = new RooRealVar("nSig","Fraction of signal events",1.0);
  nBkgComb = new RooRealVar("nBkgComb","Number of comb. background events",1.0);
  RooFormulaVar* nMisTag;
  nMisTagFrac = new RooRealVar("nMisTagFrac","Fraction of mistag",0.0,0.0,1.0);
  nMisTag     = new RooFormulaVar("nMisTag","nSig * nMisTagFrac / (1 - nMisTagFrac)", RooArgSet(*nSig,*nMisTagFrac));
  SignalT     = new RooAddPdf("SignalT","Signal + nMisTagFrac/(1-nMisTagFrac) * (MassAngleMisTag)",RooArgSet(*MassAngleMisTag,*Signal),*nMisTagFrac) ;

  nBkgPeak = new RooRealVar("nBkgPeak","Number of peaking background events",0.0,0.0,1.0);


  if ((useSignal != 0) && (usePeakB == 0) && (useCombB == 0) && (useMisTag == 0))
    *TotalPDF = new RooProdPdf("Signal","Signal Mass*Angle",RooArgSet(*MassSignal,*AngleS));
    
  else if ((useSignal == 0) && (usePeakB != 0) && (useCombB == 0) && (useMisTag == 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total pdf",RooArgSet(*BkgMassAnglePeak),RooArgSet(*nBkgPeak));
    
  else if ((useSignal != 0) && (usePeakB == 0) && (useCombB != 0) && (useMisTag == 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total pdf",RooArgSet(*Signal,*BkgMassAngleComb),*nSig);

  else if ((useSignal != 0) && (usePeakB == 0) && (useCombB == 0) && (useMisTag != 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total pdf",RooArgSet(*Signal,*MassAngleMisTag),*nMisTagFrac);

  else if ((useSignal != 0) && (usePeakB != 0) && (useCombB != 0) && (useMisTag == 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total pdf",RooArgSet(*Signal,*BkgMassAngleComb,*BkgMassAnglePeak),RooArgList(*nSig,*nBkgPeak));

  else if ((useSignal != 0) && (usePeakB == 0) && (useCombB != 0) && (useMisTag != 0)) {
    // *TotalPDF = new RooAddPdf(fitName.c_str(),"Total pdf",RooArgSet(*SignalT,*BkgMassAngleComb),*nSig);
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total extended pdf",RooArgSet(*Signal,*MassAngleMisTag,*BkgMassAngleComb),RooArgSet(*nSig,*nMisTag,*nBkgComb));
  }

  else if ((useSignal != 0) && (usePeakB != 0) && (useCombB != 0) && (useMisTag != 0))
    *TotalPDF = new RooAddPdf(fitName.c_str(),"Total pdf",RooArgSet(*SignalT,*BkgMassAngleComb,*BkgMassAnglePeak),RooArgList(*nSig,*nMisTag,*nBkgPeak));

  else
  {
    cout << "[ExtractYield::InstantiateMass3AnglesFit]\tIncorrect configuration sequence : useSignal = " << useSignal << "\tusePeakB = " << usePeakB << "\tuseCombB = " << useCombB << "\tuseMisTag = " << useMisTag << endl;
    exit (EXIT_FAILURE);
  }
}

//#################
//# DATA fitting  #
//#################

RooFitResult* MakeMass3AnglesFit (RooDataSet* dataSet, RooAbsPdf** TotalPDF, RooRealVar* x, RooRealVar* y, RooRealVar* z,RooRealVar* p,  unsigned int FitType, RooArgSet* vecConstr,  TPaveText* extText, unsigned int ID)
{
  // ###################
  // # Local variables #
  // ###################
  RooFitResult* fitResult = NULL;
  stringstream myString;
  unsigned int nElements   = 0;
  unsigned int it          = 0;
  TString legNames[6];

  const unsigned int nCanv = 4;
  TCanvas* localCanv[nCanv];
  for (unsigned int i = 0; i < nCanv; i++)
  {
    myString.clear(); myString.str("");
    myString << "localCanv" << i;
    localCanv[i] = new TCanvas(myString.str().c_str(),myString.str().c_str(),20,20,500,500);
  }

  TLegend*   legX             = NULL;
  TLegend*   legY             = NULL;
  TLegend*   legZ             = NULL;
  TLegend*   legP             = NULL;

  if ((FitType == 1) ||(FitType == 6))
  {
    // ##################
    // # Make sideband  #
    // ##################
    if (GetVar(*TotalPDF,"nSig") != NULL)
    {
      cout << "[ExtractYield::MakeMass3AnglesFit]\t@@@ Making comb. angular background sideband  @@@" << endl;
      RooDataSet* sideBands = NULL;
      RooArgSet constrSidebads;
      // #############
      // # Sidebands #
      // #############
      myString.clear(); myString.str("");
      myString << "B0MassArb < " << GetVar(*TotalPDF,"meanS")->getVal() - atof(Utility->GetGenericParam("NSigmaB0").c_str())*Utility->GetB0Width();
      myString << " || B0MassArb > " << GetVar(*TotalPDF,"meanS")->getVal() + atof(Utility->GetGenericParam("NSigmaB0").c_str())*Utility->GetB0Width();
      cout << "[ExtractYield::MakeMass3AnglesFit]\tCut for B0 sidebands : " << myString.str().c_str() << endl;
      sideBands = (RooDataSet*)dataSet->reduce(myString.str().c_str());
      delete sideBands;
      ClearVars(&constrSidebads);
    }
    // ###################
    // # Make actual fit #
    // ###################

    cout << "About to fit " << Utility->GetGenericParam("ApplyConstr").c_str() << " " << Utility->GetGenericParam("UseMINOS").c_str() << " " << MINIMIZER << endl;
    if (atoi(Utility->GetGenericParam("ApplyConstr").c_str()) == true) fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),ExternalConstraints(*vecConstr),Save(true),Minos(atoi(Utility->GetGenericParam("UseMINOS").c_str())),Minimizer(MINIMIZER));
    else                                                               fitResult = (*TotalPDF)->fitTo(*dataSet,Extended(true),Save(true),Minos(atoi(Utility->GetGenericParam("UseMINOS").c_str())),Minimizer(MINIMIZER));
    
    // RooAbsReal* nll = (*TotalPDF)->createNLL(*dataSet, Extended(true));
    
    // // Manual minimization
    // RooMinimizer minimizer(*nll);
    
    // minimizer.setVerbose(false);
    // minimizer.setMinimizerType("Minuit2");
    // // minimizer.optimizeConst(true);
    // // minimizer.setEps(1e-8)
    // // minimizer.setStrategy(2);
    // // minimizer.setPrintEvalErrors(1);
    // // minimizer.setOffsetting(True);
    // minimizer.migrad();
    // minimizer.hesse();
    // GetVar(*TotalPDF,"As5S")->setError(0.1);
    // GetVar(*TotalPDF,"P1S" )->setError(0.3);
    // GetVar(*TotalPDF,"P5pS")->setError(0.3);
    // minimizer.minos();
    // fitResult = minimizer.save();

    // ###################################################
    // # Set p.d.f. independent variables to known point #
    // ###################################################
    if (GetVar(*TotalPDF,x->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(x->getPlotLabel(),Utility->B0Mass);
    if (GetVar(*TotalPDF,y->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(y->getPlotLabel(),0.0);
    if (GetVar(*TotalPDF,z->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(z->getPlotLabel(),0.0);
    if (GetVar(*TotalPDF,p->getPlotLabel()) != NULL) (*TotalPDF)->getVariables()->setRealValue(p->getPlotLabel(),0.0);
    if (fitResult != NULL) fitResult->Print("v");

    // #####################
    // # Mass plot results #
    // #####################
    if (PLOT)
    {
      cout << "[ExtractYield::MakeMass3AnglesFit]\tPlot mass" << endl;
      TCanvas* Canv = new TCanvas("Canv","Canv",1200,900);
      Canv->Divide(2,2);
      Canv->cd(1);
      RooPlot* myFrameX = x->frame(NBINS);

      dataSet->plotOn(myFrameX, Name(MakeName(dataSet,ID).c_str()));
      legNames[nElements++] = "Data";

      if (FUNCERRBAND == true)
      {
        (*TotalPDF)->plotOn(myFrameX, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), VisualizeError(*fitResult,1,true), FillColor(kGreen-7), Project(RooArgSet(*y,*z,*p)), VLines());
        dataSet->plotOn(myFrameX, Name(MakeName(dataSet,ID).c_str()));
      }
      else (*TotalPDF)->plotOn(myFrameX, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), Project(RooArgSet(*y,*z,*p)));
      legNames[nElements++] = "Total p.d.f.";

      if (GetVar(*TotalPDF,"nSig") != NULL)
      {
        (*TotalPDF)->plotOn(myFrameX, Components(*Signal), FillStyle(3345), FillColor(kBlue), Project(RooArgSet(*y,*z,*p)), DrawOption("F"));
        (*TotalPDF)->plotOn(myFrameX, Components(*Signal), LineStyle(7),    LineColor(kBlue), Project(RooArgSet(*y,*z,*p)), DrawOption("L"));
        legNames[nElements++] = "Right-tag sig";

        (*TotalPDF)->plotOn(myFrameX, Components(*BkgMassAngleComb), LineStyle(4), LineColor(kRed), Project(RooArgSet(*y,*z,*p)));
        legNames[nElements++] = "Comb. bkg";
      }

      if (GetVar(*TotalPDF,"nMisTagFrac") != NULL)
      {
        (*TotalPDF)->plotOn(myFrameX, Components(*MassAngleMisTag), FillStyle(3354), FillColor(kGreen+1), Project(RooArgSet(*y,*z,*p)), DrawOption("F"));
        (*TotalPDF)->plotOn(myFrameX, Components(*MassAngleMisTag), LineStyle(8),    LineColor(kGreen+1), Project(RooArgSet(*y,*z,*p)), DrawOption("L"));
        legNames[nElements++] = "Mis-tag sig";
      }

      if (GetVar(*TotalPDF,"nBkgPeak") != NULL)
      {
        (*TotalPDF)->plotOn(myFrameX, Components(*BkgMassAnglePeak), LineStyle(3), LineColor(kViolet), Project(RooArgSet(*y,*z,*p)));
        legNames[nElements++] = "Peak. bkg";
      }

      //  (*TotalPDF)->paramOn(myFrameX,Format("NEU",AutoPrecision(1)),Layout(0.11,0.42,0.88),ShowConstants(true));
      TPaveText* paveTextX =new TPaveText(0.12,0.82,0.4,0.86,"NDC"); 
      paveTextX->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameX->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
      CheckGoodFit(fitResult,paveTextX);
      paveTextX->SetTextAlign(11);
      paveTextX->SetBorderSize(0.0);
      paveTextX->SetFillStyle(0);
      paveTextX->SetTextSize(0.03); 
      DrawString(LUMI,myFrameX);
      paveTextX->Paint();
      myFrameX->addObject(paveTextX);
      legX = new TLegend(0.78, 0.65, 0.97, 0.88, "");
      it = 0;
      for (unsigned int i = 0; i < 2*nElements; i++)
      {
        TString objName = myFrameX->nameOf(i);
        if ((objName == "") || (objName.Contains("paramBox") == true) || (objName.Contains("TPave") == true) || ((i > 0) && (objName == myFrameX->nameOf(i-1)))) continue;
        TObject* obj = myFrameX->findObject(objName.Data());
        legX->AddEntry(obj,legNames[it++],"PL");
        legX->SetTextFont(42);
      }
      legX->SetFillStyle(0);
      legX->SetFillColor(0);
      legX->SetTextSize(0.04);
      legX->SetBorderSize(0);
      myFrameX->Draw();
      legX->Draw("same");
      // ##########################
      // # Angular-1 plot results #
      // ##########################
      cout << "[ExtractYield::MakeMass3AnglesFit]\tPlot angle 1" << endl;
      Canv->cd(2);
      RooPlot* myFrameY = y->frame(NBINS);

      dataSet->plotOn(myFrameY, Name(MakeName(dataSet,ID).c_str()));

      if (FUNCERRBAND == true)
      {
        (*TotalPDF)->plotOn(myFrameY, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), VisualizeError(*fitResult,1,true), FillColor(kGreen-7), Project(RooArgSet(*x,*z,*p)), VLines());
        dataSet->plotOn(myFrameY, Name(MakeName(dataSet,ID).c_str()));
      }
      else (*TotalPDF)->plotOn(myFrameY, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), Project(RooArgSet(*x,*z,*p)));

      if (GetVar(*TotalPDF,"nSig")        != NULL)
      {
        (*TotalPDF)->plotOn(myFrameY, Components(*Signal), FillStyle(3345), FillColor(kBlue), Project(RooArgSet(*x,*z,*p)), DrawOption("F"));
        (*TotalPDF)->plotOn(myFrameY, Components(*Signal), LineStyle(7),    LineColor(kBlue), Project(RooArgSet(*x,*z,*p)), DrawOption("L"));
        (*TotalPDF)->plotOn(myFrameY, Components(*BkgMassAngleComb), LineStyle(4), LineColor(kRed),     Project(RooArgSet(*x,*z,*p)));
      }
      if (GetVar(*TotalPDF,"nMisTagFrac") != NULL)
      {
        (*TotalPDF)->plotOn(myFrameY, Components(*MassAngleMisTag), FillStyle(3354), FillColor(kGreen+1), Project(RooArgSet(*x,*z,*p)), DrawOption("F"));
        (*TotalPDF)->plotOn(myFrameY, Components(*MassAngleMisTag), LineStyle(8),    LineColor(kGreen+1), Project(RooArgSet(*x,*z,*p)), DrawOption("L"));
      }

      if (GetVar(*TotalPDF,"nBkgPeak")    != NULL) (*TotalPDF)->plotOn(myFrameY, Components(*BkgMassAnglePeak), LineStyle(3), LineColor(kViolet),  Project(RooArgSet(*x,*z,*p)));

      TPaveText* paveTextY = new TPaveText(0.12,0.82,0.4,0.86,"NDC");
      paveTextY->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameY->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
      paveTextY->SetTextAlign(11);
      paveTextY->SetBorderSize(0.0);
      paveTextY->SetFillStyle(0);
      paveTextY->SetTextSize(0.04);
      paveTextY->Paint();
      myFrameY->addObject(paveTextY);
      DrawString(LUMI,myFrameY);
      if ((extText != NULL) && (SETBATCH == false))
      {
        extText->Paint();
        myFrameY->addObject(extText);
      }
      myFrameY->Draw();

      legY = new TLegend(0.78, 0.65, 0.97, 0.88, "");
      it = 0;
      for (unsigned int i = 0; i < 2*nElements; i++)
      {
        TString objName = myFrameY->nameOf(i);
        if ((objName == "") || (objName.Contains("paramBox") == true) || (objName.Contains("TPave") == true) || ((i > 0) && (objName == myFrameY->nameOf(i-1)))) continue;
        TObject* obj = myFrameY->findObject(objName.Data());
        legY->AddEntry(obj,legNames[it++],"PL");
        legY->SetTextFont(42);
      }
      legY->SetFillStyle(0);
      legY->SetFillColor(0);
      legY->SetTextSize(0.04);
      legY->SetBorderSize(0);
      legY->Draw("same");


      // ##########################
      // # Angular-2 plot results #
      // ##########################
      cout << "[ExtractYield::MakeMass3AnglesFit]\tPlot angle 2" << endl;
      Canv->cd(3);
      RooPlot* myFrameZ = z->frame(NBINS);

      dataSet->plotOn(myFrameZ, Name(MakeName(dataSet,ID).c_str()));

      if (FUNCERRBAND == true)
      {
        (*TotalPDF)->plotOn(myFrameZ, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), VisualizeError(*fitResult,1,true), FillColor(kGreen-7), Project(RooArgSet(*x,*y,*p)), VLines());
        dataSet->plotOn(myFrameZ, Name(MakeName(dataSet,ID).c_str()));
      }
      else (*TotalPDF)->plotOn(myFrameZ, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), Project(RooArgSet(*x,*y,*p)));

      if (GetVar(*TotalPDF,"nSig")        != NULL)
      {
        (*TotalPDF)->plotOn(myFrameZ, Components(*Signal), FillStyle(3345), FillColor(kBlue), Project(RooArgSet(*x,*y,*p)), DrawOption("F"));
        (*TotalPDF)->plotOn(myFrameZ, Components(*Signal), LineStyle(7),    LineColor(kBlue), Project(RooArgSet(*x,*y,*p)), DrawOption("L"));
        (*TotalPDF)->plotOn(myFrameZ, Components(*BkgMassAngleComb), LineStyle(4), LineColor(kRed),     Project(RooArgSet(*x,*y,*p)));
      }
      if (GetVar(*TotalPDF,"nMisTagFrac") != NULL)
      {
        (*TotalPDF)->plotOn(myFrameZ, Components(*MassAngleMisTag), FillStyle(3354), FillColor(kGreen+1), Project(RooArgSet(*x,*y,*p)), DrawOption("F"));
        (*TotalPDF)->plotOn(myFrameZ, Components(*MassAngleMisTag),  LineStyle(8),   LineColor(kGreen+1), Project(RooArgSet(*x,*y,*p)), DrawOption("L"));
      }
      if (GetVar(*TotalPDF,"nBkgPeak")    != NULL) (*TotalPDF)->plotOn(myFrameZ, Components(*BkgMassAnglePeak), LineStyle(3), LineColor(kViolet),  Project(RooArgSet(*x,*y,*p)));


      TPaveText* paveTextZ = new TPaveText(0.12,0.82,0.4,0.86,"NDC");
      paveTextZ->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameZ->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
      paveTextZ->SetTextAlign(11);
      paveTextZ->SetBorderSize(0.0);
      paveTextZ->SetFillStyle(0);
      paveTextZ->SetTextSize(0.04);
      paveTextZ->Paint();
      myFrameZ->addObject(paveTextZ);
      DrawString(LUMI,myFrameZ);
      if ((extText != NULL) && (SETBATCH == false))
      {
        extText->Paint();
        myFrameZ->addObject(extText);
      }
      myFrameZ->Draw();

      legZ = new TLegend(0.78, 0.65, 0.97, 0.88, "");
      it = 0;
      for (unsigned int i = 0; i < 2*nElements; i++)
      {
        TString objName = myFrameZ->nameOf(i);
        if ((objName == "") || (objName.Contains("paramBox") == true) || (objName.Contains("TPave") == true) || ((i > 0) && (objName == myFrameZ->nameOf(i-1)))) continue;
        TObject* obj = myFrameZ->findObject(objName.Data());
        legZ->AddEntry(obj,legNames[it++],"PL");
        legZ->SetTextFont(42);
      }
      legZ->SetFillStyle(0);
      legZ->SetFillColor(0);
      legZ->SetTextSize(0.04);
      legZ->SetBorderSize(0);
      legZ->Draw("same");

      // ##########################
      // # Angular-3 plot results #
      // ##########################
      cout << "[ExtractYield::MakeMass3AnglesFit]\tPlot angle 3" << endl;
      Canv->cd(4);
      RooPlot* myFrameP = p->frame(NBINS);

      dataSet->plotOn(myFrameP, Name(MakeName(dataSet,ID).c_str()));

      if (FUNCERRBAND == true)
      {
        (*TotalPDF)->plotOn(myFrameP, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), VisualizeError(*fitResult,1,true), FillColor(kGreen-7), Project(RooArgSet(*x,*y,*z)), VLines());
        dataSet->plotOn(myFrameP, Name(MakeName(dataSet,ID).c_str()));
      }
      else (*TotalPDF)->plotOn(myFrameP, Name((*TotalPDF)->getPlotLabel()), LineColor(kBlack), Project(RooArgSet(*x,*y,*z)));

      if (GetVar(*TotalPDF,"nSig")        != NULL)
      {
        (*TotalPDF)->plotOn(myFrameP, Components(*Signal), FillStyle(3345), FillColor(kBlue), Project(RooArgSet(*x,*y,*z)), DrawOption("F"));
        (*TotalPDF)->plotOn(myFrameP, Components(*Signal), LineStyle(7),    LineColor(kBlue), Project(RooArgSet(*x,*y,*z)), DrawOption("L"));
        (*TotalPDF)->plotOn(myFrameP, Components(*BkgMassAngleComb), LineStyle(4), LineColor(kRed),     Project(RooArgSet(*x,*y,*z)));
      }
      if (GetVar(*TotalPDF,"nMisTagFrac") != NULL)
      {
        (*TotalPDF)->plotOn(myFrameP, Components(*MassAngleMisTag), FillStyle(3354), FillColor(kGreen+1), Project(RooArgSet(*x,*y,*z)), DrawOption("F"));
        (*TotalPDF)->plotOn(myFrameP, Components(*MassAngleMisTag),  LineStyle(8),   LineColor(kGreen+1), Project(RooArgSet(*x,*y,*z)), DrawOption("L"));
      }
      if (GetVar(*TotalPDF,"nBkgPeak")    != NULL) (*TotalPDF)->plotOn(myFrameP, Components(*BkgMassAnglePeak), LineStyle(3), LineColor(kViolet),  Project(RooArgSet(*x,*y,*z)));


      TPaveText* paveTextP = new TPaveText(0.12,0.82,0.4,0.86,"NDC");
      paveTextP->AddText(Form("%s%.2f","#chi#lower[0.4]{^{2}}/DoF = ",myFrameP->chiSquare((*TotalPDF)->getPlotLabel(),MakeName(dataSet,ID).c_str())));
      paveTextP->SetTextAlign(11);
      paveTextP->SetBorderSize(0.0);
      paveTextP->SetFillStyle(0);
      paveTextP->SetTextSize(0.04);
      paveTextP->Paint();
      myFrameP->addObject(paveTextP);
      DrawString(LUMI,myFrameP);
      if ((extText != NULL) && (SETBATCH == false))
      {
        extText->Paint();
        myFrameP->addObject(extText);
      }
      myFrameP->Draw();

      legP = new TLegend(0.78, 0.65, 0.97, 0.88, "");
      it = 0;
      for (unsigned int i = 0; i < 2*nElements; i++)
      {
        TString objName = myFrameP->nameOf(i);
        if ((objName == "") || (objName.Contains("paramBox") == true) || (objName.Contains("TPave") == true) || ((i > 0) && (objName == myFrameP->nameOf(i-1)))) continue;
        TObject* obj = myFrameP->findObject(objName.Data());
        legP->AddEntry(obj,legNames[it++],"PL");
        legP->SetTextFont(42);
      }
      legP->SetFillStyle(0);
      legP->SetFillColor(0);
      legP->SetTextSize(0.04);
      legP->SetBorderSize(0);
      legP->Draw("same");

      if (SETBATCH == true)
      {
        delete legX;
        delete legY;
        delete legZ;
        delete legP;
      }

      Canv->Modified();
      Canv->Update();
      if (SAVEPLOT == true)
	{
	  
	  myString.clear(); myString.str("");
	  myString << (*TotalPDF)->getPlotLabel() << "_Canv" << ID << ".pdf";
	  Canv->Print(myString.str().c_str());
	  
	  myString.clear(); myString.str("");
	  myString << (*TotalPDF)->getPlotLabel() << "_Canv" << ID << ".root";
	  Canv->Print(myString.str().c_str());
	}

    }
    // ###################
    // # Profiling scan  # (old style)
    // ###################
    if (PROFILENLL == true && false)
    {

      string varName;
      RooAbsReal* NLL;
      if (atoi(Utility->GetGenericParam("ApplyConstr").c_str()) == true) NLL = (*TotalPDF)->createNLL(*dataSet,Extended(true),ExternalConstraints(*vecConstr));
      else                                                               NLL = (*TotalPDF)->createNLL(*dataSet,Extended(true));
      
      // RooMinimizer RooMin(*NLL);
      // RooMin.minimize(...);
      
      
      //       RooMinuit(*NLL).migrad();
      //       RooMinuit RooMin(*NLL); 
      varName = "As5S";
      // GetVar(*TotalPDF,varName.c_str())->setRange(GetVar(*TotalPDF,varName.c_str())->getVal() + GetVar(*TotalPDF,varName.c_str())->getErrorLo() * atof(Utility->GetGenericParam("NSigmaB0").c_str()),
      //                                             GetVar(*TotalPDF,varName.c_str())->getVal() + GetVar(*TotalPDF,varName.c_str())->getErrorHi() * atof(Utility->GetGenericParam("NSigmaB0").c_str()));
      GetVar(*TotalPDF,varName.c_str())->setRange(GetVar(*TotalPDF,varName.c_str())->getVal()-0.1, GetVar(*TotalPDF,varName.c_str())->getVal()+0.1);
      localCanv[0]->cd();
      RooPlot* myFrameNLLVar1 = GetVar(*TotalPDF,varName.c_str())->frame();
      NLL->plotOn(myFrameNLLVar1, ShiftToZero());
      RooAbsReal* var1Profile = NLL->createProfile(*(GetVar(*TotalPDF,varName.c_str())));
      var1Profile->plotOn(myFrameNLLVar1, LineColor(kRed));
      myFrameNLLVar1->SetMinimum(0);
      myFrameNLLVar1->SetMaximum(5);
      DrawString(LUMI,myFrameNLLVar1);
      myFrameNLLVar1->DrawClone();

      varName = "P5pS";
      // GetVar(*TotalPDF,varName.c_str())->setRange(GetVar(*TotalPDF,varName.c_str())->getVal() + GetVar(*TotalPDF,varName.c_str())->getErrorLo() * atof(Utility->GetGenericParam("NSigmaB0").c_str()),
      //                                             GetVar(*TotalPDF,varName.c_str())->getVal() + GetVar(*TotalPDF,varName.c_str())->getErrorHi() * atof(Utility->GetGenericParam("NSigmaB0").c_str()));
      
      GetVar(*TotalPDF,varName.c_str())->setRange(GetVar(*TotalPDF,varName.c_str())->getVal()-0.3, GetVar(*TotalPDF,varName.c_str())->getVal()+0.3);
      localCanv[1]->cd();
      RooPlot* myFrameNLLVar2 = GetVar(*TotalPDF,varName.c_str())->frame();
      NLL->plotOn(myFrameNLLVar2, ShiftToZero());
      RooAbsReal* var2Profile = NLL->createProfile(*(GetVar(*TotalPDF,varName.c_str())));
      var2Profile->plotOn(myFrameNLLVar2, LineColor(kRed));
      myFrameNLLVar2->SetMinimum(0);
      myFrameNLLVar2->SetMaximum(5);
      DrawString(LUMI,myFrameNLLVar2);
      myFrameNLLVar2->DrawClone();

      varName = "P1S";
      // GetVar(*TotalPDF,varName.c_str())->setRange(GetVar(*TotalPDF,varName.c_str())->getVal() + GetVar(*TotalPDF,varName.c_str())->getErrorLo() * atof(Utility->GetGenericParam("NSigmaB0").c_str()),
      //                                             GetVar(*TotalPDF,varName.c_str())->getVal() + GetVar(*TotalPDF,varName.c_str())->getErrorHi() * atof(Utility->GetGenericParam("NSigmaB0").c_str()));
      GetVar(*TotalPDF,varName.c_str())->setRange(GetVar(*TotalPDF,varName.c_str())->getVal()-0.3, GetVar(*TotalPDF,varName.c_str())->getVal()+0.3);
      localCanv[2]->cd();
      RooPlot* myFrameNLLVar3 = GetVar(*TotalPDF,varName.c_str())->frame();
      NLL->plotOn(myFrameNLLVar3, ShiftToZero());
      RooAbsReal* var3Profile = NLL->createProfile(*(GetVar(*TotalPDF,varName.c_str())));
      var3Profile->plotOn(myFrameNLLVar3, LineColor(kRed));
      myFrameNLLVar3->SetMinimum(0);
      myFrameNLLVar3->SetMaximum(5);
      DrawString(LUMI,myFrameNLLVar3);
      myFrameNLLVar3->DrawClone();


      delete NLL;
    }

    // ##############
    // # Save plots #
    // ##############
    if (PROFILENLL) {
      for (unsigned int i = 0; i < nCanv; i++)
	{
	  myString.clear(); myString.str("");
	  myString << (*TotalPDF)->getPlotLabel() << "_localCanv" << i << "_" << ID << ".pdf";
	  localCanv[i]->Print(myString.str().c_str());
	  
	  myString.clear(); myString.str("");
	  myString << (*TotalPDF)->getPlotLabel() << "_localCanv" << i << "_" << ID << ".root";
	  localCanv[i]->Print(myString.str().c_str());
	
	}
    }
  }
  
  return fitResult;
}


void IterativeMass3AnglesFitq2Bins (RooDataSet* dataSet,
                                    bool useEffPDF,
                                    RooRealVar* x, RooRealVar* y, RooRealVar* z,RooRealVar* p,
                                    int specBin,
                                    unsigned int FitType,
                                    vector<double>* q2Bins,
                                    vector<vector<unsigned int>*>* configParam, vector<vector<string>*>* fitParam,
                                    RooArgSet* vecConstr,
                                    unsigned int ID)
{
  // ###################
  // # Local variables #
  // ###################
  stringstream myString;


  RooDataSet* dataSet_q2Bins[q2Bins->size()-1];
  RooAbsPdf*  TotalPDFq2Bins[q2Bins->size()-1];

  TPaveText*  extText[q2Bins->size()-1];

  RooFitResult* fitResult=0;


  for (unsigned int i = (specBin == -1 ? 0 : specBin); i < (specBin == -1 ? q2Bins->size()-1 : specBin+1); i++)
  {
    // ######################
    // # Make external text #
    // ######################
    extText[i] = new TPaveText(0.65,0.5,0.97,0.63,"NDC");
    extText[i]->AddText(Form("%s%.2f%s%.2f%s","q#lower[0.4]{^{2}}: ",q2Bins->operator[](i)," #font[122]{\55} ",q2Bins->operator[](i+1)," GeV#lower[0.4]{^{2}}"));
    extText[i]->SetTextAlign(11);
    extText[i]->SetBorderSize(0.0);
    extText[i]->SetFillStyle(0);
    extText[i]->SetTextSize(0.04);

    myString.clear(); myString.str("");
    myString << "(mumuMass*mumuMass) > " << q2Bins->operator[](i) << " && (mumuMass*mumuMass) <= " << q2Bins->operator[](i+1);
    if (!useToyDataset) {
      cout << "\n[ExtractYield::IterativeMass3AnglesFitq2Bins]\tCut string: " << myString.str() << endl;
      dataSet_q2Bins[i] = (RooDataSet*)dataSet->reduce(myString.str().c_str());
    } else {
      dataSet_q2Bins[i] = (RooDataSet*)dataSet;
    }
    cout << "[ExtractYield::IterativeMass3AnglesFitq2Bins]\tNumber of events : " << dataSet_q2Bins[i]->sumEntries() << endl;

    unsigned int countMisTag  = 0;
    unsigned int countGoodTag = 0;
    for (int j = 0; j < static_cast<int>(dataSet_q2Bins[i]->sumEntries()); j++)
    {
      if (dataSet_q2Bins[i]->get(j)->getRealValue("truthMatchSignal") == true)
      {
        if (dataSet_q2Bins[i]->get(j)->getRealValue("rightFlavorTag") == 0.0) countMisTag++;
        else                                                                  countGoodTag++;
      }
    }
    cout << "[ExtractYield::IterativeMassFitq2Bins]\tDynamic mis-tag fraction : " << static_cast<double>(countMisTag) / static_cast<double>(countMisTag + countGoodTag) << " = (" << countMisTag << "/(" << countMisTag << "+" << countGoodTag << "))" << endl;


    myString.clear(); myString.str("");
    myString << "TotalPDFq2Bin_" << i;
    InstantiateMass3AnglesFit(&TotalPDFq2Bins[i],useEffPDF,x,y,z,p,myString.str(),FitType,configParam,fitParam,i);

    if (FitType ==1)
    {
      //################
      // # scan values #
      //################
      for (int scan = 0; scan < 1; scan++)
      {   
        // #####################
        // # Initialize p.d.f. #
        // #####################
        CopyFitResults(TotalPDFq2Bins[i],i,fitParam);

        // #####################
        // # Apply constraints #
        // #####################
        ClearVars(vecConstr);
        BuildMassConstraints(vecConstr,TotalPDFq2Bins[i],"sign");
        BuildMassConstraints(vecConstr,TotalPDFq2Bins[i],"peak");
        BuildMassConstraints(vecConstr,TotalPDFq2Bins[i],"mistag");
        BuildAngularConstraints(vecConstr,TotalPDFq2Bins[i],"peak");

        // ###################
        // # Perform the fit #
        // ##################
        cout << "I am scanning the fitting values"<<endl;
        fitResult = MakeMass3AnglesFit(dataSet_q2Bins[i],&TotalPDFq2Bins[i],x,y,z,p,FitType,vecConstr,extText[i],ID);
        double nll= fitResult->minNll();
        cout << "nll=" <<nll <<endl;
        fileFitResults << "nll=" << nll << endl;
        if (CheckGoodFit(fitResult) == true) 
        {
          cout << "\n[ExtractYield::IterativeMass3AnglesFitq2Bins]\t@@@ Fit converged ! @@@" << endl;
          fileFitResults << "Fit converged " << endl;
        }
        else    
        {
          cout << "\n[ExtractYield::IterativeMass3AnglesFitq2Bins]\t@@@ Fit didn't converge ! @@@" << endl;
          fileFitResults << "Fit didn't converge " << endl;
        }
      }
    }
    //###########
    //# fit data#
    //###########
    else if (FitType ==6)
    {
      // #####################
      // # Initialize p.d.f. #
      // #####################
      CopyFitResults(TotalPDFq2Bins[i],i,fitParam);


      // #####################
      // # Apply constraints #
      // #####################
      ClearVars(vecConstr);
      BuildMassConstraints(vecConstr,TotalPDFq2Bins[i],"sign");
      BuildMassConstraints(vecConstr,TotalPDFq2Bins[i],"peak");
      BuildMassConstraints(vecConstr,TotalPDFq2Bins[i],"mistag");
      BuildAngularConstraints(vecConstr,TotalPDFq2Bins[i],"peak");

      // ###################
      // # Perform the fit #
      // ##################
      fitResult = MakeMass3AnglesFit(dataSet_q2Bins[i],&TotalPDFq2Bins[i],x,y,z,p,FitType,vecConstr,extText[i],ID);
      if (CheckGoodFit(fitResult) == true) cout << "\n[ExtractYield::IterativeMass3AnglesFitq2Bins]\t@@@ Fit converged ! @@@" << endl;
      else                                 cout << "\n[ExtractYield::IterativeMass3AnglesFitq2Bins]\t@@@ Fit didn't converge ! @@@" << endl;
    }

    // Open new ROOT file save save result 
    string outname = "Fitresult_std_2.root";
    if (!useToyDataset) {
      if (FitType==6)   outname = "Fitresult4_2.root";
      if (FitType==1)   outname = "Fitresult6_2.root";
    } else {
      if (FitType==6)   outname = "Fitresult5_2.root";
      if (FitType==1)   outname = "Fitresult7_2.root";
    }      
    TFile f(outname.c_str(),"UPDATE") ;
    f.cd();
    if (fitResult /*&& CheckGoodFit(fitResult) == true*/) fitResult->Write(Form("fitResult_Bin%d",i)) ;
    f.Close() ;
  }

}


int main(int argc, char** argv)
{
  if (argc >= 4)
  {
    // ##################
    // # Main variables #
    // ##################
    stringstream myString;

    string fileName           = "";
    string correct4Efficiency = "";
    string tmpFileName        = "";

    int specBin               = -1;
    int toyIndx     = 0;
    unsigned int fileIndx     = 0;
    unsigned int FitType      = atoi(argv[1]);

    scanInitVal               = false;
    scanIndx                  = 0;

    bool useEffPDF            = false;

    TFile* NtplFile           = NULL;

    useToyDataset        = false;

    if (FitType == 116)
    {
      FitType = 106;
      useToyDataset = true;
    }
    if (FitType == 16)
    {
      FitType = 6;
      useToyDataset = true;
    }


    if (((FitType == 1) ||(FitType == 6) || (FitType == 106) || (FitType== 206)) && (argc >= 4))

    {
      ParameterFILE = Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str();


      // ###################
      // # Read parameters #
      // ###################
      Utility = new Utils(false);
      Utility->ReadAllBins(ParameterFILE,&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins);


      // #################################
      // # Check that FitType is correct #
      // #################################
      if ((FitType == 1)      ||
          (FitType == 6)      ||
          (FitType == 106)    ||
          (FitType == 206))
      {
        fileName           = argv[2];
        correct4Efficiency = argv[3];
      }


      // ###################################
      // # Check that FitOption is correct #
      // ###################################
      if ((correct4Efficiency != "noEffCorr") && (correct4Efficiency != "yesEffCorr"))
      {
        cout << "[ExtractYield::main]\tIncorrect option parameter " << correct4Efficiency << endl;
        exit (EXIT_FAILURE);
      }


      // #################################
      // # Read the q^2 bin and the rest #
      // #################################
      if (argc >= 5) specBin = atoi(argv[4]);
      if ((correct4Efficiency == "yesEffCorr")) useEffPDF = true;
      if ((correct4Efficiency == "noEffCorr") || (correct4Efficiency == "yesEffCorr"))
      {
        if (argc >= 6)
        {
	  toyIndx = atoi(argv[5]);
	  cout<<"aaaaaaaaaaaaaaaaaaaaaa "<<toyIndx<<endl;
          fileIndx = 0;
          if (argc >= 7) tmpFileName = argv[6];
        }
      }

      if (argc>=8)
      {
	scanInitVal = true;
	scanIndx = atoi(argv[7]);
      }

      cout << "\n[ExtractYield::main]\t@@@ Input variables from command line @@@" << endl;
      cout << "- input/outputFile.root = " << fileName.c_str() << endl;
      cout << "- correct4Efficiency = "    << correct4Efficiency << endl;
      cout << "- tmpFileName = "           << tmpFileName.c_str() << endl;
      cout << "- specBin = "               << specBin << endl;
      cout << "- fileIndx = "              << fileIndx << endl;
      cout << "- FitType = "               << FitType << endl;
      cout << "- useEffPDF = "             << useEffPDF << endl;
      cout << "- ParameterFILE = "         << ParameterFILE.c_str() << endl;

      cout << "\n[ExtractYield::main]\t@@@ Internal settings @@@" << endl;
      cout << "NBINS = "         << NBINS << endl;
      cout << "MULTYIELD = "     << MULTYIELD << endl;
      cout << "NCOEFFPOLYBKG = " << NCOEFFPOLYBKG << endl;

      cout << "\nMAKEmumuPLOTS = " << MAKEmumuPLOTS << endl;
      cout << "SETBATCH  = "       << SETBATCH << endl;
      cout << "PLOT  = "           << PLOT << endl;
      cout << "SAVEPOLY = "        << SAVEPOLY << endl;
      cout << "SAVEPLOT = "        << SAVEPLOT << endl;
      cout << "RESETsigANG = "     << RESETsigANG << endl;
      cout << "RESETcomANG = "     << RESETcomANG << endl;
      cout << "FULLTOYS = "        << FULLTOYS << endl;
      cout << "FUNCERRBAND = "     << FUNCERRBAND << endl;
      cout << "MINIMIZER = "       << MINIMIZER << endl;
      cout << "GENPARAMS = "       << GENPARAMS << endl;

      cout << "\nPARAMETERFILEIN = " << PARAMETERFILEIN << endl;
      cout << "PARAMETERFILEOUT = "  << PARAMETERFILEOUT << endl;


      if (SETBATCH == true)
      {
        cout << "\n[ExtractYield::main]\t@@@ Setting batch mode @@@" << endl;
        gROOT->SetBatch(true);

        // #############################
        // # Turn off all the printout #
        // #############################
        RooMsgService::instance().setGlobalKillBelow(RooFit::FATAL);
      }
      TApplication* theApp = new TApplication("Applications", &argc, argv);


      // ##########################
      // # Set histo layout style #
      // ##########################
      SetStyle();

      // ##############################
      // # Initialize fit output file #
      // ##############################
      myString.clear(); myString.str("");
      if (specBin != -1) myString << "fileFitResults_" << specBin << ".txt";
      else               myString << "fileFitResults.txt";
      fileFitResults.open(myString.str().c_str(),ios_base::app);
      if (fileFitResults.good() == false)
      {
        cout << "[ExtractYield::main]\tError opening file : " << myString.str().c_str() << endl;
        CloseAllAndQuit(theApp,NtplFile);
      }

      // ###################
      // # Read parameters #
      // ###################
      Utility->ReadGenericParam(ParameterFILE);
      Utility->ReadSelectionCuts(ParameterFILE);
      Utility->ReadFitStartingValues(ParameterFILE,&fitParam,&configParam,Utility->ParFileBlockN("fitValBins"));

      CTRLfitWRKflow = Utility->GetGenericParam("CtrlFitWrkFlow");
      cout << "Use MINOS: "                          << Utility->GetGenericParam("UseMINOS").c_str() << " (0 = false; 1 = true)" << endl;
      cout << "Apply constraints: "                  << Utility->GetGenericParam("ApplyConstr").c_str() << " (0 = false; 1 = true)" << endl;
      cout << "Control fit workflow: "               << CTRLfitWRKflow.c_str() << endl;
      cout << "Control mis-tag fraction workflow: "  << Utility->GetGenericParam("CtrlMisTagWrkFlow").c_str() << endl;
      cout << "Save mis-tag fraction: "              << Utility->GetGenericParam("SaveMisTagFrac").c_str() << " (0 = false; 1 = true)" << endl;

      // #############################################################
      // # Make q^2 bins for histograms and fill efficiency matrices #
      // #############################################################
      q2BinsHisto = new double[q2Bins.size()];
      for (unsigned int i = 0; i < q2Bins.size(); i++) q2BinsHisto[i] = q2Bins[i];

      // ###############################################################################################
      // # Read other parameters : this also allow to understand if the parameter file is well written #
      // ###############################################################################################
      LUMI = Utility->ReadLumi(ParameterFILE);
      if (Utility->WhatIsThis(ParameterFILE) == 0) cout << "\n[ExtractYield::main]\t@@@ I recognize that this is a DATA file @@@" << endl;
      else                                         cout << "\n[ExtractYield::main]\t@@@ I recognize that this is a Monte Carlo file @@@" << endl;


      // ###################
      // # Select fit type #
      // ###################
      if ((FitType == 1) ||
          (FitType == 6)  ||
          (FitType == 106) ||
          (FitType == 206))
      {
        NtplFile = TFile::Open(fileName.c_str(),"READ");

	if (useToyDataset)
	{
	  // #################
	  // # Load datasets #
	  // #################
	  LoadDatasets(NtplFile, specBin, FitType, toyIndx);
	} else {
	  // #################
	  // # Make datasets #
	  // #################
        theTree  = (TTree*) NtplFile->Get("B0KstMuMu/B0KstMuMuNTuple");
        NTuple   = new B0KstMuMuSingleCandTreeContent();
        NTuple->Init();

        cout << "\n[ExtractYield::main]\t@@@ Making datasets @@@" << endl;
        MakeDatasets(NTuple,FitType);
	}

        // ##############################
        // # Select the proper fit type #
        // ##############################

        // #############################
        // # 3D-fit P5P-Fl per q^2 bin #
        // #############################
        cout << "\n[ExtractYield::main]\t@@@ Now fit invariant mass, cos(theta_K) and cos(theta_l) per mumu q^2 bins @@@" << endl;
        if ((FitType == 1) ||(FitType == 6))   IterativeMass3AnglesFitq2Bins(SingleCandNTuple_RejectPsi,
                                                                             useEffPDF,
                                                                             B0MassArb,
                                                                             ctL,
                                                                             ctK,
                                                                             phi,
                                                                             specBin,
                                                                             FitType,
                                                                             &q2Bins,
                                                                             &configParam,&fitParam,
                                                                             &vecConstr,
                                                                             fileIndx);

        if ((FitType == 106))                     Iterative3AnglesFitq2Bins(SingleCandNTuple_RejectPsi,
                                                                            useEffPDF,
                                                                            ctL,
                                                                            ctK,
                                                                            phi,
                                                                            specBin,
                                                                            FitType,
                                                                            &q2Bins,
                                                                            &configParam,&fitParam,
                                                                            &vecConstr,
                                                                            fileIndx);

        if ((FitType == 206))                      Iterative3AnglesFitq2Bins(SingleCandNTuple,
                                                                             useEffPDF,
                                                                             ctL,
                                                                             ctK,
                                                                             phi,
                                                                             specBin,
                                                                             FitType,
                                                                             &q2Bins,
                                                                             &configParam,&fitParam,
                                                                             &vecConstr,
                                                                             fileIndx);

      }
      fileFitResults.close();           

      if (SETBATCH == true)
      {
        cout << "Bye bye !" << endl;
        CloseAllAndQuit(theApp,NtplFile);
      }
      else
      {
        system("echo \" Let's rock and roll ! \"");
        theApp->Run (); 
      }
    }
    else
    {
      cout << "Wrong parameter: " << endl;
      cout << "./ExtractYield [FitType] [input/output[if toy-MC]File.root] [noEffCorr yesEffCorr]" << endl;
      cout << "               [q^2 bin to fit (0 - ...)]" << endl;

      cout << "\n --> noEffCorr     = no eff. correction" << endl;
      cout << " --> yesEffCorr    = use eff. correction" << endl;

      return EXIT_FAILURE;
    }
  }
  else
  {
    cout << "Parameter missing: " << endl;
    cout << "./ExtractYield [FitType] [input/output[if toy-MC]File.root] [noEffCorr yesEffCorr]" << endl;
    cout << "               [q^2 bin to fit (0 - ...)]" << endl;

    cout << "\n --> noEffCorr     = no eff. correction" << endl;
    cout << " --> yesEffCorr    = use eff. correction" << endl;
    cout << "\n@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  Signa  @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;
    cout << "FitType = 1: scanning the fitting values per q^2 bin on Data" << endl;
    cout << "FitType = 6: 4D P5p-Fl (B0Mass, cos(theta_K), cos(theta_l), phi) per q^2 bin on Data" << endl;
    cout << "FitType = 206: 3D P5p-Fl (B0Mass, cos(theta_K), cos(theta_l)) per q^2 bin on Gen-level" << endl;
    cout << "FitType = 106: 3D P5p-Fl (B0Mass, cos(theta_K), cos(theta_l)) per q^2 bin on Reco-level" << endl;
    cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << endl;

    return EXIT_FAILURE;
  }

  return 0;
}
