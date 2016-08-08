// #########################################################################
// # Program to compute the efficiency for the B0 --> K*0 mu+ mu- analysis #
// # in bins of dimuon q^2, cos(theta_K), cos(theta_l), and phi            #
// #########################################################################
// # Author: Stefano Lacaprara - Alessio Boletti
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
#include <Math/Functor.h>

#include <RooNumIntConfig.h>
#include <RooGenericPdf.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"
#include "RooProdPdf.h"
#include <RooFunctorBinding.h>
#include "RooBinning.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"

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

using namespace RooFit ;

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


const string fileNameGenCandidatesNoFilter = "dcap://t2-srm-02.lnl.infn.it//pnfs/lnl.infn.it/data/cms/store/user/slacapra/B0KsMuMu/Data2012B0KstMuMuResults/MonteCarlo2012/ForEfficiency/B0ToKstMuMu_GEN_NoFilter_MC_NTuples_addGENvars.root";
const string fileNameRecoCandidates        = "dcap://t2-srm-02.lnl.infn.it////pnfs/lnl.infn.it/data/cms/store/user/slacapra/B0KsMuMu/Data2012B0KstMuMuResults/MonteCarlo2012/ForEfficiency/B0ToKstMuMu_MC_NTuple_addGENvars.root";
const string fileNameSingleCand            = "dcap://t2-srm-02.lnl.infn.it////pnfs/lnl.infn.it/data/cms/store/user/slacapra/B0KsMuMu/Data2012B0KstMuMuResults/MonteCarlo2012/SingleCand/singleCand_B0ToKstMuMu_MC_NTuple.root";
const string SignalType                    = "1";

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

int Counter[4];
double Vector[4];

std::vector<RooNDKeysPdf*>* ele_KEpdf;
RooFunctorPdfBinding* EffGenPDF;
RooFunctorPdfBinding* EffRecoPDF;
RooFunctorPdfBinding* EffPDF;

vector<double> q2Bins;
double* q2Bins_ ;
double* cosThetaKBins_ ;
double* cosThetaLBins_ ;
double* phiBins_ ;
RooBinning cosThetaLBinning;
RooBinning cosThetaKBinning;
RooBinning phiBinning;
RooRealVar CosThetaK("CosThetaK","cos#theta_{K}",-1,1);
RooRealVar CosThetaMu("CosThetaMu","cos#theta_{#mu}",-1,1) ;
RooRealVar AbsCosThetaMu("AbsCosThetaMu","abs(cos#theta_{#mu})",0,1) ;
RooRealVar ThetaK("ThetaK","#theta_{K}",0,TMath::Pi());
RooRealVar ThetaMu("ThetaMu","#theta_{#mu}",0,TMath::Pi()) ;
RooRealVar PhiKstMuMuPlane("PhiKstMuMuPlane","#phi_{K*#mu#mu}",-TMath::Pi(),TMath::Pi()) ;
RooRealVar AbsPhiKstMuMuPlane("AbsPhiKstMuMuPlane","abs(#phi_{K*#mu#mu})",0,TMath::Pi()) ;

unsigned int termIdx=0;

// #######################
// # Function Definition #
// #######################
void SetStyle          ();
void ComputeEfficiency (TTree* theTree, B0KstMuMuSingleCandTreeContent* NTuple, double q2BinIdx, unsigned int type, int SignalType);

// ###########################
// # Function Implementation #
// ###########################
struct MyRatioPdf
{
  public:
    MyRatioPdf (RooAbsPdf& numPdf, RooAbsPdf& denPdf) : _numPdf(numPdf), _denPdf(denPdf)
  { 
    const RooArgSet* allvar1 = numPdf.getVariables();
    const RooArgSet* allvar2 = denPdf.getVariables();

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

      return _numPdf.getVal() / _denPdf.getVal();
    }


  private:
    RooAbsPdf& _numPdf; 
    RooAbsPdf& _denPdf; 
    RooArgList _vars; 
};

struct MyProdPdf
{
  public:
    MyProdPdf (RooAbsPdf& pdf1, RooAbsPdf& pdf2) : _pdf1(pdf1), _pdf2(pdf2)
  { 
    const RooArgSet* allvar1 = _pdf1.getVariables();
    const RooArgSet* allvar2 = _pdf2.getVariables();

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

struct MyEffPdf
{
  public:
    MyEffPdf (RooAbsPdf& pdf1, RooAbsPdf& pdf2, RooAbsPdf& pdf3, RooAbsPdf& pdf4) : _pdf1(pdf1), _pdf2(pdf2), _pdf3(pdf3), _pdf4(pdf4)
  { 
    const RooArgSet* allvar1 = pdf1.getVariables();
    const RooArgSet* allvar2 = pdf2.getVariables();
    const RooArgSet* allvar3 = pdf3.getVariables();
    const RooArgSet* allvar4 = pdf4.getVariables();

    _vars.add(*allvar1);
    _vars.add(*allvar2,false);
    _vars.add(*allvar3,false);
    _vars.add(*allvar4,false);

    delete allvar1;
    delete allvar2;
    delete allvar3;
    delete allvar4;
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

      return _pdf2.getVal() / _pdf1.getVal() * _pdf4.getVal() / _pdf3.getVal();
    }


  private:
    RooAbsPdf& _pdf1; 
    RooAbsPdf& _pdf2; 
    RooAbsPdf& _pdf3; 
    RooAbsPdf& _pdf4; 
    RooArgList _vars; 
};

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


void loadBinning() {
  vector<double> cosThetaKBins;
  vector<double> cosThetaLBins;
  vector<double> phiBins;
  Utility->ReadAllBins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"goodTag");
  q2Bins_ = Utility->MakeBinning(&q2Bins);
  cosThetaKBins_ = Utility->MakeBinning(&cosThetaKBins);
  cosThetaLBins_ = Utility->MakeBinning(&cosThetaLBins);
  phiBins_       = Utility->MakeBinning(&phiBins);
  cosThetaKBinning = RooBinning(cosThetaKBins.size()-1,cosThetaKBins_);
  cosThetaLBinning = RooBinning(cosThetaLBins.size()-1,cosThetaLBins_);
  phiBinning = RooBinning(phiBins.size()-1,phiBins_);
}

void ComputeEfficiency (TTree* theTree, B0KstMuMuSingleCandTreeContent* NTuple, double q2BinIdx, unsigned int type, int SignalType, float width)
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

  double mumuq2;
  double cosThetaKArb;
  double cosThetaMuArb;
  double phiKstMuMuPlaneArb;
  // ###################

  nEntries = theTree->GetEntries();
  cout<<nEntries<<endl;
  // ###################################
  // # Initialize efficiency structure #
  // ###################################
  NTuple->ClearNTuple();
  NTuple->SetBranchAddresses(theTree);
  cout << "\n[ComputeEfficiency::ComputeEfficiency]\t@@@ Computing efficiency type " << type;
  if      (type == 1) cout << " (before filter) @@@" << endl;
  else if (type == 2) cout << " (after filter) @@@" << endl;
  else if (type == 3) cout << " (reco events) @@@" << endl;
  else if (type == 4) cout << " (single candidate events) @@@" << endl;
  cout << "[ComputeEfficiency::ComputeEfficiency]\t@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;

  RooDataSet KEpdf_data("KEpdf_data", "KEpdf_data", RooArgSet(CosThetaK,AbsCosThetaMu,AbsPhiKstMuMuPlane));
  //RooDataSet KEpdf_data("KEpdf_data", "KEpdf_data", RooArgSet(ThetaK,ThetaMu, PhiKstMuMuPlane));

  for (int entry = 0; entry < nEntries; entry++)
  {
    theTree->GetEntry(entry);

    if ((NTuple->B0pT > Utility->GetSeleCut("B0pT")) && (fabs(NTuple->B0Eta) < Utility->GetSeleCut("B0Eta")) &&

        ((NTuple->genSignal == SignalType || NTuple->genSignal == SignalType+1)) &&

        ((type == 1 || type == 3) ||

         ((type == 2) &&
          (sqrt(NTuple->genMumPx*NTuple->genMumPx + NTuple->genMumPy*NTuple->genMumPy)   > Utility->GetPreCut("MinMupT")) &&
          (sqrt(NTuple->genMupPx*NTuple->genMupPx + NTuple->genMupPy*NTuple->genMupPy)   > Utility->GetPreCut("MinMupT")) &&
          (fabs(Utility->computeEta(NTuple->genMumPx,NTuple->genMumPy,NTuple->genMumPz)) < Utility->GetPreCut("MuEta"))   &&
          (fabs(Utility->computeEta(NTuple->genMupPx,NTuple->genMupPy,NTuple->genMupPz)) < Utility->GetPreCut("MuEta")))  ||

         ((type == 4) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag != Utility->RIGHTflavorTAG) &&
          (NTuple->B0MassArb > Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str()))            &&
          (NTuple->B0MassArb < Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str()))           &&

          ((((SignalType == Utility->B0ToKstMuMu)  || (SignalType == Utility->B0ToKstMuMu+1))  && (Utility->PsiRejection(NTuple->B0MassArb,NTuple->mumuMass->at(0),NTuple->mumuMassE->at(0),"rejectPsi",true) == true)) ||
           (((SignalType == Utility->B0ToJPsiKst)  || (SignalType == Utility->B0ToJPsiKst+1))  && (Utility->PsiRejection(NTuple->B0MassArb,NTuple->mumuMass->at(0),NTuple->mumuMassE->at(0),"keepJpsi")       == true)) ||
           (((SignalType == Utility->B0ToPsi2SKst) || (SignalType == Utility->B0ToPsi2SKst+1)) && (Utility->PsiRejection(NTuple->B0MassArb,NTuple->mumuMass->at(0),NTuple->mumuMassE->at(0),"keepPsiP")       == true))))))
    {
      mumuq2 = NTuple->mumuMass->at(0)*NTuple->mumuMass->at(0);

      if (mumuq2 < q2Bins[q2BinIdx-1] || mumuq2 > q2Bins[q2BinIdx]) continue;

      if ((type == 4) && (NTuple->rightFlavorTag == false))
      {
        cosThetaKArb  = - NTuple->CosThetaKArb;
        cosThetaMuArb = - NTuple->CosThetaMuArb;
      }
      else
      {
        cosThetaKArb  = NTuple->CosThetaKArb;
        cosThetaMuArb = NTuple->CosThetaMuArb;
      }
      phiKstMuMuPlaneArb = NTuple->PhiKstMuMuPlaneArb;

      Counter[type-1]++;

      CosThetaK.setVal(cosThetaKArb);
      CosThetaMu.setVal(cosThetaMuArb);
      AbsCosThetaMu.setVal(TMath::Abs(cosThetaMuArb));
      ThetaK.setVal(TMath::ACos(cosThetaKArb));
      ThetaMu.setVal(TMath::ACos(cosThetaMuArb));
      PhiKstMuMuPlane.setVal(phiKstMuMuPlaneArb);
      AbsPhiKstMuMuPlane.setVal(TMath::Abs(phiKstMuMuPlaneArb));
      KEpdf_data.add(RooArgSet(CosThetaK,AbsCosThetaMu,AbsPhiKstMuMuPlane));
      //KEpdf_data.add(RooArgSet(ThetaK,ThetaMu, PhiKstMuMuPlane));

      // ##################################################
      // # No event-weight for muon acceptance efficiency #
      // ##################################################
      if ((type == 1) || (type == 2)) Vector[type-1]++;

      // ########################################
      // # Add event-weight for reco efficiency #
      // ########################################
      else if ((type == 3) || (type == 4)) Vector[type-1] += NTuple->evWeight;
    }
  }

  cout << "[ComputeEfficiency::ComputeEfficiency]\t@@@ Total number of passing events: " << Counter[type-1] << " @@@" << endl;
  ele_KEpdf->at(type-1) = new RooNDKeysPdf(Form("ele_KEpdf%i",type-1), Form("ele_KEpdf%i",type-1), RooArgList(CosThetaK,AbsCosThetaMu,AbsPhiKstMuMuPlane), KEpdf_data, "AM", width);
  //ele_KEpdf->at(type-1) = new RooNDKeysPdf(Form("ele_KEpdf%i",type-1), Form("ele_KEpdf%i",type-1), RooArgList(ThetaK,ThetaMu, PhiKstMuMuPlane), KEpdf_data, "m");
}

void plot3D(unsigned int q2bin, const TH3* xyzEff, TString suffix="") {
  // TH3* h3dfull = (TH3*)xyzEff.createHistogram("hKE",CosThetaK,Binning(20),YVar(CosThetaMu,Binning(20)),ZVar(PhiKstMuMuPlane,Binning(20)));
  // TH2F* h2d1full = new TH2F("h"+(TString)h3dfull->GetName()+"1"+suffix,"histo "+(TString)h3dfull->GetName()+" bin_{#phi}=1 full:cos#theta_{K}::cos#theta_{L}"+suffix,20,-1.,1.,20,-1.,1.);
  // h3dfull->GetZaxis()->SetRange(1,5);
  // h2d1full=(TH2F*) h3dfull->Project3D("z");

  // TCanvas* c3Dfull = new TCanvas("c3Dfull","c3Dfull",800,800) ;
  // c3Dfull->Divide(2,2);
  // c3Dfull->cd(1);
  // TH2* hhcKcLpdf = (TH2*)xyzEff.createHistogram("hhcKcLpdf",CosThetaK,Binning(cosThetaKBinning),YVar(CosThetaMu,Binning(cosThetaLBinning)));
  // hhcKcLpdf->DrawCopy("lego2 fp");
  TH1* hh_pdf3xy = xyzEff->Project3D("xy");
  TH1* hh_pdf3xz = xyzEff->Project3D("xz");
  TH1* hh_pdf3yz = xyzEff->Project3D("yz");
  hh_pdf3xy->SetLineColor(kBlue) ;
  hh_pdf3xz->SetLineColor(kBlue) ;
  hh_pdf3yz->SetLineColor(kBlue) ;

  TCanvas* c = new TCanvas("Plot3d",Form("Plot3D_%s_%s_q2bin%i",xyzEff->GetName(),suffix.Data(),q2bin),800,800) ;
  c->Divide(2,2);
  c->cd(1);
  gPad->SetLeftMargin(0.15) ;
  hh_pdf3xy->GetXaxis()->SetTitleOffset(1.2) ;
  hh_pdf3xy->GetYaxis()->SetTitleOffset(1.2) ;
  hh_pdf3xy->DrawCopy("surf fp") ;
  c->cd(2);
  gPad->SetLeftMargin(0.15) ;
  hh_pdf3xz->GetXaxis()->SetTitleOffset(1.2) ;
  hh_pdf3xz->GetZaxis()->SetTitleOffset(1.2) ;
  hh_pdf3xz->DrawCopy("surf fp") ;
  c->cd(3);
  gPad->SetLeftMargin(0.15) ;
  hh_pdf3yz->GetYaxis()->SetTitleOffset(1.2) ;
  hh_pdf3yz->GetZaxis()->SetTitleOffset(1.2) ;
  hh_pdf3yz->DrawCopy("surf fp") ;
  c->Draw();
  c->Print("c"+(TString)c->GetTitle()+".pdf");
}

void plot3D(unsigned int q2bin, const RooAbsPdf& xyzEff, TString suffix="") {
  // TH3* h3dfull = (TH3*)xyzEff.createHistogram("hKE",CosThetaK,Binning(20),YVar(CosThetaMu,Binning(20)),ZVar(PhiKstMuMuPlane,Binning(20)));
  // TH2F* h2d1full = new TH2F("h"+(TString)h3dfull->GetName()+"1"+suffix,"histo "+(TString)h3dfull->GetName()+" bin_{#phi}=1 full:cos#theta_{K}::cos#theta_{L}"+suffix,20,-1.,1.,20,-1.,1.);
  // h3dfull->GetZaxis()->SetRange(1,5);
  // h2d1full=(TH2F*) h3dfull->Project3D("z");

  xyzEff.Print("v");

  // TCanvas* c3Dfull = new TCanvas("c3Dfull","c3Dfull",800,800) ;
  // c3Dfull->Divide(2,2);
  // c3Dfull->cd(1);
  // TH2* hhcKcLpdf = (TH2*)xyzEff.createHistogram("hhcKcLpdf",CosThetaK,Binning(cosThetaKBinning),YVar(CosThetaMu,Binning(cosThetaLBinning)));
  // hhcKcLpdf->DrawCopy("lego2 fp");
  //TH1* hh_pdf3xy = xyzEff.createHistogram("hh_pdfxy",CosThetaK,Binning(10),YVar(CosThetaMu,Binning(10))) ;
  //TH1* hh_pdf3xz = xyzEff.createHistogram("hh_pdfxz",CosThetaK,Binning(10),YVar(PhiKstMuMuPlane,Binning(10))) ;
  //TH1* hh_pdf3yz = xyzEff.createHistogram("hh_pdfyz",CosThetaMu,Binning(10),YVar(PhiKstMuMuPlane,Binning(10))) ;
  TH1* hh_pdf3xy = xyzEff.createHistogram("hh_pdfxy",ThetaK,Binning(10),YVar(ThetaMu,Binning(10))) ;
  TH1* hh_pdf3xz = xyzEff.createHistogram("hh_pdfxz",ThetaK,Binning(10),YVar(PhiKstMuMuPlane,Binning(10))) ;
  TH1* hh_pdf3yz = xyzEff.createHistogram("hh_pdfyz",ThetaMu,Binning(10),YVar(PhiKstMuMuPlane,Binning(10))) ;
  hh_pdf3xy->SetLineColor(kBlue) ;
  hh_pdf3xz->SetLineColor(kBlue) ;
  hh_pdf3yz->SetLineColor(kBlue) ;

  TCanvas* c = new TCanvas("Test_RooNDKeysPdf","Test_RooNDKeysPdf",800,800) ;
  c->Divide(2,2);
  c->cd(1);
  gPad->SetLeftMargin(0.15) ;
  hh_pdf3xy->GetXaxis()->SetTitleOffset(1.2) ;
  hh_pdf3xy->GetYaxis()->SetTitleOffset(1.2) ;
  hh_pdf3xy->DrawCopy("surf fp") ;
  c->cd(2);
  gPad->SetLeftMargin(0.15) ;
  hh_pdf3xz->GetXaxis()->SetTitleOffset(1.2) ;
  hh_pdf3xz->GetZaxis()->SetTitleOffset(1.2) ;
  hh_pdf3xz->DrawCopy("surf fp") ;
  c->cd(3);
  gPad->SetLeftMargin(0.15) ;
  hh_pdf3yz->GetYaxis()->SetTitleOffset(1.2) ;
  hh_pdf3yz->GetZaxis()->SetTitleOffset(1.2) ;
  hh_pdf3yz->DrawCopy("surf fp") ;
  c->Draw();
  c->Print("c"+(TString)xyzEff.GetName()+Form("_%i",q2bin)+suffix+".pdf");


  // c3Dfull->cd(2);
  // h2d1full->Draw("lego3");
  // c3Dfull->cd(3);
  // h2d1full->Draw("surf");
  // c3Dfull->cd(4);
  // h3dfull->Draw("box");


  // TH3* h3d = (TH3*)xyzEff.createHistogram("hKE",CosThetaK,Binning(cosThetaKBinning),YVar(CosThetaMu,Binning(cosThetaLBinning)),ZVar(PhiKstMuMuPlane,Binning(phiBinning)));

  // TCanvas* c3D = new TCanvas("c"+(TString)h3d->GetName()+Form("_%i",q2bin)+suffix,(TString)h3d->GetName()+Form(" bin%i",q2bin)+" "+suffix,800,800) ;
  // c3D->Divide(2,2);
  // TH2F* h2d1 = new TH2F("h"+(TString)h3d->GetName()+"1"+suffix,"h"+(TString)h3d->GetName()+"1"+suffix,cosThetaKBinning.numBoundaries()-1,cosThetaKBins_,cosThetaLBinning.numBoundaries()-1,cosThetaLBins_);
  // TH2F* h2d2 = new TH2F("h"+(TString)h3d->GetName()+"2"+suffix,"h"+(TString)h3d->GetName()+"2"+suffix,cosThetaKBinning.numBoundaries()-1,cosThetaKBins_,cosThetaLBinning.numBoundaries()-1,cosThetaLBins_);
  // TH2F* h2d3 = new TH2F("h"+(TString)h3d->GetName()+"3"+suffix,"h"+(TString)h3d->GetName()+"3"+suffix,cosThetaKBinning.numBoundaries()-1,cosThetaKBins_,cosThetaLBinning.numBoundaries()-1,cosThetaLBins_);
  // TH2F* h2d4 = new TH2F("h"+(TString)h3d->GetName()+"4"+suffix,"h"+(TString)h3d->GetName()+"4"+suffix,cosThetaKBinning.numBoundaries()-1,cosThetaKBins_,cosThetaLBinning.numBoundaries()-1,cosThetaLBins_);
  // for (int i=1; i<cosThetaKBinning.numBoundaries(); i++) for (int j=1; j<cosThetaLBinning.numBoundaries(); j++) {
  //   h2d1->SetBinContent(i, j, h3d->GetBinContent(i, j, 1) );
  //   h2d2->SetBinContent(i, j, h3d->GetBinContent(i, j, 2) );
  //   h2d3->SetBinContent(i, j, h3d->GetBinContent(i, j, 3) );
  //   h2d4->SetBinContent(i, j, h3d->GetBinContent(i, j, 4) );
  //   h2d1->SetBinError(i, j, h3d->GetBinError(i, j, 1) );
  //   h2d2->SetBinError(i, j, h3d->GetBinError(i, j, 2) );
  //   h2d3->SetBinError(i, j, h3d->GetBinError(i, j, 3) );
  //   h2d4->SetBinError(i, j, h3d->GetBinError(i, j, 4) );
  // }
  // c3D->cd(1);
  // h2d1->DrawCopy("lego2 fp");
  // c3D->cd(2);
  // h2d2->DrawCopy("lego2 fp");
  // c3D->cd(3);
  // h2d3->DrawCopy("lego2 fp");
  // c3D->cd(4);
  // h2d4->DrawCopy("lego2 fp");
  // c3D->Print((TString)c3D->GetName()+".pdf");
}


void computeEffForOneBin(int q2BinIndx, bool doPlot, float width) {
  cout << "computeEffForOneBin " << q2BinIndx << endl;
  // Open input file
  ele_KEpdf = new std::vector<RooNDKeysPdf*> (4);
  string fileNameOutput = Form("effKEpdf_b%i.root",q2BinIndx);
  TFile* fileOut = new TFile(fileNameOutput.c_str(), "UPDATE");

  int nbin=40;
  fileOut->cd();

  if (termIdx == 1 || termIdx == 0) {
    cout << "Gen den " ;
    ComputeEfficiency(theTreeGenCandidatesNoFilter,NTupleGenCandidatesNoFilter,q2BinIndx, 1 ,atoi(SignalType.c_str()),width);
    TH3* h3genDen = (TH3*)(ele_KEpdf->at(0))->createHistogram("h3genDen",CosThetaK,Binning(nbin),YVar(AbsCosThetaMu,Binning(nbin)),ZVar(AbsPhiKstMuMuPlane,Binning(nbin)));
    //TH3* h3genDen = (TH3*)(ele_KEpdf->at(0))->createHistogram("h3genDen",ThetaK,Binning(nbin),YVar(ThetaMu,Binning(nbin)),ZVar(PhiKstMuMuPlane,Binning(nbin)));
    h3genDen->Scale(Counter[0]*1./h3genDen->Integral());
    h3genDen->Write("h3genDen");
  }

  // // Create HistPdf from this TH3
  // RooArgList ctKctLPhi(CosThetaK,CosThetaMu, PhiKstMuMuPlane);
  // RooDataHist rh3genDen("rh3genDen","rh3genDen", ctKctLPhi, Import(*h3genDen,kTRUE));
  // RooHistPdf pdf_genDen(Form("pdf_genDen_q2bin%d",q2BinIndx),Form("pdf_genDen q2bin=%d",q2BinIndx), ctKctLPhi, rh3genDen, 0);
  // pdf_genDen.Write(Form("pdf_genDen_q2bin%d",q2BinIndx));

  if (termIdx == 2 || termIdx == 0) {
    cout << "Gen num " ;
    ComputeEfficiency(theTreeGenCandidatesNoFilter,NTupleGenCandidatesNoFilter,q2BinIndx, 2 ,atoi(SignalType.c_str()),width);
    TH3* h3genNum = (TH3*)(*ele_KEpdf->at(1)).createHistogram("h3genNum",CosThetaK,Binning(nbin),YVar(AbsCosThetaMu,Binning(nbin)),ZVar(AbsPhiKstMuMuPlane,Binning(nbin)));
    //TH3* h3genNum = (TH3*)(*ele_KEpdf->at(1)).createHistogram("h3genNum",ThetaK,Binning(nbin),YVar(ThetaMu,Binning(nbin)),ZVar(PhiKstMuMuPlane,Binning(nbin)));
    h3genNum->Scale(Counter[1]*1./h3genNum->Integral());
    h3genNum->Write("h3genNum");
  }
  // RooDataHist rh3genNum("rh3genNum","rh3genNum", ctKctLPhi, Import(*h3genNum,kTRUE));
  // RooHistPdf pdf_genNum(Form("pdf_genNum_q2bin%d",q2BinIndx),Form("pdf_genNum q2bin=%d",q2BinIndx), ctKctLPhi, rh3genNum, 0);
  // pdf_genNum.Write(Form("pdf_genNum_q2bin%d",q2BinIndx));

  //TH3* h3genEff = (TH3*)h3genNum->Clone("h3genEff");
  //h3genEff->Divide(h3genDen);
  //h3genEff->Scale((Counter[1]*1./Counter[0])/h3genEff->Integral());
  //cout << "Gen: N=" << Counter[1] << "/D=" << Counter[0] << "=" << h3genEff->Integral() << endl;
  //h3genEff->Write(Form("h3genEff_q2bin%i",q2BinIndx));

  // RooDataHist rh3genEff("rh3genEff","rh3genEff", ctKctLPhi, Import(*h3genEff,kTRUE));
  // RooHistPdf pdf_genEff(Form("pdf_genEff_q2bin%d",q2BinIndx),Form("pdf_genEff q2bin=%d",q2BinIndx), ctKctLPhi, rh3genEff, 0);
  // pdf_genEff.Write(Form("pdf_genEff_q2bin%d",q2BinIndx));

  //if (doPlot)  plot3D(q2BinIndx, h3genEff, "_effGen");

  // // Compute ratio for Gen efficiency
  // MyRatioPdf* ratioGen  = new MyRatioPdf(*ele_KEpdf->at(1),*ele_KEpdf->at(0));
  // ROOT::Math::Functor* effGenFunctor = new ROOT::Math::Functor(*ratioGen,ratioGen->ndim());
  // EffGenPDF             = new RooFunctorPdfBinding("EffGenPDF","EffGenPDF",*effGenFunctor,ratioGen->vars());
  // if (doPlot)  plot3D(q2BinIndx, *EffGenPDF, "_effGen");
  //cout << " Done" << endl;

  // numerator and denominator for RECO ratio
  if (termIdx == 3 || termIdx == 0) {
    cout << "Reco den " ;
    ComputeEfficiency(theTreeRecoCandidates,       NTupleRecoCandidates,       q2BinIndx, 3 ,atoi(SignalType.c_str()),width);
    TH3* h3recoDen = (TH3*)(ele_KEpdf->at(2))->createHistogram("h3recoDen",CosThetaK,Binning(nbin),YVar(AbsCosThetaMu,Binning(nbin)),ZVar(AbsPhiKstMuMuPlane,Binning(nbin)));
    //TH3* h3recoDen = (TH3*)(ele_KEpdf->at(2))->createHistogram("h3recoDen",ThetaK,Binning(nbin),YVar(ThetaMu,Binning(nbin)),ZVar(PhiKstMuMuPlane,Binning(nbin)));
    h3recoDen->Scale(Counter[2]*1./h3recoDen->Integral());
    h3recoDen->Write("h3recoDen");
  }

  // RooDataHist rh3recoDen("rh3recoDen","rh3recoDen", ctKctLPhi, Import(*h3recoDen,kTRUE));
  // RooHistPdf pdf_recoDen(Form("pdf_recoDen_q2bin%d",q2BinIndx),Form("pdf_recoDen q2bin=%d",q2BinIndx), ctKctLPhi, rh3recoDen, 0);
  // pdf_recoDen.Write(Form("pdf_recoDen_q2bin%d",q2BinIndx));

  if (termIdx == 4 || termIdx == 0) {
    cout << "Reco num " ;
    ComputeEfficiency(theTreeSingleCand,           NTupleSingleCand,           q2BinIndx, 4 ,atoi(SignalType.c_str()),width);
    TH3* h3recoNum = (TH3*)(ele_KEpdf->at(3))->createHistogram("h3recoNum",CosThetaK,Binning(nbin),YVar(AbsCosThetaMu,Binning(nbin)),ZVar(AbsPhiKstMuMuPlane,Binning(nbin)));
    //TH3* h3recoNum = (TH3*)(ele_KEpdf->at(3))->createHistogram("h3recoNum",ThetaK,Binning(nbin),YVar(ThetaMu,Binning(nbin)),ZVar(PhiKstMuMuPlane,Binning(nbin)));
    h3recoNum->Scale(Counter[3]*1./h3recoNum->Integral());
    h3recoNum->Write("h3recoNum");
  }

  // RooDataHist rh3recoNum("rh3recoNum","rh3recoNum", ctKctLPhi, Import(*h3recoNum,kTRUE));
  // RooHistPdf pdf_recoNum(Form("pdf_recoNum_q2bin%d",q2BinIndx),Form("pdf_recoNum q2bin=%d",q2BinIndx), ctKctLPhi, rh3recoNum, 0);
  // pdf_recoNum.Write(Form("pdf_recoNum_q2bin%d",q2BinIndx));

  //TH3* h3recoEff = (TH3*)h3recoNum->Clone("h3recoEff");
  //h3recoEff->Divide(h3recoDen);
  //h3recoEff->Scale((Counter[3]*1./Counter[2])/h3recoEff->Integral());
  //cout << "Reco: N=" << Counter[3] << "/D=" << Counter[2] << "=" << h3recoEff->Integral()<< endl;
  //fileOut->cd();
  //h3recoEff->Write(Form("h3recoEff_q2bin%i",q2BinIndx));

  // RooDataHist rh3recoEff("rh3recoEff","rh3recoEff", ctKctLPhi, Import(*h3recoEff,kTRUE));
  // RooHistPdf pdf_recoEff(Form("pdf_recoEff_q2bin%d",q2BinIndx),Form("pdf_recoEff q2bin=%d",q2BinIndx), ctKctLPhi, rh3recoEff, 0);
  // pdf_recoEff.Write(Form("pdf_recoEff_q2bin%d",q2BinIndx));

  //if (doPlot)  plot3D(q2BinIndx, h3recoEff, "_effReco");

  // // Compute ratio for RECO efficiency
  // MyRatioPdf* ratioReco  = new MyRatioPdf(*ele_KEpdf->at(3),*ele_KEpdf->at(2));
  // ROOT::Math::Functor* effRecoFunctor = new ROOT::Math::Functor(*ratioReco,ratioReco->ndim());
  // EffRecoPDF             = new RooFunctorPdfBinding("EffRecoPDF","EffRecoPDF",*effRecoFunctor,ratioReco->vars());
  // if (doPlot)  plot3D(q2BinIndx, *EffRecoPDF, "_effReco");
  // cout << "Done " << endl;

  // // Create a new empty workspace
  // RooWorkspace *w = new RooWorkspace("w","workspace") ;

  // // Import pdf
  // for (int i=0; i<4; ++i) w->import(*ele_KEpdf->at(i));
  // w->import(*EffGenPDF);
  // w->import(*EffRecoPDF);

  // // Print workspace contents
  // w->Print() ;

  // w->writeToFile("RooWorkspace.root") ;

  // TFile* fileOut = new TFile(fileNameOutput.c_str(), "UPDATE");
  // fileOut->cd();
  // for (int i=0; i<4; i++) ele_KEpdf->at(i)->Write(Form("term%iPdf_q2bin%i",i,q2BinIndx));
  // EffGenPDF->Write(Form("effGenPdf_q2bin%i",q2BinIndx));
  // EffRecoPDF->Write(Form("effRecoPdf_q2bin%i",q2BinIndx));
  // fileOut->Close();
  //if (doPlot) for (int i=0; i<4; i++) plot3D(q2BinIndx, *ele_KEpdf->at(i), Form("_term%i",i));

  //TH3* h3Eff = (TH3*)h3genEff->Clone("h3Eff");
  //h3Eff->Multiply(h3recoEff);
  //h3Eff->Scale((Counter[1]*1./Counter[0])*(Counter[3]*1./Counter[2])/h3Eff->Integral());
  //cout << "h3Eff " << h3Eff->Integral() << endl;
  //fileOut->cd();
  //h3Eff->Write(Form("h3Eff_q2bin%i",q2BinIndx));

  // RooDataHist rh3Eff("rh3Eff","rh3Eff", ctKctLPhi, Import(*h3Eff,kTRUE));
  // RooHistPdf pdf_Eff(Form("pdf_Eff_q2bin%d",q2BinIndx),Form("pdf_Eff q2bin=%d",q2BinIndx), ctKctLPhi, rh3Eff, 0);
  // pdf_Eff.Write(Form("pdf_Eff_q2bin%d",q2BinIndx));

  //if (doPlot)  plot3D(q2BinIndx, h3Eff, "_eff");
  fileOut->Close();


  // cout << "Doing product " << endl;
  // // Try do do double ratio, which is a product
  // MyProdPdf* prodPdf  = new MyProdPdf(*EffGenPDF, *EffRecoPDF);
  // ROOT::Math::Functor* effFunctor = new ROOT::Math::Functor(*prodPdf,prodPdf->ndim());
  // EffPDF             = new RooFunctorPdfBinding("EffPDF","EffPDF",*effFunctor,prodPdf->vars());
  // cout << "Done" << endl;

  // MyEffPdf* myeffpdf  = new MyEffPdf(*ele_KEpdf->at(0),*ele_KEpdf->at(1),*ele_KEpdf->at(2),*ele_KEpdf->at(3));
  // ROOT::Math::Functor* effFunctor = new ROOT::Math::Functor(*myeffpdf,myeffpdf->ndim());
  // EffPDF             = new RooFunctorPdfBinding("EffPDF","EffPDF",*effFunctor,myeffpdf->vars());

  // fileOut = new TFile(fileNameOutput.c_str(), "UPDATE");
  // fileOut->cd();
  // EffPDF->Write(Form("effPdf_q2bin%i",q2BinIndx));
  // fileOut->Close();
  //if (doPlot)  plot3D(q2BinIndx, *EffPDF, "_eff");
  // EffPDF->Delete();
  delete ele_KEpdf;
}

int main(int argc, char** argv)
{
  if (argc > 0)
  {
    // ##################
    // # Main variables #
    // ##################

    TApplication theApp ("Applications", &argc, argv);
    Utility = new Utils(false);
    // CosThetaK.defaultIntegratorConfig()->method2D().setLabel("RooMCIntegrator");
    // CosThetaMu.defaultIntegratorConfig()->method2D().setLabel("RooMCIntegrator");
    // PhiKstMuMuPlane.defaultIntegratorConfig()->method2D().setLabel("RooMCIntegrator");

    unsigned int q2BinIndx=0;
    bool doPlot=false;
    bool doBatch=false;
    float width=1.;
    if ( argc > 1 ) q2BinIndx = atoi(argv[1]);
    if ( argc > 2 ) doPlot = atoi(argv[2]);
    if ( argc > 3 ) doBatch = atoi(argv[3]);
    if ( argc > 4 ) termIdx = atoi(argv[4]);
    if ( argc > 5 ) width = atof(argv[5]);

    if (doBatch) {
      cout << "\n[" << argv[0] << "::main]\t@@@ Setting batch mode @@@" << endl;
      gROOT->SetBatch(doBatch);
    }

    
    TFile* NtplFileGenCandidatesNoFilter =  TFile::Open(fileNameGenCandidatesNoFilter.c_str(), "READ");
    theTreeGenCandidatesNoFilter         = (TTree*) NtplFileGenCandidatesNoFilter->Get("B0KstMuMu/B0KstMuMuNTuple");

    //cout<<"Tree entries: "<<theTreeGenCandidatesNoFilter->GetEntries()<<endl;

    NTupleGenCandidatesNoFilter          = new B0KstMuMuSingleCandTreeContent();
    NTupleGenCandidatesNoFilter->Init();

    TFile* NtplFileRecoCandidates =  TFile::Open(fileNameRecoCandidates.c_str(), "READ");
    theTreeRecoCandidates         = (TTree*) NtplFileRecoCandidates->Get("B0KstMuMu/B0KstMuMuNTuple");
    NTupleRecoCandidates          = new B0KstMuMuSingleCandTreeContent();
    NTupleRecoCandidates->Init();

    TFile* NtplFileSingleCand =  TFile::Open(fileNameSingleCand.c_str(), "READ");
    theTreeSingleCand         = (TTree*) NtplFileSingleCand->Get("B0KstMuMu/B0KstMuMuNTuple");
    NTupleSingleCand          = new B0KstMuMuSingleCandTreeContent();
    NTupleSingleCand->Init();

    Utility = new Utils(RIGHTtag);
    Utility->ReadPreselectionCut(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str());
    Utility->ReadSelectionCuts(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str());
    Utility->ReadGenericParam(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str());
    loadBinning();

    // Create a new file and close it
    //TFile* fileOut = new TFile(fileNameOutput.c_str(), "RECREATE");
    //fileOut->Close();

    // do all q2Bins at once
    if (q2BinIndx==0)
      for (q2BinIndx=1; q2BinIndx<10; ++q2BinIndx) {
        computeEffForOneBin(q2BinIndx, doPlot, width);
    }
    else if (q2BinIndx>0 && q2BinIndx<10) {
      computeEffForOneBin(q2BinIndx, doPlot, width);
    }
    else {
      cout<<"q2Bin must be greater than 0 and smaller than 10. FAILURE!"<<endl;
      delete Utility;
      if (doBatch == false) theApp.Run (); // Eventloop on air
      return EXIT_FAILURE;
    }

    cout<<Vector [0]<<"\t"<<Vector [1]<<"\t"<<Vector [2]<<"\t"<<Vector [3]<<endl;
    cout<<Counter[0]<<"\t"<<Counter[1]<<"\t"<<Counter[2]<<"\t"<<Counter[3]<<endl;

    delete Utility;
    if (doBatch == false) theApp.Run (); // Eventloop on air
    return EXIT_SUCCESS;
  }
  else
  {
    cout << "Parameter missing: " << endl;
    cout << "./" << argv[0] << " <q2Bin (0 for all)> " << " <plot histo (decfault=0/no)> " << " <batch running (default=0/no) " << endl;

    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
