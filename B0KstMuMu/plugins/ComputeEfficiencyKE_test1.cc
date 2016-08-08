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


string fileNameGenCandidatesNoFilter;
string fileNameRecoCandidates;
string fileNameSingleCand;

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

void ComputeEfficiency (TTree* theTree, B0KstMuMuSingleCandTreeContent* NTuple, double q2BinIdx, unsigned int type, float width, int SignalType, bool mist, int test)
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
  if (SignalType==3 && type != 4 && q2BinIdx==5) nEntries = nEntries/5;
  //nEntries = nEntries/1000;
  cout << "[ComputeEfficiency::ComputeEfficiency]\t@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;

  RooDataSet KEpdf_data("KEpdf_data", "KEpdf_data", RooArgSet(CosThetaK,AbsCosThetaMu,AbsPhiKstMuMuPlane));

  int entryMin = 0;

  //if (type<3) {
  entryMin = (test * nEntries)/4;
  nEntries = ((test+1) * nEntries)/4;
  //}

  for (int entry = entryMin; entry < nEntries; entry++) {
    theTree->GetEntry(entry);

    //if (type>2 && NTuple->eventN%4 != test) continue;

    if ((NTuple->B0pT > Utility->GetSeleCut("B0pT")) && (fabs(NTuple->B0Eta) < Utility->GetSeleCut("B0Eta")) &&

        ((NTuple->genSignal == SignalType || NTuple->genSignal == SignalType+1)) &&

        ((type == 1 || type == 3) ||

         ((type == 2) &&
          (sqrt(NTuple->genMumPx*NTuple->genMumPx + NTuple->genMumPy*NTuple->genMumPy)   > Utility->GetPreCut("MinMupT")) &&
          (sqrt(NTuple->genMupPx*NTuple->genMupPx + NTuple->genMupPy*NTuple->genMupPy)   > Utility->GetPreCut("MinMupT")) &&
          (fabs(Utility->computeEta(NTuple->genMumPx,NTuple->genMumPy,NTuple->genMumPz)) < Utility->GetPreCut("MuEta"))   &&
          (fabs(Utility->computeEta(NTuple->genMupPx,NTuple->genMupPy,NTuple->genMupPz)) < Utility->GetPreCut("MuEta")))  ||

         ((type == 4) && (NTuple->truthMatchSignal->at(0) == true) && (NTuple->rightFlavorTag != mist)          &&
          (NTuple->B0MassArb > Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str()))  &&
          (NTuple->B0MassArb < Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str())) &&

          ((((SignalType == Utility->B0ToKstMuMu)  || (SignalType == Utility->B0ToKstMuMu+1))  && (Utility->PsiRejection(NTuple->B0MassArb,NTuple->mumuMass->at(0),NTuple->mumuMassE->at(0),"rejectPsi",true) == true)) ||
           (((SignalType == Utility->B0ToJPsiKst)  || (SignalType == Utility->B0ToJPsiKst+1))  && (Utility->PsiRejection(NTuple->B0MassArb,NTuple->mumuMass->at(0),NTuple->mumuMassE->at(0),"keepJpsi")       == true)) ||
           (((SignalType == Utility->B0ToPsi2SKst) || (SignalType == Utility->B0ToPsi2SKst+1)) && (Utility->PsiRejection(NTuple->B0MassArb,NTuple->mumuMass->at(0),NTuple->mumuMassE->at(0),"keepPsiP")       == true))))))
    {
      mumuq2 = NTuple->mumuMass->at(0)*NTuple->mumuMass->at(0);

      if (mumuq2 < q2Bins[q2BinIdx-1] || mumuq2 > q2Bins[q2BinIdx]) continue;

      if ((type == 4) && mist)
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
  ele_KEpdf->at(type-1) = new RooNDKeysPdf(Form("ele_KEpdf%i",type-1), Form("ele_KEpdf%i",type-1), RooArgList(CosThetaK,AbsCosThetaMu,AbsPhiKstMuMuPlane), KEpdf_data, "M", width);
  //ele_KEpdf->at(type-1) = new RooNDKeysPdf(Form("ele_KEpdf%i",type-1), Form("ele_KEpdf%i",type-1), RooArgList(ThetaK,ThetaMu, PhiKstMuMuPlane), KEpdf_data, "m");
}

void computeEffForOneBin(int q2BinIndx, float width, int SignalType, bool mist, int test) {
  cout << "computeEffForOneBin " << q2BinIndx << endl;
  // Open input file
  ele_KEpdf = new std::vector<RooNDKeysPdf*> (4);
  string fileNameOutput = Form("effKEpdf_b%i.root",q2BinIndx);
  TFile* fileOut = new TFile(fileNameOutput.c_str(), "RECREATE");

  float fac = 1.;
  if (SignalType==3 && q2BinIndx==5) fac = 5.;

  int nbin=40;
  fileOut->cd();

  // denominator GEN
  if (termIdx == 1 || termIdx == 0) {
    cout << "Gen den " ;
    ComputeEfficiency(theTreeGenCandidatesNoFilter,NTupleGenCandidatesNoFilter,q2BinIndx, 1, width, SignalType, mist, test);
    TH3* h3genDen = (TH3*)(ele_KEpdf->at(0))->createHistogram("h3genDen",CosThetaK,Binning(nbin),YVar(AbsCosThetaMu,Binning(nbin)),ZVar(AbsPhiKstMuMuPlane,Binning(nbin)));
    h3genDen->Scale(Counter[0]*fac/h3genDen->Integral());
    h3genDen->Write("h3genDen");
  }

  // numerator GEN
  if (termIdx == 2 || termIdx == 0) {
    cout << "Gen num " ;
    ComputeEfficiency(theTreeGenCandidatesNoFilter,NTupleGenCandidatesNoFilter,q2BinIndx, 2, width, SignalType, mist, test);
    TH3* h3genNum = (TH3*)(*ele_KEpdf->at(1)).createHistogram("h3genNum",CosThetaK,Binning(nbin),YVar(AbsCosThetaMu,Binning(nbin)),ZVar(AbsPhiKstMuMuPlane,Binning(nbin)));
    h3genNum->Scale(Counter[1]*fac/h3genNum->Integral());
    h3genNum->Write("h3genNum");
  }

  // denominator RECO
  if (termIdx == 3 || termIdx == 0) {
    cout << "Reco den " ;
    ComputeEfficiency(theTreeRecoCandidates,       NTupleRecoCandidates,       q2BinIndx, 3, width, SignalType, mist, test);
    TH3* h3recoDen = (TH3*)(ele_KEpdf->at(2))->createHistogram("h3recoDen",CosThetaK,Binning(nbin),YVar(AbsCosThetaMu,Binning(nbin)),ZVar(AbsPhiKstMuMuPlane,Binning(nbin)));
    h3recoDen->Scale(Counter[2]*fac/h3recoDen->Integral());
    h3recoDen->Write("h3recoDen");
  }

  // numerator RECO
  if (termIdx == 4 || termIdx == 0) {
    cout << "Reco num " ;
    ComputeEfficiency(theTreeSingleCand,           NTupleSingleCand,           q2BinIndx, 4, width, SignalType, mist, test);
    TH3* h3recoNum = (TH3*)(ele_KEpdf->at(3))->createHistogram("h3recoNum",CosThetaK,Binning(nbin),YVar(AbsCosThetaMu,Binning(nbin)),ZVar(AbsPhiKstMuMuPlane,Binning(nbin)));
    h3recoNum->Scale(Counter[3]*1./h3recoNum->Integral());
    h3recoNum->Write("h3recoNum");
  }

  fileOut->Close();
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
    bool doBatch=false;
    int SignalType=1;
    bool mist=false;
    int test=0;
    if ( argc > 1 ) q2BinIndx = atoi(argv[1]);
    if ( argc > 2 ) doBatch = atoi(argv[2]);
    if ( argc > 3 ) termIdx = atoi(argv[3]);
    if ( argc > 4 ) SignalType = atoi(argv[4]);
    if ( argc > 5 ) mist = atoi(argv[5]);
    if ( argc > 6 ) test = atoi(argv[6]);

    if (doBatch) {
      cout << "\n[" << argv[0] << "::main]\t@@@ Setting batch mode @@@" << endl;
      gROOT->SetBatch(doBatch);
    }

    float width=0.5;
    if (mist) {
      width=1.0;
      termIdx=4;
    }

    cout << "computing: bin " << q2BinIndx << " term " << termIdx << " mist " << mist << " test " << test << endl;

    TString sample_name;
    if (SignalType<3) sample_name = "B0ToKstMuMu";
    else if (SignalType<5) sample_name = "B0ToJPsiKst";
    else sample_name = "B0ToPsi2SKst";

    fileNameGenCandidatesNoFilter = "dcap://t2-srm-02.lnl.infn.it//pnfs/lnl.infn.it/data/cms/store/user/slacapra/B0KsMuMu/Data2012B0KstMuMuResults/MonteCarlo2012/ForEfficiency/"+sample_name+"_GEN_NoFilter_MC_NTuples_addGENvars.root";
    fileNameRecoCandidates        = "dcap://t2-srm-02.lnl.infn.it////pnfs/lnl.infn.it/data/cms/store/user/slacapra/B0KsMuMu/Data2012B0KstMuMuResults/MonteCarlo2012/ForEfficiency/"+sample_name+"_MC_NTuple_addGENvars.root";
    fileNameSingleCand            = "dcap://t2-srm-02.lnl.infn.it////pnfs/lnl.infn.it/data/cms/store/user/slacapra/B0KsMuMu/Data2012B0KstMuMuResults/MonteCarlo2012/SingleCand/singleCand_"+sample_name+"_MC_NTuple.root";

    TFile* NtplFileGenCandidatesNoFilter =  TFile::Open(fileNameGenCandidatesNoFilter.c_str(), "READ");
    theTreeGenCandidatesNoFilter         = (TTree*) NtplFileGenCandidatesNoFilter->Get("B0KstMuMu/B0KstMuMuNTuple");
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

    // do all q2Bins at once
    if (q2BinIndx==0)
      for (q2BinIndx=1; q2BinIndx<10; ++q2BinIndx) {
        computeEffForOneBin(q2BinIndx, width, SignalType, mist, test);
    }
    else if (q2BinIndx>0 && q2BinIndx<10) {
      computeEffForOneBin(q2BinIndx, width, SignalType, mist, test);
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
