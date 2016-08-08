#include <RooGenericPdf.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooProdPdf.h"
#include "RooWorkspace.h"
#include "RooPolynomial.h"
#include "RooExponential.h"
#include "RooRandom.h"

using namespace RooFit ;

#include <TROOT.h>
#include "TFile.h"
#include "TMath.h"

#include "Utils.h"
#include "B0KstMuMuSingleCandTreeContent.h"

#define PARAMETERFILEIN "/python/ParameterFile.txt"
Utils* Utility;

string fileNameOutput = "toyDatasets.root";

float coeff_ctl_p1[9] = {-0.94024, -0.999251, -0.746712, 1.31074, 0, 1.53294, 0, 1.3881, 0.283786};
float coeff_ctl_p3[9] = {0, 0, 0, -2.11336, 0, -1.82117, 0, -1.535, 0};

float coeff_ctk_p0[9] = {0.462527, -0.885392, -1.3381, -0.550072, 0, -1.06026, 0, -0.145155, 0.218614};
float coeff_ctk_p1[9] = {1.58702, 0.243571, 1.41569, 0.589272, 0, 1.10164, 0, -0.854784, -0.777978};
float coeff_ctk_p2[9] = {-2.07989, 0.148541, 1.21027, 0, 0, 0.597256, 0, 0, 0};
float coeff_ctk_p3[9] = {0, 0.866815, -0.812571, 0, 0, -1.35426, 0, 0, 0};

float coeff_phi_p0[9] = {-1.4662e-01, -1.4662e-01, -1.4662e-01, -1.4662e-01, 0, -1.4662e-01, 0, -1.4662e-01, -1.4662e-01};

float coeff_mass_tau[9] = {0.210102,0.291326,0.248466,0.469477,0,0.224728,0,0.193928,0.512508};

RooRealVar ctK("ctK","cos#theta_{K}",-1,1);
RooRealVar ctL("ctL","cos#theta_{L}",0,1);
RooRealVar phi("phi","#phi",0,TMath::Pi());

RooRealVar B0mass("B0mass","B_{0} mass",5.27958-0.280,5.27958+0.280);

int sigYield[9] = {84,145,119,225,0,361,0,222,239};
int bkgYield[9] = {91,289,216,343,0,567,0,178,82};

vector<double> q2Bins;

void loadBinning() {
  Utility->Readq2Bins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins);
}

void GenerateSig ( int toyIndx, unsigned int q2BinIndx )
{
  TFile* NtplFileIn = TFile::Open("dcap://t2-srm-02.lnl.infn.it//pnfs/lnl.infn.it/data/cms/store/user/slacapra/B0KsMuMu/Data2012B0KstMuMuResults/MonteCarlo2012/SingleCand/singleCand_B0ToKstMuMu_MC_NTuple.root", "READ");
  TTree* theTreeIn  = (TTree*) NtplFileIn->Get("B0KstMuMu/B0KstMuMuNTuple");
  B0KstMuMuSingleCandTreeContent* NTupleIn = new B0KstMuMuSingleCandTreeContent();
  NTupleIn->Init();
  NTupleIn->ClearNTuple();
  NTupleIn->SetBranchAddresses(theTreeIn);
  int nEntries = theTreeIn->GetEntries();

  //RooDataSet* dataSe  = new RooDataSet("a","a",RooArgSet(ctK));
  RooDataSet* dataSet = new RooDataSet(Form("toy%i_sig%i",toyIndx,q2BinIndx), Form("Signal dataset for toy %i (bin %i)",toyIndx,q2BinIndx), RooArgSet(ctK,ctL,phi,B0mass));

  int cntPart = 0;
  int evPart = 0;
  bool filled = false;
  for (int entry = 0; entry < nEntries; entry++) {
    theTreeIn->GetEntry(entry);

    double mumuq2 = NTupleIn->mumuMass->at(0)*NTupleIn->mumuMass->at(0);
    if (mumuq2 < q2Bins[q2BinIndx-1] || mumuq2 > q2Bins[q2BinIndx]) continue;

    cntPart++;
    if (cntPart > sigYield[q2BinIndx-1]) {
      cntPart = 1;
      evPart++;
    }      

    if ( evPart < toyIndx ) continue;
    if ( evPart > toyIndx ) break;

    if (NTupleIn->rightFlavorTag == false) ctK.setVal( - NTupleIn->CosThetaKArb );
    else ctK.setVal( NTupleIn->CosThetaKArb );
    ctL.setVal( fabs(NTupleIn->CosThetaMuArb) );
    phi.setVal( fabs(NTupleIn->PhiKstMuMuPlaneArb) );

    B0mass.setVal( NTupleIn->B0MassArb );

    dataSet->add(RooArgSet(ctK,ctL,phi,B0mass));

    if (cntPart == sigYield[q2BinIndx-1]) filled = true;
  }

  if (!filled) return;
  
  TFile* fout = new TFile(fileNameOutput.c_str(), "UPDATE");

  // RooWorkspace* ws_s = new RooWorkspace(Form("ws%i_s%i",toyIndx,q2BinIndx));
  // ws_s->import(*dataSet);
  // ws_s->Write(0,TObject::kOverwrite);
  
  RooWorkspace* ws = (RooWorkspace*)fout->Get("ws");
  if ( ws==0 ) ws = new RooWorkspace("ws");
  ws->import(*dataSet);
  ws->Write(0,TObject::kOverwrite);
  
  fout->Close();
  delete fout;
  delete dataSet;
}

void GenerateBkg ( int toyIndx, unsigned int q2BinIndx )
{
  RooRealVar ctk_p0(Form("ctk%i_p0",q2BinIndx),"ctk_p0",coeff_ctk_p0[q2BinIndx-1]);
  RooRealVar ctk_p1(Form("ctk%i_p1",q2BinIndx),"ctk_p1",coeff_ctk_p1[q2BinIndx-1]);
  RooRealVar ctk_p2(Form("ctk%i_p2",q2BinIndx),"ctk_p2",coeff_ctk_p2[q2BinIndx-1]);
  RooRealVar ctk_p3(Form("ctk%i_p3",q2BinIndx),"ctk_p3",coeff_ctk_p3[q2BinIndx-1]);
  
  RooRealVar ctl_p0(Form("ctl%i_p0",q2BinIndx),"ctl_p0",0);
  RooRealVar ctl_p1(Form("ctl%i_p1",q2BinIndx),"ctl_p1",coeff_ctl_p1[q2BinIndx-1]);
  RooRealVar ctl_p2(Form("ctl%i_p2",q2BinIndx),"ctl_p2",0);
  RooRealVar ctl_p3(Form("ctl%i_p3",q2BinIndx),"ctl_p3",coeff_ctl_p3[q2BinIndx-1]);
  
  RooRealVar phi_p0(Form("phi%i_p0",q2BinIndx),"phi_p0",coeff_phi_p0[q2BinIndx-1]);

  RooRealVar mass_c(Form("mass%i_c",q2BinIndx),"mass_c",-1./coeff_mass_tau[q2BinIndx-1]);

  RooPolynomial ctk_pdf(Form("ctk_pdf%i",q2BinIndx),"ctk_pdf",ctK,RooArgList(ctk_p0,ctk_p1,ctk_p2,ctk_p3));
  RooPolynomial ctl_pdf(Form("ctl_pdf%i",q2BinIndx),"ctl_pdf",ctL,RooArgList(ctl_p0,ctl_p1,ctl_p2,ctl_p3));
  RooPolynomial phi_pdf(Form("phi_pdf%i",q2BinIndx),"phi_pdf",phi,RooArgList(phi_p0));
  RooExponential m_pdf(Form("m_pdf%i",q2BinIndx),"m_pdf",B0mass,mass_c);

  RooProdPdf bkg_pdf(Form("pdf_bkg%i",q2BinIndx),Form("Background PDF (bin %i)",q2BinIndx),RooArgSet(ctk_pdf,ctl_pdf,phi_pdf,m_pdf));

  RooRandom::randomGenerator()->SetSeed(toyIndx+1);

  RooDataSet* dataSet = bkg_pdf.generate(RooArgSet(ctK,ctL,phi,B0mass),bkgYield[q2BinIndx-1]);
  dataSet->SetName(Form("toy%i_bkg%i",toyIndx,q2BinIndx));
  dataSet->SetTitle(Form("Background dataset for toy %i (bin %i)",toyIndx,q2BinIndx));
  
  TFile* fout = new TFile(fileNameOutput.c_str(), "UPDATE");
  
  // RooWorkspace* ws_b = new RooWorkspace(Form("ws%i_b%i",toyIndx,q2BinIndx)) ;
  // ws_b->import(*dataSet);
  // ws_b->Write(0,TObject::kOverwrite);
  
  RooWorkspace* ws = (RooWorkspace*)fout->Get("ws");
  if ( ws==0 ) ws = new RooWorkspace("ws");
  ws->import(bkg_pdf);
  ws->import(*dataSet);
  ws->Write(0,TObject::kOverwrite);
  
  fout->Close();
  delete fout;
  delete dataSet;
}

int main(int argc, char** argv)
{
  if (argc > 0)
  {
    int toyIndx=0;
    int q2BinIndx=0;
    int samp=0;
    if ( argc > 1 ) toyIndx = atoi(argv[1]);
    if ( argc > 2 ) q2BinIndx = atoi(argv[2]);
    if ( argc > 3 ) samp = atoi(argv[3]);

    Utility = new Utils(false);
    loadBinning();

    if (toyIndx<0) {
      cout<<"toy index must be greater than 0. FAILURE!"<<endl;
      return EXIT_FAILURE;
    }

    // do all q2Bins at once
    if (q2BinIndx==0)
      for (q2BinIndx=1; q2BinIndx<10; ++q2BinIndx) {
	if (q2BinIndx==5 || q2BinIndx==7) continue;
        if (samp==0) {
	  GenerateSig(toyIndx, q2BinIndx);
	  GenerateBkg(toyIndx, q2BinIndx);
	}
	else if (samp==1)
	  GenerateSig(toyIndx, q2BinIndx);
	else if (samp==2)
	  GenerateBkg(toyIndx, q2BinIndx);
	else {
	  cout<<"Sample Index must be 1:sig, 2:bkg or 0:all"<<endl;
	  return EXIT_FAILURE;
	}
      }
    else if (q2BinIndx>0 && q2BinIndx<10 && q2BinIndx!=5 && q2BinIndx!=7) {
        if (samp==0) {
	  GenerateSig(toyIndx, q2BinIndx);
	  GenerateBkg(toyIndx, q2BinIndx);
	}
	else if (samp==1)
	  GenerateSig(toyIndx, q2BinIndx);
	else if (samp==2)
	  GenerateBkg(toyIndx, q2BinIndx);
	else {
	  cout<<"Sample Index must be 1:sig, 2:bkg or 0:all"<<endl;
	  return EXIT_FAILURE;
	}
    }
    else {
      cout<<"q2Bin must be greater than 0 and smaller than 10 and no 5 nor 7. FAILURE!"<<endl;
      return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
  }
  else
  {
    cout << "Parameter missing: " << endl;
    cout << "./" << argv[0] << " <toy index> " << " <q2Bin, 0:all (default 0)> " << " <1:sig 2:bkg 0:all (default 0)>" << endl;

    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
