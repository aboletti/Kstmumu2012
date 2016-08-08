// ################################################################################
// # Program to perform the a closure test for efficiency for angular analysis
// of the decay B0 --> K*0 mu+ mu- 
// Author: Stefano Lacaprara  <lacaprara@pd.infn.it>  INFN Padova
// #
// ################################################################################
// # Search for @TMP@ to look for temporary code options                          #
// ################################################################################
#include "Utils.h"

#include <RooGenericPdf.h>
#include "RooHistPdf.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooUniform.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "RooBinning.h"
using namespace RooFit ;

#include <TROOT.h>
#include <TApplication.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"
#include <TLorentzVector.h>

// ####################
// # Global constants #
// ####################
#define SETBATCH        false
#define PARAMETERFILEIN "/python/ParameterFile.txt"
#define ordinateRange   1e-2

#define doInterpolation 3

// ####################
// # Global variables #
// ####################
Utils* Utility;

vector<double> q2Bins;
double* q2Bins_ ;
double* cosThetaKBins_ ;
double* cosThetaLBins_ ;
double* phiBins_ ;
RooBinning cosThetaLBinning;
RooBinning cosThetaKBinning;
RooBinning phiBinning;
RooRealVar ctK("ctK","cos#theta_{K}",-1.,1.);
RooRealVar ctL("ctL","cos#theta_{L}",-1.,1.);
RooRealVar tK("tK","#theta_{K}",0,TMath::Pi());
RooRealVar tL("tL","#theta_{#mu}",0,TMath::Pi()) ;
RooRealVar phi("phi","#phi",-1.*TMath::Pi(),1.*TMath::Pi());
RooArgSet ctKctL(ctK,ctL);
RooArgSet ctKphi(ctK,phi);
RooArgSet ctLphi(ctL,phi);
RooArgSet ctKctLphi(ctK,ctL,phi);
RooArgSet tKtL(tK,tL);
RooArgSet tKphi(tK,phi);
RooArgSet tLphi(tL,phi);
RooArgSet tKtLphi(tK,tL,phi);

TFile* pdfEffHistoFile3d=0;
TFile* histEffHistoFile3d=0;
TH2D* H2Deff_ctL_phi;
TH2D* H2Deff_ctK_phi;
TH3D* H3Deff_ctK_ctL_phi;

B0KstMuMuSingleCandTreeContent* NTupleIn_reco;
B0KstMuMuSingleCandTreeContent* NTupleIn_gen;

static const int qbins=8;
TH3F* hist_ctKctLphi_q2bin2[qbins];
RooHistPdf* pdf_ctKctLphi_q2bin2[qbins];
TH2F* h_reco_ctk_ctl_phi1[qbins];
TH2F* h_reco_ctk_ctl_phi2[qbins];
TH2F* h_reco_ctk_ctl_phi3[qbins];
TH2F* h_reco_ctk_ctl_phi4[qbins];
TH2F* h_gen_ctk_ctl_phi1[qbins];
TH2F* h_gen_ctk_ctl_phi2[qbins];
TH2F* h_gen_ctk_ctl_phi3[qbins];
TH2F* h_gen_ctk_ctl_phi4[qbins];
TH2F* h_diff_ctk_ctl_phi1[qbins];
TH2F* h_diff_ctk_ctl_phi2[qbins];
TH2F* h_diff_ctk_ctl_phi3[qbins];
TH2F* h_diff_ctk_ctl_phi4[qbins];
TH1F* h_reco_ctk[qbins];
TH1F* h_gen_ctk[qbins];
TH1F* h_reco_ctl[qbins];
TH1F* h_gen_ctl[qbins];
TH1F* h_reco_phi[qbins];
TH1F* h_gen_phi[qbins];
TCanvas* c_ctest[qbins];
TCanvas* c1_ctest[qbins];

double nweight=0;
int gen_cnt = 0;
int reco_cnt = 0;

int ser_num = 8;

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

void loadEffHisto(unsigned int q2bin) {
  if (!pdfEffHistoFile3d) { 
    cout << "Loading effKEpdf_out.root" << endl;
    pdfEffHistoFile3d=TFile::Open("effKEpdf_out.root","READ");
    pdfEffHistoFile3d->ls();
  }
  cout << "Try to get " << Form("pdf_ctKctLphi_q2bin%i",q2bin) ;
  pdf_ctKctLphi_q2bin2[q2bin-1]=(RooHistPdf*) pdfEffHistoFile3d->Get(Form("pdf_ctKctLphi_q2bin%d",q2bin));
  cout << " :pdf_ctKctLphi_q2bin2[" << (q2bin-1) << "]=" << pdf_ctKctLphi_q2bin2[q2bin-1] << endl;

  /* if (!histEffHistoFile3d) { 
    cout << "Loading eff3DOutputFile.root" << endl;
    histEffHistoFile3d=TFile::Open("effKEpdf.root","READ");
    histEffHistoFile3d->ls();
  }
  cout << "Try to get " << Form("h3Eff_q2bin%d",q2bin) ;
  hist_ctKctLphi_q2bin2[q2bin-1]=(TH3F*) histEffHistoFile3d->Get(Form("h3Eff_q2bin%d",q2bin));
  cout << " :hist_ctKctLphi_q2bin2[" << (q2bin-1) << "]=" << hist_ctKctLphi_q2bin2[q2bin-1] << endl;*/
}

void loadBinning() {
  vector<double> cosThetaKBins;
  vector<double> cosThetaLBins;
  vector<double> phiBins;
  Utility->ReadAllBins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"goodTag");
  q2Bins_        = Utility->MakeBinning(&q2Bins);
  cosThetaKBins_ = Utility->MakeBinning(&cosThetaKBins);
  cosThetaLBins_ = Utility->MakeBinning(&cosThetaLBins);
  phiBins_       = Utility->MakeBinning(&phiBins);
  cosThetaLBinning=RooBinning (cosThetaLBins.size()-1,cosThetaLBins_);
  cosThetaLBinning.Print();
  cosThetaKBinning=RooBinning (cosThetaKBins.size()-1,cosThetaKBins_);
  cosThetaKBinning.Print();
  phiBinning=RooBinning(phiBins.size()-1,phiBins_);
  phiBinning.Print();

}

void initHisto(int q2bin) {
  h_reco_ctk_ctl_phi1[q2bin-1] = new TH2F(Form("h_reco_ctk_ctl_phi1_%i",q2bin-1),Form("h_reco_ctk_ctl_phi1_%i",q2bin-1),10,-1,1,10,0,1);
  h_reco_ctk_ctl_phi2[q2bin-1] = new TH2F(Form("h_reco_ctk_ctl_phi2_%i",q2bin-1),Form("h_reco_ctk_ctl_phi2_%i",q2bin-1),10,-1,1,10,0,1);
  h_reco_ctk_ctl_phi3[q2bin-1] = new TH2F(Form("h_reco_ctk_ctl_phi3_%i",q2bin-1),Form("h_reco_ctk_ctl_phi3_%i",q2bin-1),10,-1,1,10,0,1);
  h_reco_ctk_ctl_phi4[q2bin-1] = new TH2F(Form("h_reco_ctk_ctl_phi4_%i",q2bin-1),Form("h_reco_ctk_ctl_phi4_%i",q2bin-1),10,-1,1,10,0,1);
  h_gen_ctk_ctl_phi1[q2bin-1] = new TH2F(Form("h_gen_ctk_ctl_phi1_%i",q2bin-1),Form("h_gen_ctk_ctl_phi1_%i",q2bin-1),10,-1,1,10,0,1);
  h_gen_ctk_ctl_phi2[q2bin-1] = new TH2F(Form("h_gen_ctk_ctl_phi2_%i",q2bin-1),Form("h_gen_ctk_ctl_phi2_%i",q2bin-1),10,-1,1,10,0,1);
  h_gen_ctk_ctl_phi3[q2bin-1] = new TH2F(Form("h_gen_ctk_ctl_phi3_%i",q2bin-1),Form("h_gen_ctk_ctl_phi3_%i",q2bin-1),10,-1,1,10,0,1);
  h_gen_ctk_ctl_phi4[q2bin-1] = new TH2F(Form("h_gen_ctk_ctl_phi4_%i",q2bin-1),Form("h_gen_ctk_ctl_phi4_%i",q2bin-1),10,-1,1,10,0,1);
  h_diff_ctk_ctl_phi1[q2bin-1] = new TH2F(Form("h_diff_ctk_ctl_phi1_%i",q2bin-1),Form("h_diff_ctk_ctl_phi1_%i;#theta_{K};#theta_{L}",q2bin-1),10,-1,1,10,0,1);
  h_diff_ctk_ctl_phi2[q2bin-1] = new TH2F(Form("h_diff_ctk_ctl_phi2_%i",q2bin-1),Form("h_diff_ctk_ctl_phi2_%i;#theta_{K};#theta_{L}",q2bin-1),10,-1,1,10,0,1);
  h_diff_ctk_ctl_phi3[q2bin-1] = new TH2F(Form("h_diff_ctk_ctl_phi3_%i",q2bin-1),Form("h_diff_ctk_ctl_phi3_%i;#theta_{K};#theta_{L}",q2bin-1),10,-1,1,10,0,1);
  h_diff_ctk_ctl_phi4[q2bin-1] = new TH2F(Form("h_diff_ctk_ctl_phi4_%i",q2bin-1),Form("h_diff_ctk_ctl_phi4_%i;#theta_{K};#theta_{L}",q2bin-1),10,-1,1,10,0,1);
  h_gen_ctk [q2bin-1] = new TH1F(Form("h_gen_ctk%i",q2bin-1) ,Form("h_gen_ctk%i;#theta_{K}",q2bin-1) ,100,-1,1);
  h_reco_ctk[q2bin-1] = new TH1F(Form("h_reco_ctk%i",q2bin-1),Form("h_reco_ctk%i;#theta_{K}",q2bin-1),100,-1,1);
  h_gen_ctl [q2bin-1] = new TH1F(Form("h_gen_ctl%i",q2bin-1) ,Form("h_gen_ctl%i;#theta_{L}",q2bin-1) ,100,0,1);
  h_reco_ctl[q2bin-1] = new TH1F(Form("h_reco_ctl%i",q2bin-1),Form("h_reco_ctl%i;#theta_{L}",q2bin-1),100,0,1);
  h_gen_phi [q2bin-1] = new TH1F(Form("h_gen_phi%i",q2bin-1) ,Form("h_gen_phi%i;#phi",q2bin-1)       ,100,0,TMath::Pi());
  h_reco_phi[q2bin-1] = new TH1F(Form("h_reco_phi%i",q2bin-1),Form("h_reco_phi%i;#phi",q2bin-1)      ,100,0,TMath::Pi());
}

void fillHistoGen(int q2bin) {
  //cout << NTupleIn_reco->CosThetaKArb << " " <<  NTupleIn_reco->CosThetaMuArb << " " << NTupleIn_reco->PhiKstMuMuPlaneArb << endl;

  if ( (NTupleIn_gen->B0pT < Utility->GetSeleCut("B0pT")) || (fabs(NTupleIn_gen->B0Eta) > Utility->GetSeleCut("B0Eta")) ) return;

  RooArgSet* arg=pdf_ctKctLphi_q2bin2[q2bin-1]->getVariables();
  arg->setRealValue("ctK",NTupleIn_gen->CosThetaKArb);
  arg->setRealValue("ctL",TMath::Abs(NTupleIn_gen->CosThetaMuArb));
  arg->setRealValue("phi",TMath::Abs(NTupleIn_gen->PhiKstMuMuPlaneArb));

  // // cout << "pdf_ctKctLphi_q2bin2[" << (q2bin-1) << "]=" << pdf_ctKctLphi_q2bin2[q2bin-1] << endl;
  // // PrintVariables(pdf_ctKctLphi_q2bin2[q2bin-1]->getVariables(),"vars");

  double weight = pdf_ctKctLphi_q2bin2[q2bin-1]->getVal();

  // using the histogramA
  //int binX=hist_ctKctLphi_q2bin2[q2bin-1]->GetXaxis()->FindBin(TMath::ACos(NTupleIn_gen->CosThetaKArb));
  //int binY=hist_ctKctLphi_q2bin2[q2bin-1]->GetYaxis()->FindBin(TMath::ACos(NTupleIn_gen->CosThetaMuArb));
  //int binZ=hist_ctKctLphi_q2bin2[q2bin-1]->GetZaxis()->FindBin(NTupleIn_gen->PhiKstMuMuPlaneArb);
  //cout << "bins: " << binX << " " << binY << " " << binZ << endl;
  //weight = hist_ctKctLphi_q2bin2[q2bin-1]->GetBinContent(binX,binY,binZ);
  //cout << "weight " << weight<< endl;

  //double weight=1;
  // cout << "weight " << weight << endl;
  if      (TMath::Abs(NTupleIn_gen->PhiKstMuMuPlaneArb) < 0.25*TMath::Pi()) h_gen_ctk_ctl_phi1[q2bin-1]->Fill(NTupleIn_gen->CosThetaKArb,TMath::Abs(NTupleIn_gen->CosThetaMuArb), weight*nweight);
  else if (TMath::Abs(NTupleIn_gen->PhiKstMuMuPlaneArb) < 0.5*TMath::Pi())  h_gen_ctk_ctl_phi2[q2bin-1]->Fill(NTupleIn_gen->CosThetaKArb,TMath::Abs(NTupleIn_gen->CosThetaMuArb), weight*nweight);
  else if (TMath::Abs(NTupleIn_gen->PhiKstMuMuPlaneArb) < 0.75*TMath::Pi()) h_gen_ctk_ctl_phi3[q2bin-1]->Fill(NTupleIn_gen->CosThetaKArb,TMath::Abs(NTupleIn_gen->CosThetaMuArb), weight*nweight);
  else                                                                      h_gen_ctk_ctl_phi4[q2bin-1]->Fill(NTupleIn_gen->CosThetaKArb,TMath::Abs(NTupleIn_gen->CosThetaMuArb), weight*nweight);
  h_gen_ctk[q2bin-1]->Fill(NTupleIn_gen->CosThetaKArb,                   weight*nweight);
  h_gen_ctl[q2bin-1]->Fill(TMath::Abs(NTupleIn_gen->CosThetaMuArb),      weight*nweight);
  h_gen_phi[q2bin-1]->Fill(TMath::Abs(NTupleIn_gen->PhiKstMuMuPlaneArb), weight*nweight);
  gen_cnt++;
}

void fillHistoReco(int q2bin, bool mist ) {
  if ( (NTupleIn_reco->B0pT < Utility->GetSeleCut("B0pT")) || (fabs(NTupleIn_reco->B0Eta) > Utility->GetSeleCut("B0Eta")) ) return;

  float sign = 1.;
  if (mist) sign = -1.;

  if      (TMath::Abs(NTupleIn_reco->PhiKstMuMuPlaneArb) < 0.25*TMath::Pi()) h_reco_ctk_ctl_phi1[q2bin-1]->Fill(sign*NTupleIn_reco->CosThetaKArb,TMath::Abs(NTupleIn_reco->CosThetaMuArb));
  else if (TMath::Abs(NTupleIn_reco->PhiKstMuMuPlaneArb) < 0.5*TMath::Pi())  h_reco_ctk_ctl_phi2[q2bin-1]->Fill(sign*NTupleIn_reco->CosThetaKArb,TMath::Abs(NTupleIn_reco->CosThetaMuArb));
  else if (TMath::Abs(NTupleIn_reco->PhiKstMuMuPlaneArb) < 0.75*TMath::Pi()) h_reco_ctk_ctl_phi3[q2bin-1]->Fill(sign*NTupleIn_reco->CosThetaKArb,TMath::Abs(NTupleIn_reco->CosThetaMuArb));
  else                                                                       h_reco_ctk_ctl_phi4[q2bin-1]->Fill(sign*NTupleIn_reco->CosThetaKArb,TMath::Abs(NTupleIn_reco->CosThetaMuArb));
  h_reco_ctk[q2bin-1]->Fill(sign*NTupleIn_reco->CosThetaKArb);
  h_reco_ctl[q2bin-1]->Fill(TMath::Abs(NTupleIn_reco->CosThetaMuArb));
  h_reco_phi[q2bin-1]->Fill(TMath::Abs(NTupleIn_reco->PhiKstMuMuPlaneArb));
  reco_cnt++;
}

void plotTestHisto(int q2bin) {

  TFile* fout = new TFile ( Form("closureTest%i.root",q2bin), "RECREATE" );

  float scale_reco = h_reco_ctk[q2bin-1]->Integral();
  h_reco_ctk_ctl_phi1[q2bin-1]->Scale( 1/scale_reco );
  h_reco_ctk_ctl_phi2[q2bin-1]->Scale( 1/scale_reco );
  h_reco_ctk_ctl_phi3[q2bin-1]->Scale( 1/scale_reco );
  h_reco_ctk_ctl_phi4[q2bin-1]->Scale( 1/scale_reco );
  float scale_gen  = h_gen_ctk [q2bin-1]->Integral();
  h_gen_ctk_ctl_phi1[q2bin-1]->Scale( 1/scale_gen );
  h_gen_ctk_ctl_phi2[q2bin-1]->Scale( 1/scale_gen );
  h_gen_ctk_ctl_phi3[q2bin-1]->Scale( 1/scale_gen );
  h_gen_ctk_ctl_phi4[q2bin-1]->Scale( 1/scale_gen );
  
  for (int i=1; i<=h_diff_ctk_ctl_phi1[q2bin-1]->GetXaxis()->GetNbins(); i++) for (int j=1; j<=h_diff_ctk_ctl_phi1[q2bin-1]->GetYaxis()->GetNbins(); j++) {
      h_diff_ctk_ctl_phi1[q2bin-1]->SetBinContent(i, j, (h_gen_ctk_ctl_phi1[q2bin-1]->GetBinContent(i,j)-h_reco_ctk_ctl_phi1[q2bin-1]->GetBinContent(i,j))/((h_gen_ctk_ctl_phi1[q2bin-1]->GetBinContent(i,j)+h_reco_ctk_ctl_phi1[q2bin-1]->GetBinContent(i,j))>0.0001?(h_gen_ctk_ctl_phi1[q2bin-1]->GetBinContent(i,j)+h_reco_ctk_ctl_phi1[q2bin-1]->GetBinContent(i,j)):1));
      h_diff_ctk_ctl_phi2[q2bin-1]->SetBinContent(i, j, (h_gen_ctk_ctl_phi2[q2bin-1]->GetBinContent(i,j)-h_reco_ctk_ctl_phi2[q2bin-1]->GetBinContent(i,j))/((h_gen_ctk_ctl_phi2[q2bin-1]->GetBinContent(i,j)+h_reco_ctk_ctl_phi2[q2bin-1]->GetBinContent(i,j))>0.0001?(h_gen_ctk_ctl_phi2[q2bin-1]->GetBinContent(i,j)+h_reco_ctk_ctl_phi2[q2bin-1]->GetBinContent(i,j)):1));
      h_diff_ctk_ctl_phi3[q2bin-1]->SetBinContent(i, j, (h_gen_ctk_ctl_phi3[q2bin-1]->GetBinContent(i,j)-h_reco_ctk_ctl_phi3[q2bin-1]->GetBinContent(i,j))/((h_gen_ctk_ctl_phi3[q2bin-1]->GetBinContent(i,j)+h_reco_ctk_ctl_phi3[q2bin-1]->GetBinContent(i,j))>0.0001?(h_gen_ctk_ctl_phi3[q2bin-1]->GetBinContent(i,j)+h_reco_ctk_ctl_phi3[q2bin-1]->GetBinContent(i,j)):1));
      h_diff_ctk_ctl_phi4[q2bin-1]->SetBinContent(i, j, (h_gen_ctk_ctl_phi4[q2bin-1]->GetBinContent(i,j)-h_reco_ctk_ctl_phi4[q2bin-1]->GetBinContent(i,j))/((h_gen_ctk_ctl_phi4[q2bin-1]->GetBinContent(i,j)+h_reco_ctk_ctl_phi4[q2bin-1]->GetBinContent(i,j))>0.0001?(h_gen_ctk_ctl_phi4[q2bin-1]->GetBinContent(i,j)+h_reco_ctk_ctl_phi4[q2bin-1]->GetBinContent(i,j)):1));
    }

  c_ctest[q2bin-1] = new TCanvas(Form("c_ctest%i",q2bin-1),Form("c_ctest%i",q2bin-1),800,800);
  c_ctest[q2bin-1]->Divide(2,2);
  c1_ctest[q2bin-1] = new TCanvas(Form("c1_ctest%i",q2bin-1),Form("c1_ctest%i",q2bin-1),1200,400);
  c1_ctest[q2bin-1]->Divide(3,1);

  h_diff_ctk_ctl_phi1[q2bin-1]->GetXaxis()->SetTitleOffset(2.0);
  h_diff_ctk_ctl_phi1[q2bin-1]->GetYaxis()->SetTitleOffset(2.0);
  h_diff_ctk_ctl_phi2[q2bin-1]->GetXaxis()->SetTitleOffset(2.0);
  h_diff_ctk_ctl_phi2[q2bin-1]->GetYaxis()->SetTitleOffset(2.0);
  h_diff_ctk_ctl_phi3[q2bin-1]->GetXaxis()->SetTitleOffset(2.0);
  h_diff_ctk_ctl_phi3[q2bin-1]->GetYaxis()->SetTitleOffset(2.0);
  h_diff_ctk_ctl_phi4[q2bin-1]->GetXaxis()->SetTitleOffset(2.0);
  h_diff_ctk_ctl_phi4[q2bin-1]->GetYaxis()->SetTitleOffset(2.0);
  double max, min;
  max = h_diff_ctk_ctl_phi1[q2bin-1]->GetMaximum();
  min = h_diff_ctk_ctl_phi1[q2bin-1]->GetMinimum();
  if (h_diff_ctk_ctl_phi2[q2bin-1]->GetMaximum() > max) max = h_diff_ctk_ctl_phi2[q2bin-1]->GetMaximum();
  if (h_diff_ctk_ctl_phi2[q2bin-1]->GetMinimum() < min) min = h_diff_ctk_ctl_phi2[q2bin-1]->GetMinimum();
  if (h_diff_ctk_ctl_phi3[q2bin-1]->GetMaximum() > max) max = h_diff_ctk_ctl_phi3[q2bin-1]->GetMaximum();
  if (h_diff_ctk_ctl_phi3[q2bin-1]->GetMinimum() < min) min = h_diff_ctk_ctl_phi3[q2bin-1]->GetMinimum();
  if (h_diff_ctk_ctl_phi4[q2bin-1]->GetMaximum() > max) max = h_diff_ctk_ctl_phi4[q2bin-1]->GetMaximum();
  if (h_diff_ctk_ctl_phi4[q2bin-1]->GetMinimum() < min) min = h_diff_ctk_ctl_phi4[q2bin-1]->GetMinimum();
  max *= 1.1;
  min *= 0.9;
  h_diff_ctk_ctl_phi1[q2bin-1]->SetMaximum(max);
  h_diff_ctk_ctl_phi2[q2bin-1]->SetMaximum(max);
  h_diff_ctk_ctl_phi3[q2bin-1]->SetMaximum(max);
  h_diff_ctk_ctl_phi4[q2bin-1]->SetMaximum(max);
  h_diff_ctk_ctl_phi1[q2bin-1]->SetMinimum(min);
  h_diff_ctk_ctl_phi2[q2bin-1]->SetMinimum(min);
  h_diff_ctk_ctl_phi3[q2bin-1]->SetMinimum(min);
  h_diff_ctk_ctl_phi4[q2bin-1]->SetMinimum(min);

  c_ctest[q2bin-1]->cd(1);
  h_diff_ctk_ctl_phi1[q2bin-1]->DrawCopy("lego2 fp");
  c_ctest[q2bin-1]->cd(2);
  h_diff_ctk_ctl_phi2[q2bin-1]->DrawCopy("lego2 fp");
  c_ctest[q2bin-1]->cd(3);
  h_diff_ctk_ctl_phi3[q2bin-1]->DrawCopy("lego2 fp");
  c_ctest[q2bin-1]->cd(4);
  h_diff_ctk_ctl_phi4[q2bin-1]->DrawCopy("lego2 fp");
  c_ctest[q2bin-1]->SaveAs(Form("closureTest%i_plot_bin%i.pdf",ser_num,q2bin));
  //c_ctest[q2bin-1]->Print();

  h_gen_ctk[q2bin-1]->SetMarkerStyle(7);
  h_gen_ctl[q2bin-1]->SetMarkerStyle(7);
  h_gen_phi[q2bin-1]->SetMarkerStyle(7);
  h_reco_ctk[q2bin-1]->SetMinimum(0);
  h_reco_ctl[q2bin-1]->SetMinimum(0);
  h_reco_phi[q2bin-1]->SetMinimum(0);
  //if (h_gen_ctk[q2bin-1]->GetMaximum() > h_reco_ctk[q2bin-1]->GetMaximum()) h_reco_ctk[q2bin-1]->SetMaximum(1.1 * h_gen_ctk[q2bin-1]->GetMaximum());
  //if (h_gen_ctl[q2bin-1]->GetMaximum() > h_reco_ctl[q2bin-1]->GetMaximum()) h_reco_ctl[q2bin-1]->SetMaximum(1.1 * h_gen_ctl[q2bin-1]->GetMaximum());
  //if (h_gen_phi[q2bin-1]->GetMaximum() > h_reco_phi[q2bin-1]->GetMaximum()) h_reco_phi[q2bin-1]->SetMaximum(1.1 * h_gen_phi[q2bin-1]->GetMaximum());

  /*c1_ctest[q2bin-1]->cd(1);
  h_reco_ctk[q2bin-1]->DrawCopy("");
  h_gen_ctk [q2bin-1]->DrawCopy("same P");
  c1_ctest[q2bin-1]->cd(2);
  h_reco_ctl[q2bin-1]->DrawCopy("");
  h_gen_ctl [q2bin-1]->DrawCopy("same P");
  c1_ctest[q2bin-1]->cd(3);
  h_reco_phi[q2bin-1]->DrawCopy("");
  h_gen_phi [q2bin-1]->DrawCopy("same P");
  c1_ctest[q2bin-1]->SaveAs(Form("closureTest%i_plot1_bin%i.pdf",ser_num,q2bin));*/
  cout << "doid]pwo erifweojhifgo;'12t o[9" << endl;

  c1_ctest[q2bin-1]->cd(1);
  h_reco_ctk[q2bin-1]->DrawNormalized("");
  h_gen_ctk [q2bin-1]->DrawNormalized("same P");
  c1_ctest[q2bin-1]->cd(2);
  h_reco_ctl[q2bin-1]->DrawNormalized("");
  h_gen_ctl [q2bin-1]->DrawNormalized("same P");
  c1_ctest[q2bin-1]->cd(3);
  h_reco_phi[q2bin-1]->DrawNormalized("");
  h_gen_phi [q2bin-1]->DrawNormalized("same P");
  c1_ctest[q2bin-1]->SaveAs(Form("closureTest%i_plot1_bin%i.pdf",ser_num,q2bin));

  fout->cd();
  c1_ctest[q2bin-1]->Write();
  c_ctest[q2bin-1]->Write();
  h_reco_ctk[q2bin-1]->Write();
  h_gen_ctk [q2bin-1]->Write();
  h_reco_ctl[q2bin-1]->Write();
  h_gen_ctl [q2bin-1]->Write();
  h_reco_phi[q2bin-1]->Write();
  h_gen_phi [q2bin-1]->Write();
  fout->Close();
}

void performTest(int q2bin, int nev_reco=0, int nev_gen=0, bool doPlot=true, bool mist=false, int evt_cat=0, bool crossCheck=0) {
  initHisto(q2bin);

  TString name_cat = "B0ToKstMuMu";
  if (evt_cat==1) name_cat = "B0ToJPsiKst";
  else if (evt_cat==2) name_cat = "B0ToPsi2SKst";

  TFile* NtplFileIn_reco = TFile::Open("dcap://t2-srm-02.lnl.infn.it//pnfs/lnl.infn.it/data/cms/store/user/slacapra/B0KsMuMu/Data2012B0KstMuMuResults/MonteCarlo2012/SingleCand/singleCand_"+name_cat+"_MC_NTuple.root", "READ");
  TTree* theTreeIn_reco  = (TTree*) NtplFileIn_reco->Get("B0KstMuMu/B0KstMuMuNTuple");
  NTupleIn_reco = new B0KstMuMuSingleCandTreeContent();
  NTupleIn_reco->Init();
  NTupleIn_reco->ClearNTuple();
  NTupleIn_reco->SetBranchAddresses(theTreeIn_reco);
  int nEntries_reco = theTreeIn_reco->GetEntries();

  cout << "\n[ClosureTest::performTest]\t@@@ Total number of reco events in the tree: " << nEntries_reco << " @@@" << endl;
  if (nev_reco>0 && nev_reco < nEntries_reco) nEntries_reco = nev_reco;
  cout << "\n[ClosureTest::performTest]\t@@@ Number of reco events to be used: " << nEntries_reco << " @@@" << endl;
  
  TFile* NtplFileIn_gen = TFile::Open("dcap://t2-srm-02.lnl.infn.it//pnfs/lnl.infn.it/data/cms/store/user/slacapra/B0KsMuMu/Data2012B0KstMuMuResults/MonteCarlo2012/ForEfficiency/"+name_cat+"_GEN_NoFilter_MC_NTuples_addGENvars.root","READ");
  TTree* theTreeIn_gen  = (TTree*) NtplFileIn_gen->Get("B0KstMuMu/B0KstMuMuNTuple");
  NTupleIn_gen = new B0KstMuMuSingleCandTreeContent();
  NTupleIn_gen->Init();
  NTupleIn_gen->ClearNTuple();
  NTupleIn_gen->SetBranchAddresses(theTreeIn_gen);
  int nEntries_gen = theTreeIn_gen->GetEntries();

  cout << "\n[ClosureTest::performTest]\t@@@ Total number of gen events in the tree: " << nEntries_gen << " @@@" << endl;
  if (nev_gen>0 && nev_gen < nEntries_gen) nEntries_gen = nev_gen;
  cout << "\n[ClosureTest::performTest]\t@@@ Number of gen events to be used: " << nEntries_gen << " @@@" << endl;

  nweight = 1.0*theTreeIn_gen->GetEntries()/theTreeIn_reco->GetEntries() * nEntries_reco/nEntries_gen;

  cout << "\n[ClosureTest::performTest]\t@@@ RECO processing... @@@" << endl;
  for (int entry = 0; entry < nEntries_reco; entry++) {
    theTreeIn_reco->GetEntry(entry);
    if (crossCheck && NTupleIn_reco->eventN%2 == 0) continue;
    double mumuq2 = NTupleIn_reco->mumuMass->at(0)*NTupleIn_reco->mumuMass->at(0);
    if (mumuq2<q2Bins_[q2bin-1] || mumuq2>q2Bins_[q2bin]) continue;
    if ((!NTupleIn_reco->rightFlavorTag && !mist) || (NTupleIn_reco->rightFlavorTag && mist)) continue;
    if (!( (NTupleIn_reco->truthMatchSignal->at(0) == true) &&
	   (NTupleIn_reco->B0MassArb > Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str())) &&
	   (NTupleIn_reco->B0MassArb < Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str())) &&
	   ( ((evt_cat == 0) && (Utility->PsiRejection(NTupleIn_reco->B0MassArb,NTupleIn_reco->mumuMass->at(0),NTupleIn_reco->mumuMassE->at(0),"rejectPsi",true) == true)) ||
	     ((evt_cat == 1) && (Utility->PsiRejection(NTupleIn_reco->B0MassArb,NTupleIn_reco->mumuMass->at(0),NTupleIn_reco->mumuMassE->at(0),"keepJpsi")       == true)) ||
	     ((evt_cat == 2) && (Utility->PsiRejection(NTupleIn_reco->B0MassArb,NTupleIn_reco->mumuMass->at(0),NTupleIn_reco->mumuMassE->at(0),"keepPsiP")       == true))   ) ) ) continue;
    fillHistoReco(q2bin,mist);
  }
  
  cout << "\n[ClosureTest::performTest]\t@@@ GEN processing... @@@" << endl;
  for (int entry = 0; entry < nEntries_gen; entry++) {
    if (crossCheck && entry%2 == 0) continue;
    theTreeIn_gen->GetEntry(entry);
    double mumuq2 = NTupleIn_gen->mumuMass->at(0)*NTupleIn_gen->mumuMass->at(0);
    if (mumuq2<q2Bins_[q2bin-1] || mumuq2>q2Bins_[q2bin]) continue;
    fillHistoGen(q2bin);
  }
  cout << "gen  ev: " << gen_cnt << endl << "reco ev: " << reco_cnt << endl;
  if (doPlot) plotTestHisto(q2bin);
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
      Utility->ReadPreselectionCut(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str());
      Utility->ReadSelectionCuts(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str());
      Utility->ReadGenericParam(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str());

      loadBinning();

      unsigned int q2BinIndx=2;
      bool doPlot=false;
      bool doBatch=false;
      bool mist = 0;
      int evt_cat=0;
      int nev_gen=0;
      int nev_reco=0;
      bool crossCheck=false;
      if ( argc > 1 ) q2BinIndx = atoi(argv[1]);
      if ( argc > 2 ) doPlot = atoi(argv[2]);
      if ( argc > 3 ) doBatch = atoi(argv[3]);
      if ( argc > 4 ) mist = atoi(argv[4]);
      if ( argc > 5 ) evt_cat = atoi(argv[5]);
      if ( argc > 6 ) nev_gen = atoi(argv[6]);
      if ( argc > 7 ) nev_reco = atoi(argv[7]);
      if ( argc > 8 ) crossCheck = atoi(argv[8]);

      // createHistPdfCtKCtL();
      // createHistPdfCtKPhi();
      // createHistPdfCtLPhi();

      if (doBatch) {
        cout << "\n[" << argv[0] << "::main]\t@@@ Setting batch mode @@@" << endl;
        gROOT->SetBatch(doBatch);
      }


      // do all q2Bins at once
      if (q2BinIndx==0) {
        for (q2BinIndx=1; q2BinIndx<10; ++q2BinIndx) {
          // Open input file
          loadEffHisto(q2BinIndx);
          performTest(q2BinIndx, nev_reco, nev_gen, doPlot, mist, evt_cat, crossCheck);
        }
      }
      else {
        // do only one q2 bin
        loadEffHisto(q2BinIndx);
        performTest(q2BinIndx, nev_reco, nev_gen, doPlot, mist, evt_cat, crossCheck);
      }

      if (doBatch == false) theApp.Run();
    }
  else
    {
      cout << "Parameter missing: " << endl;
      cout << "./" << argv[0] << " <q2Bin (0 for all)> " << " <plot histo (decfault=0/no)> " << " <batch running (default=0/no)> " << " <GEN events (default=all)> " << " <RECO events (default=all)> " << endl;

      
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}
