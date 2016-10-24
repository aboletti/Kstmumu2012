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
// #include <RooRealVar.h>
// #include <RooGaussian.h>
// #include <RooAbsPdf.h>
// #include <RooHistPdf.h>
// #include <RooAddPdf.h>
// #include <RooDataSet.h>
// #include <RooGenericPdf.h>
// #include <RooPlot.h>
// #include <RooArgSet.h>
// #include <RooFitResult.h>
// #include <RooPolynomial.h>
// #include <RooMCStudy.h>
// #include <RooMinuit.h>
// #include <RooWorkspace.h>
// #include <RooConstVar.h>
// #include <RooRandom.h>
// #include <RooDataHist.h>
// #include <RooFunctorBinding.h>
// #include <RooStats/RooStatsUtils.h>
// #include <RooMinimizer.h>

#include <ctime>
#include <iostream>
#include <utility>
#include <sstream>

using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::vector;
using std::ios_base;
using std::pair;
using std::make_pair;
using namespace RooFit;

RooRealVar* v1;
RooRealVar* v2;
RooArgSet* vset;

RooRealVar* p0;
RooRealVar* p1;
RooRealVar* p2;

RooAbsPdf* pdf2d;
RooAbsPdf* pdf1d;

RooDataSet* systRes;

double trueP0 = 0.5;
double trueP1 = 0.8;
double trueP2 = -0.2;

double startP0 = 0.5;
double startP1 = 0.8;
double startP2 = -0.2;

double fixP0;
double fixP1;
double cov00;
double cov01;
double cov11;

TH1F* hp0_2df;
TH1F* hp1_2df;
TH1F* hp2_2df;
TH1F* hp0_1df;
TH1F* hp1_1df;
TH1F* hp2_fixf;
TH1F* hp2err_2df;
TH1F* hp2err_fixf;
TH1F* hp2err_syst;
TH1F* hp2comp_squ;
TH1F* hp2comp_lin;

TH2F* hp01_2df;
TH2F* hp02_2df;
TH2F* hp12_2df;
TH2F* hp01_1df;

TCanvas* can0;
TCanvas* can1;
TCanvas* can2;
TCanvas* can01_2d;
TCanvas* can02_2d;
TCanvas* can12_2d;
TCanvas* can01;

RooDataSet* generate(int seed) {
  p0->setVal(trueP0);
  p1->setVal(trueP1);
  p2->setVal(trueP2);

  RooRandom::randomGenerator()->SetSeed(seed+1);
  RooDataSet* ds = pdf2d->generate(*vset,1000,Name(Form("ds%i",seed)));

  return ds;
}

int systRun(RooDataSet& ds, RooArgList& ral) {
  // cout<<((RooRealVar*)ral.find("p0"))->getVal()<<" "<<((RooRealVar*)ral.find("p1"))->getVal()<<endl;
  RooGenericPdf* pdfSyst = new RooGenericPdf("pdfSyst","pdfSyst","1 + p0*v1*v1*(v1+p1) + p1*p2*v1*sin(v2)",ral);
  ((RooRealVar*)ral.find("p0"))->setConstant(true);
  ((RooRealVar*)ral.find("p1"))->setConstant(true);
  p2->setConstant(false);
  p2->setVal(startP2);

  RooFitResult* fitSyst = pdfSyst->fitTo(ds, Verbose(false), PrintLevel(-1), Save(true), PrintEvalErrors(-1));
  if (fitSyst->status()==0 && fitSyst->covQual()==3) systRes->add(RooArgSet(*p2));

  return 0;
}  
  
int fitAndFill (int seed) {

  cout<<"start generation "<<seed<<endl;
  RooDataSet* ds = generate(seed);
  // cout<<"end generation"<<endl;
  
  p0->setConstant(false);
  p1->setConstant(false);
  p2->setConstant(false);
  p0->setVal(startP0);
  p1->setVal(startP1);
  p2->setVal(startP2);

  // cout<<"start fit1"<<endl;
  RooFitResult* fit2d = pdf2d->fitTo(*ds, Verbose(false), PrintLevel(-1), Save(true), PrintEvalErrors(-1));

  float floatError = p2->getError();
  if (fit2d->status()==0 && fit2d->covQual()==3) {
    hp0_2df->Fill(p0->getVal());
    hp1_2df->Fill(p1->getVal());
    hp2_2df->Fill(p2->getVal());
    hp2err_2df->Fill(floatError);
    hp01_2df->Fill(p0->getVal(),p1->getVal());
    hp02_2df->Fill(p0->getVal(),p2->getVal());
    hp12_2df->Fill(p1->getVal(),p2->getVal());
  }

  p0->setConstant(false);
  p1->setConstant(false);
  p0->setVal(startP0);
  p1->setVal(startP1);
  
  // cout<<"start fit2"<<endl;
  RooFitResult* fit1d = pdf1d->fitTo(*ds, Verbose(false), PrintLevel(-1), Save(true), PrintEvalErrors(-1));
  
  // cout<<endl;
  // p0->Print();
  // p1->Print();
  // return 0;
  
  if (fit1d->status()==0 && fit1d->covQual()==3) {
    hp0_1df->Fill(p0->getVal());
    hp1_1df->Fill(p1->getVal());
    hp01_1df->Fill(p0->getVal(),p1->getVal());
    
    p0->setConstant(true);
    p1->setConstant(true);
    p2->setConstant(false);
    p2->setVal(startP2);
    
    // cout<<"start fit3"<<endl;
    RooFitResult* fixfit = pdf2d->fitTo(*ds, Verbose(false), PrintLevel(-1), Save(true), PrintEvalErrors(-1));
    
    if (fixfit->status()==0 && fixfit->covQual()==3) {
      hp2_fixf->Fill(p2->getVal());
      float fixError = p2->getError();
      hp2err_fixf->Fill(fixError);
      
      systRes = new RooDataSet("systRes","systRes",RooArgSet(*p2));
      for (int systToy=0; systToy<100; systToy++) {
	RooArgSet* ral = new RooArgSet(fit1d->randomizePars(),*p2);	 
	systRun(*ds,RooArgSet(*ral,*vset));
      }
      float systError = systRes->rmsVar(*p2)->getVal();
      hp2err_syst->Fill(systError);

      if (fit2d->status()==0 && fit2d->covQual()==3) {
	hp2comp_squ->Fill(sqrt(pow(systError,2)+pow(fixError,2)) - floatError);
	hp2comp_lin->Fill(systError + fixError - floatError);
      }
    }
  }
  
  // cout<<"end generation"<<endl;
  
  return 0;
}

void simpleToy() {
  RooFit::RooMsgService::instance().getStream(1).removeTopic(NumIntegration) ;
  
  v1 = new RooRealVar("v1","v1",-1,1);
  v2 = new RooRealVar("v2","v2",-1*TMath::Pi(),TMath::Pi());
  vset = new RooArgSet(*v1,*v2,"vset");

  p0 = new RooRealVar("p0","p0",-1,1);
  p1 = new RooRealVar("p1","p1",-2,2);
  p2 = new RooRealVar("p2","p2",-1,1);

  vpset2d = new RooArgSet(*v1,*v2,*p0,*p1,*p2,"vpset2d");
  vpset1d = new RooArgSet(*v1,*p0,*p1,"vpset1d");

  pdf2d = new RooGenericPdf("pdf2d","pdf2d","1 + p0*v1*v1*(v1+p1) + p1*p2*v1*sin(v2)",*vpset2d);
  pdf1d = new RooGenericPdf("pdf1d","pdf2d","1 + p0*v1*v1*(v1+p1)",*vpset1d);

  hp0_2df = new TH1F("hp0_2df","hp0_2df",40,-1,1);
  hp1_2df = new TH1F("hp1_2df","hp1_2df",40,-2,2);
  hp2_2df = new TH1F("hp2_2df","hp2_2df",40,-1,1);
  hp0_1df = new TH1F("hp0_1df","hp0_1df",40,-1,1);
  hp1_1df = new TH1F("hp1_1df","hp1_1df",40,-2,2);
  hp2_fixf = new TH1F("hp2_fixf","hp2_fixf",40,-1,1);
  hp2err_2df = new TH1F("hp2err_2df","hp2err_2df",40,0,1);
  hp2err_fixf = new TH1F("hp2err_fixf","hp2err_fixf",40,0,1);
  hp2err_syst = new TH1F("hp2err_syst","hp2err_syst",40,0,1);
  hp2comp_squ = new TH1F("hp2comp_squ","hp2comp_squ",40,-0.2,0.2);
  hp2comp_lin = new TH1F("hp2comp_lin","hp2comp_lin",40,-0.2,0.2);

  hp01_2df = new TH2F("hp01_2df","hp01_2df",50,-1,1,50,-2,2);
  hp02_2df = new TH2F("hp02_2df","hp02_2df",50,-1,1,50,-1,1);
  hp12_2df = new TH2F("hp12_2df","hp12_2df",50,-2,2,50,-1,1);
  hp01_1df = new TH2F("hp01_1df","hp01_1df",50,-1,1,50,-2,2);

  for (int seed=0; seed<1000; seed++) fitAndFill(seed);

  can0 = new TCanvas("can0");
  hp0_2df->Draw();
  hp0_1df->Draw("same");

  can1 = new TCanvas("can1");
  hp1_2df->Draw();
  hp1_1df->Draw("same");

  can2 = new TCanvas("can2");
  hp2_2df->Draw();
  hp2_fixf->Draw("same");

  can2err = new TCanvas("can2err");
  hp2err_2df->Draw();
  hp2err_fixf->Draw("same");
  hp2err_syst->Draw("same");

  can2comp = new TCanvas("can2comp");
  hp2comp_squ->Draw();
  hp2comp_lin->Draw("same");

  can01_2d = new TCanvas("can01_2d");
  hp01_2df->Draw();

  can02_2d = new TCanvas("can02_2d");
  hp02_2df->Draw();

  can12_2d = new TCanvas("can12_2d");
  hp12_2df->Draw();

  can01 = new TCanvas("can01");
  hp01_1df->Draw();

  hp0_1df->SetLineColor(2);
  hp1_1df->SetLineColor(2);
  hp2_fixf->SetLineColor(2);
  hp2err_fixf->SetLineColor(2);
  hp2err_syst->SetLineColor(3);
  hp2comp_lin->SetLineColor(2);

  can0->SaveAs("simpleToy2_p0.pdf");
  can1->SaveAs("simpleToy2_p1.pdf");
  can2->SaveAs("simpleToy2_p2.pdf");
  can01->SaveAs("simpleToy2_p01.pdf");
  can01_2d->SaveAs("simpleToy2_p01_2d.pdf");
  can02_2d->SaveAs("simpleToy2_p02_2d.pdf");
  can12_2d->SaveAs("simpleToy2_p12_2d.pdf");
  can2err->SaveAs("simpleToy2_p2err.pdf");
  can2comp->SaveAs("simpleToy2_p2comp.pdf");

  can0->SaveAs("simpleToy2_p0.root");
  can1->SaveAs("simpleToy2_p1.root");
  can2->SaveAs("simpleToy2_p2.root");
  can01->SaveAs("simpleToy2_p01.root");
  can01_2d->SaveAs("simpleToy2_p01_2d.root");
  can02_2d->SaveAs("simpleToy2_p02_2d.root");
  can12_2d->SaveAs("simpleToy2_p12_2d.root");
  can2err->SaveAs("simpleToy2_p2err.root");
  can2comp->SaveAs("simpleToy2_p2comp.root");

  return;
}
