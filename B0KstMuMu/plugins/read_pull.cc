#include <iostream>

#include "TH1.h"
#include "TLine.h"
#include "TCanvas.h"
// #include ".h"
// #include ".h"
#include "TFile.h"

// #include "RooRealVar.h"
// #include "RooFitResult.h"
// #include "RooArgList.h"
// #include ".h"
// #include ".h"
// #include ".h"

using namespace std;

double minNll;
int best;
double edm;
float best_P5p;
float best_P1;
float best_As5;

int counter = 0;

void open (int q2BinIndx, int scanIndx, bool print=false)
{
  TString filename = Form("%i/Fitresult5_2.root",scanIndx);
  if ( gSystem->AccessPathName(filename) ) return;
  TFile* f = new TFile(filename);
  if (f) {
    RooFitResult* fr = (RooFitResult*)f->Get(Form("fitResult_Bin%i",q2BinIndx)); 
    if (!fr) return;

    if ( !print ) {
      if ( fr->status()==0 && fr->covQual()==3 ) {

	  if (minNll>fr->minNll() || (minNll==fr->minNll() && edm>fr->edm())) {
	    minNll = fr->minNll();
	    best = scanIndx;
	    edm = fr->edm();
	    best_As5 = ((RooRealVar*)fr->floatParsFinal().at(0))->getVal();
	    best_P1  = ((RooRealVar*)fr->floatParsFinal().at(1))->getVal();
	    best_P5p = ((RooRealVar*)fr->floatParsFinal().at(2))->getVal();
	  }
	}
    } else {
      cout<<"(fnc="<<fr->minNll()<<",edm="<<fr->edm()<<",stat="<<fr->status()<<",covQ="<<fr->covQual()<<")"<<endl;
      cout<<"A5s\t"<<((RooRealVar*)fr->floatParsInit().at(0))->getVal()<<"  \t"<<((RooRealVar*)fr->floatParsFinal().at(0))->getVal()
	  <<"\nP1\t"<<((RooRealVar*)fr->floatParsInit().at(1))->getVal()<<"  \t"<<((RooRealVar*)fr->floatParsFinal().at(1))->getVal()
	  <<"\nP5p\t"<<((RooRealVar*)fr->floatParsInit().at(2))->getVal()<<"  \t"<<((RooRealVar*)fr->floatParsFinal().at(2))->getVal()<<endl;
    }
  } else return;

  f->Close();
  delete f;
}

void read (int q2BinIndx)
{
  if (q2BinIndx==4 || q2BinIndx==6) return;

  counter = 0;

  minNll = 9999999;
  for (int cnt=100; cnt<140; cnt++) open(q2BinIndx, cnt);
  if ( minNll == 9999999 ) {
    // cout<<q2BinIndx<<" 0 0 0 0"<<endl;
    cout<<0<<endl;
    return;
  }

  // cout<<q2BinIndx<<" "<<best<<" "<<best_As5<<" "<<best_P1<<" "<<best_P5p<<endl;
  cout<<1<<endl;
  return;
}

void read_pull (int q2BinIndx = 0)
{
  if (q2BinIndx==0) for (q2BinIndx=1;q2BinIndx<10;q2BinIndx++) read(q2BinIndx-1);
  else if (q2BinIndx>0 && q2BinIndx<10) read(q2BinIndx-1);
  
  return;
}
