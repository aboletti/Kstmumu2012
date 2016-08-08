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
using namespace RooFit;

double minNll;
int best;
double edm;

void open (int q2BinIndx, int scanIndx, bool print=false)
{
  TFile* f = new TFile(Form("Scan/%i/Fitresult5_2.root",scanIndx));
  if (f) {
    RooFitResult* fr = (RooFitResult*)f->Get(Form("fitResult_Bin%i",q2BinIndx-1)); 
    if (!fr) {cout<<q2BinIndx<<":\t"<<scanIndx<<endl; return;}

    if ( !print ) {
      if (fr->status()==0 && fr->covQual()==3 && (minNll>fr->minNll() || (minNll==fr->minNll() && edm>fr->edm()))) {
	minNll = fr->minNll();
	best = scanIndx;
	edm = fr->edm();
      }
    } else {
      cout<<scanIndx<<": "<<endl;
      cout<<"(fnc="<<fr->minNll()<<",edm="<<fr->edm()<<",stat="<<fr->status()<<",covQ="<<fr->covQual()<<")"<<endl;
      cout<<"Fl\t"<<((RooRealVar*)fr->floatParsInit().at(0))->getVal()<<"  \t"<<((RooRealVar*)fr->floatParsFinal().at(0))->getVal()
	  <<"\nP1\t"<<((RooRealVar*)fr->floatParsInit().at(1))->getVal()<<"  \t"<<((RooRealVar*)fr->floatParsFinal().at(1))->getVal()
	  <<"\nP5p\t"<<((RooRealVar*)fr->floatParsInit().at(2))->getVal()<<"  \t"<<((RooRealVar*)fr->floatParsFinal().at(2))->getVal()<<endl;
    }
  } else {
    cout<<"No file for result "<<scanIndx<<endl;
    return;
  }
  f->Close();
  delete f;
}

void read (int q2BinIndx)
{
  if (q2BinIndx==5 || q2BinIndx==7) return;

  minNll = 9999999;
  for (int cnt=0; cnt<200; cnt++) open(q2BinIndx, cnt);

  cout<<"Best "<<q2BinIndx<<": "<<best<<" ("<<minNll<<")"<<endl;

  open(q2BinIndx, best, true);

  return;
}

void read_scan (int q2BinIndx = 0)
{
  if (q2BinIndx==0) for (q2BinIndx=1;q2BinIndx<10;q2BinIndx++) read(q2BinIndx);
  else if (q2BinIndx>0 && q2BinIndx<10) read(q2BinIndx);
  
  return;
}
