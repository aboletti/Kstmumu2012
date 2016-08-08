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

TGraph* gP5p[9];
TGraph* gP1 [9];
TGraph* gA5s[9];
TGraph* gSig[9];
TGraph* gBkg[9];

TCanvas* can[9];

vector<double> vP5p (0);
vector<double> vP1  (0);
vector<double> vA5s (0);
vector<double> vSig (0);
vector<double> vBkg (0);
vector<double> vLik (0);

void open (int q2BinIndx, int scanIndx, bool print=false)
{
  TString filename = Form("Data_scan_test2/%i/Fitresult4_2.root",scanIndx);
  if ( gSystem->AccessPathName(filename) ) return;
  TFile* f = new TFile(filename);
  if (f) {
    RooFitResult* fr = (RooFitResult*)f->Get(Form("fitResult_Bin%i",q2BinIndx-1)); 
    if (!fr) {/*cout<<q2BinIndx<<":\t"<<scanIndx<<endl;*/ return;}
    // cout<<scanIndx<<" (fnc="<<fr->minNll()<<",edm="<<fr->edm()<<",stat="<<fr->status()<<",covQ="<<fr->covQual()<<")"<<endl;
    // fr->Print();
    // cout<<endl;

    if ( !print ) {
      if ( fr->status()==0 && fr->covQual()==3 ) if ( true || fr->status()>-1 ) {
	  if (q2BinIndx>1 || fr->minNll()<-622.040) cout<<q2BinIndx-1<<" "<<scanIndx<<" "<<fr->minNll()<<endl;
	  vP5p.push_back(((RooRealVar*)fr->floatParsFinal().at(2))->getVal());
	  vP1.push_back(((RooRealVar*)fr->floatParsFinal().at(1))->getVal());
	  vA5s.push_back(((RooRealVar*)fr->floatParsFinal().at(0))->getVal());
	  // vSig.push_back(((RooRealVar*)fr->floatParsFinal().at(4))->getVal());
	  // vBkg.push_back(((RooRealVar*)fr->floatParsFinal().at(3))->getVal());
	  vLik.push_back(fr->minNll());
	  // cout<<scanIndx<<endl;
	  if (minNll>fr->minNll() || (minNll==fr->minNll() && edm>fr->edm())) {
	    minNll = fr->minNll();
	    best = scanIndx;
	    edm = fr->edm();
	  }
	} //else cout<<scanIndx<<" "<<fr->status()<<" "<<fr->covQual()<<endl;
    } else {
      cout<<scanIndx<<": "<<endl;
      cout<<"(fnc="<<fr->minNll()<<",edm="<<fr->edm()<<",stat="<<fr->status()<<",covQ="<<fr->covQual()<<")"<<endl;
      cout<<"A5s\t"<<((RooRealVar*)fr->floatParsInit().at(0))->getVal()<<"  \t"<<((RooRealVar*)fr->floatParsFinal().at(0))->getVal()
	  <<"\nP1\t"<<((RooRealVar*)fr->floatParsInit().at(1))->getVal()<<"  \t"<<((RooRealVar*)fr->floatParsFinal().at(1))->getVal()
	  <<"\nP5p\t"<<((RooRealVar*)fr->floatParsInit().at(2))->getVal()<<"  \t"<<((RooRealVar*)fr->floatParsFinal().at(2))->getVal()<<endl;
    }
  } else {
    // cout<<"No file for result "<<scanIndx<<endl;
    return;
  }
  f->Close();
  delete f;
}

void read (int q2BinIndx)
{
  if (q2BinIndx==5 || q2BinIndx==7) return;

  vP5p.clear();
  vP1.clear();
  vA5s.clear();
  vSig.clear();
  vBkg.clear();
  vLik.clear();

  minNll = 9999999;
  for (int cnt=600; cnt<800; cnt++) open(q2BinIndx, cnt);
  if ( minNll == 9999999 ) {
    // cout<<"Best "<<q2BinIndx-1<<": no good results"<<endl;
    return;
  }

  // cout<<"Best "<<q2BinIndx-1<<": "<<best<<" ("<<minNll<<"), "<<vLik.size()<<" good fits"<<endl;

  // open(q2BinIndx, best, true);

  // double* aP5p = new double[vLik.size()];
  // double* aP1  = new double[vLik.size()];
  // double* aA5s = new double[vLik.size()];
  // double* aSig = new double[vLik.size()];
  // double* aBkg = new double[vLik.size()];
  // double* aLik = new double[vLik.size()];
  // for (int i=0; i<vLik.size(); i++) {
  //   aP5p[i]=vP5p[i];
  //   aP1 [i]=vP1 [i];
  //   aA5s[i]=vA5s[i];
  //   // aSig[i]=vSig[i];
  //   // aBkg[i]=vBkg[i];
  //   aLik[i]=vLik[i];
  // }
  // gP5p[q2BinIndx] = new TGraph( vLik.size(), aP5p, aLik );
  // gP1 [q2BinIndx] = new TGraph( vLik.size(), aP1 , aLik );
  // gA5s[q2BinIndx] = new TGraph( vLik.size(), aA5s, aLik );
  // // gSig[q2BinIndx] = new TGraph( vLik.size(), aSig, aLik );
  // // gBkg[q2BinIndx] = new TGraph( vLik.size(), aBkg, aLik );

  // can[q2BinIndx] = new TCanvas (Form("can%i",q2BinIndx),Form("can%i",q2BinIndx));
  // can[q2BinIndx]->Divide(3,1);
  // can[q2BinIndx]->cd(1);
  // gP5p[q2BinIndx]->Draw("AP");
  // gP5p[q2BinIndx]->SetMarkerStyle(7);
  // gP5p[q2BinIndx]->SetTitle("P5p");
  // can[q2BinIndx]->cd(2);
  // gP1[q2BinIndx]->Draw("AP");
  // gP1[q2BinIndx]->SetMarkerStyle(7);
  // gP1[q2BinIndx]->SetTitle("P1");
  // can[q2BinIndx]->cd(3);
  // gA5s[q2BinIndx]->Draw("AP");
  // gA5s[q2BinIndx]->SetMarkerStyle(7);
  // gA5s[q2BinIndx]->SetTitle("As5");
  // // can[q2BinIndx]->cd(4);
  // // gSig[q2BinIndx]->Draw("AP");
  // // gSig[q2BinIndx]->SetMarkerStyle(7);
  // // can[q2BinIndx]->cd(5);
  // // gBkg[q2BinIndx]->Draw("AP");
  // // gBkg[q2BinIndx]->SetMarkerStyle(7);

  // can[q2BinIndx]->SaveAs(Form("Data_scan_test2/like_b%i.pdf",q2BinIndx-1));
  return;
}

void read_data (int q2BinIndx = 0)
{
  if (q2BinIndx==0) for (q2BinIndx=1;q2BinIndx<10;q2BinIndx++) read(q2BinIndx);
  else if (q2BinIndx>0 && q2BinIndx<10) read(q2BinIndx);
  
  return;
}
