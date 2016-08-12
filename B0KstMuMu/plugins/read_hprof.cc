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
TGraph* g2d1 [9];
TGraph* g2d2 [9];
TGraph* gLim[9];

TCanvas* can[9];
TCanvas* can1[9];

vector<double> vP5p (0);
vector<double> vP1  (0);
vector<double> vA5s (0);
vector<double> vSig (0);
vector<double> vBkg (0);
vector<double> vLikP5p (0);
vector<double> vLikP1 (0);
vector<double> vP5pf(0);
vector<double> vP1f (0);

TF1* lim[9];
double fl[9] = {0.641004,0.799186,0.619384,0.503676,0,0.392124,0,0.476826,0.377081};

int nLim[9] = {199,199,199,199,0,199,0,199,199};

void open (int q2BinIndx, int scanIndx, bool print=false)
{
  TString filename = Form("Data_scan_hprof2/%i/Fitresult4_2.root",scanIndx);
  if ( gSystem->AccessPathName(filename) ) return;
  TFile* f = new TFile(filename);
  if (f) {
    RooFitResult* fr = (RooFitResult*)f->Get(Form("fitResult_Bin%i",q2BinIndx)); 
    if (!fr) {/*cout<<q2BinIndx<<":\t"<<scanIndx<<endl;*/ return;}
    // cout<<scanIndx<<" (fnc="<<fr->minNll()<<",edm="<<fr->edm()<<",stat="<<fr->status()<<",covQ="<<fr->covQual()<<")"<<endl;
    // fr->Print();
    // cout<<endl;

    if ( !print ) {
      if ( fr->status()==0 && fr->covQual()==3 ) if ( true || fr->status()>-1 ) {
	  // if (q2BinIndx>0 || fr->minNll()<-622.040) cout<<q2BinIndx<<" "<<scanIndx<<" "<<fr->minNll()<<endl;
	  // vA5s.push_back(((RooRealVar*)fr->floatParsFinal().at(0))->getVal());
	  // vSig.push_back(((RooRealVar*)fr->floatParsFinal().at(4))->getVal());
	  // vBkg.push_back(((RooRealVar*)fr->floatParsFinal().at(3))->getVal());

	  if (scanIndx%2==0) {
	    vP1.push_back(((RooRealVar*)fr->constPars().at(3))->getVal());
	    vP5pf.push_back(((RooRealVar*)fr->floatParsFinal().at(1))->getVal());
	    vLikP1.push_back(fr->minNll());
	  } else {
	    vP5p.push_back(((RooRealVar*)fr->constPars().at(3))->getVal());
	    vP1f.push_back(((RooRealVar*)fr->floatParsFinal().at(1))->getVal());
	    vLikP5p.push_back(fr->minNll());
	  }
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
      if (scanIndx%2==0) {
	cout<<"A5s\t"<<((RooRealVar*)fr->floatParsInit().at(0))->getVal()<<"  \t"<<((RooRealVar*)fr->floatParsFinal().at(0))->getVal()
	    <<"\nP1\t         \t"<<((RooRealVar*)fr->constPars().at(3))->getVal()
	    <<"\nP5p\t"<<((RooRealVar*)fr->floatParsInit().at(1))->getVal()<<"  \t"<<((RooRealVar*)fr->floatParsFinal().at(1))->getVal()<<endl;
      } else {
	cout<<"A5s\t"<<((RooRealVar*)fr->floatParsInit().at(0))->getVal()<<"  \t"<<((RooRealVar*)fr->floatParsFinal().at(0))->getVal()
            <<"\nP1\t"<<((RooRealVar*)fr->floatParsInit().at(1))->getVal()<<"  \t"<<((RooRealVar*)fr->floatParsFinal().at(1))->getVal()
            <<"\nP5p\t         \t"<<((RooRealVar*)fr->constPars().at(3))->getVal()<<endl;
      }
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
  if (q2BinIndx==4 || q2BinIndx==6) return;

  vP5p.clear();
  vP1.clear();
  vLikP5p.clear();
  vLikP1.clear();
  vP5pf.clear();
  vP1f.clear();

  minNll = 9999999;
  for (int cnt=0; cnt<800; cnt++) open(q2BinIndx, cnt);
  if ( minNll == 9999999 ) {
    // cout<<"Best "<<q2BinIndx<<": no good results"<<endl;
    return;
  }

  cout<<"Best "<<q2BinIndx<<": "<<best<<" ("<<minNll<<"), "<<vLikP1.size()<<", "<<vLikP5p.size()<<" good fits"<<endl;

  open(q2BinIndx, best, true);

  double* aP5pLim = new double[nLim[q2BinIndx]];
  double* aP1Lim = new double[nLim[q2BinIndx]];
  fstream fin (Form("waveP_lim%i.list",q2BinIndx),fstream::in);
  fin.ignore(100,'\n');
  fin.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) fin>>aP5pLim[i]>>aP1Lim[i];

  double* aP5p = new double[vP5p.size()];
  double* aP1  = new double[vP1.size()];
  // double* aA5s = new double[vLik.size()];
  // double* aSig = new double[vLik.size()];
  // double* aBkg = new double[vLik.size()];
  double* aLikP5p = new double[vLikP5p.size()];
  double* aLikP1 = new double[vLikP1.size()];
  double* aP5pf = new double[vP5pf.size()];
  double* aP1f  = new double[vP1f.size()];

  for (int i=0; i<vP5p.size(); i++) {
    aP5p[i]=vP5p[i];
    aP1f [i]=vP1f [i];
    aLikP5p[i]=vLikP5p[i]; 
  }
  for (int i=0; i<vP1.size(); i++) {
    aP1 [i]=vP1 [i];
    aP5pf[i]=vP5pf[i];
    aLikP1[i]=vLikP1[i];
  }

  can1[q2BinIndx] = new TCanvas (Form("can1%i",q2BinIndx),Form("can1%i",q2BinIndx));

  g2d1[q2BinIndx] = new TGraph( vP5p.size(), aP1f, aP5p);
  g2d2[q2BinIndx] = new TGraph(vP1.size(), aP1, aP5pf);
  // lim[q2BinIndx] = new TF1(Form("lim%i",q2BinIndx),"pol1",-1,1);
  // lim[q2BinIndx]->SetParameter(0,-1.*(1-(1-fl[q2BinIndx])*3./4.)/sqrt(fl[q2BinIndx]*(1-fl[q2BinIndx])));
  // lim[q2BinIndx]->SetParameter(1,-1.*(1-fl[q2BinIndx])/4./sqrt(fl[q2BinIndx]*(1-fl[q2BinIndx])));
  gLim[q2BinIndx] = new TGraph(nLim[q2BinIndx], aP1Lim, aP5pLim);

  g2d1[q2BinIndx]->GetXaxis()->SetLimits(-1,1);
  g2d1[q2BinIndx]->GetYaxis()->SetLimits(-1.2,.8);
  g2d1[q2BinIndx]->Draw("AP");
  g2d1[q2BinIndx]->SetMarkerStyle(7);
  g2d1[q2BinIndx]->SetMarkerColor(3);
  g2d2[q2BinIndx]->Draw("Psame");
  g2d2[q2BinIndx]->SetMarkerStyle(7);
  g2d2[q2BinIndx]->SetMarkerColor(2);
  // lim[q2BinIndx]->Draw("same");
  gLim[q2BinIndx]->Draw("sameP");
  gLim[q2BinIndx]->SetMarkerStyle(7);

  can1[q2BinIndx]->SaveAs(Form("Data_scan_hprof2/scanLimit_%i.pdf",q2BinIndx));

  can[q2BinIndx] = new TCanvas (Form("can%i",q2BinIndx),Form("can%i",q2BinIndx));

  gP5p[q2BinIndx] = new TGraph( vP5p.size(), aP5p, aLikP5p );
  gP1 [q2BinIndx] = new TGraph( vP1.size(), aP1 , aLikP1 );
  // gA5s[q2BinIndx] = new TGraph( vLik.size(), aA5s, aLik );
  // gSig[q2BinIndx] = new TGraph( vLik.size(), aSig, aLik );
  // gBkg[q2BinIndx] = new TGraph( vLik.size(), aBkg, aLik );

  can[q2BinIndx]->Divide(2,1);
  can[q2BinIndx]->cd(1);
  gP5p[q2BinIndx]->Draw("AP");
  gP5p[q2BinIndx]->SetMarkerStyle(7);
  gP5p[q2BinIndx]->SetTitle("P5p");
  can[q2BinIndx]->cd(2);
  gP1[q2BinIndx]->Draw("AP");
  gP1[q2BinIndx]->SetMarkerStyle(7);
  gP1[q2BinIndx]->SetTitle("P1");
  // can[q2BinIndx]->cd(3);
  // gA5s[q2BinIndx]->Draw("AP");
  // gA5s[q2BinIndx]->SetMarkerStyle(7);
  // gA5s[q2BinIndx]->SetTitle("As5");
  // can[q2BinIndx]->cd(4);
  // gSig[q2BinIndx]->Draw("AP");
  // gSig[q2BinIndx]->SetMarkerStyle(7);
  // can[q2BinIndx]->cd(5);
  // gBkg[q2BinIndx]->Draw("AP");
  // gBkg[q2BinIndx]->SetMarkerStyle(7);

  can[q2BinIndx]->SaveAs(Form("Data_scan_hprof2/like_b%i.pdf",q2BinIndx));
  return;
}

void read_hprof (int q2BinIndx = 0)
{
  if (q2BinIndx==0) for (q2BinIndx=1;q2BinIndx<10;q2BinIndx++) read(q2BinIndx-1);
  else if (q2BinIndx>0 && q2BinIndx<10) read(q2BinIndx-1);
  
  return;
}
