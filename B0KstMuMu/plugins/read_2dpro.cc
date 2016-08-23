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
double bestP1[1];
double bestP5p[1];

TGraph* gP5p[9];
TGraph* gP1 [9];

TGraph* g2d1s [9];
TGraph* g2d2s [9];
TGraph* g2d3s [9];

TGraph* gBest[9];
TGraph* gLim[9];

TCanvas* can[9];


vector<double> vP5p (0);
vector<double> vP1  (0);
vector<double> vA5s (0);
vector<double> vSig (0);
vector<double> vBkg (0);
vector<double> vLik (0);

vector<double> v1sP5p (0);
vector<double> v1sP1  (0);
vector<double> v1sLik (0);
vector<double> v2sP5p (0);
vector<double> v2sP1  (0);
vector<double> v2sLik (0);
vector<double> v3sP5p (0);
vector<double> v3sP1  (0);
vector<double> v3sLik (0);


TF1* lim[9];
double fl[9] = {0.641004,0.799186,0.619384,0.503676,0,0.392124,0,0.476826,0.377081};

int nLim[9] = {199,199,199,199,0,199,0,199,199};

void open (int q2BinIndx, int scanIndx, bool print=false)
{
  TString filename = Form("Data_scan_2dpro/%i/Fitresult4_2.root",scanIndx);
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

	  vA5s.push_back(((RooRealVar*)fr->floatParsFinal().at(0))->getVal());
	  vSig.push_back(((RooRealVar*)fr->floatParsFinal().at(2))->getVal());
	  vBkg.push_back(((RooRealVar*)fr->floatParsFinal().at(1))->getVal());

	  vP1.push_back(((RooRealVar*)fr->constPars().at(3))->getVal());
	  vP5p.push_back(((RooRealVar*)fr->constPars().at(4))->getVal());
	  // cout<<((RooRealVar*)fr->constPars().at(4))->getVal()<<endl;
	  vLik.push_back(fr->minNll());
	  if ( !(q2BinIndx==7 && fr->minNll()<-2000) ) {
	    // cout<<scanIndx<<endl;
	    if (minNll>fr->minNll() || (minNll==fr->minNll() && edm>fr->edm())) {
	      minNll = fr->minNll();
	      best = scanIndx;
	      edm = fr->edm();
	      bestP1[0] = ((RooRealVar*)fr->constPars().at(3))->getVal();
	      bestP5p[0] = ((RooRealVar*)fr->constPars().at(4))->getVal();
	    }
	  }
	} //else cout<<scanIndx<<" "<<fr->status()<<" "<<fr->covQual()<<endl;
    } else {
      cout<<"(fnc="<<fr->minNll()<<",edm="<<fr->edm()<<",stat="<<fr->status()<<",covQ="<<fr->covQual()<<")"<<endl;
	cout<<"A5s\t"<<((RooRealVar*)fr->floatParsInit().at(0))->getVal()<<"  \t"<<((RooRealVar*)fr->floatParsFinal().at(0))->getVal()
	    <<"\nP1\t         \t"<<((RooRealVar*)fr->constPars().at(3))->getVal()
            <<"\nP5p\t         \t"<<((RooRealVar*)fr->constPars().at(4))->getVal()<<endl;
    
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
  vLik.clear();
  v1sP5p.clear();
  v1sP1.clear();
  v1sLik.clear();
  v2sP5p.clear();
  v2sP1.clear();
  v2sLik.clear();
  v3sP5p.clear();
  v3sP1.clear();
  v3sLik.clear();

  minNll = 9999999;
  for (int cnt=1; cnt<=8100; cnt++) {
    if (q2BinIndx==0 && (cnt==1627 || cnt==1652 || cnt==2644 || cnt==4357 || cnt==5575 || cnt==6710) ) continue;
    if (q2BinIndx==1 && (cnt==1479 || cnt==1735) ) continue;
    if (q2BinIndx==2 && (cnt==944  || cnt==1358 || cnt==3252 || cnt==6318 || cnt==7622) ) continue;
    if (q2BinIndx==3 && (cnt==2297 || cnt==6318 || cnt==6710) ) continue;
    if (q2BinIndx==5 && (cnt==1358 || cnt==1735 || cnt==2297 || cnt==2801) ) continue;
    if (q2BinIndx==7 && (cnt==1215 || cnt==6741) ) continue;
    if (q2BinIndx==8 && (cnt==1227 || cnt==2462 || cnt==8069) ) continue;

    open(q2BinIndx, cnt);
  }
  if ( minNll == 9999999 ) {
    // cout<<"Best "<<q2BinIndx<<": no good results"<<endl;
    return;
  }

  cout<<endl<<"Best "<<q2BinIndx<<": "<<best<<" ("<<minNll<<"), "<<vLik.size()<<" good fits"<<endl;

  open(q2BinIndx, best, true);

  double* aP5pLim = new double[nLim[q2BinIndx]];
  double* aP1Lim = new double[nLim[q2BinIndx]];
  fstream fin (Form("waveP_lim%i.list",q2BinIndx),fstream::in);
  fin.ignore(100,'\n');
  fin.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) fin>>aP5pLim[i]>>aP1Lim[i];

  for (int i=0; i<vLik.size(); i++) {
    if (vLik[i] < minNll+0.5) {
      v1sLik.push_back(vLik[i]);
      v1sP1.push_back(vP1[i]);
      v1sP5p.push_back(vP5p[i]);
    } else if (vLik[i] < minNll+2) {
      v2sLik.push_back(vLik[i]);
      v2sP1.push_back(vP1[i]);
      v2sP5p.push_back(vP5p[i]);
    } else {
      v3sLik.push_back(vLik[i]);
      v3sP1.push_back(vP1[i]);
      v3sP5p.push_back(vP5p[i]);
    }
  }

  double* a1sP5p = new double[v1sP1.size()];
  double* a1sP1  = new double[v1sP1.size()];
  // double* a1sLik = new double[v1sLik.size()];
  double* a2sP5p = new double[v2sP1.size()];
  double* a2sP1  = new double[v2sP1.size()];
  // double* a2sLik = new double[v2sLik.size()];
  double* a3sP5p = new double[v3sP1.size()];
  double* a3sP1  = new double[v3sP1.size()];
  // double* a3sLik = new double[v3sLik.size()];

  // double* aA5s = new double[vLik.size()];
  // double* aSig = new double[vLik.size()];
  // double* aBkg = new double[vLik.size()];

  for (int i=0; i<v1sP1.size(); i++) {
    a1sP5p[i]=v1sP5p[i];
    a1sP1 [i]=v1sP1 [i];
  }
  for (int i=0; i<v2sP1.size(); i++) {
    a2sP5p[i]=v2sP5p[i];
    a2sP1 [i]=v2sP1 [i];
  }
  for (int i=0; i<v3sP1.size(); i++) {
    a3sP5p[i]=v3sP5p[i];
    a3sP1 [i]=v3sP1 [i];
  }

  can[q2BinIndx] = new TCanvas (Form("can%i",q2BinIndx),Form("can%i",q2BinIndx));

  // cout<<nLim[q2BinIndx]<<" "<<v2sP1.size()<<" "<<v1sP1.size()<<" "<<v3sP1.size()<<endl;


  gLim[q2BinIndx] = new TGraph(nLim[q2BinIndx], aP1Lim, aP5pLim);
  gBest[q2BinIndx] = new TGraph(1, bestP1, bestP5p);

  g2d2s[q2BinIndx] = new TGraph( v2sP1.size(), a2sP1, a2sP5p);
  g2d2s[q2BinIndx]->GetXaxis()->SetLimits(-1,1);
  // g2d1[q2BinIndx]->GetYaxis()->SetLimits(-1.2,.8);
  g2d2s[q2BinIndx]->Draw("AP");
  g2d2s[q2BinIndx]->SetMarkerStyle(7);
  g2d2s[q2BinIndx]->SetMarkerColor(3);
  
  if (v1sP1.size()>0) {
    g2d1s[q2BinIndx] = new TGraph( v1sP1.size(), a1sP1, a1sP5p);
    g2d1s[q2BinIndx]->Draw("Psame");
    g2d1s[q2BinIndx]->SetMarkerStyle(7);
    g2d1s[q2BinIndx]->SetMarkerColor(4);
  }

  if (v3sP1.size()>0) {
    g2d3s[q2BinIndx] = new TGraph( v3sP1.size(), a3sP1, a3sP5p);
    g2d3s[q2BinIndx]->Draw("Psame");
    g2d3s[q2BinIndx]->SetMarkerStyle(7);
    g2d3s[q2BinIndx]->SetMarkerColor(2);
  }

  gBest[q2BinIndx]->Draw("sameP");
  gBest[q2BinIndx]->SetMarkerStyle(34);
  gBest[q2BinIndx]->SetMarkerColor(2);

  gLim[q2BinIndx]->Draw("sameP");
  gLim[q2BinIndx]->SetMarkerStyle(7);

  can[q2BinIndx]->SaveAs(Form("Data_scan_2dpro/scan2d_b%i.pdf",q2BinIndx));

  // can[q2BinIndx] = new TCanvas (Form("can%i",q2BinIndx),Form("can%i",q2BinIndx));

  // gP5p[q2BinIndx] = new TGraph( vP5p.size(), aP5p, aLikP5p );
  // gP1 [q2BinIndx] = new TGraph( vP1.size(), aP1 , aLikP1 );
  // // gA5s[q2BinIndx] = new TGraph( vLik.size(), aA5s, aLik );
  // // gSig[q2BinIndx] = new TGraph( vLik.size(), aSig, aLik );
  // // gBkg[q2BinIndx] = new TGraph( vLik.size(), aBkg, aLik );

  // can[q2BinIndx]->Divide(2,1);
  // can[q2BinIndx]->cd(1);
  // gP5p[q2BinIndx]->Draw("AP");
  // gP5p[q2BinIndx]->SetMarkerStyle(7);
  // gP5p[q2BinIndx]->SetTitle("P5p");
  // can[q2BinIndx]->cd(2);
  // gP1[q2BinIndx]->Draw("AP");
  // gP1[q2BinIndx]->SetMarkerStyle(7);
  // gP1[q2BinIndx]->SetTitle("P1");
  // // can[q2BinIndx]->cd(3);
  // // gA5s[q2BinIndx]->Draw("AP");
  // // gA5s[q2BinIndx]->SetMarkerStyle(7);
  // // gA5s[q2BinIndx]->SetTitle("As5");
  // // can[q2BinIndx]->cd(4);
  // // gSig[q2BinIndx]->Draw("AP");
  // // gSig[q2BinIndx]->SetMarkerStyle(7);
  // // can[q2BinIndx]->cd(5);
  // // gBkg[q2BinIndx]->Draw("AP");
  // // gBkg[q2BinIndx]->SetMarkerStyle(7);

  // can[q2BinIndx]->SaveAs(Form("Data_scan_2dpro/like_b%i.pdf",q2BinIndx));

  return;
}

void read_2dpro (int q2BinIndx = 0)
{
  if (q2BinIndx==0) for (q2BinIndx=1;q2BinIndx<10;q2BinIndx++) read(q2BinIndx-1);
  else if (q2BinIndx>0 && q2BinIndx<10) read(q2BinIndx-1);
  
  return;
}
