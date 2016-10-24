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
double legaP1[1];
double legaP5p[1];

TGraph* gP5p[9];
TGraph* gP1 [9];

TGraph* g2d1s [9];
TGraph* g2d2s [9];
TGraph* g2d3s [9];

TGraph* gBest[9];
TGraph* gLega[9];
TGraph* gLim[9];
TGraph* gLimt[9];
TGraph* gLiml[9];
TGraph* gLim0[9];
TGraph* gLim1[9];
TGraph* gLim2[9];
TGraph* gLim3[9];
TGraph* gLim4[9];
TGraph* gLim5[9];
TGraph* gLim6[9];


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

vector<int> good (8100,0);

TF1* lim[9];
double fl[9] = {0.641004,0.799186,0.619384,0.503676,0,0.392124,0,0.476826,0.377081};

int nLim[9] = {199,198,199,199,0,199,0,199,199};

double maxP1 [9] = { 0.85, 0.35, 1.00, 0.00,0,-0.25,0, 0.10,-0.17};
double minP1 [9] = {-0.65,-1.00,-0.15,-0.90,0,-0.80,0,-0.75,-0.85};
double maxP5p[9] = { 0.65, 0.50,-0.60,-0.30,0,-0.45,0,-0.35,-0.32};
double minP5p[9] = {-0.45,-1.10,-1.30,-0.95,0,-0.85,0,-0.95,-0.82};

int counter = 0;

float As5BF[9] = {       1, 0.999115,-0.434302, 0.530473,0,-0.933958,0,0.0524524,-0.0542194};
float P1BF [9] = {0.150675,-0.668496, 0.498243,-0.476409,0,-0.453944,0,-0.368363, -0.547864};
float P5BF [9] = {0.105585,-0.562285,-0.953227,-0.643975,0, -0.73848,0,-0.646982, -0.549104};

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
      if ( ( fr->status()==0 && fr->covQual()==3 ) ) if ( true || ( fr->status()>-1 && fr->minNll()<0 && fr->minNll()>-1e5 ) ) {
	  good[scanIndx-1] = 1;
	  if ( ((RooRealVar*)fr->floatParsInit().at(0))->getVal() != 0 ) counter++;

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

  counter=0;

  minNll = 9999999;
  for (int cnt=1; cnt<=8100; cnt++) {
    good[cnt-1] = 0;

    if (q2BinIndx==0 && (cnt==2644 || cnt==6710) ) continue;
    if (q2BinIndx==1 && (cnt==1735) ) continue;
    if (q2BinIndx==2 && (cnt==1358 || cnt==1627 || cnt==3252) ) continue;
    if (q2BinIndx==3 && (cnt==2297) ) continue;
    if (q2BinIndx==5 && (cnt==1227 || cnt==1358 || cnt==1735 || cnt==2297 || cnt==2801) ) continue;
    if (q2BinIndx==7 && (cnt==1215) ) continue;
    if (q2BinIndx==8 && (cnt==2462 || cnt==2548 || cnt==3447) ) continue;

    if (q2BinIndx==2 && (cnt==944  /*|| cnt==4791*/ || cnt==6318/* || cnt==7261*/) ) continue;
    if (q2BinIndx==3 && (/*cnt==1605 ||*/ cnt==6318 || cnt==6710) ) continue;
    if (q2BinIndx==7 && (/*cnt==4254 ||*/ cnt==6741) ) continue;
    if (q2BinIndx==8 && (cnt==4615) ) continue; 

    // cout<<cnt<<endl;
    open(q2BinIndx, cnt);
  }
  if ( minNll == 9999999 ) {
    // cout<<"Best "<<q2BinIndx<<": no good results"<<endl;
    // return;
  }

  cout<<endl<<"Best "<<q2BinIndx<<": "<<best<<" ("<<minNll<<"), "<<vLik.size()<<" ("<<counter<<") good fits"<<endl;

  open(q2BinIndx, best, true);

  gLim [q2BinIndx] = new TGraph();
  gLiml[q2BinIndx] = new TGraph();
  gLimt[q2BinIndx] = new TGraph();
  gLim0[q2BinIndx] = new TGraph();
  gLim1[q2BinIndx] = new TGraph();
  gLim2[q2BinIndx] = new TGraph();
  gLim3[q2BinIndx] = new TGraph();
  gLim4[q2BinIndx] = new TGraph();
  gLim5[q2BinIndx] = new TGraph();
  gLim6[q2BinIndx] = new TGraph();

  double temp1, temp2;
  int gr_cnt = 0;
  fstream fin (Form("waveP_lim%i.list",q2BinIndx),fstream::in);
  fin.ignore(100,'\n');
  fin.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) {
    fin>>temp2>>temp1;
    if ( temp2!=0 || temp1!=0 ) {
      gLim[q2BinIndx]->SetPoint(gr_cnt,temp1,temp2);
      gr_cnt++;
    }
  }

  gr_cnt=0;
  fstream finl (Form("pdf_llim%i.list",q2BinIndx),fstream::in);
  finl.ignore(100,'\n');
  finl.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) {
    finl>>temp2>>temp1;
    if ( temp2!=0 || temp1!=0 ) {
      gLiml[q2BinIndx]->SetPoint(gr_cnt,temp1,temp2);
      gr_cnt++;
    }
  }

  gr_cnt=0;
  fstream fint (Form("pdf_tlim%i.list",q2BinIndx),fstream::in);
  fint.ignore(100,'\n');
  fint.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) {
    fint>>temp2>>temp1;
    if ( temp2!=0 || temp1!=0 ) {
      gLimt[q2BinIndx]->SetPoint(gr_cnt,temp1,temp2);
      gr_cnt++;
    }
  }

  gr_cnt=0;
  fstream fin0 (Form("pdf_lim%i_%i.list",q2BinIndx,-3),fstream::in);
  fin0.ignore(100,'\n');
  fin0.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) {
    fin0>>temp2>>temp1;
    if ( temp2!=0 || temp1!=0 ) {
      gLim0[q2BinIndx]->SetPoint(gr_cnt,temp1,temp2);
      gr_cnt++;
    }
  }

  gr_cnt=0;
  fstream fin1 (Form("pdf_lim%i_%i.list",q2BinIndx,-2),fstream::in);
  fin1.ignore(100,'\n');
  fin1.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) {
    fin1>>temp2>>temp1;
    if ( temp2!=0 || temp1!=0 ) {
      gLim1[q2BinIndx]->SetPoint(gr_cnt,temp1,temp2);
      gr_cnt++;
    }
  }

  gr_cnt=0;
  fstream fin2 (Form("pdf_lim%i_%i.list",q2BinIndx,-1),fstream::in);
  fin2.ignore(100,'\n');
  fin2.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) {
    fin2>>temp2>>temp1;
    if ( temp2!=0 || temp1!=0 ) {
      gLim2[q2BinIndx]->SetPoint(gr_cnt,temp1,temp2);
      gr_cnt++;
    }
  }

  gr_cnt=0;
  fstream fin3 (Form("pdf_lim%i_%i.list",q2BinIndx,0),fstream::in);
  fin3.ignore(100,'\n');
  fin3.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) {
    fin3>>temp2>>temp1;
    if ( temp2!=0 || temp1!=0 ) {
      gLim3[q2BinIndx]->SetPoint(gr_cnt,temp1,temp2);
      gr_cnt++;
    }
  }

  gr_cnt=0;
  fstream fin4 (Form("pdf_lim%i_%i.list",q2BinIndx,1),fstream::in);
  fin4.ignore(100,'\n');
  fin4.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) {
    fin4>>temp2>>temp1;
    if ( temp2!=0 || temp1!=0 ) {
      gLim4[q2BinIndx]->SetPoint(gr_cnt,temp1,temp2);
      gr_cnt++;
    }
  }

  gr_cnt=0;
  fstream fin5 (Form("pdf_lim%i_%i.list",q2BinIndx,2),fstream::in);
  fin5.ignore(100,'\n');
  fin5.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) {
    fin5>>temp2>>temp1;
    if ( temp2!=0 || temp1!=0 ) {
      gLim5[q2BinIndx]->SetPoint(gr_cnt,temp1,temp2);
      gr_cnt++;
    }
  }

  gr_cnt=0;
  fstream fin6 (Form("pdf_lim%i_%i.list",q2BinIndx,3),fstream::in);
  fin6.ignore(100,'\n');
  fin6.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) {
    fin6>>temp2>>temp1;
    if ( temp2!=0 || temp1!=0 ) {
      gLim6[q2BinIndx]->SetPoint(gr_cnt,temp1,temp2);
      gr_cnt++;
    }
  }

  // int cnt1 = 0;
  // int cnt2 = 0;
  // for (int i=0; i<8100; i++) {
  //   double xP5 = minP5p[q2BinIndx] + (i%90)/(89.0)*(maxP5p[q2BinIndx]-minP5p[q2BinIndx]);
  //   double xP1 = minP1[q2BinIndx] + (i/90)/(89.0)*(maxP1[q2BinIndx]-minP1[q2BinIndx]);
  //   if (good[i]) continue;
  //   // cnt1++;
  //   int j;
  //   for (j=0; j<nLim[q2BinIndx]; j++) if (aP1Lim[j]>xP1) break;
  //   if ( fabs(xP5) > fabs( aP5pLim[j] - (aP1Lim[j]-xP1) * (aP5pLim[j]-aP5pLim[j-1]) / (aP1Lim[j]-aP1Lim[j-1]) ) ) continue;
  //   cout<<q2BinIndx<<" "<<i+1<<endl;
  //   // cnt2++;
  // }
  // cout<<cnt1<<" "<<cnt2<<endl;

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

  legaP1[0] = P1BF[q2BinIndx];
  legaP5p[0] = P5BF[q2BinIndx];

  gBest[q2BinIndx] = new TGraph(1, bestP1, bestP5p);
  gLega[q2BinIndx] = new TGraph(1, legaP1, legaP5p);

  fstream finGEN (Form("GENpoints%i.list",q2BinIndx),fstream::in);
  for (int i=0; i<8; i++) {
    finGEN>>temp2>>temp2>>temp1;
    gBest[q2BinIndx]->SetPoint(i+1,temp2,temp1);
  }

  g2d2s[q2BinIndx] = new TGraph( v2sP1.size(), a2sP1, a2sP5p);
  g2d2s[q2BinIndx]->GetXaxis()->SetLimits(-1,1);
  // g2d1[q2BinIndx]->GetYaxis()->SetLimits(-1.2,.8);
  g2d2s[q2BinIndx]->GetXaxis()->SetTitle("P1");
  g2d2s[q2BinIndx]->GetYaxis()->SetTitle("P'5");
  g2d2s[q2BinIndx]->SetTitle("");
  g2d2s[q2BinIndx]->Draw("AP");
  g2d2s[q2BinIndx]->SetMarkerStyle(7);
  g2d2s[q2BinIndx]->SetMarkerColor(3);
  
  if (v1sP1.size()>0) {
    g2d1s[q2BinIndx] = new TGraph( v1sP1.size(), a1sP1, a1sP5p);
    g2d1s[q2BinIndx]->Draw("Psame");
    g2d1s[q2BinIndx]->SetMarkerStyle(7);
    g2d1s[q2BinIndx]->SetMarkerColor(4);
  }

  /*if (v3sP1.size()>0) {
    g2d3s[q2BinIndx] = new TGraph( v3sP1.size(), a3sP1, a3sP5p);
    g2d3s[q2BinIndx]->Draw("Psame");
    g2d3s[q2BinIndx]->SetMarkerStyle(7);
    g2d3s[q2BinIndx]->SetMarkerColor(2);
  }
  */
  gBest[q2BinIndx]->Draw("sameP");
  gBest[q2BinIndx]->SetMarkerStyle(34);
  gBest[q2BinIndx]->SetMarkerColor(2);
  
  // gLega[q2BinIndx]->Draw("sameP");
  // gLega[q2BinIndx]->SetMarkerStyle(34);
  // gLega[q2BinIndx]->SetMarkerColor(2);

  gLim [q2BinIndx]->Draw("samePL");
  gLim [q2BinIndx]->SetMarkerStyle(7);
  gLim [q2BinIndx]->SetMarkerColor(6);
  gLim [q2BinIndx]->SetLineColor(6);
  gLim [q2BinIndx]->SetLineWidth(2);
  // gLiml[q2BinIndx]->Draw("sameP");
  // gLiml[q2BinIndx]->SetMarkerStyle(7);
  // gLiml[q2BinIndx]->SetMarkerColor(7);
  // gLimt[q2BinIndx]->Draw("sameP");
  // gLimt[q2BinIndx]->SetMarkerStyle(7);
  // gLimt[q2BinIndx]->SetMarkerColor(6);

  gLim0[q2BinIndx]->Draw("sameP");
  gLim0[q2BinIndx]->SetMarkerStyle(7);
  gLim0[q2BinIndx]->SetMarkerColor(1);
  gLim1[q2BinIndx]->Draw("sameP");
  gLim1[q2BinIndx]->SetMarkerStyle(7);
  gLim1[q2BinIndx]->SetMarkerColor(12);
  gLim2[q2BinIndx]->Draw("sameP");
  gLim2[q2BinIndx]->SetMarkerStyle(7);
  gLim2[q2BinIndx]->SetMarkerColor(13);
  gLim3[q2BinIndx]->Draw("sameP");
  gLim3[q2BinIndx]->SetMarkerStyle(7);
  gLim3[q2BinIndx]->SetMarkerColor(14);
  gLim4[q2BinIndx]->Draw("sameP");
  gLim4[q2BinIndx]->SetMarkerStyle(7);
  gLim4[q2BinIndx]->SetMarkerColor(15);
  gLim5[q2BinIndx]->Draw("sameP");
  gLim5[q2BinIndx]->SetMarkerStyle(7);
  gLim5[q2BinIndx]->SetMarkerColor(16);
  gLim6[q2BinIndx]->Draw("sameP");
  gLim6[q2BinIndx]->SetMarkerStyle(7);
  gLim6[q2BinIndx]->SetMarkerColor(17);

  // can[q2BinIndx]->SaveAs(Form("Data_scan_2dpro/bound_b%i.pdf",q2BinIndx));
  // can[q2BinIndx]->SaveAs(Form("Data_scan_2dpro/scan2d_b%i_v6.pdf",q2BinIndx));

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
