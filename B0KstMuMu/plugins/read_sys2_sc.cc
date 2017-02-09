#include <iostream>

#include "TH1.h"
#include "TLine.h"
#include "TCanvas.h"
// #include ".h"
// #include ".h"
#include "TFile.h"
#include "TError.h"

// #include "RooRealVar.h"
// #include "RooFitResult.h"
// #include "RooArgList.h"
// #include ".h"
// #include ".h"
// #include ".h"

using namespace std;
using namespace RooFit;

TH1F* h_As5[9];
TH1F* h_P1[9];
TH1F* h_P5p[9];

TH1F* h_As5b[9];
TH1F* h_P1b[9];
TH1F* h_P5pb[9];

TH1F* h_As5e[9];
TH1F* h_P1e[9];
TH1F* h_P5pe[9];

TCanvas* c[9];
TLine* l1[9];
TLine* l2[9];
TLine* l3[9];
TLine* l4[9];
TLine* l5[9];
TLine* l6[9];
TLine* l1b[9];
TLine* l2b[9];
TLine* l3b[9];

float err_cut = 0.1;

// float As5_sv[9] = {0.756,0.856,0.651,0.530,0,0.413,0,0.324,0.293};
// float P1_sv[9] = {0.690,-0.216,0.110,-0.260,0,-0.235,0,-0.780,-0.613};
// float P5p_sv[9] = {0.680,-0.215,-0.236,-0.860,0,-0.875,0,-0.470,-0.682};

float As5_sv[9] = {0.755113,0.858203,0.658011,0.527329,0,0.419112,0,0.322128,0.267017};
float P1_sv[9] = {0.67814,-0.202957,0.119352,-0.315382,0,-0.27662,0,-0.692714,-0.622906};
float P5p_sv[9] = {0.672541,-0.223946,-0.225956,-0.827075,0,-0.85071,0,-0.553885,-0.614661};

double minNll;
int best;
double edm;
int nfound;
bool isDone;

int max_scan = 300;

void plot (int q2BinIndx, int i, int scanIndx)
{
    TFile* f = new TFile(Form("Sys2_data/%05i/%i/Fitresult4_2.root",i,scanIndx));
    if (f) {
      RooFitResult* fr = (RooFitResult*)f->Get(Form("fitResult_Bin%i",q2BinIndx-1)); if (!fr) {cout<<q2BinIndx<<":\t"<<i<<" something went wrong!!"<<endl; f->Close(); return;}

      h_As5[q2BinIndx-1]->Fill( ((RooRealVar*)fr->floatParsFinal().at(0))->getVal() );
      h_P1[q2BinIndx-1]->Fill( ((RooRealVar*)fr->floatParsFinal().at(1))->getVal() );
      h_P5p[q2BinIndx-1]->Fill( ((RooRealVar*)fr->floatParsFinal().at(2))->getVal() );
      h_As5e[q2BinIndx-1]->Fill( ((RooRealVar*)fr->floatParsFinal().at(0))->getError() );
      h_P1e[q2BinIndx-1]->Fill( ((RooRealVar*)fr->floatParsFinal().at(1))->getError() );
      h_P5pe[q2BinIndx-1]->Fill( ((RooRealVar*)fr->floatParsFinal().at(2))->getError() );

      if (fr->status()==0 && fr->covQual()==3) {
	h_As5b[q2BinIndx-1]->Fill( ((RooRealVar*)fr->floatParsFinal().at(0))->getVal() );
	h_P1b[q2BinIndx-1]->Fill( ((RooRealVar*)fr->floatParsFinal().at(1))->getVal() );
	h_P5pb[q2BinIndx-1]->Fill( ((RooRealVar*)fr->floatParsFinal().at(2))->getVal() );
      }

      if ( ((RooRealVar*)fr->floatParsFinal().at(1))->getVal() == ((RooRealVar*)fr->floatParsInit().at(1))->getVal() || ((RooRealVar*)fr->floatParsFinal().at(2))->getVal() == ((RooRealVar*)fr->floatParsInit().at(2))->getVal() ) cout<<q2BinIndx<<":\t"<<i<<"\tParameters unchanged"<<endl;
    } else {
      cout<<"No file for result "<<i<<" something went wrong!"<<endl;
      delete f;
      return;
    }
    f->Close();
    delete f;

  return;
}

void open (int q2BinIndx, int i, int scanIndx, bool print=false)
{
  gErrorIgnoreLevel=4000;
  
  TString filename = Form("Sys2_data/%05i/%i/Fitresult4_2.root",i,scanIndx);
  if (gSystem->AccessPathName(filename)) { cout<<"No file"<<endl; return;}
    
  // TFile* f = new TFile(Form("Closure_sys2_rt_0.5/%05i/%i/Fitresult5_2.root",i,scanIndx));
  TFile* f = new TFile(filename);
  // TFile* f = new TFile(Form("Scan/%i/Fitresult5_2.root",scanIndx));
  if (!f->IsZombie()) {
    RooFitResult* fr = (RooFitResult*)f->Get(Form("fitResult_Bin%i",q2BinIndx-1));
    if (!fr) {
      cout<<"!fr"<<endl;
      f->Close();
      delete f;
      return;
    }
    if (fr->status()==0 && fr->covQual()==3) nfound++;
    isDone = true;

    if ( !print ) {
      if (fr->status()==0 && fr->covQual()==3 && (minNll>fr->minNll() || (minNll==fr->minNll() && edm>fr->edm()))) {
	minNll = fr->minNll();
	best = scanIndx;
	edm = fr->edm();
      }
    } else {
      cout<<scanIndx<<": "<<endl;
      cout<<"(fnc="<<fr->minNll()<<",edm="<<fr->edm()<<",stat="<<fr->status()<<",covQ="<<fr->covQual()<<")"<<endl;
      cout<<"As5\t"<<((RooRealVar*)fr->floatParsInit().at(0))->getVal()<<"  \t"<<((RooRealVar*)fr->floatParsFinal().at(0))->getVal()
	  <<"\nP1\t"<<((RooRealVar*)fr->floatParsInit().at(1))->getVal()<<"  \t"<<((RooRealVar*)fr->floatParsFinal().at(1))->getVal()
	  <<"\nP5p\t"<<((RooRealVar*)fr->floatParsInit().at(2))->getVal()<<"  \t"<<((RooRealVar*)fr->floatParsFinal().at(2))->getVal()<<endl;
    }
  } else {
    cout<<"IsZombie"<<endl;
    delete f;
    return;
  }
  f->Close();
  delete f;

  return;
}

void read (int q2BinIndx)
{
  if (q2BinIndx==5 || q2BinIndx==7) return;

  h_As5[q2BinIndx-1] = new TH1F(Form("h_As5%i",q2BinIndx-1),Form("h_As5%i",q2BinIndx-1),100,-1,1);
  h_P1[q2BinIndx-1] = new TH1F(Form("h_P1%i",q2BinIndx-1),Form("h_P1%i",q2BinIndx-1),100,-1,1);
  h_P5p[q2BinIndx-1] = new TH1F(Form("h_P5p%i",q2BinIndx-1),Form("h_P5p%i",q2BinIndx-1),100,-1.5,1.5);
  h_As5e[q2BinIndx-1] = new TH1F(Form("h_As5e%i",q2BinIndx-1),Form("h_As5e%i",q2BinIndx-1),100,0,0.2);
  h_P1e[q2BinIndx-1] = new TH1F(Form("h_P1e%i",q2BinIndx-1),Form("h_P1e%i",q2BinIndx-1),100,0,1);
  h_P5pe[q2BinIndx-1] = new TH1F(Form("h_P5pe%i",q2BinIndx-1),Form("h_P5pe%i",q2BinIndx-1),100,0,1);
  h_As5b[q2BinIndx-1] = new TH1F(Form("h_As5b%i",q2BinIndx-1),Form("h_As5b%i",q2BinIndx-1),100,-1,1);
  h_P1b[q2BinIndx-1] = new TH1F(Form("h_P1b%i",q2BinIndx-1),Form("h_P1b%i",q2BinIndx-1),100,-1.5,1.5);
  h_P5pb[q2BinIndx-1] = new TH1F(Form("h_P5pb%i",q2BinIndx-1),Form("h_P5pb%i",q2BinIndx-1),100,-1,1);

  float As5_0=0, P1_0=0, P5p_0=0;
  float As5e_0=0, P1e_0=0, P5pe_0=0;
  // TFile* f_0 = new TFile("../../../B0Analysis/B0KstMuMu/plugins/Fitresult.root");
  // if (f_0) {
  //   const RooArgList& aList= ((RooFitResult*)f_0->Get(Form("fitResult_Bin%i",q2BinIndx-1)))->floatParsFinal();
  //   As5_0 = ((RooRealVar*)aList.at(0))->getVal();
  //   P1_0 = ((RooRealVar*)aList.at(1))->getVal();
  //   P5p_0 = ((RooRealVar*)aList.at(2))->getVal();
  //   As5e_0 = ((RooRealVar*)aList.at(0))->getError();
  //   P1e_0 = ((RooRealVar*)aList.at(1))->getError();
  //   P5pe_0 = ((RooRealVar*)aList.at(2))->getError();
  // } else {
  //   cout<<"No file for result main"<<endl;
  //   return;
  // }
  // f_0->Close();
  // delete f_0;

  float As5_0b=0, P1_0b=0, P5p_0b=0;
  // TFile* f_0 = new TFile("Fitresult.root");
  // if (f_0) {
  //   const RooArgList& aList= ((RooFitResult*)f_0->Get(Form("fitResult_Bin%i",q2BinIndx-1)))->floatParsFinal();
  //   As5_0b = ((RooRealVar*)aList.at(0))->getVal();
  //   P1_0b = ((RooRealVar*)aList.at(1))->getVal();
  //   P5p_0b = ((RooRealVar*)aList.at(2))->getVal();
  // } else {
  //   cout<<"No file for result main"<<endl;
  //   return;
  // }
  // f_0->Close();
  // delete f_0;

  // cout<<1<<endl;
  // for (int i=0; i<=100; i++) { 
  //   TFile* f = new TFile(Form("../../../B0Analysis/B0KstMuMu/plugins/Closure_sys2_rt_0.5/%05i/Fitresult5.root",i));
  //   if (f) {
  //     RooFitResult* fr = (RooFitResult*)f->Get(Form("fitResult_Bin%i",q2BinIndx-1)); if (!fr) {cout<<q2BinIndx<<":\t"<<i<<endl; continue;}
  //     // if (fr->status()!=0 || fr->covQual()!=3) continue;
  //     //cout<<fr->covQual()<<" ";

  //     if ( ((RooRealVar*)fr->floatParsFinal().at(1))->getVal() == ((RooRealVar*)fr->floatParsInit().at(1))->getVal() || ((RooRealVar*)fr->floatParsFinal().at(2))->getVal() == ((RooRealVar*)fr->floatParsInit().at(2))->getVal() ) cout<<q2BinIndx<<":\t"<<i<<"\tParameters unchanged"<<endl;
  //   } else {
  //     cout<<"No file for result "<<i<<endl;
  //     return;
  //   }
  //   f->Close();
  //   delete f;
  // }
  // cout<<2<<endl;
  
  int nDone = 0;
  for (int i=0; i<=100; i++) { //cout<<i<<endl;

    minNll = 9999999;
    nfound=0;
    isDone = false;
    for (int cnt=200; cnt<=249; cnt++) {
      cout<<i<<" "<<cnt<<endl;
      open(q2BinIndx, i, cnt);
    }
    for (int cnt=275; cnt<=max_scan; cnt++) {
      cout<<i<<" "<<cnt<<endl;
      open(q2BinIndx, i, cnt);
    }
    if (isDone) nDone++;
   
    // cout<<nfound<<endl;

    if (minNll==9999999 || nfound<3) {
      // cout<<"No good fits at bin "<<q2BinIndx<<" toy "<<i<<" (found "<<nfound<<"/10)"<<endl;
      continue;	    
    }
    // cout<<"Best "<<q2BinIndx<<": "<<best<<" ("<<minNll<<")"<<endl;
    // cout<<q2BinIndx<<" "<<i<<endl;
    plot(q2BinIndx, i, best);

  }

  // cout<<q2BinIndx<<" As5 :\t"<<h_As5[q2BinIndx-1]->GetMean()<<"\t+-"<<h_As5[q2BinIndx-1]->GetRMS()<<endl;
  // cout<<q2BinIndx<<" P1 :\t"<<h_P1[q2BinIndx-1]->GetMean()<<"\t+-"<<h_P1[q2BinIndx-1]->GetRMS()<<endl;
  // cout<<q2BinIndx<<" P5p:\t"<<h_P5p[q2BinIndx-1]->GetMean()<<"\t+-"<<h_P5p[q2BinIndx-1]->GetRMS()<<endl;

  cout<<"Converging "<<h_As5b[q2BinIndx-1]->GetEntries()<<" / "<<nDone<<endl;
  cout<<q2BinIndx<<" As5:\t"<<h_As5b[q2BinIndx-1]->GetMean()<<"\t+-"<<h_As5b[q2BinIndx-1]->GetRMS()<<endl;
  cout<<q2BinIndx<<" P1 :\t"<<h_P1b[q2BinIndx-1]->GetMean()<<"\t+-"<<h_P1b[q2BinIndx-1]->GetRMS()<<endl;
  cout<<q2BinIndx<<" P5p:\t"<<h_P5pb[q2BinIndx-1]->GetMean()<<"\t+-"<<h_P5pb[q2BinIndx-1]->GetRMS()<<endl;
   
  c[q2BinIndx-1] = new TCanvas (Form("c%i",q2BinIndx-1),Form("c%i",q2BinIndx-1));
  c[q2BinIndx-1]->Divide(2,3);

  c[q2BinIndx-1]->cd(1);
  h_As5[q2BinIndx-1]->Draw();
  h_As5[q2BinIndx-1]->SetLineWidth(1);
  // h_As5b[q2BinIndx-1]->Draw("same");
  h_As5b[q2BinIndx-1]->SetLineColor(2);
  h_As5b[q2BinIndx-1]->SetLineWidth(1);
  // l1[q2BinIndx-1] = new TLine(As5_sv[q2BinIndx-1],0,As5_sv[q2BinIndx-1], h_As5[q2BinIndx-1]->GetMaximum());
  // l1[q2BinIndx-1]->Draw();
  // l1[q2BinIndx-1]->SetLineColor(8);
  // l1[q2BinIndx-1]->SetLineWidth(1);
  // l1b[q2BinIndx-1] = new TLine(As5_0b,0,As5_0b, h_As5[q2BinIndx-1]->GetMaximum());
  // l1b[q2BinIndx-1]->Draw();
  // l1b[q2BinIndx-1]->SetLineColor(2);
  // l1b[q2BinIndx-1]->SetLineWidth(1);

  c[q2BinIndx-1]->cd(3);
  h_P1[q2BinIndx-1]->Draw();
  h_P1[q2BinIndx-1]->SetLineWidth(1);
  // h_P1b[q2BinIndx-1]->Draw("same");
  h_P1b[q2BinIndx-1]->SetLineColor(2);
  h_P1b[q2BinIndx-1]->SetLineWidth(1);
  // l2[q2BinIndx-1] = new TLine(P1_sv[q2BinIndx-1],0,P1_sv[q2BinIndx-1], h_P1[q2BinIndx-1]->GetMaximum());
  // l2[q2BinIndx-1]->Draw();
  // l2[q2BinIndx-1]->SetLineColor(8);
  // l2[q2BinIndx-1]->SetLineWidth(1);
  // l2b[q2BinIndx-1] = new TLine(P1_0b,0,P1_0b, h_P1[q2BinIndx-1]->GetMaximum());
  // l2b[q2BinIndx-1]->Draw();
  // l2b[q2BinIndx-1]->SetLineColor(2);
  // l2b[q2BinIndx-1]->SetLineWidth(1);

  c[q2BinIndx-1]->cd(5);
  h_P5p[q2BinIndx-1]->Draw();
  h_P5p[q2BinIndx-1]->SetLineWidth(1);
  // h_P5pb[q2BinIndx-1]->Draw("same");
  h_P5pb[q2BinIndx-1]->SetLineColor(2);
  h_P5pb[q2BinIndx-1]->SetLineWidth(1);
  // l3[q2BinIndx-1] = new TLine(P5p_sv[q2BinIndx-1],0,P5p_sv[q2BinIndx-1], h_P5p[q2BinIndx-1]->GetMaximum());
  // l3[q2BinIndx-1]->Draw();
  // l3[q2BinIndx-1]->SetLineColor(8);
  // l3[q2BinIndx-1]->SetLineWidth(1);
  // l3b[q2BinIndx-1] = new TLine(P5p_0b,0,P5p_0b, h_P5p[q2BinIndx-1]->GetMaximum());
  // l3b[q2BinIndx-1]->Draw();
  // l3b[q2BinIndx-1]->SetLineColor(2);
  // l3b[q2BinIndx-1]->SetLineWidth(1);

  c[q2BinIndx-1]->cd(2);
  h_As5e[q2BinIndx-1]->Draw();
  // l4[q2BinIndx-1] = new TLine(As5e_0,0,As5e_0, h_As5e[q2BinIndx-1]->GetMaximum());
  // l4[q2BinIndx-1]->Draw();
  // l4[q2BinIndx-1]->SetLineWidth(1);

  c[q2BinIndx-1]->cd(4);
  h_P1e[q2BinIndx-1]->Draw();
  // l5[q2BinIndx-1] = new TLine(P1e_0,0,P1e_0, h_P1e[q2BinIndx-1]->GetMaximum());
  // l5[q2BinIndx-1]->Draw();
  // l5[q2BinIndx-1]->SetLineWidth(1);

  c[q2BinIndx-1]->cd(6);
  h_P5pe[q2BinIndx-1]->Draw();
  // l6[q2BinIndx-1] = new TLine(P5pe_0,0,P5pe_0, h_P5pe[q2BinIndx-1]->GetMaximum());
  // l6[q2BinIndx-1]->Draw();
  // l6[q2BinIndx-1]->SetLineWidth(1);

  c[q2BinIndx-1]->SaveAs(Form("plots/Sys2_data_b%i.pdf",q2BinIndx-1));

  return;
}

void read_sys2_sc (int q2BinIndx = 0)
{
  if (q2BinIndx==0) for (q2BinIndx=1;q2BinIndx<10;q2BinIndx++) read(q2BinIndx);
  else if (q2BinIndx>0 && q2BinIndx<10) read(q2BinIndx);
  
  return;
}
