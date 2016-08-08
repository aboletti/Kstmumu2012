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

TH1F* h_Fl[9];
TH1F* h_P1[9];
TH1F* h_P5p[9];

TH1F* h_Flb[9];
TH1F* h_P1b[9];
TH1F* h_P5pb[9];

TH1F* h_Fle[9];
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

// float Fl_sv[9] = {0.756,0.856,0.651,0.530,0,0.413,0,0.324,0.293};
// float P1_sv[9] = {0.690,-0.216,0.110,-0.260,0,-0.235,0,-0.780,-0.613};
// float P5p_sv[9] = {0.680,-0.215,-0.236,-0.860,0,-0.875,0,-0.470,-0.682};

float Fl_sv[9] = {0.755113,0.858203,0.658011,0.527329,0,0.419112,0,0.322128,0.267017};
float P1_sv[9] = {0.67814,-0.202957,0.119352,-0.315382,0,-0.27662,0,-0.692714,-0.622906};
float P5p_sv[9] = {0.672541,-0.223946,-0.225956,-0.827075,0,-0.85071,0,-0.553885,-0.614661};

double minNll;
int best;
double edm;
int nfound;

void plot (int q2BinIndx, int i, int scanIndx)
{
    TFile* f = new TFile(Form("LHCb_val/%i/%i/Fitresult5_2.root",i,scanIndx));
    if (f) {
      RooFitResult* fr = (RooFitResult*)f->Get(Form("fitResult_Bin%i",q2BinIndx-1)); if (!fr) {cout<<q2BinIndx<<":\t"<<i<<" something went wrong!!"<<endl; continue;}

      h_Fl[q2BinIndx-1]->Fill( ((RooRealVar*)fr->floatParsFinal().at(0))->getVal() );
      h_P1[q2BinIndx-1]->Fill( ((RooRealVar*)fr->floatParsFinal().at(1))->getVal() );
      h_P5p[q2BinIndx-1]->Fill( ((RooRealVar*)fr->floatParsFinal().at(2))->getVal() );
      h_Fle[q2BinIndx-1]->Fill( ((RooRealVar*)fr->floatParsFinal().at(0))->getError() );
      h_P1e[q2BinIndx-1]->Fill( ((RooRealVar*)fr->floatParsFinal().at(1))->getError() );
      h_P5pe[q2BinIndx-1]->Fill( ((RooRealVar*)fr->floatParsFinal().at(2))->getError() );

      if (fr->status()==0 && fr->covQual()==3) {
	h_Flb[q2BinIndx-1]->Fill( ((RooRealVar*)fr->floatParsFinal().at(0))->getVal() );
	h_P1b[q2BinIndx-1]->Fill( ((RooRealVar*)fr->floatParsFinal().at(1))->getVal() );
	h_P5pb[q2BinIndx-1]->Fill( ((RooRealVar*)fr->floatParsFinal().at(2))->getVal() );
      }

      if ( ((RooRealVar*)fr->floatParsFinal().at(1))->getVal() == ((RooRealVar*)fr->floatParsInit().at(1))->getVal() || ((RooRealVar*)fr->floatParsFinal().at(2))->getVal() == ((RooRealVar*)fr->floatParsInit().at(2))->getVal() ) cout<<q2BinIndx<<":\t"<<i<<"\tParameters unchanged"<<endl;
    } else {
      cout<<"No file for result "<<i<<" something went wrong!"<<endl;
      return;
    }
    f->Close();
    delete f;

  return;
}

void open (int q2BinIndx, int i, int scanIndx, bool print=false)
{
  gErrorIgnoreLevel=4000;

  TFile* f = new TFile(Form("LHCb_val/%i/%i/Fitresult5_2.root",i,scanIndx));
  // TFile* f = new TFile(Form("Scan/%i/Fitresult5_2.root",scanIndx));
  if (!f->IsZombie()) {
    RooFitResult* fr = (RooFitResult*)f->Get(Form("fitResult_Bin%i",q2BinIndx-1));
    if (!fr) return;
    nfound++;

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
    // cout<<"No file for result "<<scanIndx<<endl;
    delete f;
    return;
  }
  f->Close();
  delete f;
}

void read (int q2BinIndx)
{
  if (q2BinIndx==5 || q2BinIndx==7) return;

  h_Fl[q2BinIndx-1] = new TH1F(Form("h_Fl%i",q2BinIndx-1),Form("h_Fl%i",q2BinIndx-1),100,0,1);
  h_P1[q2BinIndx-1] = new TH1F(Form("h_P1%i",q2BinIndx-1),Form("h_P1%i",q2BinIndx-1),100,-1,1);
  h_P5p[q2BinIndx-1] = new TH1F(Form("h_P5p%i",q2BinIndx-1),Form("h_P5p%i",q2BinIndx-1),100,-1,1);
  h_Fle[q2BinIndx-1] = new TH1F(Form("h_Fle%i",q2BinIndx-1),Form("h_Fle%i",q2BinIndx-1),100,0,0.2);
  h_P1e[q2BinIndx-1] = new TH1F(Form("h_P1e%i",q2BinIndx-1),Form("h_P1e%i",q2BinIndx-1),100,0,1);
  h_P5pe[q2BinIndx-1] = new TH1F(Form("h_P5pe%i",q2BinIndx-1),Form("h_P5pe%i",q2BinIndx-1),100,0,1);
  h_Flb[q2BinIndx-1] = new TH1F(Form("h_Flb%i",q2BinIndx-1),Form("h_Flb%i",q2BinIndx-1),100,0,1);
  h_P1b[q2BinIndx-1] = new TH1F(Form("h_P1b%i",q2BinIndx-1),Form("h_P1b%i",q2BinIndx-1),100,-1,1);
  h_P5pb[q2BinIndx-1] = new TH1F(Form("h_P5pb%i",q2BinIndx-1),Form("h_P5pb%i",q2BinIndx-1),100,-1,1);

  float Fl_0=0, P1_0=0, P5p_0=0;
  float Fle_0=0, P1e_0=0, P5pe_0=0;
  // TFile* f_0 = new TFile("../../../B0Analysis/B0KstMuMu/plugins/Fitresult.root");
  // if (f_0) {
  //   const RooArgList& aList= ((RooFitResult*)f_0->Get(Form("fitResult_Bin%i",q2BinIndx-1)))->floatParsFinal();
  //   Fl_0 = ((RooRealVar*)aList.at(0))->getVal();
  //   P1_0 = ((RooRealVar*)aList.at(1))->getVal();
  //   P5p_0 = ((RooRealVar*)aList.at(2))->getVal();
  //   Fle_0 = ((RooRealVar*)aList.at(0))->getError();
  //   P1e_0 = ((RooRealVar*)aList.at(1))->getError();
  //   P5pe_0 = ((RooRealVar*)aList.at(2))->getError();
  // } else {
  //   cout<<"No file for result main"<<endl;
  //   return;
  // }
  // f_0->Close();
  // delete f_0;

  float Fl_0b, P1_0b, P5p_0b;
  TFile* f_0 = new TFile("Fitresult.root");
  if (f_0) {
    const RooArgList& aList= ((RooFitResult*)f_0->Get(Form("fitResult_Bin%i",q2BinIndx-1)))->floatParsFinal();
    Fl_0b = ((RooRealVar*)aList.at(0))->getVal();
    P1_0b = ((RooRealVar*)aList.at(1))->getVal();
    P5p_0b = ((RooRealVar*)aList.at(2))->getVal();
  } else {
    cout<<"No file for result main"<<endl;
    return;
  }
  f_0->Close();
  delete f_0;

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
  
  for (int i=1; i<=80; i++) {

    minNll = 9999999;
    nfound=0;
    for (int cnt=100; cnt<110; cnt++) open(q2BinIndx, i, cnt);
    
    if (minNll==9999999) {
      // cout<<"No good fits at bin "<<q2BinIndx<<" toy "<<i<<" (found "<<nfound<<"/10)"<<endl;
      continue;	    
    }
    // cout<<"Best "<<q2BinIndx<<": "<<best<<" ("<<minNll<<")"<<endl;
    // cout<<q2BinIndx<<" "<<i<<endl;
    plot(q2BinIndx, i, best);

  }

  // cout<<q2BinIndx<<" Fl :\t"<<h_Fl[q2BinIndx-1]->GetMean()<<"\t+-"<<h_Fl[q2BinIndx-1]->GetRMS()<<endl;
  // cout<<q2BinIndx<<" P1 :\t"<<h_P1[q2BinIndx-1]->GetMean()<<"\t+-"<<h_P1[q2BinIndx-1]->GetRMS()<<endl;
  // cout<<q2BinIndx<<" P5p:\t"<<h_P5p[q2BinIndx-1]->GetMean()<<"\t+-"<<h_P5p[q2BinIndx-1]->GetRMS()<<endl;

  cout<<"Converging "<<h_Flb[q2BinIndx-1]->GetEntries()<<" / "<<h_Fl[q2BinIndx-1]->GetEntries()<<endl;
  cout<<q2BinIndx<<" Fl :\t"<<h_Flb[q2BinIndx-1]->GetMean()<<"\t+-"<<h_Flb[q2BinIndx-1]->GetRMS()<<endl;
  cout<<q2BinIndx<<" P1 :\t"<<h_P1b[q2BinIndx-1]->GetMean()<<"\t+-"<<h_P1b[q2BinIndx-1]->GetRMS()<<endl;
  cout<<q2BinIndx<<" P5p:\t"<<h_P5pb[q2BinIndx-1]->GetMean()<<"\t+-"<<h_P5pb[q2BinIndx-1]->GetRMS()<<endl;
   
  c[q2BinIndx-1] = new TCanvas (Form("c%i",q2BinIndx-1),Form("c%i",q2BinIndx-1));
  c[q2BinIndx-1]->Divide(2,3);

  c[q2BinIndx-1]->cd(1);
  h_Fl[q2BinIndx-1]->Draw();
  h_Fl[q2BinIndx-1]->SetLineWidth(1);
  // h_Flb[q2BinIndx-1]->Draw("same");
  h_Flb[q2BinIndx-1]->SetLineColor(2);
  h_Flb[q2BinIndx-1]->SetLineWidth(1);
  l1[q2BinIndx-1] = new TLine(Fl_sv[q2BinIndx-1],0,Fl_sv[q2BinIndx-1], h_Fl[q2BinIndx-1]->GetMaximum());
  l1[q2BinIndx-1]->Draw();
  l1[q2BinIndx-1]->SetLineColor(8);
  l1[q2BinIndx-1]->SetLineWidth(1);
  l1b[q2BinIndx-1] = new TLine(Fl_0b,0,Fl_0b, h_Fl[q2BinIndx-1]->GetMaximum());
  l1b[q2BinIndx-1]->Draw();
  l1b[q2BinIndx-1]->SetLineColor(2);
  l1b[q2BinIndx-1]->SetLineWidth(1);

  c[q2BinIndx-1]->cd(3);
  h_P1[q2BinIndx-1]->Draw();
  h_P1[q2BinIndx-1]->SetLineWidth(1);
  // h_P1b[q2BinIndx-1]->Draw("same");
  h_P1b[q2BinIndx-1]->SetLineColor(2);
  h_P1b[q2BinIndx-1]->SetLineWidth(1);
  l2[q2BinIndx-1] = new TLine(P1_sv[q2BinIndx-1],0,P1_sv[q2BinIndx-1], h_P1[q2BinIndx-1]->GetMaximum());
  l2[q2BinIndx-1]->Draw();
  l2[q2BinIndx-1]->SetLineColor(8);
  l2[q2BinIndx-1]->SetLineWidth(1);
  l2b[q2BinIndx-1] = new TLine(P1_0b,0,P1_0b, h_P1[q2BinIndx-1]->GetMaximum());
  l2b[q2BinIndx-1]->Draw();
  l2b[q2BinIndx-1]->SetLineColor(2);
  l2b[q2BinIndx-1]->SetLineWidth(1);

  c[q2BinIndx-1]->cd(5);
  h_P5p[q2BinIndx-1]->Draw();
  h_P5p[q2BinIndx-1]->SetLineWidth(1);
  // h_P5pb[q2BinIndx-1]->Draw("same");
  h_P5pb[q2BinIndx-1]->SetLineColor(2);
  h_P5pb[q2BinIndx-1]->SetLineWidth(1);
  l3[q2BinIndx-1] = new TLine(P5p_sv[q2BinIndx-1],0,P5p_sv[q2BinIndx-1], h_P5p[q2BinIndx-1]->GetMaximum());
  l3[q2BinIndx-1]->Draw();
  l3[q2BinIndx-1]->SetLineColor(8);
  l3[q2BinIndx-1]->SetLineWidth(1);
  l3b[q2BinIndx-1] = new TLine(P5p_0b,0,P5p_0b, h_P5p[q2BinIndx-1]->GetMaximum());
  l3b[q2BinIndx-1]->Draw();
  l3b[q2BinIndx-1]->SetLineColor(2);
  l3b[q2BinIndx-1]->SetLineWidth(1);

  c[q2BinIndx-1]->cd(2);
  h_Fle[q2BinIndx-1]->Draw();
  // l4[q2BinIndx-1] = new TLine(Fle_0,0,Fle_0, h_Fle[q2BinIndx-1]->GetMaximum());
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

  c[q2BinIndx-1]->SaveAs(Form("plots/LHCb_bin%i.pdf",q2BinIndx));

  return;
}

void read_LHCb (int q2BinIndx = 0)
{
  if (q2BinIndx==0) for (q2BinIndx=1;q2BinIndx<10;q2BinIndx++) read(q2BinIndx);
  else if (q2BinIndx>0 && q2BinIndx<10) read(q2BinIndx);
  
  return;
}
