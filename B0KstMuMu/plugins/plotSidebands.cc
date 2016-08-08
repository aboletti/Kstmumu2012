#include <TROOT.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TF12.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TText.h>
#include <TFitResult.h>
#include <Math/Functor.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

using std::cout;
using std::endl;
using std::string;
using std::stringstream;
using std::vector;
using std::ios_base;
using std::make_pair;

double sigMean[9] = {5.26619,5.27389,5.27604,5.27553,5.27831,5.27451,5.27826,5.27621,5.27964};
double sigSig1[9] = {0.0304184,0.0308886,0.0311103,0.0301922,0.0302283,0.0311809,0.032646,0.0316837,0.031461};
double sigSig2[9] = {0.0921757,0.0979804,0.104148,0.0890352,0.0638594,0.0942737,0.062256,0.094884,0.0766426};
double sigSigf[9] = {0.831529,0.844362,0.857988,0.85582,0.829155,0.872235,0.910166,0.863067,0.828301};
double q2Bins[10] = {1,2,4.3,6,8.68,10.09,12.86,14.18,16,19};

TH1F* plotK_L0 [11];
TH1F* plotK_L1 [11];
TH1F* plotK_P0 [11];
TH1F* plotK_P1 [11];
TH1F* plotK_M0 [11];
TH1F* plotK_M1 [11];
TH1F* plotL_K0 [11];
TH1F* plotL_K1 [11];
TH1F* plotL_P0 [11];
TH1F* plotL_P1 [11];
TH1F* plotL_M0 [11];
TH1F* plotL_M1 [11];
TH1F* plotP_L0 [11];
TH1F* plotP_L1 [11];
TH1F* plotP_K0 [11];
TH1F* plotP_K1 [11];
TH1F* plotP_M0 [11];
TH1F* plotP_M1 [11];
TH1F* plotM_L0 [11];
TH1F* plotM_L1 [11];
TH1F* plotM_P0 [11];
TH1F* plotM_P1 [11];
TH1F* plotM_K0 [11];
TH1F* plotM_K1 [11];

TCanvas* can [11];

bool PsiRejection1 (double M, double Q, double QE) {
  if (fabs(Q - 3.097) < 3 * QE)
    return true;
  return false;
}

bool PsiRejection2 (double M, double Q, double QE) {
  if (fabs(Q - 3.686) < 3 * QE)
    return true;
  return false;
}

bool PsiRejection3 (double M, double Q, double QE) {
  if (Q < 3.097 &&
      fabs((M - 5.28) - (Q - 3.097)) > 0.16 &&
      fabs((M - 5.28) - (Q - 3.686)) > 0.06)
    return true;
	  
  if (Q > 3.686 &&
      fabs((M - 5.28) - (Q - 3.097)) > 0.06 &&
      fabs((M - 5.28) - (Q - 3.686)) > 0.03)
    return true;

  if (Q > 3.097 && Q < 3.686 &&
      fabs((M - 5.28) - (Q - 3.097)) > 0.06 &&
      fabs((M - 5.28) - (Q - 3.686)) > 0.06)
    return true;

  return false;
}

void fillHistos (double K, double L, double P, double M, int indx) {
  if (K<0) {
    plotL_K0[indx]->Fill(L);
    plotP_K0[indx]->Fill(P);
    plotM_K0[indx]->Fill(M);
  } else {
    plotL_K1[indx]->Fill(L);
    plotP_K1[indx]->Fill(P);
    plotM_K1[indx]->Fill(M);
  }
  if (L<0.5) {
    plotK_L0[indx]->Fill(K);
    plotP_L0[indx]->Fill(P);
    plotM_L0[indx]->Fill(M);
  } else {
    plotK_L1[indx]->Fill(K);
    plotP_L1[indx]->Fill(P);
    plotM_L1[indx]->Fill(M);
  }
  if (P<TMath::Pi()/2.) {
    plotL_P0[indx]->Fill(L);
    plotK_P0[indx]->Fill(K);
    plotM_P0[indx]->Fill(M);
  } else {
    plotL_P1[indx]->Fill(L);
    plotK_P1[indx]->Fill(K);
    plotM_P1[indx]->Fill(M);
  }
  if (M<5.275) {
    plotL_M0[indx]->Fill(L);
    plotP_M0[indx]->Fill(P);
    plotK_M0[indx]->Fill(K);
  } else {
    plotL_M1[indx]->Fill(L);
    plotP_M1[indx]->Fill(P);
    plotK_M1[indx]->Fill(K);
  } //cout<<indx<<" ";
}  

void plotSidebands () {

  TFile* f = TFile::Open("dcap://t2-srm-02.lnl.infn.it//pnfs/lnl.infn.it/data/cms/store/user/slacapra/B0KsMuMu/Data2012B0KstMuMuResults/singleCand_B0ToKstMuMu_Data2012ABCD_NTuples.root");
        
  TTree* t = (TTree*)f->Get("B0KstMuMu/B0KstMuMuNTuple");

  // gStyle->SetLineWidth(3);
  // TH1::SetDefaultSumw2();
  for (int i=0; i<11; i++) {
    plotK_L0 [i] = new TH1F(Form("plotK_L0%i",i),Form("plotK_L0%i",i),6,-1,1);
    plotK_L1 [i] = new TH1F(Form("plotK_L1%i",i),Form("plotK_L1%i",i),6,-1,1);
    plotK_P0 [i] = new TH1F(Form("plotK_P0%i",i),Form("plotK_P0%i",i),6,-1,1);
    plotK_P1 [i] = new TH1F(Form("plotK_P1%i",i),Form("plotK_P1%i",i),6,-1,1);
    plotK_M0 [i] = new TH1F(Form("plotK_M0%i",i),Form("plotK_M0%i",i),6,-1,1);
    plotK_M1 [i] = new TH1F(Form("plotK_M1%i",i),Form("plotK_M1%i",i),6,-1,1);
    plotL_K0 [i] = new TH1F(Form("plotL_K0%i",i),Form("plotL_K0%i",i),6,0,1);
    plotL_K1 [i] = new TH1F(Form("plotL_K1%i",i),Form("plotL_K1%i",i),6,0,1);
    plotL_P0 [i] = new TH1F(Form("plotL_P0%i",i),Form("plotL_P0%i",i),6,0,1);
    plotL_P1 [i] = new TH1F(Form("plotL_P1%i",i),Form("plotL_P1%i",i),6,0,1);
    plotL_M0 [i] = new TH1F(Form("plotL_M0%i",i),Form("plotL_M0%i",i),6,0,1);
    plotL_M1 [i] = new TH1F(Form("plotL_M1%i",i),Form("plotL_M1%i",i),6,0,1);
    plotP_L0 [i] = new TH1F(Form("plotP_L0%i",i),Form("plotP_L0%i",i),6,0,TMath::Pi());
    plotP_L1 [i] = new TH1F(Form("plotP_L1%i",i),Form("plotP_L1%i",i),6,0,TMath::Pi());
    plotP_K0 [i] = new TH1F(Form("plotP_K0%i",i),Form("plotP_K0%i",i),6,0,TMath::Pi());
    plotP_K1 [i] = new TH1F(Form("plotP_K1%i",i),Form("plotP_K1%i",i),6,0,TMath::Pi());
    plotP_M0 [i] = new TH1F(Form("plotP_M0%i",i),Form("plotP_M0%i",i),6,0,TMath::Pi());
    plotP_M1 [i] = new TH1F(Form("plotP_M1%i",i),Form("plotP_M1%i",i),6,0,TMath::Pi());
    plotM_L0 [i] = new TH1F(Form("plotM_L0%i",i),Form("plotM_L0%i",i),10,5,5.56);
    plotM_L1 [i] = new TH1F(Form("plotM_L1%i",i),Form("plotM_L1%i",i),10,5,5.56);
    plotM_P0 [i] = new TH1F(Form("plotM_P0%i",i),Form("plotM_P0%i",i),10,5,5.56);
    plotM_P1 [i] = new TH1F(Form("plotM_P1%i",i),Form("plotM_P1%i",i),10,5,5.56);
    plotM_K0 [i] = new TH1F(Form("plotM_K0%i",i),Form("plotM_K0%i",i),10,5,5.56);
    plotM_K1 [i] = new TH1F(Form("plotM_K1%i",i),Form("plotM_K1%i",i),10,5,5.56);
    plotK_L0 [i]->SetLineWidth(3);
    plotK_L1 [i]->SetLineWidth(3);
    plotK_P0 [i]->SetLineWidth(3);
    plotK_P1 [i]->SetLineWidth(3);
    plotK_M0 [i]->SetLineWidth(3);
    plotK_M1 [i]->SetLineWidth(3);
    plotL_K0 [i]->SetLineWidth(3);
    plotL_K1 [i]->SetLineWidth(3);
    plotL_P0 [i]->SetLineWidth(3);
    plotL_P1 [i]->SetLineWidth(3);
    plotL_M0 [i]->SetLineWidth(3);
    plotL_M1 [i]->SetLineWidth(3);
    plotP_L0 [i]->SetLineWidth(3);
    plotP_L1 [i]->SetLineWidth(3);
    plotP_K0 [i]->SetLineWidth(3);
    plotP_K1 [i]->SetLineWidth(3);
    plotP_M0 [i]->SetLineWidth(3);
    plotP_M1 [i]->SetLineWidth(3);
    plotM_L0 [i]->SetLineWidth(3);
    plotM_L1 [i]->SetLineWidth(3);
    plotM_P0 [i]->SetLineWidth(3);
    plotM_P1 [i]->SetLineWidth(3);
    plotM_K0 [i]->SetLineWidth(3);
    plotM_K1 [i]->SetLineWidth(3);
  }

  double K,L,P,M,Q,QE;
  std::vector<double> *QV;
  std::vector<double> *QEV;
  t->SetBranchAddress("CosThetaKArb"      ,&K); 
  t->SetBranchAddress("CosThetaMuArb"     ,&L);
  t->SetBranchAddress("PhiKstMuMuPlaneArb",&P);
  t->SetBranchAddress("B0MassArb"         ,&M);
  t->SetBranchAddress("mumuMass"          ,&QV); 
  t->SetBranchAddress("mumuMassE"         ,&QEV);
  int Q2indx;
  double mumuq2;

  int n_ev = t->GetEntries();
  for (int i=0; i<n_ev; i++) {
    t->GetEntry(i);

    if (M<5 || M>5.56) continue;

    Q = QV->at(0);
    QE = QEV->at(0);
    // if (Q>0.5) cout<<K<<" "<<L<<" "<<P<<" "<<M<<" "<<Q<<endl;
    mumuq2 = Q*Q;
    if (mumuq2<q2Bins[0] || mumuq2>q2Bins[9]) continue;
    for (Q2indx=0; Q2indx<9; Q2indx++) if (mumuq2<q2Bins[Q2indx+1]) break;
    // cout<<Q2indx<<endl;

    if (M<sigMean[Q2indx]+3*sqrt(sigSigf[Q2indx]*sigSig1[Q2indx]*sigSig1[Q2indx]+(1-sigSigf[Q2indx])*sigSig2[Q2indx]*sigSig2[Q2indx]) &&
	M>sigMean[Q2indx]-3*sqrt(sigSigf[Q2indx]*sigSig1[Q2indx]*sigSig1[Q2indx]+(1-sigSigf[Q2indx])*sigSig2[Q2indx]*sigSig2[Q2indx])) continue;
    
    if (Q2indx!=4 && Q2indx!=6 && PsiRejection3(M,Q,QE) == false) continue;
    if (Q2indx==4 && PsiRejection1(M,Q,QE) == false) continue;
    if (Q2indx==6 && PsiRejection2(M,Q,QE) == false) continue;

    if (L<0) L=-1*L;
    if (P<0) P=-1*P;

    fillHistos(K,L,P,M,Q2indx);
    if (Q2indx<4) fillHistos(K,L,P,M,9);
    if (Q2indx>6) fillHistos(K,L,P,M,10); //cout<<endl;
  }
  
   for (int i=0; i<11; i++) {
    can[i] = new TCanvas(Form("can%i",i),Form("can%i",i));
    can[i]->Divide(2,2);
    can[i]->cd(1);
    plotK_L0[i]->SetLineColor(1);
    plotK_L1[i]->SetLineColor(2);
    plotK_P0[i]->SetLineColor(3);
    plotK_P1[i]->SetLineColor(4);
    plotK_M0[i]->SetLineColor(6);
    plotK_M1[i]->SetLineColor(7);
    plotK_L0[i]->Scale(1./plotK_L0[i]->Integral());
    plotK_L1[i]->Scale(1./plotK_L1[i]->Integral());
    plotK_P0[i]->Scale(1./plotK_P0[i]->Integral());
    plotK_P1[i]->Scale(1./plotK_P1[i]->Integral());
    plotK_M0[i]->Scale(1./plotK_M0[i]->Integral());
    plotK_M1[i]->Scale(1./plotK_M1[i]->Integral());
    plotK_L0[i]->SetMinimum(0);
    plotK_L0[i]->SetMaximum(0.6);
    plotK_L0[i]->Draw("L");
    plotK_L1[i]->Draw("sameL");
    plotK_P0[i]->Draw("sameL");
    plotK_P1[i]->Draw("sameL");
    plotK_M0[i]->Draw("sameL");
    plotK_M1[i]->Draw("sameL");


    can[i]->cd(2);
    plotL_K0[i]->SetLineColor(30);
    plotL_K1[i]->SetLineColor(40);
    plotL_P0[i]->SetLineColor(3);
    plotL_P1[i]->SetLineColor(4);
    plotL_M0[i]->SetLineColor(6);
    plotL_M1[i]->SetLineColor(7);
    plotL_K0[i]->Scale(1./plotL_K0[i]->Integral());
    plotL_K1[i]->Scale(1./plotL_K1[i]->Integral());
    plotL_P0[i]->Scale(1./plotL_P0[i]->Integral());
    plotL_P1[i]->Scale(1./plotL_P1[i]->Integral());
    plotL_M0[i]->Scale(1./plotL_M0[i]->Integral());
    plotL_M1[i]->Scale(1./plotL_M1[i]->Integral());
    plotL_K0[i]->SetMinimum(0);
    plotL_K0[i]->SetMaximum(0.6);
    plotL_K0[i]->Draw("L");
    plotL_K1[i]->Draw("sameL");
    plotL_P0[i]->Draw("sameL");
    plotL_P1[i]->Draw("sameL");
    plotL_M0[i]->Draw("sameL");
    plotL_M1[i]->Draw("sameL");

    can[i]->cd(3);
    plotP_L0[i]->SetLineColor(1);
    plotP_L1[i]->SetLineColor(2);
    plotP_K0[i]->SetLineColor(30);
    plotP_K1[i]->SetLineColor(40);
    plotP_M0[i]->SetLineColor(6);
    plotP_M1[i]->SetLineColor(7);
    plotP_L0[i]->Scale(1./plotP_L0[i]->Integral());
    plotP_L1[i]->Scale(1./plotP_L1[i]->Integral());
    plotP_K0[i]->Scale(1./plotP_K0[i]->Integral());
    plotP_K1[i]->Scale(1./plotP_K1[i]->Integral());
    plotP_M0[i]->Scale(1./plotP_M0[i]->Integral());
    plotP_M1[i]->Scale(1./plotP_M1[i]->Integral());
    plotP_L0[i]->SetMinimum(0);
    plotP_L0[i]->SetMaximum(0.6);
    plotP_L0[i]->Draw("L");
    plotP_L1[i]->Draw("sameL");
    plotP_K0[i]->Draw("sameL");
    plotP_K1[i]->Draw("sameL");
    plotP_M0[i]->Draw("sameL");
    plotP_M1[i]->Draw("sameL");

    can[i]->cd(4);
    plotM_L0[i]->SetLineColor(1);
    plotM_L1[i]->SetLineColor(2);
    plotM_P0[i]->SetLineColor(3);
    plotM_P1[i]->SetLineColor(4);
    plotM_K0[i]->SetLineColor(30);
    plotM_K1[i]->SetLineColor(40);
    plotM_L0[i]->Scale(1./plotM_L0[i]->Integral());
    plotM_L1[i]->Scale(1./plotM_L1[i]->Integral());
    plotM_P0[i]->Scale(1./plotM_P0[i]->Integral());
    plotM_P1[i]->Scale(1./plotM_P1[i]->Integral());
    plotM_K0[i]->Scale(1./plotM_K0[i]->Integral());
    plotM_K1[i]->Scale(1./plotM_K1[i]->Integral());
    plotM_L0[i]->SetMaximum(0.5);
    plotM_L0[i]->Draw("L");
    plotM_L1[i]->Draw("sameL");
    plotM_P0[i]->Draw("sameL");
    plotM_P1[i]->Draw("sameL");
    plotM_K0[i]->Draw("sameL");
    plotM_K1[i]->Draw("sameL");

    can[i]->SaveAs(Form("sidebands%i.pdf",i));
  }
}



