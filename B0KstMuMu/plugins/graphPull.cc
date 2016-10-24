#include <iostream>
#include <string>

#include "TH1.h"
#include "TLine.h"
#include "TCanvas.h"

#include "TFile.h"

using namespace std;

TGraphAsymmErrors*  graph1_P1[9];
TGraphAsymmErrors*  graph2_P1[9];
TGraphAsymmErrors*  graph3_P1[9];
TGraphAsymmErrors*  graph1_P5[9];
TGraphAsymmErrors*  graph2_P5[9];
TGraphAsymmErrors*  graph3_P5[9];

TH1F* width1_P1[9];
TH1F* width2_P1[9];
TH1F* width3_P1[9];
TH1F* width1_P5[9];
TH1F* width2_P5[9];
TH1F* width3_P5[9];

TLine* line_P1[9];
TLine* line_P5[9];
TLine* mean_P1[9];
TLine* mean_P5[9];

TCanvas* cs[9];
TCanvas* cw[9];

TFile* fout;

float As5BF[9] = {       1, 0.999115,-0.434302, 0.530473,0,-0.933958,0,0.0524524,-0.0542194};
float P1BF [9] = {0.150675,-0.668496, 0.498243,-0.476409,0,-0.453944,0,-0.368363, -0.547864};
float P5BF [9] = {0.105585,-0.562285,-0.953227,-0.643975,0, -0.73848,0,-0.646982, -0.549104};

float P1mean [9] = {0.150675,-0.668496, 0.498243,-0.476409,0,-0.453944,0,-0.368363, -0.547864};
float P5mean [9] = {0.105585,-0.562285,-0.953227,-0.643975,0, -0.73848,0,-0.646982, -0.549104};

void read (int q2BinIndx)
{
  if (q2BinIndx==4 || q2BinIndx==6) return;

  fout = new TFile("Data_pull/graphOut.root","UPDATE");
  fout->cd();
  graph1_P1[q2BinIndx] = (TGraphAsymmErrors*)fout->Get(Form("graph%i_P1_m1",q2BinIndx));
  graph2_P1[q2BinIndx] = (TGraphAsymmErrors*)fout->Get(Form("graph%i_P1_m2",q2BinIndx));
  graph3_P1[q2BinIndx] = (TGraphAsymmErrors*)fout->Get(Form("graph%i_P1_m3",q2BinIndx));
  graph1_P5[q2BinIndx] = (TGraphAsymmErrors*)fout->Get(Form("graph%i_P5_m1",q2BinIndx));
  graph2_P5[q2BinIndx] = (TGraphAsymmErrors*)fout->Get(Form("graph%i_P5_m2",q2BinIndx));
  graph3_P5[q2BinIndx] = (TGraphAsymmErrors*)fout->Get(Form("graph%i_P5_m3",q2BinIndx));

  width1_P1[q2BinIndx] = new TH1F(Form("width1_P1%i",q2BinIndx),Form("width1_P1%i",q2BinIndx),40,0,2);
  width2_P1[q2BinIndx] = new TH1F(Form("width2_P1%i",q2BinIndx),Form("width2_P1%i",q2BinIndx),40,0,2);
  width3_P1[q2BinIndx] = new TH1F(Form("width3_P1%i",q2BinIndx),Form("width3_P1%i",q2BinIndx),40,0,2);
  width1_P5[q2BinIndx] = new TH1F(Form("width1_P5%i",q2BinIndx),Form("width1_P5%i",q2BinIndx),40,0,2);
  width2_P5[q2BinIndx] = new TH1F(Form("width2_P5%i",q2BinIndx),Form("width2_P5%i",q2BinIndx),40,0,2);
  width3_P5[q2BinIndx] = new TH1F(Form("width3_P5%i",q2BinIndx),Form("width3_P5%i",q2BinIndx),40,0,2);

  for (int i=0; i<graph1_P1[q2BinIndx]->GetN(); i++) {
    width1_P1[q2BinIndx]->Fill(graph1_P1[q2BinIndx]->GetErrorXhigh(i)+graph1_P1[q2BinIndx]->GetErrorXlow(i));
    width2_P1[q2BinIndx]->Fill(graph2_P1[q2BinIndx]->GetErrorXhigh(i)+graph2_P1[q2BinIndx]->GetErrorXlow(i));
    width3_P1[q2BinIndx]->Fill(graph3_P1[q2BinIndx]->GetErrorXhigh(i)+graph3_P1[q2BinIndx]->GetErrorXlow(i));
    width1_P5[q2BinIndx]->Fill(graph1_P5[q2BinIndx]->GetErrorXhigh(i)+graph1_P5[q2BinIndx]->GetErrorXlow(i));
    width2_P5[q2BinIndx]->Fill(graph2_P5[q2BinIndx]->GetErrorXhigh(i)+graph2_P5[q2BinIndx]->GetErrorXlow(i));
    width3_P5[q2BinIndx]->Fill(graph3_P5[q2BinIndx]->GetErrorXhigh(i)+graph3_P5[q2BinIndx]->GetErrorXlow(i));
  }


  /*  mean_P1[q2BinIndx] = new TLine(P1mean[q2BinIndx],0,P1mean[q2BinIndx],101);
  mean_P5[q2BinIndx] = new TLine(P5mean[q2BinIndx],0,P5mean[q2BinIndx],101);
  line_P1[q2BinIndx] = new TLine(P1BF[q2BinIndx],0,P1BF[q2BinIndx],101);
  line_P5[q2BinIndx] = new TLine(P5BF[q2BinIndx],0,P5BF[q2BinIndx],101);

  cs[q2BinIndx] = new TCanvas(Form("cs%i",q2BinIndx));
  cs[q2BinIndx]->Divide(2,1);
  TVirtualPad* p1 = cs[q2BinIndx]->cd(1);
  p1->SetMargin(0,0,0.05,0);
  graph1_P1[q2BinIndx]->Draw("PA");
  graph2_P1[q2BinIndx]->Draw("P");
  graph3_P1[q2BinIndx]->Draw("P");
  mean_P1[q2BinIndx]->Draw("same");
  line_P1[q2BinIndx]->Draw("same");
  TVirtualPad* p2 = cs[q2BinIndx]->cd(2);
  p2->SetMargin(0,0,0.05,0);
  graph1_P5[q2BinIndx]->Draw("PA");
  graph2_P5[q2BinIndx]->Draw("P");
  graph3_P5[q2BinIndx]->Draw("P");
  mean_P5[q2BinIndx]->Draw("same");
  line_P5[q2BinIndx]->Draw("same");

  graph1_P1[q2BinIndx]->SetLineColor(4);
  graph1_P1[q2BinIndx]->SetMarkerStyle(6);
  graph2_P1[q2BinIndx]->SetLineColor(2);
  graph2_P1[q2BinIndx]->SetMarkerStyle(6);
  graph3_P1[q2BinIndx]->SetLineColor(3);
  graph3_P1[q2BinIndx]->SetMarkerStyle(6);
  graph1_P5[q2BinIndx]->SetLineColor(4);
  graph1_P5[q2BinIndx]->SetMarkerStyle(6);
  graph2_P5[q2BinIndx]->SetLineColor(2);
  graph2_P5[q2BinIndx]->SetMarkerStyle(6);
  graph3_P5[q2BinIndx]->SetLineColor(3);
  graph3_P5[q2BinIndx]->SetMarkerStyle(6);

  mean_P1[q2BinIndx]->SetLineStyle(2);
  mean_P5[q2BinIndx]->SetLineStyle(2);

  graph1_P1[q2BinIndx]->GetYaxis()->SetLabelSize(0);
  graph1_P1[q2BinIndx]->GetYaxis()->SetTickLength(0);
  graph1_P5[q2BinIndx]->GetYaxis()->SetLabelSize(0);
  graph1_P5[q2BinIndx]->GetYaxis()->SetTickLength(0);
  graph1_P1[q2BinIndx]->GetXaxis()->SetLimits(-1.01,1.01);
  graph1_P1[q2BinIndx]->GetYaxis()->SetRangeUser(-1,100);
  graph1_P5[q2BinIndx]->GetYaxis()->SetRangeUser(-1,100);
*/
  cw[q2BinIndx] = new TCanvas(Form("cw%i",q2BinIndx));
  cw[q2BinIndx]->Divide(2,1);
  cw[q2BinIndx]->cd(1);
  width3_P1[q2BinIndx]->Draw();
  width2_P1[q2BinIndx]->Draw("same");
  width1_P1[q2BinIndx]->Draw("same");
  cw[q2BinIndx]->cd(2);
  width3_P5[q2BinIndx]->Draw();
  width2_P5[q2BinIndx]->Draw("same");
  width1_P5[q2BinIndx]->Draw("same");

  width2_P1[q2BinIndx]->SetLineColor(2);
  width3_P1[q2BinIndx]->SetLineColor(3);
  width2_P5[q2BinIndx]->SetLineColor(2);
  width3_P5[q2BinIndx]->SetLineColor(3);

  // cs[q2BinIndx]->SaveAs(Form("Data_pull/result_graph%i.pdf",q2BinIndx));
  // cs[q2BinIndx]->SaveAs(Form("Data_pull/result_graph%i.root",q2BinIndx));

  // delete fout;
  
  return;
}

void graphPull (int q2BinIndx = 0)
{
  if (q2BinIndx==0) for (q2BinIndx=1;q2BinIndx<10;q2BinIndx++) read(q2BinIndx-1);
  else if (q2BinIndx>0 && q2BinIndx<10) read(q2BinIndx-1);

  return;
}
