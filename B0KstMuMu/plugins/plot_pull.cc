#include <iostream>
#include <string>

#include "TH1.h"
#include "TLine.h"
#include "TCanvas.h"

#include "TFile.h"

using namespace std;

TH1F* hcent_P1[9];
TH1F* hlow_P1[9];
TH1F* hhig_P1[9];
TH1F* herhi_P1[9];
TH1F* herlo_P1[9];
TH1F* hpull_P1[9];
TH1F* hcent_P5[9];
TH1F* hlow_P5[9];
TH1F* hhig_P5[9];
TH1F* herhi_P5[9];
TH1F* herlo_P5[9];
TH1F* hpull_P5[9];
TH1F* hs;
TH1F* hdist_P1[9];
TH1F* hdist_P5[9];
TH1F* hdist[9];
TH1F* hlastval[9];

TGraphAsymmErrors* graph_P1[9];
TGraphAsymmErrors* graph_P5[9];

TLine* linehi_P1;
TLine* linelo_P1;
TLine* linehi_P5;
TLine* linelo_P5;
TLine* line_P1;
TLine* line_P5;

TCanvas* c1c;
TCanvas* c1e;
TCanvas* c1p;
TCanvas* c5c;
TCanvas* c5e;
TCanvas* c5p;

TCanvas* cs;

TFile* fout;

double MINhi_P1[9] = {0.457,0.582,0.342,0.264,0,0.159,0,0.163,0.222};
double MINlo_P1[9] = {0.439,0.422,0.374,0.245,0,0.440,0,0.709,0.205};
double MINhi_P5[9] = {0.318,0.369,0.179,0.214,0,0.072,0,-0.20,0.186};
double MINlo_P5[9] = {0.326,0.300,0.392,0.190,0,0.277,0,-0.20,0.173};

float P1BF [9] = {0.150675,-0.668496, 0.498243,-0.476409,0,-0.453944,0,-0.368363, -0.547864};
float P5BF [9] = {0.105585,-0.562285,-0.953227,-0.643975,0, -0.73848,0,-0.646982, -0.549104};

int meth = 3;

void read (int q2BinIndx, int genIndx, bool plot=true)
{
  if (q2BinIndx==4 || q2BinIndx==6) return;

  if (plot) {
    hcent_P1[q2BinIndx] = new TH1F(Form("hcent_P1%i",q2BinIndx),Form("hcent_P1%i",q2BinIndx),24,-1,1);
    hcent_P5[q2BinIndx] = new TH1F(Form("hcent_P5%i",q2BinIndx),Form("hcent_P5%i",q2BinIndx),24,-1.5,1.5);
    hhig_P1[q2BinIndx] = new TH1F(Form("hhig_P1%i",q2BinIndx),Form("hhig_P1%i",q2BinIndx),24,-1,1);
    hlow_P1[q2BinIndx] = new TH1F(Form("hlow_P1%i",q2BinIndx),Form("hlow_P1%i",q2BinIndx),24,-1,1);
    hhig_P5[q2BinIndx] = new TH1F(Form("hhig_P5%i",q2BinIndx),Form("hhig_P5%i",q2BinIndx),24,-1.5,1.5);
    hlow_P5[q2BinIndx] = new TH1F(Form("hlow_P5%i",q2BinIndx),Form("hlow_P5%i",q2BinIndx),24,-1.5,1.5);
    herhi_P1[q2BinIndx] = new TH1F(Form("herhi_P1%i",q2BinIndx),Form("herhi_P1%i",q2BinIndx),24,-0.2,1);
    herhi_P5[q2BinIndx] = new TH1F(Form("herhi_P5%i",q2BinIndx),Form("herhi_P5%i",q2BinIndx),24,-0.2,1);
    herlo_P1[q2BinIndx] = new TH1F(Form("herlo_P1%i",q2BinIndx),Form("herlo_P1%i",q2BinIndx),24,-0.2,1);
    herlo_P5[q2BinIndx] = new TH1F(Form("herlo_P5%i",q2BinIndx),Form("herlo_P5%i",q2BinIndx),24,-0.2,1);
    hlastval[q2BinIndx] = new TH1F(Form("h%ilastval%i",meth,q2BinIndx),Form("h%ilastval%i",meth,q2BinIndx),40,0,2);

    hservice = new TH1F("hs","hs",20,0,20);
    
    graph_P1[q2BinIndx] = new TGraphAsymmErrors();
    graph_P5[q2BinIndx] = new TGraphAsymmErrors();
    graph_P1[q2BinIndx]->SetName(Form("graph%i_P1_m%i",q2BinIndx,meth));
    graph_P5[q2BinIndx]->SetName(Form("graph%i_P5_m%i",q2BinIndx,meth));
  }

  hpull_P1[q2BinIndx] = new TH1F(Form("hpull_P1%i",q2BinIndx),Form("hpull_P1%i",q2BinIndx),24,-4,4);
  hpull_P5[q2BinIndx] = new TH1F(Form("hpull_P5%i",q2BinIndx),Form("hpull_P5%i",q2BinIndx),24,-4,4);
  hdist_P1[q2BinIndx] = new TH1F(Form("hdist_P1%i",q2BinIndx),Form("hdist_P1%i",q2BinIndx),200,0,1);
  hdist_P5[q2BinIndx] = new TH1F(Form("hdist_P5%i",q2BinIndx),Form("hdist_P5%i",q2BinIndx),200,0,1);
  hdist[q2BinIndx] = new TH1F(Form("hdist%i",q2BinIndx),Form("hdist%i",q2BinIndx),200,0,1);

  float Chi2,NDf,Edm,NCalls,p0,p1,p2,p3,p4,p5,freeP1,freeP5,lastval,area;
  float cent_P1=0;
  float cent_P5=0;
  float erhi_P1=0;
  float erhi_P5=0;
  float erlo_P1=0;
  float erlo_P5=0;
  float pull_P1=0;
  float pull_P5=0;

  int err1_cnt=0;
  int err2_cnt=0;
  int err3_cnt=0;

  vector<float> vcent_P1 (0);
  vector<float> verlo_P1 (0);
  vector<float> verhi_P1 (0);
  vector<float> vcent_P5 (0);
  vector<float> verlo_P5 (0);
  vector<float> verhi_P5 (0);

  double aP5pLim [198];
  double aP1Lim [198];
  fstream finl (Form("waveP_lim%i.list",q2BinIndx),fstream::in);
  finl.ignore(100,'\n');
  finl.ignore(100,'\n');
  for (int indx=0; indx<198; indx++) {
    finl>>aP5pLim[indx]>>aP1Lim[indx];
  }
  finl.close();

  if (genIndx>0) {
    fstream fingen (Form("GENpoints%i.list",q2BinIndx),fstream::in);
    int nEntry=0;
    do {
      fingen>>nEntry>>P1BF[q2BinIndx]>>P5BF[q2BinIndx];
    } while (nEntry<genIndx);
  }

  float ylow = -5;
  float xlow = -5;
  float xhigh = 1;
  for (int indx=1; indx<198; indx++) {
    if (ylow==-5 && aP1Lim[indx]>P1BF[q2BinIndx]) ylow = aP5pLim[indx-1];
    if (xlow==-5 && fabs(aP5pLim[indx])>fabs(P5BF[q2BinIndx])) xlow = aP1Lim[indx];
    if (xlow!=-5 && xhigh==1 && fabs(aP5pLim[indx])<fabs(P5BF[q2BinIndx])) xhigh = aP1Lim[indx-1];
  }
  // cout<<ylow<<" "<<xlow<<" "<<xhigh<<endl;

  
  TF2* fgaus = new TF2("f2","exp([0] - (x-[1])*(x-[1])/2/[3]/[3] - (y-[2])*(y-[2])/2/[4]/[4] - (x-[1])/[3]*(y-[2])/[4]*[5])",-1.,1.,-1.5,1.5);

  for (int toyIndx=0; toyIndx<100; toyIndx++) { //if (toyIndx==42) continue;
    // cout<<toyIndx<<endl;
    TString filename = Form("Data_pull/%i/pull_m%i_result%i_v4.log",toyIndx,meth,q2BinIndx+1); 
    if (genIndx>0) filename = Form("Data_pull%i/%i/pull_m%i_result%i_v4.log",genIndx,toyIndx,meth,q2BinIndx+1); 
    if ( gSystem->AccessPathName(filename) ) continue;
    fstream fin (filename,fstream::in);

    // char s[100];
    string s;
    bool out=false;
    do {
      fin>>s;
      if (s=="No") {
	err1_cnt++;
	out=true;
	break;
      }
      if (s=="Invalid") {
	err2_cnt++;
	out=true;
	break;
      }
      // if (s=='EOF') {
      //   out=true;
      //   break;
      // }
    }
    while (s!="=");
    if (out) continue;
    fin>>Chi2;
    do fin>>s;
    while (s!="=");
    fin>>NDf;
    {
      if (NDf<3) {
	err2_cnt++;
        continue;
      }
    }
    do fin>>s;
    while (s!="=");
    fin>>Edm;
    do fin>>s;
    while (s!="=");
    fin>>NCalls;
    do fin>>s;
    while (s!="=");
    fin>>p0;
    do fin>>s;
    while (s!="=");
    fin>>p1;
    do fin>>s;
    while (s!="=");
    fin>>p2;
    do fin>>s;
    while (s!="=");
    fin>>p3;
    do fin>>s;
    while (s!="=");
    fin>>p4;
    do fin>>s;
    while (s!="=");
    fin>>p5;
    fin.ignore(100,'\n');
    fin>>s;
    if (s=="Limit") {
      err3_cnt++;
      continue;
    }
    freeP1 = atof(s.c_str());
    fin>>freeP5;
    fin>>cent_P1;
    fin>>erhi_P1;
    fin>>erlo_P1;
    // if (genIndx==0) fin>>pull_P1;
    fin>>cent_P5;
    fin>>erhi_P5;
    fin>>erlo_P5;
    // if (genIndx==0) fin>>pull_P5;

    fin>>lastval>>area;

    fin.close();

    // if (erhi_P1==-1) {
    //   // cout<<toyIndx<<endl;
    //   continue;
    // }
    // float err_P1, err_P5;
    // if (cent_P1-P1BF[q2BinIndx]>0) err_P1 = erhi_P1;
    // else err_P1 = erlo_P1;
    // if (cent_P5-P5BF[q2BinIndx]>0) err_P5 = erhi_P5;
    // else err_P5 = erlo_P5;

    if (plot) {
      hcent_P1[q2BinIndx]->Fill(cent_P1);
      hhig_P1[q2BinIndx]->Fill(erhi_P1);
      hlow_P1[q2BinIndx]->Fill(erlo_P1);
      herhi_P1[q2BinIndx]->Fill(erhi_P1-cent_P1);
      herlo_P1[q2BinIndx]->Fill(cent_P1-erlo_P1);
      // hpull_P1[q2BinIndx]->Fill(pull_P1);
      hcent_P5[q2BinIndx]->Fill(cent_P5);
      hhig_P5[q2BinIndx]->Fill(erhi_P5);
      hlow_P5[q2BinIndx]->Fill(erlo_P5);
      herhi_P5[q2BinIndx]->Fill(erhi_P5-cent_P5);
      herlo_P5[q2BinIndx]->Fill(cent_P5-erlo_P5);
      // hpull_P5[q2BinIndx]->Fill(pull_P5);
    }

    vcent_P1.push_back(cent_P1);
    verhi_P1.push_back(erhi_P1);
    verlo_P1.push_back(erlo_P1);
    vcent_P5.push_back(cent_P5);
    verhi_P5.push_back(erhi_P5);
    verlo_P5.push_back(erlo_P5);

    if (plot) {
      float yshift = 0;
      if (meth==2) yshift = -0.3;
      if (meth==3) yshift = 0.3;
      graph_P1[q2BinIndx]->SetPoint(vcent_P1.size()-1,cent_P1,toyIndx+yshift);
      graph_P5[q2BinIndx]->SetPoint(vcent_P1.size()-1,cent_P5,toyIndx+yshift);
      graph_P1[q2BinIndx]->SetPointError(vcent_P1.size()-1,cent_P1-erlo_P1,erhi_P1-cent_P1,0,0);
      graph_P5[q2BinIndx]->SetPointError(vcent_P1.size()-1,cent_P5-erlo_P5,erhi_P5-cent_P5,0,0);
      // if(erhi_P1<=0 || erlo_P1<=0) cout<<toyIndx<<" (P1)"<<endl;
      // if(erhi_P5<=0 || erlo_P5<=0) cout<<toyIndx<<" (P5p)"<<endl;

      hservice->Fill(NDf);
    }
    
    fgaus->SetParameters(0,p1,p2,p3,p4,p5);
    float ymas = 0;
    float xmas = 0;
    for (float xpoint = xlow; xpoint <= xhigh; xpoint+=0.001) {
      float val = fgaus->Eval(xpoint,P5BF[q2BinIndx]);
      if (val > xmas) xmas = val;
    }
    for(float ypoint =ylow; ypoint <= fabs(ylow); ypoint+=0.001) {
      float val= fgaus->Eval(P1BF[q2BinIndx],ypoint);
      if (val > ymas) ymas = val;
    }
    float peak = fgaus->Eval(cent_P1,cent_P5);
    hdist_P1[q2BinIndx]->Fill(ymas/peak);
    hdist_P5[q2BinIndx]->Fill(xmas/peak);
    hdist[q2BinIndx]->Fill(fgaus->Eval(P1BF[q2BinIndx],P5BF[q2BinIndx])/peak);

    if (plot) hlastval[q2BinIndx]->Fill(-1*log(lastval/peak));

    // if (P5BF[q2BinIndx] < erlo_P5) cout<<toyIndx<<" "<<erlo_P5<<" "<<cent_P5<<" "<<erhi_P5<<endl; 
    // if (erhi_P5-cent_P5<0 || cent_P5-erlo_P5<0) cout<<toyIndx<<" "<<erlo_P5<<" "<<cent_P5<<" "<<erhi_P5<<" "<<-1*log(lastval/peak)<<endl;
    // if (erhi_P1-cent_P1<0 || cent_P1-erlo_P1<0) cout<<toyIndx<<" "<<erlo_P1<<" "<<cent_P1<<" "<<erhi_P1<<" "<<-1*log(lastval/peak)<<endl;
    // if ((fabs(pull_P1)>1 && ymas>exp(-0.5)) || (fabs(pull_P1)<1 && ymas<exp(-0.5))) cout<<toyIndx<<" "<<pull_P1<<" "<<ymas<<" "<<ylow<<endl;
  }

  // float refval_P1 = hcent_P1[q2BinIndx]->GetMean();
  // float refval_P5 = hcent_P5[q2BinIndx]->GetMean();
  int coverP1=0;
  int coverP5=0;
  float refval_P1 = P1BF[q2BinIndx];
  float refval_P5 = P5BF[q2BinIndx];
  float err_P1, err_P5; 
  for (size_t i=0; i<vcent_P1.size(); i++) {
    if (vcent_P1[i]-refval_P1<0) err_P1 = verhi_P1[i]-vcent_P1[i];
    else err_P1 = vcent_P1[i]-verlo_P1[i];
    if (vcent_P5[i]-refval_P5<0) err_P5 = verhi_P5[i]-vcent_P5[i];
    else err_P5 = vcent_P5[i]-verlo_P5[i];
    hpull_P1[q2BinIndx]->Fill((vcent_P1[i]-refval_P1)/err_P1);
    hpull_P5[q2BinIndx]->Fill((vcent_P5[i]-refval_P5)/err_P5);
    if ( refval_P1<verhi_P1[i] && refval_P1>verlo_P1[i] )  coverP1++;
    // cout<<verlo_P5[i]<<refval_P5<<verhi_P5[i]<<endl;
    if ( refval_P5<verhi_P5[i] && refval_P5>verlo_P5[i] )  coverP5++;
  }
  
  // double lim=-1;
  // double limup=-1;
  // double limdo=-1;
  // for (int j=hdist_P1[q2BinIndx]->GetNbinsX()-1; j>0; j--) {
  //   if ( lim==-1 && hdist_P1[q2BinIndx]->Integral(j,hdist_P1[q2BinIndx]->GetNbinsX()+1) > 0.6827*hdist_P1[q2BinIndx]->Integral(0,hdist_P1[q2BinIndx]->GetNbinsX()+1) ) lim = hdist_P1[q2BinIndx]->GetBinCenter(j);
  //   if ( limup==-1 && hdist_P1[q2BinIndx]->Integral(j,hdist_P1[q2BinIndx]->GetNbinsX()+1) > 0.6362*hdist_P1[q2BinIndx]->Integral(0,hdist_P1[q2BinIndx]->GetNbinsX()+1) ) limup = hdist_P1[q2BinIndx]->GetBinCenter(j);
  //   if ( limdo==-1 && hdist_P1[q2BinIndx]->Integral(j,hdist_P1[q2BinIndx]->GetNbinsX()+1) > 0.7292*hdist_P1[q2BinIndx]->Integral(0,hdist_P1[q2BinIndx]->GetNbinsX()+1) ) limdo = hdist_P1[q2BinIndx]->GetBinCenter(j);
  // }
  // cout<<q2BinIndx<<" like lim P1: "<<lim<<"["<<limup<<","<<limdo<<"] ("<<-1*log(lim)<<"["<<-1*log(limup)<<","<<-1*log(limdo)<<"])"<<endl;
  // lim=-1;
  // limup=-1;
  // limdo=-1;
  // for (int j=hdist_P5[q2BinIndx]->GetNbinsX()-1; j>0; j--) {
  //   if ( lim==-1 && hdist_P5[q2BinIndx]->Integral(j,hdist_P5[q2BinIndx]->GetNbinsX()+1) > 0.6827*hdist_P5[q2BinIndx]->Integral(0,hdist_P5[q2BinIndx]->GetNbinsX()+1) ) lim = hdist_P5[q2BinIndx]->GetBinCenter(j);
  //   if ( limup==-1 && hdist_P5[q2BinIndx]->Integral(j,hdist_P5[q2BinIndx]->GetNbinsX()+1) > 0.6362*hdist_P5[q2BinIndx]->Integral(0,hdist_P5[q2BinIndx]->GetNbinsX()+1) ) limup = hdist_P5[q2BinIndx]->GetBinCenter(j);
  //   if ( limdo==-1 && hdist_P5[q2BinIndx]->Integral(j,hdist_P5[q2BinIndx]->GetNbinsX()+1) > 0.7292*hdist_P5[q2BinIndx]->Integral(0,hdist_P5[q2BinIndx]->GetNbinsX()+1) ) limdo = hdist_P5[q2BinIndx]->GetBinCenter(j);
  // }
  // cout<<q2BinIndx<<" like lim P5: "<<lim<<"["<<limup<<","<<limdo<<"] ("<<-1*log(lim)<<"["<<-1*log(limup)<<","<<-1*log(limdo)<<"])"<<endl;
  // cout<<hdist_P1[q2BinIndx]->Integral(hdist_P1[q2BinIndx]->FindBin(exp(-0.5)),hdist_P1[q2BinIndx]->GetNbinsX()+1)<<" "<<hdist_P5[q2BinIndx]->Integral(hdist_P5[q2BinIndx]->FindBin(exp(-0.5)),hdist_P5[q2BinIndx]->GetNbinsX()+1)<<endl;
  // lim=-1;
  // limup=-1;
  // limdo=-1;
  // for (int j=hdist[q2BinIndx]->GetNbinsX()-1; j>0; j--) {
  //   if ( lim==-1 && hdist[q2BinIndx]->Integral(j,hdist[q2BinIndx]->GetNbinsX()+1) > 0.6827*hdist[q2BinIndx]->Integral(0,hdist[q2BinIndx]->GetNbinsX()+1) ) lim = hdist[q2BinIndx]->GetBinCenter(j);
  //   if ( limup==-1 && hdist[q2BinIndx]->Integral(j,hdist[q2BinIndx]->GetNbinsX()+1) > 0.6362*hdist[q2BinIndx]->Integral(0,hdist[q2BinIndx]->GetNbinsX()+1) ) limup = hdist[q2BinIndx]->GetBinCenter(j);
  //   if ( limdo==-1 && hdist[q2BinIndx]->Integral(j,hdist[q2BinIndx]->GetNbinsX()+1) > 0.7292*hdist[q2BinIndx]->Integral(0,hdist[q2BinIndx]->GetNbinsX()+1) ) limdo = hdist[q2BinIndx]->GetBinCenter(j);
  // }
  // cout<<q2BinIndx<<" like lim 2D: "<<lim<<"["<<limup<<","<<limdo<<"] ("<<-1*log(lim)<<"["<<-1*log(limup)<<","<<-1*log(limdo)<<"])"<<endl;
  // lim=-1;
  // limup=-1;
  // limdo=-1;
  // for (int j=hdist[q2BinIndx]->GetNbinsX()-1; j>0; j--) {
  //   if ( lim==-1 && hdist[q2BinIndx]->Integral(j,hdist[q2BinIndx]->GetNbinsX()+1) > 0.39347*hdist[q2BinIndx]->Integral(0,hdist[q2BinIndx]->GetNbinsX()+1) ) lim = hdist[q2BinIndx]->GetBinCenter(j);
  //   if ( limup==-1 && hdist[q2BinIndx]->Integral(j,hdist[q2BinIndx]->GetNbinsX()+1) > 0.34462*hdist[q2BinIndx]->Integral(0,hdist[q2BinIndx]->GetNbinsX()+1) ) limup = hdist[q2BinIndx]->GetBinCenter(j);
  //   if ( limdo==-1 && hdist[q2BinIndx]->Integral(j,hdist[q2BinIndx]->GetNbinsX()+1) > 0.44232*hdist[q2BinIndx]->Integral(0,hdist[q2BinIndx]->GetNbinsX()+1) ) limdo = hdist[q2BinIndx]->GetBinCenter(j);
  // }
  // cout<<q2BinIndx<<" like lim 39: "<<lim<<"["<<limup<<","<<limdo<<"] ("<<-1*log(lim)<<"["<<-1*log(limup)<<","<<-1*log(limdo)<<"])"<<endl;

  // for (int j=hdist[q2BinIndx]->GetNbinsX(); j>0; j--) if ( hdist[q2BinIndx]->Integral(j,hdist[q2BinIndx]->GetNbinsX()+1) > 0.39347*hdist[q2BinIndx]->Integral(0,hdist[q2BinIndx]->GetNbinsX()+1) ) {
  //     cout<<q2BinIndx<<" like lim 39: "<<hdist[q2BinIndx]->GetBinCenter(j)-hdist[q2BinIndx]->GetBinWidth(j)<<" ("<<-1*log(hdist[q2BinIndx]->GetBinCenter(j)-hdist[q2BinIndx]->GetBinWidth(j))<<")"<<endl;
  //     break;
  //   }
  
  // cout<<q2BinIndx<<" P1 cov:"<<(hpull_P1[q2BinIndx]->Integral(0,25)>0?hpull_P1[q2BinIndx]->Integral(10,15)/hpull_P1[q2BinIndx]->Integral(0,25):-1)<<" "<<1.*coverP1/vcent_P1.size()<<" ("<<hpull_P1[q2BinIndx]->Integral(10,15)<<"/"<<hpull_P1[q2BinIndx]->Integral(0,25)<<")"<<endl;
  // cout<<q2BinIndx<<" P5 cov:"<<(hpull_P5[q2BinIndx]->Integral(0,25)>0?hpull_P5[q2BinIndx]->Integral(10,15)/hpull_P5[q2BinIndx]->Integral(0,25):-1)<<" "<<1.*coverP5/vcent_P1.size()<<" ("<<hpull_P5[q2BinIndx]->Integral(10,15)<<"/"<<hpull_P5[q2BinIndx]->Integral(0,25)<<")"<<endl;
  // cout<<q2BinIndx<<" Errors:"<<err1_cnt<<" "<<err2_cnt<<" "<<err3_cnt<<endl;
  if (vcent_P1.size()>0) cout<<q2BinIndx<<" "<<genIndx<<" "<<1.*coverP1/vcent_P1.size()<<" "<<1.*coverP5/vcent_P5.size()<<" "<<coverP1<<" "<<coverP5<<" "<<vcent_P1.size()<<" "<<hpull_P1[q2BinIndx]->GetRMS()<<" "<<hpull_P5[q2BinIndx]->GetRMS()<<endl;

  if (plot) {
    linehi_P1 = new TLine(MINhi_P1[q2BinIndx],0,MINhi_P1[q2BinIndx],herhi_P1[q2BinIndx]->GetMaximum());
    linelo_P1 = new TLine(MINlo_P1[q2BinIndx],0,MINlo_P1[q2BinIndx],herhi_P1[q2BinIndx]->GetMaximum());
    linehi_P5 = new TLine(MINhi_P5[q2BinIndx],0,MINhi_P5[q2BinIndx],herhi_P5[q2BinIndx]->GetMaximum());
    linelo_P5 = new TLine(MINlo_P5[q2BinIndx],0,MINlo_P5[q2BinIndx],herhi_P5[q2BinIndx]->GetMaximum());
    line_P1 = new TLine(P1BF[q2BinIndx],0,P1BF[q2BinIndx],hcent_P1[q2BinIndx]->GetMaximum());
    line_P5 = new TLine(P5BF[q2BinIndx],0,P5BF[q2BinIndx],hcent_P5[q2BinIndx]->GetMaximum());

    c1c = new TCanvas("c1c");
    hcent_P1[q2BinIndx]->Draw();
    // hhig_P1[q2BinIndx]->Draw("same");
    // hlow_P1[q2BinIndx]->Draw("same");
    line_P1->Draw("same");
    c1e = new TCanvas("c1e");
    herhi_P1[q2BinIndx]->Draw();
    herlo_P1[q2BinIndx]->Draw("same");
    linehi_P1->Draw("same");
    linelo_P1->Draw("same");
    c1p = new TCanvas("c1p");
    hpull_P1[q2BinIndx]->Draw();
    c5c = new TCanvas("c5c");
    hcent_P5[q2BinIndx]->Draw();
    // hhig_P5[q2BinIndx]->Draw("same");
    // hlow_P5[q2BinIndx]->Draw("same");
    line_P5->Draw("same");
    c5e = new TCanvas("c5e"); 
    herhi_P5[q2BinIndx]->Draw();
    herlo_P5[q2BinIndx]->Draw("same");
    linehi_P5->Draw("same");
    linelo_P5->Draw("same");
    c5p = new TCanvas("c5p");
    hpull_P5[q2BinIndx]->Draw();
    
    herlo_P1[q2BinIndx]->SetLineColor(2);
    herlo_P5[q2BinIndx]->SetLineColor(2);
    hhig_P1[q2BinIndx]->SetLineColor(2);
    hlow_P1[q2BinIndx]->SetLineColor(3);
    hhig_P5[q2BinIndx]->SetLineColor(2);
    hlow_P5[q2BinIndx]->SetLineColor(3);
    hcent_P1[q2BinIndx]->SetLineWidth(2);
    hcent_P5[q2BinIndx]->SetLineWidth(2);
    linelo_P1->SetLineColor(2);
    linelo_P5->SetLineColor(2);
    linehi_P1->SetLineColor(4);
    linehi_P5->SetLineColor(4);
    line_P1->SetLineColor(2);
    line_P5->SetLineColor(2);
    linelo_P1->SetLineWidth(2);
    linelo_P5->SetLineWidth(2);
    linehi_P1->SetLineWidth(2);
    linehi_P5->SetLineWidth(2);
    line_P1->SetLineWidth(2);
    line_P5->SetLineWidth(2);

    string dirname = "Data_pull/";
    if (genIndx>0) dirname = Form("Data_pull%i/",genIndx);
    c1c->SaveAs((dirname+Form("P1cent%i_m%i_%i_v4.pdf",q2BinIndx,meth,genIndx)).c_str());
    c5c->SaveAs((dirname+Form("P5cent%i_m%i_%i_v4.pdf",q2BinIndx,meth,genIndx)).c_str());
    c1e->SaveAs((dirname+Form("P1err%i_m%i_%i_v4.pdf",q2BinIndx,meth,genIndx)).c_str());
    c5e->SaveAs((dirname+Form("P5err%i_m%i_%i_v4.pdf",q2BinIndx,meth,genIndx)).c_str());
    c1p->SaveAs((dirname+Form("P1pull%i_m%i_%i_v4.pdf",q2BinIndx,meth,genIndx)).c_str());
    c5p->SaveAs((dirname+Form("P5pull%i_m%i_%i_v4.pdf",q2BinIndx,meth,genIndx)).c_str());
    c1c->SaveAs((dirname+Form("P1cent%i_m%i_%i_v4.root",q2BinIndx,meth,genIndx)).c_str());
    c5c->SaveAs((dirname+Form("P5cent%i_m%i_%i_v4.root",q2BinIndx,meth,genIndx)).c_str());
    c1e->SaveAs((dirname+Form("P1err%i_m%i_%i_v4.root",q2BinIndx,meth,genIndx)).c_str());
    c5e->SaveAs((dirname+Form("P5err%i_m%i_%i_v4.root",q2BinIndx,meth,genIndx)).c_str());
    c1p->SaveAs((dirname+Form("P1pull%i_m%i_%i_v4.root",q2BinIndx,meth,genIndx)).c_str());
    c5p->SaveAs((dirname+Form("P5pull%i_m%i_%i_v4.root",q2BinIndx,meth,genIndx)).c_str());

    cs = new TCanvas("cs");
    hlastval[q2BinIndx]->Draw();
    // hdist_P1[q2BinIndx]->Draw();
    // hdist_P5[q2BinIndx]->Draw("same");
    // hdist[q2BinIndx]->Draw("same");
    // hdist_P1[q2BinIndx]->SetLineColor(1);
    // hdist_P5[q2BinIndx]->SetLineColor(2);
    // hservice->Draw();

    // fout = new TFile("graphOut.root","UPDATE");
    // fout->cd();
    // graph_P1[q2BinIndx]->Write(0,TFile::kOverwrite);
    // graph_P5[q2BinIndx]->Write(0,TFile::kOverwrite);
    // fout->Close();

    fout = new TFile("limitOut.root","UPDATE");
    fout->cd();
    hlastval[q2BinIndx]->Write(0,TFile::kOverwrite);
    fout->Close();

  } else {
    delete hpull_P1[q2BinIndx];
    delete hpull_P5[q2BinIndx];
    delete hdist_P1[q2BinIndx];
    delete hdist_P5[q2BinIndx];
    delete hdist[q2BinIndx];
  }

  return;
}

void plot_pull (int q2BinIndx = 0, int genIndx=-1, bool plot=true)
{
  if (q2BinIndx==0) for (q2BinIndx=1;q2BinIndx<10;q2BinIndx++) {
      if (genIndx==-1) {
	for (genIndx=0;genIndx<9;genIndx++) read(q2BinIndx-1, genIndx, false);
	genIndx=-1;
      }
      else if (genIndx>-1 && genIndx<9) read(q2BinIndx-1, genIndx, false);
    }
  else if (q2BinIndx>0 && q2BinIndx<10) {
    if (genIndx==-1) {
      for (genIndx=0;genIndx<9;genIndx++) read(q2BinIndx-1, genIndx, false);
      genIndx=-1;
    }
    else if (genIndx>-1 && genIndx<9) read(q2BinIndx-1, genIndx, plot);
  }
  return;
}
