#include <iostream>
#include <list>
#include <utility>

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

TGraph* g2d[9];

TGraph* g2d1s [9];
TGraph* g2d2s [9];
TGraph* g2d3s [9];

TGraph* gBest[9];
TGraph* gLim[9];
TGraph* gLim0[9];
TGraph* gLim1[9];
TGraph* gLim2[9];
TGraph* gLim3[9];
TGraph* gLim4[9];
TGraph* gLim5[9];
TGraph* gLim6[9];

TGraph2D* grlike[9];
TF2* fgaus[9];

TCanvas* can[9];

TH2F* hlike[9];
TH2S* hboun[9];

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

// vector<int> good (8100,0);

TF1* lim[9];
double fl[9] = {0.641004,0.799186,0.619384,0.503676,0,0.392124,0,0.476826,0.377081};

int nLim[9] = {199,198,199,199,0,199,0,199,199};

int counter = 0;

float As5BF[9] = {       1, 0.999115,-0.434302, 0.530473,0,-0.933958,0,0.0524524,-0.0542194};
float P1BF [9] = {0.150675,-0.668496, 0.498243,-0.476409,0,-0.453944,0,-0.368363, -0.547864};
float P5BF [9] = {0.105585,-0.562285,-0.953227,-0.643975,0, -0.73848,0,-0.646982, -0.549104};

float level = 0;

float minP1=0;
float minP5=0;

void open (int q2BinIndx, int toyIndx, int scanIndx, int genIndx=0, bool print=false)
{
  TString filename = Form("Data_pull/%i/sc_%i/Fitresult5_2.root",toyIndx,scanIndx);
  if ( genIndx>0 ) filename = Form("Data_pull%i/%i/sc_%i/Fitresult5_2.root",genIndx,toyIndx,scanIndx);
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
	  // good[scanIndx-1] = 1;

	  // cout<<q2BinIndx<<" "<<scanIndx<<" "<<fr->minNll()<<endl;

	  if (level==0) level=fr->minNll();

	  grlike[q2BinIndx]->SetPoint(counter,((RooRealVar*)fr->constPars().at(3))->getVal(),((RooRealVar*)fr->constPars().at(4))->getVal(),exp(level-fr->minNll()));
	  // grlike[q2BinIndx]->SetPointError(counter,0,0,0.01);
	  // cout<<counter<<" "<<((RooRealVar*)fr->constPars().at(3))->getVal()<<" "<<((RooRealVar*)fr->constPars().at(4))->getVal()<<" "<<fr->minNll()<<endl;
	  counter++;

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
	      minP1 = ((RooRealVar*)fr->constPars().at(3))->getVal();
	      minP5 = ((RooRealVar*)fr->constPars().at(4))->getVal();
	      // bestP1[0] = ((RooRealVar*)fr->constPars().at(3))->getVal();
	      // bestP5p[0] = ((RooRealVar*)fr->constPars().at(4))->getVal();
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

void read (int q2BinIndx, int toyIndx, int genIndx)
{
  if (q2BinIndx==4 || q2BinIndx==6) return;

  // cout<<q2BinIndx<<" "<<toyIndx<<" "<<genIndx<<endl;
  counter=0;

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

  level=0;

  grlike[q2BinIndx] = new TGraph2D ();
  grlike[q2BinIndx]->SetName(Form("grlike%i",q2BinIndx));
  grlike[q2BinIndx]->SetTitle("");

  minNll = 9999999;
  for (int cnt=1; cnt<42; cnt++) {
    // good[cnt-1] = 0;

    // cout<<cnt<<endl;
    open(q2BinIndx, toyIndx, cnt, genIndx);
  }
  if ( minNll == 9999999 ) {
    cout<<"No good scan results"<<endl;
    return;
  }

  // cout<<endl<<"Best "<<q2BinIndx<<": "<<best<<" ("<<minNll<<"), "<<vLik.size()<<" ("<<counter<<") good fits"<<endl;

  // open(q2BinIndx, toyIndx, best, genIndx, true);

  // double* aP5pLim = new double[nLim[q2BinIndx]];
  // double* aP1Lim = new double[nLim[q2BinIndx]];
  vector<double> aP5pLim (0);
  vector<double> aP1Lim (0);
  fstream fin (Form("waveP_lim%i.list",q2BinIndx),fstream::in);
  fin.ignore(100,'\n');
  fin.ignore(100,'\n');
  double uno, due;
  fin>>uno>>due;
  do {
    aP5pLim.push_back(uno);
    aP1Lim.push_back(due); 
    fin>>uno>>due;
  } while (due>aP1Lim[aP1Lim.size()-1]);
  // for (int i=0; i<nLim[q2BinIndx]; i++) {
  //   fin>>aP5pLim[i]>>aP1Lim[i];
  //   // if ( aP5pLim[i]==0 && aP1Lim[i]==0 )  nLim[q2BinIndx]=i;
  // }
  fin.close();

  double* aP5pLim0 = new double[nLim[q2BinIndx]];
  double* aP1Lim0 = new double[nLim[q2BinIndx]];
  fstream fin0 (Form("pdf_lim%i_%i.list",q2BinIndx,-3),fstream::in);
  fin0.ignore(100,'\n');
  fin0.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) fin0>>aP5pLim0[i]>>aP1Lim0[i];
  fin0.close();

  double* aP5pLim1 = new double[nLim[q2BinIndx]];
  double* aP1Lim1 = new double[nLim[q2BinIndx]];
  fstream fin1 (Form("pdf_lim%i_%i.list",q2BinIndx,-2),fstream::in);
  fin1.ignore(100,'\n');
  fin1.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) fin1>>aP5pLim1[i]>>aP1Lim1[i];
  fin1.close();

  double* aP5pLim2 = new double[nLim[q2BinIndx]];
  double* aP1Lim2 = new double[nLim[q2BinIndx]];
  fstream fin2 (Form("pdf_lim%i_%i.list",q2BinIndx,-1),fstream::in);
  fin2.ignore(100,'\n');
  fin2.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) fin2>>aP5pLim2[i]>>aP1Lim2[i];
  fin2.close();

  double* aP5pLim3 = new double[nLim[q2BinIndx]];
  double* aP1Lim3 = new double[nLim[q2BinIndx]];
  fstream fin3 (Form("pdf_lim%i_%i.list",q2BinIndx,0),fstream::in);
  fin3.ignore(100,'\n');
  fin3.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) fin3>>aP5pLim3[i]>>aP1Lim3[i];
  fin3.close();

  double* aP5pLim4 = new double[nLim[q2BinIndx]];
  double* aP1Lim4 = new double[nLim[q2BinIndx]];
  fstream fin4 (Form("pdf_lim%i_%i.list",q2BinIndx,1),fstream::in);
  fin4.ignore(100,'\n');
  fin4.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) fin4>>aP5pLim4[i]>>aP1Lim4[i];
  fin4.close();

  double* aP5pLim5 = new double[nLim[q2BinIndx]];
  double* aP1Lim5 = new double[nLim[q2BinIndx]];
  fstream fin5 (Form("pdf_lim%i_%i.list",q2BinIndx,2),fstream::in);
  fin5.ignore(100,'\n');
  fin5.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) fin5>>aP5pLim5[i]>>aP1Lim5[i];
  fin5.close();

  double* aP5pLim6 = new double[nLim[q2BinIndx]];
  double* aP1Lim6 = new double[nLim[q2BinIndx]];
  fstream fin6 (Form("pdf_lim%i_%i.list",q2BinIndx,3),fstream::in);
  fin6.ignore(100,'\n');
  fin6.ignore(100,'\n');
  for (int i=0; i<nLim[q2BinIndx]; i++) fin6>>aP5pLim6[i]>>aP1Lim6[i];
  fin6.close();

  fstream fin7 (Form("Data_pull/%i/bestfit.list",toyIndx),fstream::in);
  // for (int i=-2; i<q2BinIndx; i++) fin7.ignore(100,'\n');
  char s[100];
  fin7.getline(s,100);
  do {
    fin7.ignore(200,'\n');
    fin7>>bestP1[0];
  } while ( bestP1[0] != q2BinIndx );
  fin7>>bestP1[0];
  if (bestP1[0]==0) {
    bestP1 [0] = P1BF[q2BinIndx];
    bestP5p[0] = P5BF[q2BinIndx];
  } else fin7>>bestP1[0]>>bestP1[0]>>bestP5p[0];
  fin7.close();

  fgaus[q2BinIndx] = new TF2("f2","exp([0] - (x-[1])*(x-[1])/2/[3]/[3] - (y-[2])*(y-[2])/2/[4]/[4] - (x-[1])/[3]*(y-[2])/[4]*[5])",-1.,1.,-1.5,1.5);
  fgaus[q2BinIndx]->SetParameters(0,bestP1[0],bestP5p[0],0.2,0.2,0);
  // fgaus[q2BinIndx]->FixParameter(1,bestP1[0]);
  // fgaus[q2BinIndx]->FixParameter(2,bestP5p[0]);
  // fgaus[q2BinIndx]->FixParameter(5,0);
  fgaus[q2BinIndx]->SetParLimits(1,-2,2);
  fgaus[q2BinIndx]->SetParLimits(2,-3,3);
  fgaus[q2BinIndx]->SetParLimits(3,0.01,1);
  fgaus[q2BinIndx]->SetParLimits(4,0.01,1);
  fgaus[q2BinIndx]->SetParLimits(5,-1,1);

  TFitResultPtr fr = grlike[q2BinIndx]->Fit(fgaus[q2BinIndx],"S"/*+"Q0"*/);

  /*
  // cout<<fr->Status()<<" "<<fr->Ndf()<<endl;
  fstream fout ("failfit1.list",fstream::out|fstream::app);
  // cout<<fout.is_open()<<endl;
  if (fr->Status()!=0) fout<<q2BinIndx<<" "<<genIndx<<" " <<toyIndx<<" "<<minP1<<" "<<minP5<<endl;
  else if (fr->Ndf()<8) fout<<q2BinIndx<<" "<<genIndx<<" " <<toyIndx<<" "<<fgaus[q2BinIndx]->GetParameter(1)<<" "<<fgaus[q2BinIndx]->GetParameter(2)<<endl;
  fout.close();
  return;
  */
  if (fr->Status()!=0) {
    cout<<"Invalid fit result"<<endl;
    return;
  }

  hlike[q2BinIndx] = new TH2F(Form("hlike%i",q2BinIndx),Form("hlike%i",q2BinIndx),600,-1,1,900,-1.5,1.5);
  hboun[q2BinIndx] = new TH2S(Form("hboun%i",q2BinIndx),Form("hboun%i",q2BinIndx),600,-1,1,900,-1.5,1.5);

  double total=0;
  list<int> order_i;
  list<int> order_j;
  list<int>::iterator place_i;
  list<int>::iterator place_j;
  // list<int>::iterator pre_i;
  // list<int>::iterator pre_j;
  float max_phy=0;
  float max_P1, max_P5;

  for (int i=1; i<601; i++) for (int j=1; j<901; j++) {
      hboun[q2BinIndx]->SetBinContent(i,j,0);
      float xP1 = hlike[q2BinIndx]->GetXaxis()->GetBinCenter(i);
      float xP5 = hlike[q2BinIndx]->GetYaxis()->GetBinCenter(j);
      float P5lim=2;
      if (xP1>=aP1Lim[aP5pLim.size()-1]) P5lim = fabs(aP5pLim[aP5pLim.size()-1]);
      else if (xP1<=aP1Lim[0]) P5lim = fabs(aP5pLim[0]);
      else { 
	int indx;
	for (indx=1; indx<aP1Lim.size(); indx++) if (aP1Lim[indx]>xP1) break;
	P5lim = fabs(aP5pLim[indx] - (aP1Lim[indx]-xP1) * (aP5pLim[indx]-aP5pLim[indx-1]) / (aP1Lim[indx]-aP1Lim[indx-1]) );
      }
      if ( fabs(xP5) > P5lim ) hlike[q2BinIndx]->SetBinContent(i,j,0);  
      else {
	hlike[q2BinIndx]->SetBinContent(i,j,fgaus[q2BinIndx]->Eval(xP1,xP5)/exp(fgaus[q2BinIndx]->GetParameter(0)));
	total += hlike[q2BinIndx]->GetBinContent(i,j);
	if (hlike[q2BinIndx]->GetBinContent(i,j) > max_phy) {
	  max_phy = hlike[q2BinIndx]->GetBinContent(i,j);
	  max_P1 = xP1;
	  max_P5 = xP5;
	}
      }
    }
  cout<<0<<endl;
  for (int i=1; i<601; i++) for (int j=1; j<901; j++) if (hlike[q2BinIndx]->GetBinContent(i,j)>max_phy*exp(-1)) {
	if (order_i.end()==order_i.begin()) {
	  order_i.push_front(i);
	  order_j.push_front(j);
	} else {
	  place_i=order_i.begin();
	  place_j=order_j.begin();
	  if (hlike[q2BinIndx]->GetBinContent(*place_i,*place_j) < hlike[q2BinIndx]->GetBinContent(i,j)) {
	    order_i.push_front(i);
	    order_j.push_front(j);
	  } else {
	    do {
	      place_i++;
	      place_j++;
	      if (place_i == order_i.end()) break;
	    } while (hlike[q2BinIndx]->GetBinContent(*place_i,*place_j) > hlike[q2BinIndx]->GetBinContent(i,j));
	    order_i.insert(place_i,i);
	    order_j.insert(place_j,j);
	  }
	}
      }
  cout<<0<<endl;


      // 	  for (place_i=order_i.begin(); pos<order_i.size(); pos++) if (hlike[q2BinIndx]->GetBinContent(order_i[pos],order_j[pos])<hlike[q2BinIndx]->GetBinContent(i,j)) break;
      // 	order_i.insert(order_i.begin()+pos,i);
      // 	order_j.insert(order_j.begin()+pos,j);
      // }

  // cout<<"Histo filled"<<endl;
  double prob=0;

  // do {
  //   double newmax = 0;
  //   int imax, jmax;
  //   for (int i=1; i<601; i++) for (int j=1; j<901; j++) if (hboun[q2BinIndx]->GetBinContent(i,j)==0 && hlike[q2BinIndx]->GetBinContent(i,j)>newmax ) {
  // 	  newmax = hlike[q2BinIndx]->GetBinContent(i,j);
  // 	  imax = i;
  // 	  jmax = j;
  // 	}
  //   if (newmax==0) {
  //     cout<<"Limit not reached... please check"<<endl;
  //     break;
  //   }
  //   hboun[q2BinIndx]->SetBinContent(i,j,1);
  //   prob += newmax;
  // } while (prob < total*0.68);

  // for (int i=0; prob <= total*0.6827; i++) {

  double P1up=-1;
  double P1do=1;
  double P5up=-1.5;
  double P5do=1.5;

  float lastval;
  place_i=order_i.begin();
  place_j=order_j.begin();
  do {
    hboun[q2BinIndx]->SetBinContent(*place_i,*place_j,1);
    prob += lastval = hlike[q2BinIndx]->GetBinContent(*place_i,*place_j);
    if (hlike[q2BinIndx]->GetXaxis()->GetBinLowEdge(*place_i) < P1do) P1do = hlike[q2BinIndx]->GetXaxis()->GetBinLowEdge(*place_i);
    if (hlike[q2BinIndx]->GetXaxis()->GetBinUpEdge(*place_i) > P1up) P1up = hlike[q2BinIndx]->GetXaxis()->GetBinUpEdge(*place_i);
    if (hlike[q2BinIndx]->GetYaxis()->GetBinLowEdge(*place_j) < P5do) P5do = hlike[q2BinIndx]->GetYaxis()->GetBinLowEdge(*place_j);
    if (hlike[q2BinIndx]->GetYaxis()->GetBinUpEdge(*place_j) > P5up) P5up = hlike[q2BinIndx]->GetYaxis()->GetBinUpEdge(*place_j);
    place_i++;
    place_j++;
    if (place_i==order_i.end()) {
      cout<<"Limit not reached... please check"<<endl;
      break;
    }
  } while (prob <= total*0.39347);
  cout<<0<<endl;


  // for (int i=0; prob <= total*0.39347; i++) {
  //   if (i==order_i.size()) {
  //     cout<<"Limit not reached... please check"<<endl;
  //     break;
  //   }
  //   hboun[q2BinIndx]->SetBinContent(order_i[i],order_j[i],1);
  //   prob += hlike[q2BinIndx]->GetBinContent(order_i[i],order_j[i]);
  //   lastval = hlike[q2BinIndx]->GetBinContent(order_i[i],order_j[i]);
  // }

  // grlike[q2BinIndx]->Draw("PCOL");
  // hlike[q2BinIndx]->Draw("CONT");
  // hboun[q2BinIndx]->Draw("CONT");
  // cout<<prob<<" "<<total<<endl;
  // return;

  // cout<<"Limit found"<<endl;

  // TGraph* boundary = new TGraph();
  // short last=0;
  // int bou_cnt = 0;


  // for (int i=1; i<601; i++) for (int j=1; j<901; j++) { 
  //     if ( j>1 && abs(hboun[q2BinIndx]->GetBinContent(i,j)-last) == 1 ) {
  // 	boundary->SetPoint( bou_cnt, hlike[q2BinIndx]->GetXaxis()->GetBinCenter(i), 0.5*(hlike[q2BinIndx]->GetYaxis()->GetBinCenter(j)+hlike[q2BinIndx]->GetYaxis()->GetBinCenter(j-1)) );
  // 	bou_cnt++;
  // 	if ( P5up<0.5*(hlike[q2BinIndx]->GetYaxis()->GetBinCenter(j)+hlike[q2BinIndx]->GetYaxis()->GetBinCenter(j-1)) )
  // 	  P5up = 0.5*(hlike[q2BinIndx]->GetYaxis()->GetBinCenter(j)+hlike[q2BinIndx]->GetYaxis()->GetBinCenter(j-1));
  // 	if ( P5do>0.5*(hlike[q2BinIndx]->GetYaxis()->GetBinCenter(j)+hlike[q2BinIndx]->GetYaxis()->GetBinCenter(j-1)) )
  //         P5do = 0.5*(hlike[q2BinIndx]->GetYaxis()->GetBinCenter(j)+hlike[q2BinIndx]->GetYaxis()->GetBinCenter(j-1));
  //     }
  //     last = hboun[q2BinIndx]->GetBinContent(i,j);
  //   }
  // for (int j=1; j<901; j++) for (int i=1; i<601; i++) {
  //     if ( i>1 && abs(hboun[q2BinIndx]->GetBinContent(i,j)-last) == 1 ) {
  // 	boundary->SetPoint( bou_cnt, 0.5*(hlike[q2BinIndx]->GetXaxis()->GetBinCenter(i)+hlike[q2BinIndx]->GetXaxis()->GetBinCenter(i-1)), hlike[q2BinIndx]->GetYaxis()->GetBinCenter(j) );
  // 	bou_cnt++;
  // 	if ( P1up<0.5*(hlike[q2BinIndx]->GetXaxis()->GetBinCenter(i)+hlike[q2BinIndx]->GetXaxis()->GetBinCenter(i-1)) )
  // 	  P1up = 0.5*(hlike[q2BinIndx]->GetXaxis()->GetBinCenter(i)+hlike[q2BinIndx]->GetXaxis()->GetBinCenter(i-1));
  // 	if ( P1do>0.5*(hlike[q2BinIndx]->GetXaxis()->GetBinCenter(i)+hlike[q2BinIndx]->GetXaxis()->GetBinCenter(i-1)) )
  //         P1do = 0.5*(hlike[q2BinIndx]->GetXaxis()->GetBinCenter(i)+hlike[q2BinIndx]->GetXaxis()->GetBinCenter(i-1));
  //     }
  //     last = hboun[q2BinIndx]->GetBinContent(i,j);
  //   }

  // cout<<"Limit drawn"<<endl;

  if (bestP1[0]==P1BF[q2BinIndx]) cout<<"0 0"<<endl;
  else cout<<bestP1[0]<<" "<<bestP5p[0]<<endl;

  bestP1 [0]=max_P1;
  bestP5p[0]=max_P5;

  cout<</*"P1 result: "<<*/bestP1[0]<<" "<<P1up<<" "<<P1do<<endl; //" "<<(P1BF[q2BinIndx]>bestP1[0]?(bestP1[0]-P1BF[q2BinIndx])/(P1up-bestP1[0]):(bestP1[0]-P1BF[q2BinIndx])/(bestP1[0]-P1do))<<endl;
  cout<</*"P5 result: "<<*/bestP5p[0]<<" "<<P5up<<" "<<P5do<<endl; //" "<<(P5BF[q2BinIndx]>bestP5p[0]?(bestP5p[0]-P5BF[q2BinIndx])/(P5up-bestP5p[0]):(bestP5p[0]-P5BF[q2BinIndx])/(bestP5p[0]-P5do))<<endl;

  cout<<lastval<<" 0.39347"<<endl;

  /*for (int i=0; i<vLik.size(); i++) {
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

  double* aP5p = new double[vP1.size()];
  double* aP1  = new double[vP1.size()];
  for (int i=0; i<vP1.size(); i++) {
    aP5p[i]=vP5p[i];
    aP1 [i]=vP1 [i];
    // cout<<vLik[i]-parab[q2BinIndx]->Eval(vP1 [i],vP5p[i])<<"  ";
    // if ( i%4 == 3 ) cout<<endl;
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

  gLim [q2BinIndx] = new TGraph(nLim[q2BinIndx], aP1Lim , aP5pLim );
  gLim0[q2BinIndx] = new TGraph(nLim[q2BinIndx], aP1Lim0, aP5pLim0);
  gLim1[q2BinIndx] = new TGraph(nLim[q2BinIndx], aP1Lim1, aP5pLim1);
  gLim2[q2BinIndx] = new TGraph(nLim[q2BinIndx], aP1Lim2, aP5pLim2);
  gLim3[q2BinIndx] = new TGraph(nLim[q2BinIndx], aP1Lim3, aP5pLim3);
  gLim4[q2BinIndx] = new TGraph(nLim[q2BinIndx], aP1Lim4, aP5pLim4);
  gLim5[q2BinIndx] = new TGraph(nLim[q2BinIndx], aP1Lim5, aP5pLim5);
  gLim6[q2BinIndx] = new TGraph(nLim[q2BinIndx], aP1Lim6, aP5pLim6);
  gBest[q2BinIndx] = new TGraph(1, bestP1, bestP5p);

  // fgaus[q2BinIndx]->Draw("CONT1");
  // grlike[q2BinIndx]->Draw("CONT1");

  // g2d[q2BinIndx] = new TGraph( vP1.size(), aP1, aP5p);
  // g2d[q2BinIndx]->GetXaxis()->SetLimits(-1,1);
  // g2d[q2BinIndx]->GetXaxis()->SetTitle("P1");
  // g2d[q2BinIndx]->GetYaxis()->SetTitle("P'5");
  // g2d[q2BinIndx]->SetTitle("");
  // g2d[q2BinIndx]->Draw("AP");
  // g2d[q2BinIndx]->SetMarkerStyle(7);
  // g2d[q2BinIndx]->SetMarkerColor(4);


  boundary->GetXaxis()->SetLimits(-1,1);
  boundary->GetXaxis()->SetTitle("P1");
  boundary->GetYaxis()->SetTitle("P'5");
  boundary->SetTitle("");
  boundary->Draw("AP");
  boundary->SetMarkerStyle(7);
  boundary->SetMarkerColor(2);

  if (v2sP1.size()>0) {
    g2d2s[q2BinIndx] = new TGraph( v2sP1.size(), a2sP1, a2sP5p);
    g2d2s[q2BinIndx]->GetXaxis()->SetLimits(-1,1);
    // g2d1[q2BinIndx]->GetYaxis()->SetLimits(-1.2,.8);
    g2d2s[q2BinIndx]->GetXaxis()->SetTitle("P1");
    g2d2s[q2BinIndx]->GetYaxis()->SetTitle("P'5");
    g2d2s[q2BinIndx]->SetTitle("");
    g2d2s[q2BinIndx]->Draw("Psame");
    g2d2s[q2BinIndx]->SetMarkerStyle(7);
    g2d2s[q2BinIndx]->SetMarkerColor(3);
  }
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

  // gBest[q2BinIndx]->Draw("sameP");
  // gBest[q2BinIndx]->SetMarkerStyle(34);
  // gBest[q2BinIndx]->SetMarkerColor(0);

  gLim [q2BinIndx]->Draw("samePL");
  gLim [q2BinIndx]->SetMarkerStyle(7);
  gLim [q2BinIndx]->SetMarkerColor(6);
  gLim [q2BinIndx]->SetLineColor(6);
  gLim [q2BinIndx]->SetLineWidth(2);

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

  // can[q2BinIndx]->SaveAs(Form("Data_o/scan2d_b%i_v5.pdf",q2BinIndx));

  */  return;
}

void fit_pull (int toyIndx=0, int q2BinIndx = 0, int genIndx=0)
{
  // for (q2BinIndx=1;q2BinIndx<10;q2BinIndx++) for (genIndx=0;genIndx<9;genIndx++) for (toyIndx=0;toyIndx<100;toyIndx++) read(q2BinIndx-1,toyIndx,genIndx);

  if (q2BinIndx==0) for (q2BinIndx=1;q2BinIndx<10;q2BinIndx++) read(q2BinIndx-1,toyIndx,genIndx);
  else if (q2BinIndx>0 && q2BinIndx<10) read(q2BinIndx-1,toyIndx,genIndx);
  
  return;
}
