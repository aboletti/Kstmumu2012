#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

string STEFANO="/lustre/cmswork/lacaprar/Ana/CMSSW_5_3_28/src/B0Analysis/B0KstMuMu/plugins/LHCbTest_Frac/";

using namespace RooFit ;

void getFitResults(int ibin=0){
  double bins[7]={1.8, 3.25, 5.0, 7.0,  11.75,  16, 18};
  double binse[7]={0.7, 1.25, 1.0, 1.0, 0.75,  1.0, 1.0};

  double FlIn[7]={0.660 ,0.876 ,0.611,.579 ,0.493,.349 ,.354 };
  double P5pIn[7]={0.289,-0.066,-.300,-.505,-.655,-.662,-.676};
  double P1In[7]={-0.451,0.571 , .180,-.199,0.745,-.436,-.581};

  double FlFit[7]={-2,-2,-2,-2,-2,-2,-2,};
  double FlFitE[7]={0,0,0,0,0,0,0};

  double P5pFit[7]={-2,-2,-2,-2,-2,-2,-2};
  double P5pFitE[7]={0,0,0,0,0,0,0};

  double P1Fit[7]={-2,-2,-2,-2,-2,-2,-2};
  double P1FitE[7]={0,0,0,0,0,0,0};

  gStyle->SetStatColor(kWhite);
  gStyle->SetStatStyle(0);

  gStyle->SetStatFont(42);
  gStyle->SetStatFontSize(0.04);
  gStyle->SetStatH(0.1);
  gStyle->SetStatW(0.5);
  gStyle->SetTitleXOffset(0.9);
  gStyle->SetTitleYOffset(0.9);

  TMultiGraph* mlFl=new TMultiGraph();
  TMultiGraph* mlP1=new TMultiGraph();
  TMultiGraph* mlP5p=new TMultiGraph();


  TLegend *legend=new TLegend(0.5,0.65,0.88,0.85);
  legend->SetTextFont(72);
  legend->SetTextSize(0.05);
  legend->SetFillColor(0);

  TGraph *gFlLHCb = new TGraph(7, bins, FlIn);
  gFlLHCb->SetMarkerStyle(20);
  gFlLHCb->SetMarkerColor(kBlue);
  mlFl->Add(gFlLHCb);
  legend->AddEntry(gFlLHCb,"LHCb","p");

  TGraph *gP1LHCb = new TGraph(7, bins, P1In);
  gP1LHCb->SetMarkerStyle(20);
  gP1LHCb->SetMarkerColor(kBlue);
  mlP1->Add(gP1LHCb);

  TGraph *gP5pLHCb = new TGraph(7, bins, P5pIn);
  gP5pLHCb->SetMarkerStyle(20);
  gP5pLHCb->SetMarkerColor(kBlue);
  mlP5p->Add(gP5pLHCb);


  TH1F *hFlBest = new TH1F("hFlBest", "Fl", 50, -1., 1.);
  TH1F *hP1Best= new TH1F("hP1Best", "P1", 50, -1., 1.);
  TH1F *hP5pBest = new TH1F("hP5pBest", "P5p", 50, -1., 1.);

  TH1F *hFl = new TH1F("hFl", "Fl", 50, -1., 1.);
  TH1F *hP1= new TH1F("hP1", "P1", 50, -1., 1.);
  TH1F *hP5p = new TH1F("hP5p", "P5p", 50, -1., 1.);


  TLatex* tt=new TLatex();
  tt->SetNDC();
  tt->SetTextSize(0.08);

  TLine* tl=new TLine();
  tl->SetLineWidth(2);
  tl->SetLineColor(kRed+2);

  TLine* tlg=new TLine();
  tlg->SetLineWidth(2);
  tlg->SetLineColor(kGreen+2);

  // int bin=0;
  for (int bin=0; bin<1; bin++)
  /*for (int ibin=0; ibin<9;++ibin)*/ {
    if (ibin==4 || ibin==6) continue;
    int nToys=101; //10;
    int nJobs=10; //20;
    int goodToy=0;
    for (int toy=0; toy<nToys; ++toy) {
      /* if (toy%10==0)*/ cout << "------------------------------------"<<endl<<"Toy " << toy << endl;
      double bestEdm=1E10;
      double bestNLL=1E10;
      double bestFl=-100,bestP1=-100,bestP5p=-100;
      for (int i=0; i<nJobs;++i) {
        TString filename=STEFANO+Form("../LHCbTest/PayLoad6_Toy%i/Job%i/Fitresult5_2.root",toy,i);
	if (ibin>6) filename=STEFANO+Form("./PayLoad6_Toy%i/Job%i/Fitresult5_2.root",toy,i);
        // TString filename=STEFANO+Form("./Prod2/PayLoad106_Toy%i_%i/Fitresult1.root",toy,toy);
        // TString filename=STEFANO+Form("./Prod2/PayLoad206_Toy%i_%i/Fitresult3.root",toy,toy);
        //cout << "Opening " << filename << endl;

        if ( gSystem->AccessPathName(filename)) continue;
        TFile *_file0 = TFile::Open(filename);
        if (!_file0) {
	  delete _file0;
	  continue;
	}
        //cout << "Job " << i << endl;
        RooFitResult* res=(RooFitResult*)_file0->Get(Form("fitResult_Bin%i",ibin));

        //cout << "res " << res << endl;
        if (!res) continue;
	res->printValue(cout);
        //res->Print();
        // cout << " Toy " << toy << "/" << i << " bin "<< ibin << " res " << res->status() << " " << res->covQual()<< endl;
        if (!(res->covQual() == 3 && res->status() == 0)) continue;
        //cout << "Toy " << toy << " bin "<< ibin << " " << filename << endl;

        RooArgList l=res->floatParsFinal();
        double edm= res->edm();
        double minNll= res->minNll();
        //double Fl= ((RooRealVar*)l.find("FlS"))->getValV();
        //double FlE= ((RooRealVar*)l.find("FlS"))->getValV();
        double Fl=FlIn[bin];
        double P1= ((RooRealVar*)l.find("P1S"))->getValV();
        double P1E= ((RooRealVar*)l.find("P1S"))->getAsymErrorHi();
        double P5p= ((RooRealVar*)l.find("P5pS"))->getValV();
        // cout << i << " Edm " << edm << " Nll " << minNll << " FL " << Fl << " P1=" << P1 << "+/-" << P1E << " P5=" << P5p << endl;

        // if (edm<bestEdm) {
        //   bestEdm=edm;
        //   bestFl = Fl ;
        //   bestP1 = P1 ;
        //   bestP5p= P5p;
        // }
        if (minNll<bestNLL) {
          bestNLL=minNll;
          bestFl = Fl ;
          bestP1 = P1 ;
          bestP5p= P5p;
        }

        // cout << "Fl " << Fl << endl;
        // cout << "P1 " << P1 << endl;
        // cout << "P5p " << P5p << endl;
        hFl->Fill(Fl);
        hP1->Fill(P1);
        hP5p->Fill(P5p);

        _file0->Close();
        delete _file0;
      }
      if (bestFl > -1) {
	// cout<<"*";
	hFlBest->Fill(bestFl);
	hP1Best->Fill(bestP1);
	hP5pBest->Fill(bestP5p);
      }
    }

    FlFit[bin]=hFlBest->GetMean();
    FlFitE[bin]=hFlBest->GetRMS()/sqrt(hFlBest->Integral());
    P1Fit[bin]=hP1Best->GetMean();
    P1FitE[bin]=hP1Best->GetRMS()/sqrt(hP1Best->Integral());
    P5pFit[bin]=hP5pBest->GetMean();
    P5pFitE[bin]=hP5pBest->GetRMS()/sqrt(hP5pBest->Integral());
    // cout << "Bin " << bin << endl;
    // cout << "Fl " << hFlBest->Integral(1,50)< < " " << FlFit[bin] << " " <<  FlIn[bin] << endl;
    // cout << "P1 " << hP1Best->Integral(1,50) << " " << P1Fit[bin] << " " <<  P1In[bin] << endl;
    // cout << "P5p " << hP5pBest->Integral(1,50) << " " << P5pFit[bin] << " " << P5pIn[bin] << endl;


    if (true) {
      TString nameBest=Form("cBin%i_Best6",bin);
      // TString nameBest=Form("cBin%i_Best106",bin);
      // TString nameBest=Form("cBin%i_Best206",bin);
      TCanvas* cBinBest=new TCanvas(nameBest,nameBest,800,400);
      cBinBest->Divide(3);
      cBinBest->cd(1);
      hFlBest->Draw();
      tl->DrawLine(FlIn[bin],0,FlIn[bin],hFlBest->GetMaximum()*0.80);
      tt->DrawLatex(0.1, 0.02, Form("F_{L}^{in}=%0.3f",FlIn[bin]));

      cBinBest->cd(2);
      hP1Best->Draw();
      tl->DrawLine(P1In[bin],0,P1In[bin],hP1Best->GetMaximum()*0.80);
      tt->DrawLatex(0.1, 0.02, Form("P_{1}^{in}=%0.3f",P1In[bin]));


      cBinBest->cd(3);
      hP5pBest->Draw();
      tl->DrawLine(P5pIn[bin],0,P5pIn[bin],hP5pBest->GetMaximum()*0.80);
      tt->DrawLatex(0.1, 0.02, Form("P'_{5}^{in}=%0.3f",P5pIn[bin]));

      // cBinBest->Print(nameBest+".pdf");
      // cBinBest->Print(nameBest+".png");
    }


    if (false) {
      TString name=Form("cBin%i",bin);
      TCanvas* cBin=new TCanvas(name,name,800,400);
      cBin->Divide(3);
      cBin->cd(1);
      hFl->Draw();
      tl->DrawLine(FlIn[bin],0,FlIn[bin],hFl->GetMaximum()*0.80);
      tlg->DrawLine(bestFl,0,bestFl,hFl->GetMaximum()*0.80);

      tt->DrawLatex(0.1, 0.02, Form("F_{L}^{in}=%0.3f F_{L}^{fit}=%0.3",FlIn[bin],bestFl));

      cBin->cd(2);
      hP1->Draw();
      tl->DrawLine(P1In[bin],0,P1In[bin],hP1->GetMaximum()*0.80);
      tlg->DrawLine(bestP1,0,bestP1,hP1->GetMaximum()*0.80);

      tt->DrawLatex(0.3, 0.02, Form("P_{1}^{in}=%0.3f P_{1}^{fit}=%0.3f",P1In[bin],bestP1));

      cBin->cd(3);
      hP5p->Draw();
      tl->DrawLine(P5pIn[bin],0,P5pIn[bin],hP5p->GetMaximum()*0.80);
      tlg->DrawLine(bestP5p,0,bestP5p,hP5p->GetMaximum()*0.80);
      tt->DrawLatex(0.5, 0.02, Form("P'_{5}^{in}=%0.3f P'_{5}^{fit}=%0.3f",P5pIn[bin],bestP5p));

      // cBin->Print(name+".pdf");
    }
    bin++;
    hFlBest->Reset();
    hP1Best->Reset();
    hP5pBest->Reset();
  }

  TGraphErrors *gFl = new TGraphErrors(7, bins, FlFit, binse, FlFitE);
  gFl->SetMarkerStyle(24);                                    
  gFl->SetMarkerColor(kRed);                                  
  legend->AddEntry(gFl,"Toy","pe");

  TGraphErrors *gP1 = new TGraphErrors(7, bins, P1Fit, binse, P1FitE);
  gP1->SetMarkerStyle(24);                                    
  gP1->SetMarkerColor(kRed);                                  

  TGraphErrors *gP5p = new TGraphErrors(7, bins, P5pFit, binse, P5pFitE);
  gP5p->SetMarkerStyle(24);
  gP5p->SetMarkerColor(kRed);

  mlFl->Add(gFl);
  mlP1->Add(gP1);
  mlP5p->Add(gP5p);

  TString nameBestGr="cBin_BestGraph6";
  // TString nameBestGr="cBin_BestGraph106";
  // TString nameBestGr="cBin_BestGraph206";
  TCanvas* cBinBestGraph=new TCanvas(nameBestGr,"Graph ",800,400);
  cBinBestGraph->Divide(3);

  cBinBestGraph->cd(1);
  mlFl->Draw("AP");
  mlFl->GetYaxis()->SetRangeUser(0.,1.);
  mlFl->GetYaxis()->SetTitle("F_{L}");
  mlFl->GetXaxis()->SetTitle("q^{2}");
  legend->Draw();
  gPad->Modified();

  cBinBestGraph->cd(2);
  mlP1->Draw("AP");
  mlP1->GetYaxis()->SetRangeUser(-1,1.);
  mlP1->GetYaxis()->SetTitle("P_{1}");
  mlP1->GetXaxis()->SetTitle("q^{2}");
  gPad->Modified();

  cBinBestGraph->cd(3);
  mlP5p->Draw("AP");
  mlP5p->GetYaxis()->SetRangeUser(-1,1.);
  mlP5p->GetYaxis()->SetTitle("P'_{5}");
  mlP5p->GetXaxis()->SetTitle("q^{2}");
  gPad->Modified();

  // cBinBestGraph->Print(nameBestGr+".pdf");
  // cBinBestGraph->Print(nameBestGr+".png");

}
