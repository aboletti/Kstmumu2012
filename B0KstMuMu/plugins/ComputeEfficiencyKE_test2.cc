#include <RooNumIntConfig.h>
#include <RooGenericPdf.h>
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooPolynomial.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"
#include "RooProdPdf.h"
#include <RooFunctorBinding.h>
#include "RooBinning.h"
#include "RooWorkspace.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooRandom.h"

using namespace RooFit ;

#include <TROOT.h>
#include <TApplication.h>
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TMath.h"
#include <TLorentzVector.h>

const string fileNameOutput = "effKEpdf_out.root";

RooRealVar ctK("ctK","cos#theta_{K}",-1,1);
RooRealVar ctL("ctL","cos#theta_{L}",0,1);
RooRealVar phi("phi","#phi",0,TMath::Pi());
RooArgList ctKctLphi(ctK,ctL,phi);

string filename1;
string filename2;
string filename3;
string filename4;

int n_events1;
int n_events2;
int n_events3;
int n_events4;

RooDataSet* data1;
RooDataSet* data2;
RooDataSet* data3;
RooDataSet* data4;

TString hname;

bool Generate (int mist) {
  // TFile* fin1 = new TFile( filename1.c_str() );
  // TFile* fin2 = new TFile( filename2.c_str() );
  // TH3* h3rNm = (TH3*)fin2->Get("h3recoNum");
  // h3rNm->SetName("h3recoNum_mist");
  // TH3* h3gD  = (TH3*)fin1->Get("h3genDen" );
  // TH3* h3gN  = (TH3*)fin1->Get("h3genNum" );
  // TH3* h3rD  = (TH3*)fin1->Get("h3recoDen");
  // TH3* h3rN  = (TH3*)fin1->Get("h3recoNum");

  if ( mist ) {
    TFile* fin = new TFile( filename2.c_str() );
    TH3* h3 = (TH3*)fin->Get("h3recoNum");

    if (h3) {cout << "histo loaded" << endl;

      n_events1 = h3->Integral();
   
      RooDataHist dataHist("dataHist","dataHist", ctKctLphi, Import(*h3,kTRUE)); cout << "dataHist created" << endl;
      RooHistPdf histPdf("histPdf","histPdf", ctKctLphi, dataHist, 0); cout << "pdf created" << endl;

      data1 = histPdf.generate( ctKctLphi, n_events1, Name("data1"), AutoBinned(kTRUE), Extended()); cout << "dataset generated" << endl;
    } else {
      cout << "Cannot load histos" << endl;
      return false;
    }
    delete h3;
    fin->Close();
    delete fin;
  } else {
    TFile* fin = new TFile( filename1.c_str() );
    TH3* h3gD = (TH3*)fin->Get("h3genDen" );
    TH3* h3gN = (TH3*)fin->Get("h3genNum" );
    TH3* h3rD = (TH3*)fin->Get("h3recoDen");
    TH3* h3rN = (TH3*)fin->Get("h3recoNum");

    if (h3gN && h3gD && h3rD && h3rN) {cout << "histo loaded" << endl;

      n_events1 = h3gD->Integral();
      n_events2 = h3gN->Integral();
      n_events3 = h3rD->Integral();
      n_events4 = h3rN->Integral();
   
      RooDataHist dataHist1("dataHist1","dataHist1", ctKctLphi, Import(*h3gD,kTRUE));
      RooDataHist dataHist2("dataHist2","dataHist2", ctKctLphi, Import(*h3gN,kTRUE));
      RooDataHist dataHist3("dataHist3","dataHist3", ctKctLphi, Import(*h3rD,kTRUE));
      RooDataHist dataHist4("dataHist4","dataHist4", ctKctLphi, Import(*h3rN,kTRUE)); cout << "dataHist created" << endl;
      RooHistPdf histPdf1("histPdf1","histPdf1", ctKctLphi, dataHist1, 0);
      RooHistPdf histPdf2("histPdf2","histPdf2", ctKctLphi, dataHist2, 0);
      RooHistPdf histPdf3("histPdf3","histPdf3", ctKctLphi, dataHist3, 0);
      RooHistPdf histPdf4("histPdf4","histPdf4", ctKctLphi, dataHist4, 0); cout << "pdf created" << endl;

      data1 = histPdf1.generate( ctKctLphi, n_events1, Name("data1"), AutoBinned(kTRUE), Extended());
      data2 = histPdf2.generate( ctKctLphi, n_events2, Name("data2"), AutoBinned(kTRUE), Extended());
      data3 = histPdf3.generate( ctKctLphi, n_events3, Name("data3"), AutoBinned(kTRUE), Extended());
      data4 = histPdf4.generate( ctKctLphi, n_events4, Name("data4"), AutoBinned(kTRUE), Extended()); cout << "dataset generated" << endl;
    } else {
      cout << "Cannot load histos" << endl;
      return false;
    }
    delete h3gD;
    delete h3gN;
    delete h3rD;
    delete h3rN;
    fin->Close();
    delete fin;
  }
  // fin1->Close();
  // fin2->Close();
  return true;
}

void ComputeEfficiency (int mist) {
  int nbin = 40;

  if ( mist ) {
    int events1 = (int)data1->sumEntries();

    TFile* fout1 = new TFile( filename4.c_str(), "RECREATE" );

    RooNDKeysPdf* KEpdf1 = new RooNDKeysPdf( "KEpdf1", "KEpdf1", ctKctLphi, *data1, "M", 1.0); cout << "new KEpdf created" << endl;
    TH3* h3new = (TH3*) KEpdf1->createHistogram(hname,ctK,Binning(nbin),YVar(ctL,Binning(nbin)),ZVar(phi,Binning(nbin)));
    h3new->Scale(events1/h3new->Integral());
    cout<<n_events1<<" "<<events1<<" "<<h3new->Integral()<<endl;
    h3new->Write(hname);
    delete h3new;
    delete KEpdf1;

    fout1->Close(); cout << "new file closed" << endl;
    delete fout1;
  } else {
    int events1 = (int)data1->sumEntries(); //numEntries();
    int events2 = (int)data2->sumEntries();
    int events3 = (int)data3->sumEntries();
    int events4 = (int)data4->sumEntries();

    TFile* fout1 = new TFile( filename3.c_str(), "RECREATE" );

    RooNDKeysPdf* KEpdf1 = new RooNDKeysPdf( "KEpdf1", "KEpdf1", ctKctLphi, *data1, "M", 0.5);
    TH3* h3genDen = (TH3*) KEpdf1->createHistogram("h3genDen",ctK,Binning(nbin),YVar(ctL,Binning(nbin)),ZVar(phi,Binning(nbin)));
    h3genDen->Scale(events1/h3genDen->Integral());
    cout<<n_events1<<" "<<events1<<" "<<h3genDen->Integral()<<endl;
    h3genDen->Write("h3genDen");
    delete h3genDen;
    delete KEpdf1;

    RooNDKeysPdf* KEpdf2 = new RooNDKeysPdf( "KEpdf2", "KEpdf2", ctKctLphi, *data2, "M", 0.5);
    TH3* h3genNum = (TH3*) KEpdf2->createHistogram("h3genNum",ctK,Binning(nbin),YVar(ctL,Binning(nbin)),ZVar(phi,Binning(nbin)));
    h3genNum->Scale(events2/h3genNum->Integral());
    cout<<n_events2<<" "<<events2<<" "<<h3genNum->Integral()<<endl;
    h3genNum->Write("h3genNum");
    delete h3genNum;
    delete KEpdf2;

    RooNDKeysPdf* KEpdf3 = new RooNDKeysPdf( "KEpdf3", "KEpdf3", ctKctLphi, *data3, "M", 0.5);
    TH3* h3recoDen = (TH3*) KEpdf3->createHistogram("h3recoDen",ctK,Binning(nbin),YVar(ctL,Binning(nbin)),ZVar(phi,Binning(nbin)));
    h3recoDen->Scale(events3/h3recoDen->Integral());
    cout<<n_events3<<" "<<events3<<" "<<h3recoDen->Integral()<<endl;
    h3recoDen->Write("h3recoDen");
    delete h3recoDen;
    delete KEpdf3;

    RooNDKeysPdf* KEpdf4 = new RooNDKeysPdf( "KEpdf4", "KEpdf4", ctKctLphi, *data4, "M", 0.5); cout << "new KEpdf created" << endl;
    TH3* h3recoNum = (TH3*) KEpdf4->createHistogram("h3recoNum",ctK,Binning(nbin),YVar(ctL,Binning(nbin)),ZVar(phi,Binning(nbin))); cout << "new effHisto created" << endl;
    h3recoNum->Scale(events4/h3recoNum->Integral());
    cout<<n_events4<<" "<<events4<<" "<<h3recoNum->Integral()<<endl;
    h3recoNum->Write("h3recoNum");
    delete h3recoNum;
    delete KEpdf4;

    fout1->Close(); cout << "new file closed" << endl;
    delete fout1;
  }
}

int main(int argc, char** argv)
{
  if (argc > 0)
  {
    int q2BinIndx=0;
    string test_indx="ttttt";
    int mist=0;
    string signalname="";
    if ( argc > 1 ) q2BinIndx = atoi(argv[1]);
    if ( argc > 2 ) test_indx = argv[2];
    if ( argc > 3 ) mist = atoi(argv[3]);
    if ( argc > 4 ) signalname = argv[4];

    string folder = "/lustre/cmswork/boletti/Kstmumu/CMSSW_5_3_28/src/B0Analysis/B0KstMuMu/plugins/";
    filename1 = folder + "Closure_rt" + signalname + Form("_0.5/effKEpdf_b%i.root",q2BinIndx);
    filename2 = folder + "Closure_wt" + signalname + Form("_1.0/effKEpdf_b%i.root",q2BinIndx);
    filename3 = folder + "Closure_sys2_rt" + signalname + "_0.5/" + test_indx + Form("/effKEpdf_b%i.root",q2BinIndx);
    filename4 = folder + "Closure_sys2_wt" + signalname + "_1.0/" + test_indx + Form("/effKEpdf_b%i.root",q2BinIndx);

    if ( argc > 2 ) RooRandom::randomGenerator()->SetSeed(atoi(argv[2])+1);

    // do all q2Bins at once
    if (q2BinIndx==0)
      for (q2BinIndx=1; q2BinIndx<10; ++q2BinIndx) {
	if (Generate (mist)) ComputeEfficiency(mist);
	else return EXIT_FAILURE;
      }
    else if (q2BinIndx>0 && q2BinIndx<10) {
      if (Generate(mist)) ComputeEfficiency(mist);
      else return EXIT_FAILURE;
    }
    else {
      cout<<"q2Bin must be greater than 0 and smaller than 10. FAILURE!"<<endl;
      return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
  }
  else
  {
    cout << "Parameter missing: " << endl;
    cout << "./" << argv[0] << " <q2Bin (0 for all)> " << " <plot histo (decfault=0/no)> " << " <batch running (default=0/no) " << endl;

    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
