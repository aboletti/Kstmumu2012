#include <RooGenericPdf.h>
#include "RooHistPdf.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooUniform.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "RooBinning.h"
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
//const string folder;
TFile* fin;
TFile* fin1;

RooRealVar ctK("ctK","cos#theta_{K}",-1,1);
RooRealVar ctL("ctL","cos#theta_{L}",0,1);
RooRealVar phi("phi","#phi",0,TMath::Pi());
RooArgList ctKctLphi(ctK,ctL,phi);

void Convert (int q2BinIndx, bool mist=0) {
  fin = new TFile( Form("effKEpdf_b%i.root",q2BinIndx) );
  //if (mist) fin1 = new TFile( "../"+folder+Form("/effKEpdf_b%i.root",q2BinIndx) );
  if (mist) fin1 = new TFile( Form("../Closure_cc_rt_p2_0.5/effKEpdf_b%i.root",q2BinIndx) );
  //if (mist) fin1 = new TFile( Form("../Closure_cc_rt_jp_0.5/effKEpdf_b%i.root",q2BinIndx) );
  //if (mist) fin1 = new TFile( Form("../Closure_cc_rt_0.5/effKEpdf_b%i.root",q2BinIndx) );
  else fin1 = fin;
  TH3* h3genNum  = (TH3*)fin1->Get("h3genNum" );
  TH3* h3genDen  = (TH3*)fin1->Get("h3genDen" );
  TH3* h3recoDen = (TH3*)fin1->Get("h3recoDen");
  TH3* h3recoNum = (TH3*)fin->Get("h3recoNum");

  if (h3genNum && h3genDen && h3recoNum && h3recoDen) {

    TFile* fout = new TFile(fileNameOutput.c_str(), "UPDATE");

    TH3* h3genEff = (TH3*)h3genNum->Clone("h3genEff");
    h3genEff->Divide(h3genDen);
    TH3* h3recoEff = (TH3*)h3recoNum->Clone("h3recoEff");
    h3recoEff->Divide(h3recoDen);

    TH3* h3eff = (TH3*)h3genEff->Clone("h3eff");
    h3eff->Multiply(h3recoEff);

    RooDataHist dataHist("dataHist","dataHist", ctKctLphi, Import(*h3eff,kTRUE));
    RooHistPdf histPdf(Form("pdf_ctKctLphi_q2bin%i",q2BinIndx),Form("pdf_ctKctLphi_q2bin%i",q2BinIndx), ctKctLphi, dataHist, 0);

    histPdf.Write();
    fout->Close();
  }

  fin->Close();
  
}

int main(int argc, char** argv)
{
  if (argc > 0)
  {
    unsigned int q2BinIndx=0;
    bool mist=0;
    //folder = "";
    if ( argc > 1 ) q2BinIndx = atoi(argv[1]);
    if ( argc > 2 ) mist = atoi(argv[2]);
    //if ( argc > 3 ) folder = argv[3];

    //if (folder=="") mist=0;

    TFile* fout = new TFile(fileNameOutput.c_str(), "RECREATE");
    fout->Close();

    // do all q2Bins at once
    if (q2BinIndx==0)
      for (q2BinIndx=1; q2BinIndx<10; ++q2BinIndx) {
        Convert(q2BinIndx, mist);
    }
    else if (q2BinIndx>0 && q2BinIndx<10) {
        Convert(q2BinIndx, mist);
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
