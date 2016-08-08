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

//const string folder;
TFile* fin;
TFile* fin1;

string fileNameOutput;
string sample_str;

RooRealVar ctK("ctK","cos#theta_{K}",-1,1);
RooRealVar ctL("ctL","cos#theta_{L}",0,1);
RooRealVar phi("phi","#phi",0,TMath::Pi());
RooArgList ctKctLphi(ctK,ctL,phi);

void Convert (int q2BinIndx, string exe, bool mist=0) {

  sample_str = "";
  if (q2BinIndx==5) sample_str = "_jp";
  if (q2BinIndx==7) sample_str = "_p2";
  fileNameOutput = ( mist? "Closure_sys2_wt" : "Closure_sys2_rt" ) + sample_str + ( mist? "_1.0/" : "_0.5/" ) + exe + "/effKEpdf_out.root";

  fin = TFile::Open( ("Closure_sys2_rt" + sample_str + "_0.5/" + exe + Form("/effKEpdf_b%i.root",q2BinIndx)).c_str() );
  if (mist) fin1 = TFile::Open( ("Closure_sys2_wt" + sample_str + "_1.0/" + exe + Form("/effKEpdf_b%i.root",q2BinIndx)).c_str() );
  else fin1 = fin;

  if (fin && fin1) {

    string hname = "h3recoNum";
    if (mist && atoi(exe.c_str())!=0) hname = "__ctK_ctL_phi";
    TH3* h3genNum  = (TH3*)fin->Get("h3genNum" );
    TH3* h3genDen  = (TH3*)fin->Get("h3genDen" );
    TH3* h3recoDen = (TH3*)fin->Get("h3recoDen");
    TH3* h3recoNum = (TH3*)fin1->Get(hname.c_str());

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
    } else cout<<"Error: histo not found in bin "<<q2BinIndx<<" of run "<<exe<<endl;

    fin->Close();
    if (mist) fin1->Close();
  
  } else cout<<"Error: file not found in bin "<<q2BinIndx<<" of run "<<exe<<endl;

}

int main(int argc, char** argv)
{
  if (argc > 1)
  {
    unsigned int q2BinIndx=0;
    bool mist=0;
    string exe = "00000";
    if ( argc > 1 ) exe = argv[1];
    if ( argc > 2 ) q2BinIndx = atoi(argv[2]);
    if ( argc > 3 ) mist = atoi(argv[3]);

    if ( q2BinIndx == 0 || ( q2BinIndx != 5 && q2BinIndx != 7 ) ) {
      fileNameOutput = ( mist? "Closure_sys2_wt_1.0/" : "Closure_sys2_rt_0.5/" ) + exe + "/effKEpdf_out.root";
      TFile* fout = new TFile(fileNameOutput.c_str(), "RECREATE");
      fout->Close();
    }
    if ( q2BinIndx == 0 || q2BinIndx == 5 ) {
      fileNameOutput = ( mist? "Closure_sys2_wt_1.0/" : "Closure_sys2_rt_jp_0.5/" ) + exe + "/effKEpdf_out.root";
      TFile* fout = new TFile(fileNameOutput.c_str(), "RECREATE");
      fout->Close();
    }
    if ( q2BinIndx == 0 || q2BinIndx == 7 ) {
      fileNameOutput = ( mist? "Closure_sys2_wt_1.0/" : "Closure_sys2_rt_p2_0.5/" ) + exe + "/effKEpdf_out.root";
      TFile* fout = new TFile(fileNameOutput.c_str(), "RECREATE");
      fout->Close();
    }

    // do all q2Bins at once
    if (q2BinIndx==0)
      for (q2BinIndx=1; q2BinIndx<10; ++q2BinIndx) {
        Convert(q2BinIndx, exe, mist);
    }
    else if (q2BinIndx>0 && q2BinIndx<10) {
        Convert(q2BinIndx, exe, mist);
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
