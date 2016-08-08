// Stefano Lacaprara  <lacaprara@pd.infn.it>  INFN Padova
// 24-Sep-2015
//
// Example code to read the efficiency as RooAbsPdf
//

void readPdfExample(int q2bin=1) {
  // Open the file in READ mode
  TFile* file=TFile::Open("effKEpdf_out_RT.root","READ");

  // This is the name of the pdf "pdf_ctKctLphi_q2bin1-9"
  TString pdfName=Form("pdf_ctKctLphi_q2bin%d",q2bin);

  // read the ojbect in memory as a RooAbsPdf*
  RooAbsPdf* pdf_ctKctLphi=(RooHistPdf*)file->Get(pdfName);

  // try to change variable names
  RooArgSet* set=pdf_ctKctLphi->getVariables();
  //set->Print();
  set->find("ctK")->SetName("z");
  set->find("ctL")->SetName("y");
  set->find("phi")->SetName("p");
  pdf_ctKctLphi->getVariables()->Print("v");
  return;


  // These are the RooVariables used in the pdf
  RooRealVar ctK("ctK","cos#theta_{K}",-1.,1.);
  RooRealVar ctL("ctL","cos#theta_{L}",-1.,1.);
  RooRealVar phi("phi","#phi",-1.*TMath::Pi(),1.*TMath::Pi());

  // create a canvas just to show some plot
  TCanvas* c1=new TCanvas(Form("canProjection%d",q2bin),"",800,800);
  c1->Divide(2);

  // Now you can do whatever you want, eg plot a projection
  RooPlot* frameCtK=ctK.frame();
  pdf_ctKctLphi->plotOn(frameCtK);
  c1->cd(1);
  frameCtK->Draw();


  // Or create a 2D histogram
  c1->cd(2);
  TH2* h2_pdf_ctKctL=pdf_ctKctLphi->createHistogram("ctK,ctL");
  h2_pdf_ctKctL->Draw("surf3 fp");

}
