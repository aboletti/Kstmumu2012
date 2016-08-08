using namespace RooFit;
RooRealVar ctK("ctK","cos#theta_{K}",-1.,1.);
RooRealVar ctL("ctL","cos#theta_{L}",0.,1.);
RooRealVar phi("phi","#phi",0.,1.*TMath::Pi());
TFile* file=TFile::Open("effKEpdf_out.root","READ");
//TFile* file=TFile::Open("Closure_wt_1.0/effKEpdf_out.root","READ");
//TFile* file=TFile::Open("../plugins/Efficiency/effKEpdf_out_RT.root","READ");
//TFile* file=TFile::Open("Closure_rt_0.5/effKEpdf_out.root","READ");

void plotKEPdf(unsigned int q2bin) {
  TString pdfName=Form("pdf_ctKctLphi_q2bin%d",q2bin);
  TCanvas* c1=new TCanvas(Form("canProjection%d",q2bin),"",800,800);
  RooAbsPdf* term0Pdf = (RooHistPdf*)file->Get(pdfName);
  TH2* h2_pdf_ctKctL=term0Pdf->createHistogram("ctK,ctL");
  h2_pdf_ctKctL->Draw("surf3 fp");

}

void plotPdf(unsigned int q2bin) {
  gStyle->SetOptStat(0);
  gStyle->SetTitleOffset(1.1,"XYZ");
  gROOT->ForceStyle();

  RooAbsPdf* pdf_ctKctLphi=(RooHistPdf*)file->Get(Form("pdf_ctKctLphi_q2bin%d",q2bin));

  TH3* h3dfull = (TH3*)pdf_ctKctLphi->createHistogram("hKE",ctK,Binning(40),YVar(ctL,Binning(40)),ZVar(phi,Binning(40)));

  TString pro1D[3]={"x","y","z"};
  TString pro2D[3]={"xy","xz","yz"};

  TCanvas* c1=new TCanvas(Form("canProjection%d",q2bin),Form("canProjection%d",q2bin),800,800);
  c1->Divide(3,2);
  int ic=1;
  for (int ipro=0; ipro<3; ++ipro) {
    c1->cd(ic);
    TH2* h2=h3dfull->Project3D(pro2D[ipro]);
    h2->SetMinimum(0.);
    h2->Draw("surf3 fp");
    // c1->cd(++ic);
    // TH2* h2_pdf_ctKphi=pdf_ctKctLphi->createHistogram("ctK,phi");
    // h2_pdf_ctKphi->Draw("surf3 fp");
    // c1->cd(3);
    // TH2* h2_pdf_ctLphi=pdf_ctKctLphi->createHistogram("ctL,phi");
    // h2_pdf_ctLphi->Draw("surf3 fp");
    // // c1->cd(4);
    // // TH2* h2_pdf_ctKctLphi=pdf_ctKctLphi->createHistogram("ctK,ctL,phi");
    // // h2_pdf_ctKctLphi->Draw("box2");

    c1->cd(3+ic++);
    TH1* h1=h3dfull->Project3D(pro1D[ipro]);
    h1->SetMinimum(0.);
    h1->Draw("L");
    // c1->cd(4);
    // RooPlot* frameCtK=ctK.frame();
    // pdf_ctKctLphi->plotOn(frameCtK);
    // frameCtK->Draw("");
    // c1->cd(5);
    // RooPlot* frameCtL=ctL.frame();
    // pdf_ctKctLphi->plotOn(frameCtL);
    // frameCtL->Draw("");
    // c1->cd(6);
    // RooPlot* framePhi=phi.frame();
    // pdf_ctKctLphi->plotOn(framePhi);
    // framePhi->Draw("");
  }
  c1->Print(Form("canProjection%d_xInverted.pdf",q2bin));
  //c1->Print(Form("canProjection%d.pdf",q2bin));
}

void plotPdf() {
  for (unsigned int q2bin=1; q2bin<=10; ++q2bin) plotPdf(q2bin);
}
