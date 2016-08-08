void plotBin(int bin=1) {
  gStyle->SetTitleOffset(1.1,"XYZ");
  gROOT->ForceStyle();
  TFile* finOrig=TFile::Open("effHisto3DFile.root","READ");
  TH3* hEff3Orig = (TH3*) finOrig->Get(Form("hisFunc3D_%i",bin));

  //TFile* fin=TFile::Open(Form("bin_%i/effKEpdf.root",bin),"READ");
  TFile* fin=TFile::Open("effKEpdf.root","READ");
  TCanvas* c1=new TCanvas(Form("cEff1D_q2bin%i",bin),Form("Efficiency q2bin=%i",bin),1200,1200);
  TCanvas* c2=new TCanvas(Form("cEff2D_q2bin%i",bin),Form("Efficiency q2bin=%i",bin),1200,1200);
  c1->Divide(3,3);
  c2->Divide(3,3);
  TH3* hGenEff3 = (TH3*) fin->Get(Form("h3genEff_q2bin%i",bin));
  TH3* hRecoEff3 = (TH3*) fin->Get(Form("h3recoEff_q2bin%i",bin));
  TH3* hEff3 = hGenEff3->Clone("hGenEff3");
  hEff3->Multiply(hRecoEff3);
  //TH3* hEff3 = (TH3*) fin->Get(Form("h3Eff_q2bin%i",bin));
  TString pro1D[3]={"x","y","z"};
  TString pro2D[3]={"xy","xz","yz"};

  TLatex tl;
  tl.SetNDC();
  tl.SetTextSize(0.07);
  int ic=1;
  for (int ipro=0; ipro<3; ++ipro) {
    c1->cd(ic);
    hGenEff3->Project3D(pro1D[ipro])->SetMinimum(0.);
    hGenEff3->Project3D(pro1D[ipro])->Draw("L");
    tl.DrawLatex(0.2,0.85,"Gen #epsilon");
    c2->cd(ic++);
    hGenEff3->Project3D(pro2D[ipro])->Draw("surf");
    tl.DrawLatex(0.1,0.9,"Gen #epsilon");
    c1->cd(ic);
    hRecoEff3->Project3D(pro1D[ipro])->SetMinimum(0.0);
    hRecoEff3->Project3D(pro1D[ipro])->Draw("L");
    tl.DrawLatex(0.2,0.85,"Reco #epsilon");
    c2->cd(ic++);
    hRecoEff3->Project3D(pro2D[ipro])->Draw("surf");
    tl.DrawLatex(0.1,0.9,"Reco #epsilon");
    c1->cd(ic);
    hEff3Orig->Project3D(pro1D[ipro])->SetMinimum(0.);
    hEff3Orig->Project3D(pro1D[ipro])->Draw("hist");
    hEff3->Project3D(pro1D[ipro])->DrawNormalized("L same",3.2);
    tl.DrawLatex(0.2,0.85,"Tot #epsilon");
    c2->cd(ic++);
    hEff3->Project3D(pro2D[ipro])->Draw("surf ");
    hEff3Orig->Project3D(pro2D[ipro])->Draw("lego same");
    tl.DrawLatex(0.1,0.9,"Tot #epsilon");
  }
  c1->Print((c1->GetName())+TString(".pdf"));
  c2->Print((c2->GetName())+TString(".pdf"));


}
