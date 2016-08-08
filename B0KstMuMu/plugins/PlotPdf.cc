#include "Utils.h"

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

#define PARAMETERFILEIN "/python/ParameterFile.txt"
Utils* Utility;

B0KstMuMuSingleCandTreeContent* NTupleIn_reco;
B0KstMuMuSingleCandTreeContent* NTupleIn_reco1;
B0KstMuMuSingleCandTreeContent* NTupleIn_gen;

vector<double> q2Bins;

void loadBinning() {
  vector<double> cosThetaKBins;
  vector<double> cosThetaLBins;
  vector<double> phiBins;
  Utility->ReadAllBins(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str(),&q2Bins,&cosThetaKBins,&cosThetaLBins,&phiBins,"goodTag");
}

int getNbin ( float q2 ) {
  for (int i=0; i<10; i++) if ( q2 < q2Bins[i] ) return i;
  return 0;
}

void Fill (int signalType) {

  TString signalString = "";
  if (signalType==1) signalString = "_jp";
  if (signalType==2) signalString = "_p2";

  TString name_cat = "B0ToKstMuMu";
  if (signalType==1) name_cat = "B0ToJPsiKst";
  if (signalType==2) name_cat = "B0ToPsi2SKst";

  TFile* fout = new TFile ( "histo"+signalString+".root", "RECREATE" );

  TH3D* h3genNum_rt[9];
  TH3D* h3genDen_rt[9];
  TH3D* h3recoDen_rt[9];
  TH3D* h3recoNum_rt[9];
  TH3D* h3recoNum_wt[9];

  for (int i=0; i<9; i++) {
    h3genNum_rt [i] = new TH3D(Form("h3genNum_rt%i" ,i+1),Form("h3genNum_rt_bin%i" ,i+1),40,-1,1,40,0,1,40,0,TMath::Pi());
    h3genDen_rt [i] = new TH3D(Form("h3genDen_rt%i" ,i+1),Form("h3genDen_rt_bin%i" ,i+1),40,-1,1,40,0,1,40,0,TMath::Pi());
    h3recoDen_rt[i] = new TH3D(Form("h3recoNum_rt%i",i+1),Form("h3recoNum_rt_bin%i",i+1),40,-1,1,40,0,1,40,0,TMath::Pi());
    h3recoNum_rt[i] = new TH3D(Form("h3recoDen_rt%i",i+1),Form("h3recoDen_rt_bin%i",i+1),40,-1,1,40,0,1,40,0,TMath::Pi());
    h3recoNum_wt[i] = new TH3D(Form("h3recoDen_wt%i",i+1),Form("h3recoDen_wt_bin%i",i+1),40,-1,1,40,0,1,40,0,TMath::Pi());
  }

  int nbin;

  TFile* NtplFileIn_reco = TFile::Open("dcap://t2-srm-02.lnl.infn.it//pnfs/lnl.infn.it/data/cms/store/user/slacapra/B0KsMuMu/Data2012B0KstMuMuResults/MonteCarlo2012/SingleCand/singleCand_"+name_cat+"_MC_NTuple.root", "READ");
  TTree* theTreeIn_reco  = (TTree*) NtplFileIn_reco->Get("B0KstMuMu/B0KstMuMuNTuple");
  NTupleIn_reco = new B0KstMuMuSingleCandTreeContent();
  NTupleIn_reco->Init();
  NTupleIn_reco->ClearNTuple();
  NTupleIn_reco->SetBranchAddresses(theTreeIn_reco);
  int nEntries_reco = theTreeIn_reco->GetEntries();
  cout << "\n[ClosureTest::performTest]\t@@@ Total number of selected reco events in the tree: " << nEntries_reco << " @@@" << endl;
  for (int entry = 0; entry < nEntries_reco; entry++) {
    theTreeIn_reco->GetEntry(entry);
    nbin = getNbin( TMath::Power(NTupleIn_reco->mumuMass->at(0),2) );
    if ( (nbin==0) || (NTupleIn_reco->B0pT < Utility->GetSeleCut("B0pT")) || (fabs(NTupleIn_reco->B0Eta) > Utility->GetSeleCut("B0Eta")) ) continue;
    if ( (NTupleIn_reco->truthMatchSignal->at(0) == true) &&
	 (NTupleIn_reco->B0MassArb > Utility->B0Mass - atof(Utility->GetGenericParam("B0MassIntervalLeft").c_str()))            &&
	 (NTupleIn_reco->B0MassArb < Utility->B0Mass + atof(Utility->GetGenericParam("B0MassIntervalRight").c_str()))           &&
	 ( ((signalType == 0) && (Utility->PsiRejection(NTupleIn_reco->B0MassArb,NTupleIn_reco->mumuMass->at(0),NTupleIn_reco->mumuMassE->at(0),"rejectPsi",true) == true)) ||
	   ((signalType == 1) && (Utility->PsiRejection(NTupleIn_reco->B0MassArb,NTupleIn_reco->mumuMass->at(0),NTupleIn_reco->mumuMassE->at(0),"keepJpsi")       == true)) ||
	   ((signalType == 2) && (Utility->PsiRejection(NTupleIn_reco->B0MassArb,NTupleIn_reco->mumuMass->at(0),NTupleIn_reco->mumuMassE->at(0),"keepPsiP")       == true))   ) ) {
      if ( NTupleIn_reco->rightFlavorTag ) h3recoNum_rt[nbin-1]->Fill( NTupleIn_reco->CosThetaKArb, TMath::Abs(NTupleIn_reco->CosThetaMuArb), TMath::Abs(NTupleIn_reco->PhiKstMuMuPlaneArb) );
      else h3recoNum_wt[nbin-1]->Fill( -1*NTupleIn_reco->CosThetaKArb, TMath::Abs(NTupleIn_reco->CosThetaMuArb), TMath::Abs(NTupleIn_reco->PhiKstMuMuPlaneArb) );
    }
  }
  
  TFile* NtplFileIn_reco1 = TFile::Open("dcap://t2-srm-02.lnl.infn.it////pnfs/lnl.infn.it/data/cms/store/user/slacapra/B0KsMuMu/Data2012B0KstMuMuResults/MonteCarlo2012/ForEfficiency/"+name_cat+"_MC_NTuple_addGENvars.root", "READ");
  TTree* theTreeIn_reco1  = (TTree*) NtplFileIn_reco1->Get("B0KstMuMu/B0KstMuMuNTuple");
  NTupleIn_reco1 = new B0KstMuMuSingleCandTreeContent();
  NTupleIn_reco1->Init();
  NTupleIn_reco1->ClearNTuple();
  NTupleIn_reco1->SetBranchAddresses(theTreeIn_reco1);
  int nEntries_reco1 = theTreeIn_reco1->GetEntries();
  cout << "\n[ClosureTest::performTest]\t@@@ Total number of reco events in the tree: " << nEntries_reco1 << " @@@" << endl;
  for (int entry = 0; entry < nEntries_reco1; entry++) {
    theTreeIn_reco1->GetEntry(entry);
    nbin = getNbin( TMath::Power(NTupleIn_reco1->mumuMass->at(0),2) );
    if ( (nbin==0) || (NTupleIn_reco1->B0pT < Utility->GetSeleCut("B0pT")) || (fabs(NTupleIn_reco1->B0Eta) > Utility->GetSeleCut("B0Eta")) ) continue;
    h3recoDen_rt[nbin-1]->Fill( NTupleIn_reco1->CosThetaKArb, TMath::Abs(NTupleIn_reco1->CosThetaMuArb), TMath::Abs(NTupleIn_reco1->PhiKstMuMuPlaneArb) );
  }
  
  TFile* NtplFileIn_gen = TFile::Open("dcap://t2-srm-02.lnl.infn.it//pnfs/lnl.infn.it/data/cms/store/user/slacapra/B0KsMuMu/Data2012B0KstMuMuResults/MonteCarlo2012/ForEfficiency/"+name_cat+"_GEN_NoFilter_MC_NTuples_addGENvars.root","READ");
  TTree* theTreeIn_gen  = (TTree*) NtplFileIn_gen->Get("B0KstMuMu/B0KstMuMuNTuple");
  NTupleIn_gen = new B0KstMuMuSingleCandTreeContent();
  NTupleIn_gen->Init();
  NTupleIn_gen->ClearNTuple();
  NTupleIn_gen->SetBranchAddresses(theTreeIn_gen);
  int nEntries_gen = theTreeIn_gen->GetEntries();
  cout << "\n[ClosureTest::performTest]\t@@@ Total number of gen events in the tree: " << nEntries_gen << " @@@" << endl;
  for (int entry = 0; entry < nEntries_gen; entry++) {
    theTreeIn_gen->GetEntry(entry);
    nbin = getNbin( TMath::Power(NTupleIn_gen->mumuMass->at(0),2) );
    if ( (nbin==0) || (NTupleIn_gen->B0pT < Utility->GetSeleCut("B0pT")) || (fabs(NTupleIn_gen->B0Eta) > Utility->GetSeleCut("B0Eta")) ) continue;
    h3genDen_rt[nbin-1]->Fill( NTupleIn_gen->CosThetaKArb, TMath::Abs(NTupleIn_gen->CosThetaMuArb), TMath::Abs(NTupleIn_gen->PhiKstMuMuPlaneArb) );
    if ( (sqrt(NTupleIn_gen->genMumPx*NTupleIn_gen->genMumPx + NTupleIn_gen->genMumPy*NTupleIn_gen->genMumPy)   > Utility->GetPreCut("MinMupT")) &&
	 (sqrt(NTupleIn_gen->genMupPx*NTupleIn_gen->genMupPx + NTupleIn_gen->genMupPy*NTupleIn_gen->genMupPy)   > Utility->GetPreCut("MinMupT")) &&
	 (fabs(Utility->computeEta(NTupleIn_gen->genMumPx,NTupleIn_gen->genMumPy,NTupleIn_gen->genMumPz)) < Utility->GetPreCut("MuEta"))         &&
	 (fabs(Utility->computeEta(NTupleIn_gen->genMupPx,NTupleIn_gen->genMupPy,NTupleIn_gen->genMupPz)) < Utility->GetPreCut("MuEta"))            )
      h3genNum_rt[nbin-1]->Fill( NTupleIn_gen->CosThetaKArb, TMath::Abs(NTupleIn_gen->CosThetaMuArb), TMath::Abs(NTupleIn_gen->PhiKstMuMuPlaneArb) );
  }

  fout->cd();
  for (int i=0; i<9; i++) {
    h3genNum_rt [i]->Write(Form("h3genNum_rt%i" ,i));
    h3genDen_rt [i]->Write(Form("h3genDen_rt%i" ,i));
    h3recoDen_rt[i]->Write(Form("h3recoNum_rt%i",i));
    h3recoNum_rt[i]->Write(Form("h3recoDen_rt%i",i));
    h3recoNum_wt[i]->Write(Form("h3recoDen_wt%i",i));
  }
  fout->Close();
}

int main(int argc, char** argv)
{
  if (argc > 0)
  {
    int signalType=0;
    if ( argc > 1 ) signalType = atoi(argv[1]);

    Utility = new Utils(false);
    Utility->ReadPreselectionCut(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str());
    Utility->ReadSelectionCuts(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str());
    Utility->ReadGenericParam(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str());
    loadBinning();

    Fill(signalType);

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
