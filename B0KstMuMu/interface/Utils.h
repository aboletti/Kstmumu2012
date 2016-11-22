#ifndef UTILS_H
#define UTILS_H

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TH3.h>
#include <TMatrixTSym.h>
#include <TGraphAsymmErrors.h>

#if ROOFIT
#include <TH3.h>
#include <TTree.h>
#include <RooArgSet.h>
#include <RooRealVar.h>
#include <RooAbsPdf.h>
#include <RooHistPdf.h>
#include <RooHistFunc.h>
#include <RooFitResult.h>
#endif

#include <string>
#include <vector>

#include "B0KstMuMuSingleCandTreeContent.h"
#include "B0KstMuMuTreeContent.h"


class Utils
{
 public:
  
  Utils(bool rightFlavorTag = true);
  ~Utils() {};

  // #################################
  // # Data structure for efficiency #
  // #################################
  struct _effStruct
  {
    // ###########################################################
    // # Numerators and Denominators = number of events * weight #
    // ###########################################################
    double *Num1, *Num2;
    double *Den1, *Den2;
    
    // ###################################################
    // # Poissonian errors = number of events * weight^2 #
    // ###################################################
    double *Err2PoisNum1, *Err2PoisNum2;
    double *Err2PoisDen1, *Err2PoisDen2;
    
    // #################
    // # Weight errors #
    // #################
    double *Err2WeigNum1, *Err2WeigNum2;
    double *Err2WeigDen1, *Err2WeigDen2;

    // #########################
    // # Correspondence :      #
    // # GenFilter     <--> N1 #
    // # SingleCand    <--> N2 #
    // # GenNoFilter   <--> D1 #
    // # AllCandFilter <--> D2 #
    // #########################
  };
  typedef struct _effStruct effStruct;

  struct _effValue
  {
    double Num1, Num2;
    double Den1, Den2;
    
    // Poissonian errors
    double Err2PoisNum1, Err2PoisNum2;
    double Err2PoisDen1, Err2PoisDen2;
    
    // Weight errors
    double Err2WeigNum1, Err2WeigNum2;
    double Err2WeigDen1, Err2WeigDen2;
  };
  typedef struct _effValue effValue;


  double computeInvMass (double Px1,
                         double Py1,
                         double Pz1,
                         double mass1,
                         double Px2,
                         double Py2,
                         double Pz2,
                         double mass2);
  
  double computeEta (double Px,
                     double Py,
                     double Pz);
  
  double computePhi (double Px, double Py);
  
  double computeEtaPhiDistance (double Px1,
				double Py1,
				double Pz1,
                                double Px2,
                                double Py2,
				double Pz2);
  
  void computeLS (double Vx,
		  double Vy,
		  double Vz,
		  double Wx,
		  double Wy,
		  double Wz,
		  double VxErr2,
		  double VyErr2,
		  double VzErr2,
		  double VxyCov,
		  double VxzCov,
		  double VyzCov,
		  double WxErr2,
		  double WyErr2,
		  double WzErr2,
		  double WxyCov,
		  double WxzCov,
		  double WyzCov,
		  double* deltaD,
		  double* deltaDErr);

  void computeCosAlpha (double Vx,
			double Vy,
			double Vz,
			double Wx,
			double Wy,
			double Wz,
			double VxErr2,
			double VyErr2,
			double VzErr2,
			double VxyCov,
			double VxzCov,
			double VyzCov,
			double WxErr2,
			double WyErr2,
			double WzErr2,
			double WxyCov,
			double WxzCov,
			double WyzCov,
			double* cosAlpha,
			double* cosAlphaErr);
#if ROOFIT
  RooAbsReal* ReadRTEffPDF (unsigned int q2BinIndx,RooRealVar* z, RooRealVar* y,RooRealVar* p);
  RooAbsReal* ReadWTEffPDF (unsigned int q2BinIndx,RooRealVar* z, RooRealVar* y,RooRealVar* p);
#endif

  void ReadAllBins  (std::string fileName, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, std::string signalType = "goodTag");
  void Readq2Bins   (std::string fileName, std::vector<double>* q2Bins);
  void ReadHLTpaths (std::string fileName, std::vector<std::string>* TrigTable);
  int SearchBin  (double val2Search, std::vector<double>* bins);
  int GetJPsiBin (std::vector<double>* q2Bins);
  int GetPsiPBin (std::vector<double>* q2Bins);

  bool ValIsInPsi (std::vector<double>* q2Bins, double q2Val);

  unsigned int HLTpathForEvFraction (double evtFrac);
  unsigned int IsInTriggerTable     (B0KstMuMuTreeContent* NTupleIn, double* HLTCutVar1, double* HLTCutVar2, int index = 0, double evtFrac = -1.0);
  unsigned int GetNHLTCat           ();
  double ReadLumi                   (std::string fileName);

  void ReadTriggerPathsANDCutsANDEntries (std::string fileName);
  void ReadFitStartingValues             (std::string fileName, std::vector<std::vector<std::string>*>* vecParam, std::vector<std::vector<unsigned int>*>* configParam, const unsigned int dataBlockN);
  void ReadFitSystematics                (std::string fileName, std::vector<std::vector<double>*>* vecParam);
  void ReadParVsq2Bins                   (std::string fileName, std::string praName, std::vector<std::string>** vecParam);

  void MakeGraphVar (std::string fileName, TGraphAsymmErrors** graph, std::string varName, double offset = 0.0);

  void InitEffFuncThetaL (TF1* fitFun);
  void InitEffFuncThetaK (TF1* fitFun);
  void InitEffFuncPhi    (TF1* fitFun);

  void AddConstraint1D              (TH1D** histo, std::string constrType, double abscissaErr, double YerrRescale, double Yval, double Yerr, unsigned int ID);
  void AddConstraintThetaL          (TH1D** histo, unsigned int q2Indx, unsigned int cosThetaKBinIndx, unsigned int ID);
  void AddConstraint2D              (TH2D** histo, double abscissaErr, double ZerrRescale, unsigned int ID, std::string toBeConstr, double scaleConstr, double constrXerr, std::vector< std::pair <double,double> >* constraints = NULL, std::vector<std::string>* toBeAdded = NULL);
  void AddConstraintThetaKThetaL    (TH2D** histo, std::vector<double>* cosThetaKBins, unsigned int q2Indx, int SignalType, unsigned int ID);
  void AddConstraint3D              (TH3D** histo, double abscissaErr, double Tval, double Terr, double TerrRescale, unsigned int ID, std::vector<int> toBeAdded[]);
  void AddConstraintThetaKThetaLPhi (int SignalType);

  int WhatIsThis (std::string fileName);

  void SaveFitValues (std::string fileName, std::vector<std::string>* vecParStr, int indx, std::string howOpen, std::string str = "");

#if ROOFIT
  std::string Transformer (std::string varName, bool doIt, double& varValOut, double& varValOutELo, double& varValOutEHi, RooFitResult* fitResult = NULL, RooRealVar* varValIn1 = NULL, RooRealVar* varValIn2 = NULL, RooRealVar* varValIn3 = NULL, RooRealVar* varValIn4 = NULL);
  void AntiTransformer (std::string varName, bool doIt, double& varValOut, double& varValOutELo, double& varValOutEHi, RooRealVar* varValIn1 = NULL, RooRealVar* varValIn2 = NULL);
#endif
  
  std::string MakeAnalysisPATH (std::string relativePath);
  unsigned int ParFileBlockN   (std::string blockName);

  unsigned int GetFitParamIndx    (std::string varName);
  unsigned int GetConfigParamIndx (std::string varName);

  bool PsiRejection   (double myB0Mass, double myMuMuMass, double myMuMuMassE, std::string seleType, bool B0andPsiCut = false);
  bool ChooseBestCand (B0KstMuMuTreeContent* NTuple, unsigned int DoTrigCheck, double evFraction, int* BestCandIndx, bool* B0notB0bar, int* TrigCat, unsigned int* countCands);
  bool FlavorTagger   (B0KstMuMuTreeContent* NTuple, unsigned int i, bool* B0notB0bar);

  void   ReadSelectionCuts (std::string fileName);
  bool   SetSeleCut        (std::string cutName, double val);
  double GetSeleCut        (std::string cutName);

  void   ReadPreselectionCut (std::string fileName);
  bool   SetPreCut           (std::string cutName, double val);
  double GetPreCut           (std::string cutName);

  void   ReadGenericParam     (std::string fileName);
  bool   SetGenericParam      (std::string parName, std::string val);
  std::string GetGenericParam (std::string parName);

  double GetB0Width ();

  double* MakeBinning (std::vector<double>* STLvec);



  double muonMass;
  double pionMass;
  double kaonMass;
  double kstMass;
  double B0Mass;
  double JPsiMass;
  double PsiPMass;

  double JPsiBF;
  double JPsiKpiBF;
  double KstMuMuBF;
  double KstKpiMuMuBF;
  double PsiPBF;
  double PsiPKpiBF;

  double muonMassErr;
  double pionMassErr;
  double kaonMassErr;
  double B0MassErr;
  double kstSigma;

  double PI;

  unsigned int NcoeffThetaL;
  unsigned int NcoeffThetaK;
  unsigned int NcoeffPhi;

  bool RIGHTflavorTAG;

  int B0ToKstMuMu;
  int B0ToJPsiKst;
  int B0ToPsi2SKst;

  unsigned int nFitParam;
  unsigned int nConfigParam;


 private:

  TF1* KstMassShape;

  std::vector<std::string> HLTpath;
  std::vector<double> VecHLTCutVar1;
  std::vector<double> VecHLTCutVar2;
  std::vector<double> VecHLTentries;
  std::vector<std::string> TrigTable;

  std::vector<double> PreCuts;
  std::vector<double> SeleCuts;
  std::vector<std::string> GenericPars;

  double ProbThreshold;
  double scrambleFraction;

  unsigned int nFitObserv;

  std::string DirEfficiency;

  std::string Histo2DEffNameOkTagSig;
  std::string Histo2DEffNameOkTagJPsi;
  std::string Histo2DEffNameOkTagPsi2S;

  std::string Histo3DEffNameOkTagSig;
  std::string Histo3DEffNameOkTagJPsi;
  std::string Histo3DEffNameOkTagPsi2S;

  std::string Histo2DEffNameMisTagSig;
  std::string Histo2DEffNameMisTagJPsi;
  std::string Histo2DEffNameMisTagPsi2S;

  std::string Histo3DEffNameMisTagSig;
  std::string Histo3DEffNameMisTagJPsi;
  std::string Histo3DEffNameMisTagPsi2S;
};

#endif
