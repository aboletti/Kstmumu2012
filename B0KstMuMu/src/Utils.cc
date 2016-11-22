#include "../interface/Utils.h"
#include "../interface/ReadParameters.h"

#include <TAxis.h>
#include <TMath.h>
#include <TFile.h>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <fstream>

#include <RooFormulaVar.h>

// ####################
// # Global constants #
// ####################
#define YvalueOutsideLimits 20.0 // Value given to bins with zero error in order not to show them
#define ANALYPATH "ANALYPATH"

Utils::Utils (bool rightFlavorTag)
{
  muonMass     = 0.10565837;
  pionMass     = 0.13957018;
  kaonMass     = 0.493677;
  kstMass      = 0.896;
  B0Mass       = 5.27958;
  JPsiMass     = 3.096916;
  PsiPMass     = 3.686109;

  JPsiBF       =  7.87e-5; // B0 --> J/psi(mu+mu-) K*0          (1.32+/-0.06 * 5.96+/-0.03)
  JPsiKpiBF    =  5.25e-5; // B0 --> J/psi(mu+mu-) K*0(K+pi-)   (1.32+/-0.06 * 5.96+/-0.03 * 2/3)

  KstMuMuBF    =  1.05e-6; // B0 --> K*0 mu+mu-
  KstKpiMuMuBF =  7e-7;    // B0 --> K*0(K+pi-) mu+mu-          (1.05+/-0.1 * 2/3)

  PsiPBF       = 47.72e-7; // B0 --> psi(2S)(mu+mu-) K*0        (6.04+/-0.4 * 7.9+/-0.9)
  PsiPKpiBF    = 31.81e-7; // B0 --> psi(2S)(mu+mu-) K*0(K+pi-) (6.04+/-0.4 * 7.9+/-0.9 * 2/3)

  muonMassErr  = 3.5e-9;
  pionMassErr  = 3.5e-7;
  kaonMassErr  = 1.6e-5;
  B0MassErr    = 1.7e-4;
  kstSigma     = 0.05;

  nFitParam    = 67;
  nConfigParam = 4;
  nFitObserv   = 3; // FL --- P5p --- P1 --- P2 --- BF

  NcoeffThetaL = 6;
  NcoeffThetaK = 4;
  NcoeffPhi    = 4;

  PI = 3.141592653589793;

  ProbThreshold = 0.0; // Confidence Level for accepting the null hypothesis: "the two mass hypothesis are statistically indistinguishable"
  // if (F-test < ProbThreshold) --> accept the null hypothesis
  // if (F-test > ProbThreshold) --> reject the null hypothesis
  scrambleFraction = 0.0; // Fraction of events with random CP-tagging
  KstMassShape = new TF1("KstMassShape",
                         "2*sqrt(2)*[0]*[1]* sqrt([0]*[0] * ([0]*[0] + [1]*[1])) / (TMath::Pi()*sqrt([0]*[0] + sqrt([0]*[0] * ([0]*[0] + [1]*[1])))) / ((x*x - [0]*[0]) * (x*x - [0]*[0]) + [0]*[0]*[1]*[1])",
                         0.0,kstMass*2.);
  // Breit-Wigner distribution:
  // [0]: mass of the resonance
  // [1]: width of the resonance
  KstMassShape->SetParameter(0,kstMass);
  KstMassShape->SetParameter(1,kstSigma);
  KstMassShape->SetParNames("Mean","Width");

  // Define whether to compute the efficiency with good-tagged or mis-tagged events
  RIGHTflavorTAG = rightFlavorTag;

  // Define names of the files containing the histogram of the efficiency
  DirEfficiency  = "";
  //DirEfficiency  = "efficiency/";

  // ###############################################
  // # ===> Define codes to identify MC type <===  #
  // #                                             #
  // # 1 = B0    --> K*0(K+pi-) mu+mu-             #
  // # 2 = B0bar --> K*0bar(K-pi+) mu+mu-          #
  // #                                             #
  // # 3 = B0    --> K*0(K+pi-) J/psi(mu+mu-)      #
  // # 4 = B0bar --> K*0bar(K-pi+) J/psi(mu+mu-)   #
  // #                                             #
  // # 5 = B0    --> K*0(K-pi+) psi(2S)(mu+mu-)    #
  // # 6 = B0bar --> K*0bar(K-pi+) psi(2S)(mu+mu-) #
  // ###############################################
  B0ToKstMuMu  = 1;
  B0ToJPsiKst  = 3;
  B0ToPsi2SKst = 5;

  // ################################
  // # Print out internal variables #
  // ################################
  // std::cout << "\n@@@@@@ Utils class settings : private @@@@@@" << std::endl;
  // std::cout << "Analysis environment variable: " << ANALYPATH << std::endl;
  // std::cout << "nFitObserv: "                << nFitObserv << std::endl;
  // std::cout << "ProbThreshold: "             << ProbThreshold << std::endl;
  // std::cout << "scrambleFraction: "          << scrambleFraction << std::endl;
  // std::cout << "DirEfficiency: "             << DirEfficiency << std::endl;

  // std::cout << "\n@@@@@@ Utils class settings : public  @@@@@@" << std::endl;
  // std::cout << "NcoeffThetaL: "              << NcoeffThetaL << std::endl;
  // std::cout << "NcoeffThetaK: "              << NcoeffThetaK << std::endl;
  // std::cout << "NcoeffPhi: "                 << NcoeffPhi << std::endl;
  // std::cout << "RIGHTflavorTAG: "            << RIGHTflavorTAG << std::endl;
  // std::cout << "B0ToKstMuMu: "               << B0ToKstMuMu << std::endl;
  // std::cout << "B0ToJPsiKst: "               << B0ToJPsiKst << std::endl;
  // std::cout << "B0ToPsi2SKst: "              << B0ToPsi2SKst << std::endl;
  // std::cout << "nFitParam: "                 << nFitParam << std::endl;
  // std::cout << "nConfigParam: "              << nConfigParam << std::endl;
  // std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
  // std::cout << "@@@ Consider to double-check values for: @@@" << std::endl;
  // std::cout << "- Utils::AddConstraintThetaL" << std::endl;
  // std::cout << "- Utils::AddConstraintThetaKThetaL" << std::endl;
  // std::cout << "- Utils::AddConstraintThetaKThetaLPhi" << std::endl;
  // std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
}

double Utils::computeInvMass (double Px1,
                              double Py1,
                              double Pz1,
                              double mass1,
                              double Px2,
                              double Py2,
                              double Pz2,
                              double mass2)
{
  double Energy1 = sqrt(Px1*Px1 + Py1*Py1 + Pz1*Pz1 + mass1*mass1);
  double Energy2 = sqrt(Px2*Px2 + Py2*Py2 + Pz2*Pz2 + mass2*mass2);
  return sqrt((Energy1+Energy2) * (Energy1+Energy2) - ((Px1+Px2) * (Px1+Px2) + (Py1+Py2) * (Py1+Py2) + (Pz1+Pz2) * (Pz1+Pz2)));
}

double Utils::computeEta (double Px,
                          double Py,
                          double Pz)
{
  double P = sqrt(Px*Px + Py*Py + Pz*Pz);
  return 0.5*log((P + Pz) / (P - Pz));
}

double Utils::computePhi (double Px, double Py)
{
  double phi = atan(Py / Px);
  if (Px < 0 && Py < 0) phi = phi - PI;
  if (Px < 0 && Py > 0) phi = phi + PI;
  return phi;
}

double Utils::computeEtaPhiDistance (double Px1,
                                     double Py1,
                                     double Pz1,
                                     double Px2,
                                     double Py2,
                                     double Pz2)
{
  double phi1 = computePhi (Px1,Py1);
  double eta1 = computeEta (Px1,Py1,Pz1);
  double phi2 = computePhi (Px2,Py2);
  double eta2 = computeEta (Px2,Py2,Pz2);
  return sqrt((eta1-eta2) * (eta1-eta2) + (phi1-phi2) * (phi1-phi2));
}

void Utils::computeLS (double Vx,
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
                       double* deltaDErr)
{
  *deltaD = sqrt((Vx-Wx) * (Vx-Wx) + (Vy-Wy) * (Vy-Wy) + (Vz-Wz) * (Vz-Wz));
  if (*deltaD > 0.)
    *deltaDErr = sqrt((Vx-Wx) * (Vx-Wx) * VxErr2 +
                      (Vy-Wy) * (Vy-Wy) * VyErr2 +
                      (Vz-Wz) * (Vz-Wz) * VzErr2 +

                      (Vx-Wx) * (Vy-Wy) * 2.*VxyCov +
                      (Vx-Wx) * (Vz-Wz) * 2.*VxzCov +
                      (Vy-Wy) * (Vz-Wz) * 2.*VyzCov +

                      (Vx-Wx) * (Vx-Wx) * WxErr2 +
                      (Vy-Wy) * (Vy-Wy) * WyErr2 +
                      (Vz-Wz) * (Vz-Wz) * WzErr2 +

                      (Vx-Wx) * (Vy-Wy) * 2.*WxyCov +
                      (Vx-Wx) * (Vz-Wz) * 2.*WxzCov +
                      (Vy-Wy) * (Vz-Wz) * 2.*WyzCov) / *deltaD;
  else *deltaDErr = 0.;
}

void Utils::computeCosAlpha (double Vx,
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
                             double* cosAlphaErr)
{
  double Vnorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
  double Wnorm = sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
  double VdotW = Vx*Wx + Vy*Wy + Vz*Wz;

  if ((Vnorm > 0.) && (Wnorm > 0.))
  {
    *cosAlpha = VdotW / (Vnorm * Wnorm);
    *cosAlphaErr = sqrt( ((Vx*Wnorm - VdotW*Wx) * (Vx*Wnorm - VdotW*Wx) * WxErr2 +
                          (Vy*Wnorm - VdotW*Wy) * (Vy*Wnorm - VdotW*Wy) * WyErr2 +
                          (Vz*Wnorm - VdotW*Wz) * (Vz*Wnorm - VdotW*Wz) * WzErr2 +

                          (Vx*Wnorm - VdotW*Wx) * (Vy*Wnorm - VdotW*Wy) * 2.*WxyCov +
                          (Vx*Wnorm - VdotW*Wx) * (Vz*Wnorm - VdotW*Wz) * 2.*WxzCov +
                          (Vy*Wnorm - VdotW*Wy) * (Vz*Wnorm - VdotW*Wz) * 2.*WyzCov) /
                         (Wnorm*Wnorm*Wnorm*Wnorm) +

                         ((Wx*Vnorm - VdotW*Vx) * (Wx*Vnorm - VdotW*Vx) * VxErr2 +
                          (Wy*Vnorm - VdotW*Vy) * (Wy*Vnorm - VdotW*Vy) * VyErr2 +
                          (Wz*Vnorm - VdotW*Vz) * (Wz*Vnorm - VdotW*Vz) * VzErr2 +

                          (Wx*Vnorm - VdotW*Vx) * (Wy*Vnorm - VdotW*Vy) * 2.*VxyCov +
                          (Wx*Vnorm - VdotW*Vx) * (Wz*Vnorm - VdotW*Vz) * 2.*VxzCov +
                          (Wy*Vnorm - VdotW*Vy) * (Wz*Vnorm - VdotW*Vz) * 2.*VyzCov) /
                         (Vnorm*Vnorm*Vnorm*Vnorm) ) / (Wnorm*Vnorm);
  }
  else
  {
    *cosAlpha = 0.;
    *cosAlphaErr = 0.;
  }
}

#if ROOFIT
using namespace RooFit;
RooAbsReal* Utils::ReadRTEffPDF (unsigned int q2BinIndx, RooRealVar* z, RooRealVar* y,RooRealVar* p)
{
  cout << "\n[Utils::ReadRTEffPDF]\t@@@ 3D angular*efficiency p.d.f. @@@" << endl;

  //#######################
  //# read eff parameters #
  //#######################

  std::stringstream myString;
  myString << (DirEfficiency != "" ? MakeAnalysisPATH(DirEfficiency).c_str() : "") << "effKEpdf_out_RT.root";
  std::cout << "[Utils::ReadRTEffPDF]\tReading 2D binned efficiency file : " << myString.str().c_str() << std::endl;
  TFile* file=TFile::Open(myString.str().c_str(),"READ");

  TString pdfName=Form("pdf_ctKctLphi_q2bin%d",q2BinIndx+1);
  RooHistPdf*  EffPDFm =(RooHistPdf*)file->Get(pdfName);
  RooHistFunc* EffPDF = new RooHistFunc("EffPDF","EffPDF",RooArgSet(*z,*y,*p),EffPDFm->dataHist(),1);
  // //EffPDF->setInterpolationOrder(1);
  // std::cout << "\n[ExtractYield::ReadRTEffPDF] \t@@@ Reading EFF parameters fit signal" <<  std::endl;  

  // // change variable names
  // RooArgSet* set=EffPDFm->getVariables();
  // set->find("ctK")->SetName("z");
  // set->find("ctL")->SetName("y");
  // set->find("phi")->SetName("p");

  // // TH3* h3 = (TH3*) EffPDFm ->createHistogram("ctK,ctL,phi");
  // // RooDataHist* EffPDF = new RooDataHist("Eff","Eff",RooArgList(*z,*y,*p),h3);
  // RooAbsPdf* EffPDF = new RooHistPdf("EffPDF","EffPDF",RooArgSet(*z,*y,*p),EffPDFm->dataHist(),1);
  return EffPDF;
}

RooAbsReal* Utils::ReadWTEffPDF (unsigned int q2BinIndx, RooRealVar* z, RooRealVar* y,RooRealVar* p)
{
  cout << "\n[Utils::ReadWTEffPDF]\t@@@ 3D angular*efficiency p.d.f. @@@" << endl;

  //# read eff parameters #
  //#######################
  std::stringstream myString;
  myString << (DirEfficiency != "" ? MakeAnalysisPATH(DirEfficiency).c_str() : "") << "effKEpdf_out_WT.root";
  std::cout << "[Utils::ReadWTEffPDF]\tReading 2D binned efficiency file : " << myString.str().c_str() << std::endl;

  TFile* file=TFile::Open(myString.str().c_str(),"READ");
  TString pdfName=Form("pdf_ctKctLphi_q2bin%d",q2BinIndx+1);
  RooHistPdf* EffPDFm =(RooHistPdf*)file->Get(pdfName);
  RooHistFunc* EffPDF = new RooHistFunc("EffPDF","EffPDF",RooArgSet(*z,*y,*p),EffPDFm->dataHist(),1);
  // //EffPDF->setInterpolationOrder(1);

  // // change variable names
  // std::cout << "\n[ExtractYield::ReadWTEffPDF] \t@@@ Reading EFF parameters fit signal" << std::endl;
  // RooArgSet* set=EffPDFm->getVariables();
  // set->find("ctK")->SetName("z");
  // set->find("ctL")->SetName("y");
  // set->find("phi")->SetName("p");

  // // TH3* histo =(TH3*) EffPDFm ->createHistogram("ctK,ctL,phi");
  // // RooAbsPdf* Eff = new RooDataHist("Eff","Eff",RooArgList(*z,*y,*p),h3);

  // // RooAbsPdf* EffPDF = new RooHistPdf("EffPDF","EffPDF",RooArgSet(*z,*y,*p),EffPDFm->dataHist(),1);
  // RooFormulaVar minusY("minusY","-(@0)",RooArgList(*y));
  // RooAbsPdf* EffPDF = new RooHistPdf("EffPDF","EffPDF",RooArgSet(*z,minusY,*p),RooArgSet(*z,*y,*p),EffPDFm->dataHist(),1);

  // // RooArgSet* params = EffPDF ->getVariables();
  // // params->setRealValue("ctK",0);
  return EffPDF;
}

#endif
void Utils::ReadAllBins (std::string fileName, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, std::string signalType)
  // ##########################
  // # signalType = "goodtag" #
  // # signalType = "mistag"  #
  // ##########################
{
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // #################
  // # Read q^2 bins #
  // #################
  std::cout << "\n[Utils::ReadAllBins]\tAll bins from file : " << fileName.c_str() << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("q2"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
  {
    q2Bins->push_back(atof(ParVector[i].c_str()));
    // std::cout << "Bin " << i << "\t" << q2Bins->back() << std::endl;
  }

  if (false) {
    if (signalType == "goodTag")
      {
	// #######################
	// # Read cosThetaK bins #
	// #######################
	std::cout << "\n@@@ cos(theta_K) bins from file @@@" << std::endl;
	ParVector.clear();
	ParameterFile->ReadFromFile(ParFileBlockN("thetaKokTag"),&ParVector);
	for (unsigned int i = 0; i < ParVector.size(); i++)
	  {
	    cosThetaKBins->push_back(atof(ParVector[i].c_str()));
	    std::cout << "Bin " << i << "\t" << cosThetaKBins->back() << std::endl;
	  }


	// #######################
	// # Read cosThetaL bins #
	// #######################
	std::cout << "\n@@@ cos(theta_l) bins from file @@@" << std::endl;
	ParVector.clear();
	ParameterFile->ReadFromFile(ParFileBlockN("thetaLokTag"),&ParVector);
	for (unsigned int i = 0; i < ParVector.size(); i++)
	  {
	    cosThetaLBins->push_back(atof(ParVector[i].c_str()));
	    std::cout << "Bin " << i << "\t" << cosThetaLBins->back() << std::endl;
	  }


	// #################
	// # Read phi bins #
	// #################
	std::cout << "\n@@@ phi bins from file @@@" << std::endl;
	ParVector.clear();
	ParameterFile->ReadFromFile(ParFileBlockN("phiokTag"),&ParVector);
	for (unsigned int i = 0; i < ParVector.size(); i++)
	  {
	    phiBins->push_back(atof(ParVector[i].c_str()));
	    std::cout << "Bin " << i << "\t" << phiBins->back() << std::endl;
	  }
      }
    else if (signalType == "misTag")
      {
	// #######################
	// # Read cosThetaK bins #
	// #######################
	std::cout << "\n@@@ cos(theta_K) bins from file @@@" << std::endl;
	ParVector.clear();
	ParameterFile->ReadFromFile(ParFileBlockN("thetaKmisTag"),&ParVector);
	for (unsigned int i = 0; i < ParVector.size(); i++)
	  {
	    cosThetaKBins->push_back(atof(ParVector[i].c_str()));
	    std::cout << "Bin " << i << "\t" << cosThetaKBins->back() << std::endl;
	  }


	// #######################
	// # Read cosThetaL bins #
	// #######################
	std::cout << "\n@@@ cos(theta_l) bins from file @@@" << std::endl;
	ParVector.clear();
	ParameterFile->ReadFromFile(ParFileBlockN("thetaLmisTag"),&ParVector);
	for (unsigned int i = 0; 2*i < ParVector.size(); i++)
	  {
	    cosThetaLBins->push_back(atof(ParVector[i].c_str()));
	    std::cout << "Bin " << i << "\t" << cosThetaLBins->back() << std::endl;
	  }


	// #################
	// # Read phi bins #
	// #################
	std::cout << "\n@@@ phi bins from file @@@" << std::endl;
	ParVector.clear();
	ParameterFile->ReadFromFile(ParFileBlockN("phimisTag"),&ParVector);
	for (unsigned int i = 0; i < ParVector.size(); i++)
	  {
	    phiBins->push_back(atof(ParVector[i].c_str()));
	    std::cout << "Bin " << i << "\t" << phiBins->back() << std::endl;
	  }
      }
    else
      {
	std::cout << "[Utils::ReadAllBins]\tError wrong parameter name : " << signalType << std::endl;
	exit (EXIT_FAILURE);
      }
  }

  ParVector.clear();
  delete ParameterFile;
}

void Utils::Readq2Bins (std::string fileName, std::vector<double>* q2Bins)
{
  std::vector<std::string>* ParVector = NULL;
  ReadParVsq2Bins(fileName,"q2",&ParVector);


  std::cout << "\n[Utils::Readq2Bins]\tq^2 bins from file : " << fileName.c_str() << std::endl;
  for (unsigned int i = 0; i < ParVector->size(); i++)
  {
    q2Bins->push_back(atof(ParVector->operator[](i).c_str()));
    std::cout << "Bin " << i << "\t" << q2Bins->back() << std::endl;
  }


  ParVector->clear();
  delete ParVector;
}

void Utils::ReadHLTpaths (std::string fileName, std::vector<std::string>* TrigTable)
{
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ###############################
  // # Read HLT-trigger table bins #
  // ###############################
  std::cout << "\n[Utils::ReadHLTpaths]\tHLT-trigger table from file : " << fileName.c_str() << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("HLTpath"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
  {
    TrigTable->push_back(ParVector[i]);
    std::cout << "Trigger path from config file : " << TrigTable->operator[](i).c_str() << std::endl;
  }


  ParVector.clear();
  delete ParameterFile;
}

int Utils::SearchBin (double val2Search, std::vector<double>* bins)
{
  unsigned int i = 0;
  for (i = 0; i < bins->size() - 1; i++) if ((val2Search >= bins->operator[](i)) && (val2Search < bins->operator[](i+1))) break;

  if (i != bins->size() - 1) return i;
  return -1;
}

int Utils::GetJPsiBin (std::vector<double>* q2Bins)
{
  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    if ((q2Bins->operator[](i) < (JPsiMass*JPsiMass)) && (q2Bins->operator[](i+1) > (JPsiMass*JPsiMass)))
      return i;

  return -1;
}

int Utils::GetPsiPBin (std::vector<double>* q2Bins)
{
  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    if ((q2Bins->operator[](i) < (PsiPMass*PsiPMass)) && (q2Bins->operator[](i+1) > (PsiPMass*PsiPMass)))
      return i;

  return -1;
}

bool Utils::ValIsInPsi (std::vector<double>* q2Bins, double q2Val)
{
  unsigned int JPsibin = GetJPsiBin(q2Bins);
  unsigned int PsiPbin = GetPsiPBin(q2Bins);

  if (((q2Bins->operator[](JPsibin) < q2Val) && (q2Bins->operator[](JPsibin+1) > q2Val)) ||
      ((q2Bins->operator[](PsiPbin) < q2Val) && (q2Bins->operator[](PsiPbin+1) > q2Val)))
    return true;

  return false;
}

unsigned int Utils::HLTpathForEvFraction (double evtFrac)
{
  // #####################################################################################
  // # This method is used together with the method: "ReadTriggerPathsANDCutsANDEntries" #
  // #####################################################################################

  double totalLumi = 0.0;
  double countLumi = 0.0;
  unsigned int it  = 0;

  for (unsigned int j = 0; j < VecHLTentries.size(); j++) totalLumi = totalLumi + VecHLTentries[j];
  if (totalLumi == 0.0) totalLumi = 1.0;

  for (it = 0; it < HLTpath.size(); it++)
  {
    if (evtFrac < (countLumi + VecHLTentries[it]) / totalLumi) break;
    countLumi = countLumi + VecHLTentries[it];
  }

  if (it == HLTpath.size())
  {
    std::cout << "[Utils::HLTpathForEvFraction]\tI didn't find the HLT category corresponding to the event fraction : " << evtFrac;
    std::cout << " (total HLT luminosity from config. file = " << totalLumi << ")" << std::endl;
    exit (EXIT_FAILURE);
  }

  return it+1;
}

unsigned int Utils::IsInTriggerTable (B0KstMuMuTreeContent* NTupleIn, double* HLTCutVar1, double* HLTCutVar2, int index, double evtFrac)
  // #####################################################################
  // # if index == -2 just check if the global trigger was fired         #
  // # if index == -1 just split the sample in HLT categories            #
  // # else           check both global and muon triggers                #
  // #####################################################################
  // # if evtFrac >= 0 associate HLT category by evtFrac & index         #
  // # if evtFrac < 0  associate HLT category by index (not index == -1) #
  // #####################################################################
{
  // #####################################################################################
  // # This method is used together with the method: "ReadTriggerPathsANDCutsANDEntries" #
  // #####################################################################################

  unsigned int it;
  unsigned int HLTpathIndx;

  for (unsigned int j = 0; j < HLTpath.size(); j++)
  {
    if (evtFrac >= 0.0) HLTpathIndx = HLTpathForEvFraction(evtFrac)-1;
    else                HLTpathIndx = j;
    *HLTCutVar1 = VecHLTCutVar1[HLTpathIndx];
    *HLTCutVar2 = VecHLTCutVar2[HLTpathIndx];
    if (index == -1) return HLTpathIndx+1;

    // ########################
    // # Global trigger check #
    // ########################
    for (it = 0; it < NTupleIn->TrigTable->size(); it++) if (NTupleIn->TrigTable->operator[](it).find(HLTpath[HLTpathIndx].c_str()) != std::string::npos) break;

    // ######################
    // # Muon trigger check #
    // ######################
    if ((it < NTupleIn->TrigTable->size()) &&
        ((index == -2) ||
         ((NTupleIn->mumTrig->at(index).find(HLTpath[HLTpathIndx].c_str()) != std::string::npos) && (NTupleIn->mupTrig->at(index).find(HLTpath[HLTpathIndx].c_str()) != std::string::npos))))
      return HLTpathIndx+1;

    if (evtFrac >= 0.0) break;
  }

  return 0;
}

unsigned int Utils::GetNHLTCat ()
{
  // #####################################################################################
  // # This method is used together with the method: "ReadTriggerPathsANDCutsANDEntries" #
  // #####################################################################################

  return HLTpath.size();
}

double Utils::ReadLumi (std::string fileName)
{
  double val = 0.0;
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ##############################
  // # Read integrated luminosity #
  // ##############################
  ParameterFile->ReadFromFile(ParFileBlockN("lumi"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++) val = val + atof(ParVector[i].c_str());
  std::cout << "\n[Utils::ReadLumi]\t@@@ Recorded luminosity: " << val << " fb-1 @@@" << std::endl;


  ParVector.clear();
  delete ParameterFile;
  return val;
}

void Utils::ReadTriggerPathsANDCutsANDEntries (std::string fileName)
{
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ############################################
  // # Read HLT-trigger paths and relative cuts #
  // ############################################
  ParameterFile->ReadFromFile(ParFileBlockN("HLTcuts"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i=i+4)
  {
    HLTpath.push_back(ParVector[i]);
    VecHLTCutVar1.push_back(atof(ParVector[i+1].c_str()));
    VecHLTCutVar2.push_back(atof(ParVector[i+2].c_str()));
    VecHLTentries.push_back(atof(ParVector[i+3].c_str()));

    std::cout << "\nRead trigger path from config file : " << HLTpath.back() << std::endl;
    std::cout << "Read first cut value: "                  << VecHLTCutVar1.back() << std::endl;
    std::cout << "Read second cut value: "                 << VecHLTCutVar2.back() << std::endl;
    std::cout << "Read entries in Data: "                  << VecHLTentries.back() << std::endl;
  }


  ParVector.clear();
  delete ParameterFile;
}

void Utils::ReadFitStartingValues (std::string fileName, std::vector<std::vector<std::string>*>* vecParam, std::vector<std::vector<unsigned int>*>* configParam, const unsigned int dataBlockN)
{
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ###################
  // # Clear vector(s) #
  // ###################
  vecParam->clear();
  configParam->clear();


  // ############################
  // # Read fit-starting values #
  // ############################
  ParameterFile->ReadFromFile(dataBlockN,&ParVector);

  for (unsigned int j = 0; j < nFitParam; j++) vecParam->push_back(new std::vector<std::string>);
  for (unsigned int j = 0; j < nConfigParam; j++) configParam->push_back(new std::vector<unsigned int>);

  for (unsigned int i = 0; i < ParVector.size(); i=i+nFitParam+nConfigParam)
  {

    // std::cout << "\nRead set-" << static_cast<int>(static_cast<double>(i)/static_cast<double>(nFitParam+nConfigParam)) << " of fit-parameter starting values" << std::endl;

    for (unsigned int j = 0; j < nFitParam; j++)
    {
      std::stringstream rawString(ParVector[i+j]);
      vecParam->operator[](j)->push_back(rawString.str());
      // std::cout << "Fit parameter-" << j << " value: " << vecParam->operator[](j)->back() << std::endl;
    }
    for (unsigned int j = 0; j < nConfigParam; j++)
    {
      configParam->operator[](j)->push_back(atoi(ParVector[i+nFitParam+j].c_str()));
      // std::cout << "Config. parameter-" << j << " value: " << configParam->operator[](j)->back() << std::endl;
    }
  }


  ParVector.clear();
  delete ParameterFile;
}

void Utils::ReadFitSystematics (std::string fileName, std::vector<std::vector<double>*>* vecParam)
  // ###########################################
  // # vecParam[0] --> Fl +err                 #
  // # vecParam[1] --> Fl -err                 #
  // # vecParam[2] --> P5p +err                #
  // # vecParam[3] --> P5p -err                #
  // # vecParam[4] --> P1 +err                 #
  // # vecParam[5] --> P1 -err                 #
  // # vecParam[6] --> P2 +err                 #
  // # vecParam[7] --> P2 -err                 #
  // # vecParam[8] --> Branching-Fraction +err #
  // # vecParam[9] --> Branching-Fraction -err #
  // ###########################################
{
  double val;
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ###################
  // # Clear vector(s) #
  // ###################
  vecParam->clear();


  // #########################################
  // # Read fit-observable systematic errors #
  // #########################################
  ParameterFile->ReadFromFile(ParFileBlockN("fitSyst"),&ParVector);

  for (unsigned int j = 0; j < nFitObserv*2; j++) vecParam->push_back(new std::vector<double>);

  for (unsigned int i = 0; i < ParVector.size(); i=i+nFitObserv)
  {

    std::cout << "\nRead set-" << static_cast<int>(static_cast<double>(i)/static_cast<double>(nFitObserv)) << " of fit-observable systematic errors" << std::endl;

    for (unsigned int j = 0; j < nFitObserv; j++)
    {
      std::stringstream rawString(ParVector[i+j]);
      rawString >> val;
      vecParam->operator[](j*2)->push_back(val);
      rawString >> val;
      vecParam->operator[](j*2+1)->push_back(val);
      std::cout << "Fit observable-" << j << " systematic error: +" << vecParam->operator[](j*2)->back() << "/-" << vecParam->operator[](j*2+1)->back() << std::endl;
    }
  }


  ParVector.clear();
  delete ParameterFile;
}

void Utils::ReadParVsq2Bins (std::string fileName, std::string praName, std::vector<std::string>** vecParam)
{
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());

  if (*vecParam == NULL) *vecParam = new std::vector<std::string>;
  else                   (*vecParam)->clear();

  ParameterFile->ReadFromFile(ParFileBlockN(praName.c_str()),*vecParam);

  std::cout << "\n[Utils::ReadParVsq2Bins]\tReading parameters vs q^2 from file : " << fileName << std::endl;
  for (unsigned int i = 0; i < (*vecParam)->size(); i++)
    std::cout << "Parameter value and errors for q2 bin " << i << ": " << (*vecParam)->operator[](i) << std::endl;

  delete ParameterFile;
}

void Utils::AddConstraint1D (TH1D** histo, std::string constrType, double abscissaErr, double YerrRescale, double Yval, double Yerr, unsigned int ID)
  // #########################################################
  // # constrType = "justErrors" : simply rescale the errors #
  // # constrType = "low"        : at lower boundary         #
  // # constrType = "both"       : at both boundaries        #
  // # constrType = "high"       : at higher boundary        #
  // #########################################################
{
  std::stringstream myString;

  TH1D* newHisto;
  TH1D* tmpHisto;

  unsigned int nNewBins;
  if ((constrType == "low") || (constrType == "high")) nNewBins = (*histo)->GetNbinsX()+2;
  else if (constrType == "both")                       nNewBins = (*histo)->GetNbinsX()+3;
  else if (constrType == "justErrors")                 nNewBins = (*histo)->GetNbinsX()+1;
  else { std::cout << "[Utils::AddConstraint1D]\tError wrong parameter name : " << constrType << std::endl; exit (EXIT_FAILURE); }
  double* reBins;
  reBins = new double[nNewBins];


  if (constrType == "low")
  {
    reBins[0] = (*histo)->GetBinLowEdge(1) - abscissaErr;
    for (int i = 1; i <= (*histo)->GetNbinsX(); i++) reBins[i] = (*histo)->GetBinLowEdge(i);
    reBins[(*histo)->GetNbinsX()+1] = (*histo)->GetBinLowEdge((*histo)->GetNbinsX()) + (*histo)->GetBinWidth((*histo)->GetNbinsX());
  }
  else if (constrType == "high")
  {
    for (int i = 1; i <= (*histo)->GetNbinsX(); i++) reBins[i-1] = (*histo)->GetBinLowEdge(i);
    reBins[(*histo)->GetNbinsX()] = (*histo)->GetBinLowEdge((*histo)->GetNbinsX()) + (*histo)->GetBinWidth((*histo)->GetNbinsX());
    reBins[(*histo)->GetNbinsX()+1] = reBins[(*histo)->GetNbinsX()] + abscissaErr;
  }
  else if (constrType == "both")
  {
    reBins[0] = (*histo)->GetBinLowEdge(1) - abscissaErr;
    for (int i = 1; i <= (*histo)->GetNbinsX(); i++) reBins[i] = (*histo)->GetBinLowEdge(i);
    reBins[(*histo)->GetNbinsX()+1] = (*histo)->GetBinLowEdge((*histo)->GetNbinsX()) + (*histo)->GetBinWidth((*histo)->GetNbinsX());
    reBins[(*histo)->GetNbinsX()+2] = reBins[(*histo)->GetNbinsX()+1] + abscissaErr;
  }
  else
  {
    for (int i = 1; i <= (*histo)->GetNbinsX(); i++) reBins[i-1] = (*histo)->GetBinLowEdge(i);
    reBins[(*histo)->GetNbinsX()] = (*histo)->GetBinLowEdge((*histo)->GetNbinsX()) + (*histo)->GetBinWidth((*histo)->GetNbinsX());
  }


  myString.clear(); myString.str("");
  myString << (*histo)->GetName() << "_" << ID;
  newHisto = new TH1D(myString.str().c_str(),myString.str().c_str(),nNewBins-1,reBins);


  if (constrType == "low")
  {
    newHisto->SetBinContent(1, Yval);
    newHisto->SetBinError(1,Yerr);
    for (int i = 1; i <= (*histo)->GetNbinsX(); i++)
    {
      newHisto->SetBinContent(i+1,(*histo)->GetBinContent(i));
      newHisto->SetBinError(i+1,(*histo)->GetBinError(i));
    }
  }
  else if (constrType == "high")
  {
    for (int i = 1; i <= (*histo)->GetNbinsX(); i++)
    {
      newHisto->SetBinContent(i,(*histo)->GetBinContent(i));
      newHisto->SetBinError(i,(*histo)->GetBinError(i));
    }
    newHisto->SetBinContent((*histo)->GetNbinsX()+1,Yval);
    newHisto->SetBinError((*histo)->GetNbinsX()+1,Yerr);
  }
  else if (constrType == "both")
  {
    newHisto->SetBinContent(1,Yval);
    newHisto->SetBinError(1,Yerr);
    for (int i = 1; i <= (*histo)->GetNbinsX(); i++)
    {
      newHisto->SetBinContent(i+1,(*histo)->GetBinContent(i));
      newHisto->SetBinError(i+1,(*histo)->GetBinError(i));
    }
    newHisto->SetBinContent((*histo)->GetNbinsX()+2,Yval);
    newHisto->SetBinError((*histo)->GetNbinsX()+2,Yerr);
  }
  else
  {
    for (int i = 1; i <= (*histo)->GetNbinsX(); i++)
    {
      newHisto->SetBinContent(i,(*histo)->GetBinContent(i));
      newHisto->SetBinError(i,(*histo)->GetBinError(i) * YerrRescale);
    }
  }


  tmpHisto = *histo;
  *histo = newHisto;
  delete tmpHisto;
  delete reBins;
}

void Utils::AddConstraintThetaL (TH1D** histo, unsigned int q2Indx, unsigned int cosThetaKBinIndx, unsigned int ID)
{
  double abscissaErr = 1e-2;


  if ((q2Indx == 0) && (cosThetaKBinIndx == 0)) AddConstraint1D(histo,"both",abscissaErr,1.0,2e-4,1e-5,ID);
  if ((q2Indx == 0) && (cosThetaKBinIndx == 1)) AddConstraint1D(histo,"both",abscissaErr,1.0,3e-4,1e-5,ID);
  if ((q2Indx == 0) && (cosThetaKBinIndx == 2)) AddConstraint1D(histo,"both",abscissaErr,1.0,4e-4,1e-5,ID);
  if ((q2Indx == 0) && (cosThetaKBinIndx == 3)) AddConstraint1D(histo,"both",abscissaErr,1.0,2e-4,1e-5,ID);
  if ((q2Indx == 0) && (cosThetaKBinIndx == 4)) AddConstraint1D(histo,"both",abscissaErr,1.0,1e-5,1e-5,ID);

  if ((q2Indx == 1) && (cosThetaKBinIndx == 0)) AddConstraint1D(histo,"low",abscissaErr,1.0,1e-5,1e-5,ID);
  if ((q2Indx == 1) && (cosThetaKBinIndx == 1)) AddConstraint1D(histo,"both",abscissaErr,1.0,1e-5,1e-5,ID);
  if ((q2Indx == 1) && (cosThetaKBinIndx == 2)) AddConstraint1D(histo,"both",abscissaErr,1.0,1e-5,1e-5,ID);
  if ((q2Indx == 1) && (cosThetaKBinIndx == 3)) AddConstraint1D(histo,"both",abscissaErr,1.0,1e-5,1e-5,ID);
  if ((q2Indx == 1) && (cosThetaKBinIndx == 4)) AddConstraint1D(histo,"both",abscissaErr,1.0,1e-5,1e-5,ID);

  if ((q2Indx == 2) && (cosThetaKBinIndx == 0)) AddConstraint1D(histo,"both",abscissaErr,1.0,1e-5,1e-5,ID);
  if ((q2Indx == 2) && (cosThetaKBinIndx == 1)) AddConstraint1D(histo,"both",abscissaErr,1.0,1e-5,1e-5,ID);
  if ((q2Indx == 2) && (cosThetaKBinIndx == 2)) AddConstraint1D(histo,"both",abscissaErr,1.0,1e-5,1e-5,ID);
  if ((q2Indx == 2) && (cosThetaKBinIndx == 3)) AddConstraint1D(histo,"high",abscissaErr,1.0,1e-5,1e-5,ID);
  if ((q2Indx == 2) && (cosThetaKBinIndx == 4)) AddConstraint1D(histo,"low",abscissaErr,1.0,2e-4,1e-5,ID);

  if ((q2Indx == 3) && (cosThetaKBinIndx == 0)) AddConstraint1D(histo,"low",abscissaErr,1.0,2e-4,1e-5,ID);
  if ((q2Indx == 3) && (cosThetaKBinIndx == 1)) AddConstraint1D(histo,"high",abscissaErr,1.0,1e-4,1e-5,ID);
  if ((q2Indx == 3) && (cosThetaKBinIndx == 2)) AddConstraint1D(histo,"low",abscissaErr,1.0,1e-4,1e-5,ID);
  if ((q2Indx == 3) && (cosThetaKBinIndx == 3)) AddConstraint1D(histo,"high",abscissaErr,1.0,1e-4,1e-5,ID);
  if ((q2Indx == 3) && (cosThetaKBinIndx == 4)) AddConstraint1D(histo,"low",abscissaErr,1.0,5e-4,1e-5,ID);
}

void Utils::AddConstraint2D (TH2D** histo, double abscissaErr, double ZerrRescale, unsigned int ID, std::string toBeConstr, double scaleConstr, double constrXerr, std::vector< std::pair <double,double> >* constraints, std::vector<std::string>* toBeAdded)
  // ################################################################################################################################################
  // # toBeConstr = Y --> Add constraints to Y axes (= cosThetaL) to both sides, according to the toBeAdded variable (= X axes binning = cosThetaK) #
  // # toBeConstr = X --> Add constraints to X axes (= cosThetaK) to the whole positive or negative side, according to the toBeConstr variable      #
  // ################################################################################################################################################
  // # toBeConstr = "justErrors" : simply rescale the errors                                                                                        #
  // # toBeConstr = "Xlow"       : at lower boundary                                                                                                #
  // # toBeConstr = "Xboth"      : at both boundaries                                                                                               #
  // # toBeConstr = "Xhigh"      : at higher boundary                                                                                               #
  // # toBeConstr = "Y" :                                                                                                                           #
  // # toBeAdded = "low"         : at lower boundary                                                                                                #
  // # toBeAdded = "both"        : at both boundaries                                                                                               #
  // # toBeAdded = "high"        : at higher boundary                                                                                               #
  // ################################################################################################################################################
{
  std::stringstream myString;

  double* reBinsX = NULL;
  double* reBinsY = NULL;

  TH2D* newHisto = NULL;
  TH2D* tmpHisto;

  TAxis* XAxis;
  TAxis* YAxis;


  if (toBeConstr == "Y")
  {
    unsigned int nNewBinsX;
    nNewBinsX = (*histo)->GetNbinsX()+1;
    reBinsX = new double[nNewBinsX];
    XAxis = (*histo)->GetXaxis();

    for (int i = 1; i <= (*histo)->GetNbinsX(); i++) reBinsX[i-1] = XAxis->GetBinLowEdge(i);
    reBinsX[(*histo)->GetNbinsX()] = XAxis->GetBinLowEdge((*histo)->GetNbinsX()) + XAxis->GetBinWidth((*histo)->GetNbinsX());


    unsigned int nNewBinsY;
    nNewBinsY = (*histo)->GetNbinsY()+3;
    reBinsY = new double[nNewBinsY];
    YAxis = (*histo)->GetYaxis();

    reBinsY[0] = YAxis->GetBinLowEdge(1) - abscissaErr;
    for (int i = 1; i <= (*histo)->GetNbinsY(); i++) reBinsY[i] = YAxis->GetBinLowEdge(i);
    reBinsY[(*histo)->GetNbinsY()+1] = YAxis->GetBinLowEdge((*histo)->GetNbinsY()) + YAxis->GetBinWidth((*histo)->GetNbinsY());
    reBinsY[(*histo)->GetNbinsY()+2] = reBinsY[(*histo)->GetNbinsY()+1] + abscissaErr;


    myString.clear(); myString.str("");
    myString << (*histo)->GetName() << "_" << ID;
    newHisto = new TH2D(myString.str().c_str(),myString.str().c_str(),nNewBinsX-1,reBinsX,nNewBinsY-1,reBinsY);


    for (unsigned int i = 0; i < toBeAdded->size(); i++)
    {
      if ((toBeAdded->operator[](i) == "low") || (toBeAdded->operator[](i) == "both"))
      {
        newHisto->SetBinContent(i+1,1,constraints->operator[](i).first);
        newHisto->SetBinError(i+1,1,constraints->operator[](i).second * ZerrRescale);
      }

      for (int j = 1; j <= (*histo)->GetNbinsY(); j++)
      {
        newHisto->SetBinContent(i+1,j+1,(*histo)->GetBinContent(i+1,j));
        newHisto->SetBinError(i+1,j+1,(*histo)->GetBinError(i+1,j) * ZerrRescale);
      }

      if ((toBeAdded->operator[](i) == "high") || (toBeAdded->operator[](i) == "both"))
      {
        newHisto->SetBinContent(i+1,(*histo)->GetNbinsY()+2,constraints->operator[](i).first);
        newHisto->SetBinError(i+1,(*histo)->GetNbinsY()+2,constraints->operator[](i).second * ZerrRescale);
      }
    }
  }
  else if (toBeConstr == "Xhigh")
  {
    unsigned int nNewBinsX;
    nNewBinsX = (*histo)->GetNbinsX()+2;
    reBinsX = new double[nNewBinsX];
    XAxis = (*histo)->GetXaxis();

    for (int i = 1; i <= (*histo)->GetNbinsX(); i++) reBinsX[i-1] = XAxis->GetBinLowEdge(i);
    reBinsX[(*histo)->GetNbinsX()] = XAxis->GetBinLowEdge((*histo)->GetNbinsX()) + XAxis->GetBinWidth((*histo)->GetNbinsX());
    reBinsX[(*histo)->GetNbinsX()+1] = reBinsX[(*histo)->GetNbinsX()] + abscissaErr;


    unsigned int nNewBinsY;
    nNewBinsY = (*histo)->GetNbinsY()+1;
    reBinsY = new double[nNewBinsY];
    YAxis = (*histo)->GetYaxis();

    for (int i = 1; i <= (*histo)->GetNbinsY(); i++) reBinsY[i-1] = YAxis->GetBinLowEdge(i);
    reBinsY[(*histo)->GetNbinsY()] = YAxis->GetBinLowEdge((*histo)->GetNbinsY()) + YAxis->GetBinWidth((*histo)->GetNbinsY());


    myString.clear(); myString.str("");
    myString << (*histo)->GetName() << "_" << ID;
    newHisto = new TH2D(myString.str().c_str(),myString.str().c_str(),nNewBinsX-1,reBinsX,nNewBinsY-1,reBinsY);


    for (int i = 1; i <= (*histo)->GetNbinsX(); i++)
      for (int j = 1; j <= (*histo)->GetNbinsY(); j++)
      {
        newHisto->SetBinContent(i,j,(*histo)->GetBinContent(i,j));
        newHisto->SetBinError(i,j,(*histo)->GetBinError(i,j) * ZerrRescale);
      }

    for (int j = 1; j <= newHisto->GetNbinsY(); j++)
    {
      newHisto->SetBinContent(newHisto->GetNbinsX(),j,newHisto->GetBinContent(newHisto->GetNbinsX()-1,j) / scaleConstr);
      newHisto->SetBinError(newHisto->GetNbinsX(),j,constrXerr * ZerrRescale);
    }
  }
  else if (toBeConstr == "Xlow")
  {
    unsigned int nNewBinsX;
    nNewBinsX = (*histo)->GetNbinsX()+2;
    reBinsX = new double[nNewBinsX];
    XAxis = (*histo)->GetXaxis();

    reBinsX[0] = XAxis->GetBinLowEdge(1) - abscissaErr;
    for (int i = 1; i <= (*histo)->GetNbinsX(); i++) reBinsX[i] = XAxis->GetBinLowEdge(i);
    reBinsX[(*histo)->GetNbinsX()+1] = XAxis->GetBinLowEdge((*histo)->GetNbinsX()) + XAxis->GetBinWidth((*histo)->GetNbinsX());


    unsigned int nNewBinsY;
    nNewBinsY = (*histo)->GetNbinsY()+1;
    reBinsY = new double[nNewBinsY];
    YAxis = (*histo)->GetYaxis();

    for (int i = 1; i <= (*histo)->GetNbinsY(); i++) reBinsY[i-1] = YAxis->GetBinLowEdge(i);
    reBinsY[(*histo)->GetNbinsY()] = YAxis->GetBinLowEdge((*histo)->GetNbinsY()) + YAxis->GetBinWidth((*histo)->GetNbinsY());


    myString.clear(); myString.str("");
    myString << (*histo)->GetName() << "_" << ID;
    newHisto = new TH2D(myString.str().c_str(),myString.str().c_str(),nNewBinsX-1,reBinsX,nNewBinsY-1,reBinsY);


    for (int i = 1; i <= (*histo)->GetNbinsX(); i++)
      for (int j = 1; j <= (*histo)->GetNbinsY(); j++)
      {
        newHisto->SetBinContent(i+1,j,(*histo)->GetBinContent(i,j));
        newHisto->SetBinError(i+1,j,(*histo)->GetBinError(i,j) * ZerrRescale);
      }

    for (int j = 1; j <= newHisto->GetNbinsY(); j++)
    {
      newHisto->SetBinContent(1,j,newHisto->GetBinContent(2,j) / scaleConstr);
      newHisto->SetBinError(1,j,constrXerr * ZerrRescale);
    }
  }
  else if (toBeConstr == "Xboth")
  {
    unsigned int nNewBinsX;
    nNewBinsX = (*histo)->GetNbinsX()+3;
    reBinsX = new double[nNewBinsX];
    XAxis = (*histo)->GetXaxis();

    reBinsX[0] = XAxis->GetBinLowEdge(1) - abscissaErr;
    for (int i = 1; i <= (*histo)->GetNbinsX(); i++) reBinsX[i] = XAxis->GetBinLowEdge(i);
    reBinsX[(*histo)->GetNbinsX()+1] = XAxis->GetBinLowEdge((*histo)->GetNbinsX()) + XAxis->GetBinWidth((*histo)->GetNbinsX());
    reBinsX[(*histo)->GetNbinsX()+2] = reBinsX[(*histo)->GetNbinsX()+1] + abscissaErr;


    unsigned int nNewBinsY;
    nNewBinsY = (*histo)->GetNbinsY()+1;
    reBinsY = new double[nNewBinsY];
    YAxis = (*histo)->GetYaxis();

    for (int i = 1; i <= (*histo)->GetNbinsY(); i++) reBinsY[i-1] = YAxis->GetBinLowEdge(i);
    reBinsY[(*histo)->GetNbinsY()] = YAxis->GetBinLowEdge((*histo)->GetNbinsY()) + YAxis->GetBinWidth((*histo)->GetNbinsY());


    myString.clear(); myString.str("");
    myString << (*histo)->GetName() << "_" << ID;
    newHisto = new TH2D(myString.str().c_str(),myString.str().c_str(),nNewBinsX-1,reBinsX,nNewBinsY-1,reBinsY);


    for (int i = 1; i <= (*histo)->GetNbinsX(); i++)
      for (int j = 1; j <= (*histo)->GetNbinsY(); j++)
      {
        newHisto->SetBinContent(i+1,j,(*histo)->GetBinContent(i,j));
        newHisto->SetBinError(i+1,j,(*histo)->GetBinError(i,j) * ZerrRescale);
      }

    for (int j = 1; j <= newHisto->GetNbinsY(); j++)
    {
      newHisto->SetBinContent(newHisto->GetNbinsX(),j,newHisto->GetBinContent(newHisto->GetNbinsX()-1,j) / scaleConstr);
      newHisto->SetBinError(newHisto->GetNbinsX(),j,constrXerr * ZerrRescale);

      newHisto->SetBinContent(1,j,newHisto->GetBinContent(2,j) / scaleConstr);
      newHisto->SetBinError(1,j,constrXerr * ZerrRescale);
    }
  }
  else if (toBeConstr == "justErrors")
  {
    unsigned int nNewBinsX;
    nNewBinsX = (*histo)->GetNbinsX()+1;
    reBinsX = new double[nNewBinsX];
    XAxis = (*histo)->GetXaxis();

    for (int i = 1; i <= (*histo)->GetNbinsX(); i++) reBinsX[i-1] = XAxis->GetBinLowEdge(i);
    reBinsX[(*histo)->GetNbinsX()] = XAxis->GetBinLowEdge((*histo)->GetNbinsX()) + XAxis->GetBinWidth((*histo)->GetNbinsX());


    unsigned int nNewBinsY;
    nNewBinsY = (*histo)->GetNbinsY()+1;
    reBinsY = new double[nNewBinsY];
    YAxis = (*histo)->GetYaxis();

    for (int i = 1; i <= (*histo)->GetNbinsY(); i++) reBinsY[i-1] = YAxis->GetBinLowEdge(i);
    reBinsY[(*histo)->GetNbinsY()] = YAxis->GetBinLowEdge((*histo)->GetNbinsY()) + YAxis->GetBinWidth((*histo)->GetNbinsY());


    myString.clear(); myString.str("");
    myString << (*histo)->GetName() << "_" << ID;
    newHisto = new TH2D(myString.str().c_str(),myString.str().c_str(),nNewBinsX-1,reBinsX,nNewBinsY-1,reBinsY);


    for (int i = 1; i <= (*histo)->GetNbinsX(); i++)
      for (int j = 1; j <= (*histo)->GetNbinsY(); j++)
      {
        newHisto->SetBinContent(i,j,(*histo)->GetBinContent(i,j));
        newHisto->SetBinError(i,j,(*histo)->GetBinError(i,j) * ZerrRescale);
      }
  }
  else
  {
    std::cout << "[Utils::AddConstraint2D]\tError wrong parameter name : " << toBeConstr << std::endl;
    exit (EXIT_FAILURE);
  }


  newHisto->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
  newHisto->GetXaxis()->SetTitleOffset(1.8);
  newHisto->SetYTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
  newHisto->GetYaxis()->SetTitleOffset(1.8);
  newHisto->SetZTitle("Efficiency");


  tmpHisto = *histo;
  *histo = newHisto;
  delete tmpHisto;
  delete reBinsX;
  delete reBinsY;
}

void Utils::AddConstraintThetaKThetaL (TH2D** histo, std::vector<double>* cosThetaKBins, unsigned int q2Indx, int SignalType, unsigned int ID)
{
  double abscissaErr = 1e-2;
  std::vector<std::string> toBeAdded;
  std::vector< std::pair<double,double> > constraints;

  for (unsigned int i = 0; i < cosThetaKBins->size()-1; i++)
  {
    toBeAdded.push_back("no");
    constraints.push_back(std::make_pair(0.0,0.0));
  }


  if (q2Indx == 0)
  {
    toBeAdded[0]   = "both";
    constraints[0] = std::make_pair(3e-4,1e-5);

    toBeAdded[1]   = "both";
    constraints[1] = std::make_pair(3e-4,1e-5);

    toBeAdded[2]   = "both";
    constraints[2] = std::make_pair(5e-4,1e-5);

    toBeAdded[3]   = "both";
    constraints[3] = std::make_pair(2e-4,1e-5);

    toBeAdded[4]   = "both";
    constraints[4] = std::make_pair(2e-4,1e-5);
  }
  else if (q2Indx == 1)
  {
    toBeAdded[0]   = "low";
    constraints[0] = std::make_pair(1e-5,1e-5);

    toBeAdded[1]   = "both";
    constraints[1] = std::make_pair(1e-5,1e-5);

    toBeAdded[2]   = "both";
    constraints[2] = std::make_pair(1e-5,1e-5);

    toBeAdded[3]   = "both";
    constraints[3] = std::make_pair(1e-5,1e-5);

    toBeAdded[4]   = "both";
    constraints[4] = std::make_pair(1e-5,1e-5);
  }
  else if (q2Indx == 2)
  {
    toBeAdded[0]   = "both";
    constraints[0] = std::make_pair(1e-5,1e-5);

    toBeAdded[1]   = "both";
    constraints[1] = std::make_pair(1e-5,1e-5);

    toBeAdded[2]   = "both";
    constraints[2] = std::make_pair(1e-5,1e-5);

    toBeAdded[3]   = "high";
    constraints[3] = std::make_pair(1e-5,1e-5);

    toBeAdded[4]   = "low";
    constraints[4] = std::make_pair(2e-4,1e-5);
  }
  else if (q2Indx == 3)
  {
    toBeAdded[0]   = "both";
    constraints[0] = std::make_pair(2e-4,1e-5);

    toBeAdded[1]   = "both";
    constraints[1] = std::make_pair(1e-4,1e-5);

    toBeAdded[2]   = "low";
    constraints[2] = std::make_pair(1e-4,1e-5);

    toBeAdded[3]   = "high";
    constraints[3] = std::make_pair(1e-4,1e-5);

    toBeAdded[4]   = "low";
    constraints[4] = std::make_pair(5e-4,1e-5);
  }

  else if ((SignalType == B0ToPsi2SKst) && (q2Indx == 6))
  {
    // #####################
    // # B0 --> psi(2S) K* #
    // #####################
    toBeAdded[4]   = "high";
    constraints[4] = std::make_pair(9e-4,1e-5);
  }


  AddConstraint2D(histo,abscissaErr,1.0,ID,"Y",0.0,0.0,&constraints,&toBeAdded);


  if      (q2Indx == 0) AddConstraint2D(histo,abscissaErr,6.0,ID,"justErrors",0.0,0.0);
  else if (q2Indx == 1) AddConstraint2D(histo,abscissaErr,2.0,ID,"justErrors",0.0,0.0);
  else if (q2Indx == 2) AddConstraint2D(histo,abscissaErr,2.0,ID,"justErrors",0.0,0.0);
  else if (q2Indx == 3) AddConstraint2D(histo,abscissaErr,3.0,ID,"justErrors",0.0,0.0);

  // #####################
  // # B0 --> psi(2S) K* #
  // #####################
  else if ((SignalType == B0ToPsi2SKst) && (q2Indx == 6)) AddConstraint2D(histo,abscissaErr,2.0,ID,"justErrors",0.0,0.0);
}

void Utils::AddConstraint3D (TH3D** histo, double abscissaErr, double Tval, double Terr, double TerrRescale, unsigned int ID, std::vector<int> toBeAdded[])
  // @TMP@ : to be reviewed
{
  std::stringstream myString;

  double* reBinsX = NULL;
  double* reBinsY = NULL;
  double* reBinsZ = NULL;

  TH3D* newHisto = NULL;
  TH3D* tmpHisto;

  TAxis* XAxis;
  TAxis* YAxis;
  TAxis* ZAxis;

  unsigned int nNewBinsX;
  unsigned int nNewBinsY;
  unsigned int nNewBinsZ;

  int deltaX;
  int deltaY;
  int deltaZ;


  std::cout << "\n[Utils::AddConstraint3D]" << std::endl;
  std::cout << "Old binnig value (theta_K,theta_l,phi): ";
  std::cout << (*histo)->GetNbinsX() << "," << (*histo)->GetNbinsY() << "," << (*histo)->GetNbinsZ() << std::endl;

  if (toBeAdded[0].size() == 0)
  {
    nNewBinsX = (*histo)->GetNbinsX()+1;
    reBinsX = new double[nNewBinsX];
    XAxis = (*histo)->GetXaxis();

    for (int i = 1; i <= (*histo)->GetNbinsX(); i++) reBinsX[i-1] = XAxis->GetBinLowEdge(i);
    reBinsX[(*histo)->GetNbinsX()] = XAxis->GetBinLowEdge((*histo)->GetNbinsX()) + XAxis->GetBinWidth((*histo)->GetNbinsX());

    deltaX = 0;
  }
  else
  {
    nNewBinsX = (*histo)->GetNbinsX()+3;
    reBinsX = new double[nNewBinsX];
    XAxis = (*histo)->GetXaxis();

    reBinsX[0] = XAxis->GetBinLowEdge(1) - abscissaErr;
    for (int i = 1; i <= (*histo)->GetNbinsX(); i++) reBinsX[i] = XAxis->GetBinLowEdge(i);
    reBinsX[(*histo)->GetNbinsX()+1] = XAxis->GetBinLowEdge((*histo)->GetNbinsX()) + XAxis->GetBinWidth((*histo)->GetNbinsX());
    reBinsX[(*histo)->GetNbinsX()+2] = reBinsX[(*histo)->GetNbinsX()+1] + abscissaErr;

    deltaX = 1;
  }

  if (toBeAdded[1].size() == 0)
  {
    nNewBinsY = (*histo)->GetNbinsY()+1;
    reBinsY = new double[nNewBinsY];
    YAxis = (*histo)->GetYaxis();

    for (int i = 1; i <= (*histo)->GetNbinsY(); i++) reBinsY[i-1] = YAxis->GetBinLowEdge(i);
    reBinsY[(*histo)->GetNbinsY()] = YAxis->GetBinLowEdge((*histo)->GetNbinsY()) + YAxis->GetBinWidth((*histo)->GetNbinsY());

    deltaY = 0;
  }
  else
  {
    nNewBinsY = (*histo)->GetNbinsY()+3;
    reBinsY = new double[nNewBinsY];
    YAxis = (*histo)->GetYaxis();

    reBinsY[0] = YAxis->GetBinLowEdge(1) - abscissaErr;
    for (int i = 1; i <= (*histo)->GetNbinsY(); i++) reBinsY[i] = YAxis->GetBinLowEdge(i);
    reBinsY[(*histo)->GetNbinsY()+1] = YAxis->GetBinLowEdge((*histo)->GetNbinsY()) + YAxis->GetBinWidth((*histo)->GetNbinsY());
    reBinsY[(*histo)->GetNbinsY()+2] = reBinsY[(*histo)->GetNbinsY()+1] + abscissaErr;

    deltaY = 1;
  }

  if (toBeAdded[2].size() == 0)
  {
    nNewBinsZ = (*histo)->GetNbinsZ()+1;
    reBinsZ = new double[nNewBinsZ];
    ZAxis = (*histo)->GetZaxis();

    for (int i = 1; i <= (*histo)->GetNbinsZ(); i++) reBinsZ[i-1] = ZAxis->GetBinLowEdge(i);
    reBinsZ[(*histo)->GetNbinsZ()] = ZAxis->GetBinLowEdge((*histo)->GetNbinsZ()) + ZAxis->GetBinWidth((*histo)->GetNbinsZ());

    deltaZ = 0;
  }
  else
  {
    nNewBinsZ = (*histo)->GetNbinsZ()+3;
    reBinsZ = new double[nNewBinsZ];
    ZAxis = (*histo)->GetZaxis();

    reBinsZ[0] = ZAxis->GetBinLowEdge(1) - abscissaErr;
    for (int i = 1; i <= (*histo)->GetNbinsZ(); i++) reBinsZ[i] = ZAxis->GetBinLowEdge(i);
    reBinsZ[(*histo)->GetNbinsZ()+1] = ZAxis->GetBinLowEdge((*histo)->GetNbinsZ()) + ZAxis->GetBinWidth((*histo)->GetNbinsZ());
    reBinsZ[(*histo)->GetNbinsZ()+2] = reBinsZ[(*histo)->GetNbinsZ()+1] + abscissaErr;

    deltaZ = 1;
  }


  myString.clear(); myString.str("");
  myString << (*histo)->GetName() << "_" << ID;
  newHisto = new TH3D(myString.str().c_str(),myString.str().c_str(),nNewBinsX-1,reBinsX,nNewBinsY-1,reBinsY,nNewBinsZ-1,reBinsZ);

  std::cout << "New binnig value (theta_K,theta_l,phi): ";
  std::cout << newHisto->GetNbinsX() << "," << newHisto->GetNbinsY() << "," << newHisto->GetNbinsZ() << std::endl;


  for (int i = 1; i <= newHisto->GetNbinsX(); i++)
    for (int j = 1; j <= newHisto->GetNbinsY(); j++)
      for (int k = 1; k <= newHisto->GetNbinsZ(); k++)
      {
        if (((toBeAdded[0].size() != 0) && ((i == 1) || (i == newHisto->GetNbinsX()))) ||
            ((toBeAdded[1].size() != 0) && ((j == 1) || (j == newHisto->GetNbinsY()))) ||
            ((toBeAdded[2].size() != 0) && ((k == 1) || (k == newHisto->GetNbinsZ()))))
        {
          newHisto->SetBinContent(i,j,k,0.0);
          newHisto->SetBinError(i,j,k,0.0);
        }
        else
        {
          newHisto->SetBinContent(i,j,k,(*histo)->GetBinContent(i-deltaX,j-deltaY,k-deltaZ));
          newHisto->SetBinError(i,j,k,(*histo)->GetBinError(i-deltaX,j-deltaY,k-deltaZ) * TerrRescale);
        }

        for (unsigned int n = 0; n < toBeAdded[0].size(); n++)
          if ((toBeAdded[0][n] == i) && (toBeAdded[1][n] == j) && (toBeAdded[2][n] == k))
          {
            newHisto->SetBinContent(i,j,k,Tval);
            newHisto->SetBinError(i,j,k,Terr);

            std::cout << "Added constraint " << Tval << "+/-" << Terr << " to bin [#bin convention 1...N] (theta_K,theta_l,phi): ";
            std::cout << i << "," << j << "," << k << std::endl;
            break;
          }
      }


  newHisto->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
  newHisto->GetXaxis()->SetTitleOffset(1.8);
  newHisto->SetYTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
  newHisto->GetYaxis()->SetTitleOffset(1.8);
  newHisto->SetZTitle("#phi");
  newHisto->GetZaxis()->SetTitleOffset(1.8);


  tmpHisto = *histo;
  *histo = newHisto;
  delete tmpHisto;
  delete reBinsX;
  delete reBinsY;
  delete reBinsZ;
}

void Utils::AddConstraintThetaKThetaLPhi (int SignalType)
{
  std::cout << "[Utils::AddConstraintThetaKThetaLPhi]\t @TMP@ Not implemented yet : " << SignalType << std::endl;
  exit (EXIT_FAILURE);
}

int Utils::WhatIsThis (std::string fileName)
{
  int val = 0;
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ###############################
  // # Read data type (Data or MC) #
  // ###############################
  ParameterFile->ReadFromFile(ParFileBlockN("dtype"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++) val = atoi(ParVector[i].c_str());


  ParVector.clear();
  delete ParameterFile;
  return val;
}

void Utils::SaveFitValues (std::string fileName, std::vector<std::string>* vecParStr, int indx, std::string howOpen, std::string str)
  // #################################################
  // # If indx == -1 then use str within default str #
  // # If indx == -2 then use only str               #
  // #################################################
{
  std::stringstream myString;

  myString.clear(); myString.str("");
  if      (indx == -1) myString << "#@@@@@@@@@@@@@@@@@@@@@@@@@@@" << str.c_str()   <<  "@@@@@@@@@@@@@@@@@@@@@@@@@@@#";
  else if (indx == -2) myString << str.c_str();
  else                 myString << "#@@@@@@@@@@@@@@@@@@@@@@@@@@@ q^2 bin " << indx << " @@@@@@@@@@@@@@@@@@@@@@@@@@@#";
  vecParStr->insert(vecParStr->begin(),myString.str());

  vecParStr->insert(vecParStr->end(),"");
  vecParStr->insert(vecParStr->end(),"");

  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str(),howOpen.c_str());
  ParameterFile->WriteToFile(vecParStr);
  delete ParameterFile;
}


std::string Utils::MakeAnalysisPATH (std::string relativePath)
{
  std::stringstream myString;
  char* absolutePath = getenv(ANALYPATH);

  if (absolutePath == NULL)
  {
    std::cout << "[Utils::MakeAnalysisPATH]\tAnalysis environment variable " << ANALYPATH << " not defined" << std::endl;
    exit (EXIT_FAILURE);
  }
  myString.clear(); myString.str("");
  myString << absolutePath << relativePath;

  return myString.str();
}


unsigned int Utils::ParFileBlockN (std::string blockName)
{
  // #########################
  // # Parameter file blocks #
  // #########################
  if      (blockName == "HLTpath")        return 1;
  else if (blockName == "precuts")        return 2;
  else if (blockName == "selecuts")       return 3;

  else if (blockName == "q2")             return 4;

  else if (blockName == "thetaKokTag")    return 5;
  else if (blockName == "thetaLokTag")    return 6;
  else if (blockName == "phiokTag")       return 7;

  else if (blockName == "thetaKmisTag")   return 8;
  else if (blockName == "thetaLmisTag")   return 9;
  else if (blockName == "phimisTag")      return 10;

  else if (blockName == "HLTcuts")        return 11;
  else if (blockName == "fitValBins")     return 12;

  else if (blockName == "fitSyst")        return 13;

  else if (blockName == "genericpar")     return 14;
  else if (blockName == "lumi")           return 15;
  else if (blockName == "dtype")          return 16;

  std::cout << "[Utils::ParFileBlockN]\tError wrong index name : " << blockName << std::endl;
  exit (EXIT_FAILURE);
}

unsigned int Utils::GetFitParamIndx (std::string varName)
{
  if      (varName == "meanS")          return 0;
  else if (varName == "sigmaS1")        return 1;
  else if (varName == "sigmaS2")        return 2;
  else if (varName == "fracMassS")      return 3;

  else if (varName == "var1")           return 4;
  else if (varName == "var2")           return 5;
  else if (varName == "fracMassBExp")   return 6;

  else if (varName == "sigmaMisTag1")   return 7;
  else if (varName == "sigmaMisTag2")   return 8;
  else if (varName == "fracMisTag")     return 9;

  else if (varName == "meanR1")         return 10;
  else if (varName == "sigmaR1")        return 11;
  else if (varName == "meanR2")         return 12;
  else if (varName == "sigmaR2")        return 13;
  else if (varName == "fracMassBRPeak") return 14;

  else if (varName == "meanL1")         return 15;
  else if (varName == "sigmaL1")        return 16;
  else if (varName == "meanL2")         return 17;
  else if (varName == "sigmaL2")        return 18;
  else if (varName == "fracMassBLPeak") return 19;

  else if (varName == "fracMassBPeak")  return 20;

  else if (varName == "nMisTagFrac")    return 21;
  else if (varName == "nBkgPeak")       return 22;
  else if (varName == "nSig")           return 23;
  else if (varName == "nPolyP1")        return 24;
  else if (varName == "p1Poly0")        return 25;
  else if (varName == "p1Poly1")        return 26;
  else if (varName == "p1Poly2")        return 27;
  else if (varName == "p1Poly3")        return 28;
  else if (varName == "p1Poly4")        return 29;

  else if (varName == "nPolyC1")        return 30;
  else if (varName == "c1Poly0")        return 31;
  else if (varName == "c1Poly1")        return 32;
  else if (varName == "c1Poly2")        return 33;
  else if (varName == "c1Poly3")        return 34;
  else if (varName == "c1Poly4")        return 35;

  else if (varName == "nPolyP2")        return 36;
  else if (varName == "p2Poly0")        return 37;
  else if (varName == "p2Poly1")        return 38;
  else if (varName == "p2Poly2")        return 39;
  else if (varName == "p2Poly3")        return 40;
  else if (varName == "p2Poly4")        return 41;

  else if (varName == "nPolyC2")        return 42;
  else if (varName == "c2Poly0")        return 43;
  else if (varName == "c2Poly1")        return 44;
  else if (varName == "c2Poly2")        return 45;
  else if (varName == "c2Poly3")        return 46;
  else if (varName == "c2Poly4")        return 47;

  else if (varName == "nPolyP3")        return 48;
  else if (varName == "p3Poly0")        return 49;
  else if (varName == "p3Poly1")        return 50;
  else if (varName == "p3Poly2")        return 51;
  else if (varName == "p3Poly3")        return 52;
  else if (varName == "p3Poly4")        return 53;

  else if (varName == "nPolyC3")        return 54;
  else if (varName == "c3Poly0")        return 55;
  else if (varName == "c3Poly1")        return 56;
  else if (varName == "c3Poly2")        return 57;
  else if (varName == "c3Poly3")        return 58;
  else if (varName == "c3Poly4")        return 59;

  else if (varName == "FlS")            return 60;
  else if (varName == "P5pS")           return 61;
  else if (varName == "P1S")            return 62;
  else if (varName == "P2S")            return 63;
  else if (varName == "FsS")            return 64;
  else if (varName == "AsS")            return 65;
  else if (varName == "As5S")           return 66;

  std::cout << "[Utils::GetFitParamIndx]\tError wrong index name : " << varName << std::endl;
  exit (EXIT_FAILURE);
}

unsigned int Utils::GetConfigParamIndx (std::string varName)
{
  if      (varName == "SigType")      return 0;
  else if (varName == "PeakBkgType")  return 1;
  else if (varName == "CombBkgType")  return 2;
  else if (varName == "MistagType")   return 3;

  std::cout << "[Utils::GetConfigParamIndx]\tError wrong index name : " << varName << std::endl;
  exit (EXIT_FAILURE);
}

bool Utils::PsiRejection (double myB0Mass, double myMuMuMass, double myMuMuMassE, std::string seleType, bool B0andPsiCut)
  // ###########################
  // # seleType == "keepJpsi"  #
  // # seleType == "keepPsiP"  #
  // # seleType == "rejectPsi" #
  // # seleType == "keepPsi"   #
  // ###########################
{
  // ####################################################################
  // # This method is used together with the method: "ReadGenericParam" #
  // ####################################################################

  if (seleType == "keepJpsi")
  {
    if ((fabs(myMuMuMass - JPsiMass) < atof(GetGenericParam("NSigmaPsi").c_str()) * myMuMuMassE) && (B0andPsiCut == false))

      return true;
  }
  else if (seleType == "keepPsiP")
  {
    if ((fabs(myMuMuMass - PsiPMass) < atof(GetGenericParam("NSigmaPsi").c_str()) * myMuMuMassE) && (B0andPsiCut == false))

      return true;
  }
  else if (seleType == "rejectPsi")
  {
    if ((fabs(myMuMuMass - JPsiMass) > atof(GetGenericParam("NSigmaPsi").c_str()) * myMuMuMassE) &&
        (fabs(myMuMuMass - PsiPMass) > atof(GetGenericParam("NSigmaPsi").c_str()) * myMuMuMassE) &&

        ((B0andPsiCut == false) ||

         ((B0andPsiCut == true) &&

          (((myMuMuMass < JPsiMass) &&
            (!((fabs((myB0Mass - B0Mass) - (myMuMuMass - JPsiMass)) < atof(GetGenericParam("B&psiMassJpsiLo").c_str()))    ||
               (fabs((myB0Mass - B0Mass) - (myMuMuMass - PsiPMass)) < atof(GetGenericParam("B&psiMassPsiPLo").c_str()))))) ||

           ((myMuMuMass > PsiPMass) &&
            (!((fabs((myB0Mass - B0Mass) - (myMuMuMass - JPsiMass)) < atof(GetGenericParam("B&psiMassJpsiHi").c_str()))    ||
               (fabs((myB0Mass - B0Mass) - (myMuMuMass - PsiPMass)) < atof(GetGenericParam("B&psiMassPsiPHi").c_str()))))) ||

           ((myMuMuMass > JPsiMass) && (myMuMuMass < PsiPMass) &&
            (!((fabs((myB0Mass - B0Mass) - (myMuMuMass - JPsiMass)) < atof(GetGenericParam("B&psiMassJpsiHi").c_str()))     ||
               (fabs((myB0Mass - B0Mass) - (myMuMuMass - PsiPMass)) < atof(GetGenericParam("B&psiMassPsiPLo").c_str())))))))))

      return true;
  }
  else if (seleType == "keepPsi")
  {
    if (((fabs(myMuMuMass - JPsiMass) < atof(GetGenericParam("NSigmaPsi").c_str()) * myMuMuMassE)  ||
         (fabs(myMuMuMass - PsiPMass) < atof(GetGenericParam("NSigmaPsi").c_str()) * myMuMuMassE)) && (B0andPsiCut == false))

      return true;
  }
  else
  {
    std::cout << "[Utils::PsiRejection]\tSelection type not valid : " << seleType << std::endl;
    exit (EXIT_FAILURE);
  }

  return false;
}

bool Utils::ChooseBestCand (B0KstMuMuTreeContent* NTuple, unsigned int DoTrigCheck, double evFraction, int* BestCandIndx, bool* B0notB0bar, int* TrigCat, unsigned int* countCands)
  // ##############################################################################
  // # DoTrigCheck: to allow check on trigger requirements                        #
  // # 0: do not perform any trigger check                                        #
  // # 1: perform normal trigger check (i.e. both global and muon triggers are    #
  // #    associated to at least one trigger in configuration file)               #
  // # 2: perform partial trigger check (i.e. global trigger is associated to at  #
  // #    least one trigger in configuration file)                                #
  // # 3: perform normal trigger check and split the sample in HLT categories     #
  // # 4: do NOT perform any trigger check and split the sample in HLT categories #
  // ##############################################################################
{
  // #####################################################################
  // # This method is used together with the method: "ReadSelectionCuts" #
  // #####################################################################

  double MuMuVtxCL = GetSeleCut("MuMuVtxCL");
  double MinMupT   = GetSeleCut("MinMupT");
  double BestVal   = 0.0;
  double BestValTmp;

  *countCands   =  0;
  *BestCandIndx = -1;
  *TrigCat      =  0;
  for (unsigned int i = 0; i < NTuple->bMass->size(); i++)
  {
    // ##############################################
    // # Candidate selection through kinematic cuts #
    // ##############################################
    if (((DoTrigCheck == 0)                                                                           ||
         ((DoTrigCheck == 1) && (IsInTriggerTable(NTuple, &MuMuVtxCL, &MinMupT, i) >= 1))             ||
         ((DoTrigCheck == 2) && (IsInTriggerTable(NTuple, &MuMuVtxCL, &MinMupT, -2) >= 1))            ||
         ((DoTrigCheck == 3) && (IsInTriggerTable(NTuple, &MuMuVtxCL, &MinMupT, i, evFraction) >= 1)) ||
         ((DoTrigCheck == 4) && (IsInTriggerTable(NTuple, &MuMuVtxCL, &MinMupT, -1, evFraction) >= 1))) &&

        // #####################
        // # Choose good muons #
        // #####################
        (NTuple->mumNTrkLayers->at(i) > static_cast<int>(rint(GetSeleCut("NTrkLayers")))) &&
        (NTuple->mumNPixLayers->at(i) > static_cast<int>(rint(GetSeleCut("NPixLayers")))) &&
        (NTuple->mumNormChi2->at(i) < GetSeleCut("NormChi2"))                             &&
        (fabs(NTuple->mumdxyVtx->at(i)) < GetSeleCut("dxyVtx"))                           &&
        (fabs(NTuple->mumdzVtx->at(i)) < GetSeleCut("dzVtx"))                             &&
        (NTuple->mumCat->at(i).find("TMOneStationTight") != std::string::npos)            &&

        (NTuple->mupNTrkLayers->at(i) > static_cast<int>(rint(GetSeleCut("NTrkLayers")))) &&
        (NTuple->mupNPixLayers->at(i) > static_cast<int>(rint(GetSeleCut("NPixLayers")))) &&
        (NTuple->mupNormChi2->at(i) < GetSeleCut("NormChi2"))                             &&
        (fabs(NTuple->mupdxyVtx->at(i)) < GetSeleCut("dxyVtx"))                           &&
        (fabs(NTuple->mupdzVtx->at(i)) < GetSeleCut("dzVtx"))                             &&
        (NTuple->mupCat->at(i).find("TMOneStationTight") != std::string::npos)            &&

        // #########################################
        // # Check that hadron track is NOT a muon #
        // #########################################
        (NTuple->kstTrkmMuMatch->at(i).find("TrackerMuonArbitrated") == std::string::npos) &&
        (NTuple->kstTrkpMuMatch->at(i).find("TrackerMuonArbitrated") == std::string::npos) &&

        // #####################
        // # B0 selection cuts #
        // #####################
        (NTuple->bLBS->at(i)/NTuple->bLBSE->at(i) > GetSeleCut("B0LsBS"))                                          &&
        (NTuple->bVtxCL->at(i) > GetSeleCut("B0VtxCL"))                                                            &&
        (NTuple->bCosAlphaBS->at(i) > GetSeleCut("B0cosAlpha"))                                                    &&
        (sqrt(NTuple->bPx->at(i)*NTuple->bPx->at(i) + NTuple->bPy->at(i)*NTuple->bPy->at(i)) > GetSeleCut("B0pT")) &&
        (fabs(computeEta(NTuple->bPx->at(i),NTuple->bPy->at(i),NTuple->bPz->at(i))) < GetSeleCut("B0Eta"))         &&

        // #######################
        // # Muon selection cuts #
        // #######################
        (sqrt(NTuple->mumPx->at(i)*NTuple->mumPx->at(i) + NTuple->mumPy->at(i)*NTuple->mumPy->at(i)) > MinMupT) &&
        (sqrt(NTuple->mupPx->at(i)*NTuple->mupPx->at(i) + NTuple->mupPy->at(i)*NTuple->mupPy->at(i)) > MinMupT) &&

        // #########################
        // # Dimuon selection cuts #
        // #########################
        (NTuple->mumuVtxCL->at(i) > MuMuVtxCL) &&

        // #########################
        // # Hadron selection cuts #
        // #########################
        (NTuple->kstTrkmHighPurity->at(i) == true) &&
        (NTuple->kstTrkpHighPurity->at(i) == true) &&

        (fabs(NTuple->kstTrkmDCABS->at(i)/NTuple->kstTrkmDCABSE->at(i)) > GetSeleCut("HadDCASBS")) &&
        (fabs(NTuple->kstTrkpDCABS->at(i)/NTuple->kstTrkpDCABSE->at(i)) > GetSeleCut("HadDCASBS")) &&

        (sqrt(NTuple->kstTrkmPx->at(i)*NTuple->kstTrkmPx->at(i) + NTuple->kstTrkmPy->at(i)*NTuple->kstTrkmPy->at(i)) > GetSeleCut("HadpT")) &&
        (sqrt(NTuple->kstTrkpPx->at(i)*NTuple->kstTrkpPx->at(i) + NTuple->kstTrkpPy->at(i)*NTuple->kstTrkpPy->at(i)) > GetSeleCut("HadpT")) &&

        // #####################
        // # K* selection cuts #
        // #####################
        ((fabs(NTuple->kstMass->at(i) - kstMass) < GetSeleCut("KstMass")) || (fabs(NTuple->kstBarMass->at(i) - kstMass) < GetSeleCut("KstMass"))) &&

        // #####################
        // # phi rejection cut #
        // #####################
        (computeInvMass(NTuple->kstTrkmPx->at(i),NTuple->kstTrkmPy->at(i),NTuple->kstTrkmPz->at(i),kaonMass,
                        NTuple->kstTrkpPx->at(i),NTuple->kstTrkpPy->at(i),NTuple->kstTrkpPz->at(i),kaonMass) > GetSeleCut("KKMass")))
                        {
                          (*countCands)++;

                          BestValTmp = NTuple->bVtxCL->at(i);
                          if (BestValTmp > BestVal)
                          {
                            if      (DoTrigCheck == 1) *TrigCat = IsInTriggerTable(NTuple, &MuMuVtxCL, &MinMupT, i);
                            else if (DoTrigCheck == 2) *TrigCat = IsInTriggerTable(NTuple, &MuMuVtxCL, &MinMupT, -2);
                            else if (DoTrigCheck == 3) *TrigCat = IsInTriggerTable(NTuple, &MuMuVtxCL, &MinMupT, i, evFraction);
                            else if (DoTrigCheck == 4) *TrigCat = IsInTriggerTable(NTuple, &MuMuVtxCL, &MinMupT, -1, evFraction);

                            BestVal       = BestValTmp;
                            *BestCandIndx = i;
                          }
                        }
  }


  if ((*BestCandIndx != -1) && (FlavorTagger(NTuple, *BestCandIndx, B0notB0bar) == true)) return true;
  return false;
}

bool Utils::FlavorTagger (B0KstMuMuTreeContent* NTuple, unsigned int i, bool* B0notB0bar)
{
  // #####################################################################
  // # Computation of the probability related to the K*0 mass hypothesis #
  // #####################################################################
  double distKst    = fabs(KstMassShape->Integral(kstMass,NTuple->kstMass->at(i)));
  double distKstBar = fabs(KstMassShape->Integral(kstMass,NTuple->kstBarMass->at(i)));
  if (distKst < distKstBar) *B0notB0bar = true;
  else                      *B0notB0bar = false;


  // # ###########################################
  // # Scramble the tagging of the CP-eigenstate #
  // # ###########################################
  double rndDice = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
  if (rndDice < scrambleFraction)
  {
    rndDice = static_cast<double>(rand()) / static_cast<double>(RAND_MAX);
    if (rndDice < 0.5)
    {
      if (*B0notB0bar == false)
      {
        // #######################
        // # Swap the K* - K*bar #
        // #######################

        distKst = NTuple->kstMass->at(i);
        NTuple->kstMass->at(i) = NTuple->kstBarMass->at(i);
        NTuple->kstBarMass->at(i) = distKst;

        distKst = NTuple->kstTrkpPx->at(i);
        NTuple->kstTrkpPx->at(i) = NTuple->kstTrkmPx->at(i);
        NTuple->kstTrkmPx->at(i) = distKst;

        distKst = NTuple->kstTrkpPy->at(i);
        NTuple->kstTrkpPy->at(i) = NTuple->kstTrkmPy->at(i);
        NTuple->kstTrkmPy->at(i) = distKst;

        distKst = NTuple->kstTrkpPz->at(i);
        NTuple->kstTrkpPz->at(i) = NTuple->kstTrkmPz->at(i);
        NTuple->kstTrkmPz->at(i) = distKst;
      }
      *B0notB0bar = true;
    }
    else
    {
      if (*B0notB0bar == true)
      {
        // #######################
        // # Swap the K* - K*bar #
        // #######################

        distKst = NTuple->kstMass->at(i);
        NTuple->kstMass->at(i) = NTuple->kstBarMass->at(i);
        NTuple->kstBarMass->at(i) = distKst;

        distKst = NTuple->kstTrkpPx->at(i);
        NTuple->kstTrkpPx->at(i) = NTuple->kstTrkmPx->at(i);
        NTuple->kstTrkmPx->at(i) = distKst;

        distKst = NTuple->kstTrkpPy->at(i);
        NTuple->kstTrkpPy->at(i) = NTuple->kstTrkmPy->at(i);
        NTuple->kstTrkmPy->at(i) = distKst;

        distKst = NTuple->kstTrkpPz->at(i);
        NTuple->kstTrkpPz->at(i) = NTuple->kstTrkmPz->at(i);
        NTuple->kstTrkmPz->at(i) = distKst;
      }
      *B0notB0bar = false;
    }
  }


  // ####################################################################
  // # Compute the Fisher-test to see whether the (mK*0 - K*0 PDG mass) #
  // # is significantly different from (mK*0bar - K*0 PDG mass)         #
  // ####################################################################
  double VarianceRatio;
  double chi2Kst    = pow((NTuple->kstMass->at(i)    - kstMass) / NTuple->kstMassE->at(i),2.);
  double chi2KstBar = pow((NTuple->kstBarMass->at(i) - kstMass) / NTuple->kstBarMassE->at(i),2.);
  if (chi2Kst > chi2KstBar)
    VarianceRatio = chi2Kst    / chi2KstBar;
  else VarianceRatio = chi2KstBar / chi2Kst;

  if (TMath::FDistI(VarianceRatio,1.,1.) > ProbThreshold) return true;
  return false;
}

void Utils::ReadSelectionCuts (std::string fileName)
  // #############################
  // # MuMuVtxCL  = SeleCuts[0]  #
  // # MinMupT    = SeleCuts[1]  #
  // # B0LsBS     = SeleCuts[2]  #
  // # B0VtxCL    = SeleCuts[3]  #
  // # B0cosAlpha = SeleCuts[4]  #
  // # B0pT       = SeleCuts[5]  #
  // # B0Eta      = SeleCuts[6]  #
  // # HadDCASBS  = SeleCuts[7]  #
  // # HadpT      = SeleCuts[8]  #
  // # KstMass    = SeleCuts[9]  #
  // # KKMass     = SeleCuts[10] #
  // # NTrkLayers = SeleCuts[11] #
  // # NPixLayers = SeleCuts[12] #
  // # NormChi2   = SeleCuts[13] #
  // # dxyVtx     = SeleCuts[14] #
  // # dzVtx      = SeleCuts[15] #
  // #############################
{
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // #######################
  // # Read selection cuts #
  // #######################
  std::cout << "\n[Utils::ReadSelectionCuts]\tSelection cuts from file : " << fileName.c_str() << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("selecuts"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
  {
    SeleCuts.push_back(atof(ParVector[i].c_str()));
    // std::cout << "Selection cut #" << i << " from config file : " << SeleCuts[i] << std::endl;
  }


  ParVector.clear();
  delete ParameterFile;
}

bool Utils::SetSeleCut (std::string cutName, double val)
{
  if      (cutName == "MuMuVtxCL")  SeleCuts[0]  = val;
  else if (cutName == "MinMupT")    SeleCuts[1]  = val;
  else if (cutName == "B0LsBS")     SeleCuts[2]  = val;
  else if (cutName == "B0VtxCL")    SeleCuts[3]  = val;
  else if (cutName == "B0cosAlpha") SeleCuts[4]  = val;
  else if (cutName == "B0pT")       SeleCuts[5]  = val;
  else if (cutName == "B0Eta")      SeleCuts[6]  = val;
  else if (cutName == "HadDCASBS")  SeleCuts[7]  = val;
  else if (cutName == "HadpT")      SeleCuts[8]  = val;
  else if (cutName == "KstMass")    SeleCuts[9]  = val;
  else if (cutName == "KKMass")     SeleCuts[10] = val;
  else if (cutName == "NTrkLayers") SeleCuts[11] = val;
  else if (cutName == "NPixLayers") SeleCuts[12] = val;
  else if (cutName == "NormChi2")   SeleCuts[13] = val;
  else if (cutName == "dxyVtx")     SeleCuts[14] = val;
  else if (cutName == "dzVtx")      SeleCuts[15] = val;
  else return false;

  return true;
}

double Utils::GetSeleCut (std::string cutName)
{
  if      (cutName == "MuMuVtxCL")  return SeleCuts[0];
  else if (cutName == "MinMupT")    return SeleCuts[1];
  else if (cutName == "B0LsBS")     return SeleCuts[2];
  else if (cutName == "B0VtxCL")    return SeleCuts[3];
  else if (cutName == "B0cosAlpha") return SeleCuts[4];
  else if (cutName == "B0pT")       return SeleCuts[5];
  else if (cutName == "B0Eta")      return SeleCuts[6];
  else if (cutName == "HadDCASBS")  return SeleCuts[7];
  else if (cutName == "HadpT")      return SeleCuts[8];
  else if (cutName == "KstMass")    return SeleCuts[9];
  else if (cutName == "KKMass")     return SeleCuts[10];
  else if (cutName == "NTrkLayers") return SeleCuts[11];
  else if (cutName == "NPixLayers") return SeleCuts[12];
  else if (cutName == "NormChi2")   return SeleCuts[13];
  else if (cutName == "dxyVtx")     return SeleCuts[14];
  else if (cutName == "dzVtx")      return SeleCuts[15];
  else
  {
    std::cout << "[Utils::GetSeleCut]\tSelection cut not valid : " << cutName << std::endl;
    exit (EXIT_FAILURE);
  }
}

void Utils::ReadPreselectionCut (std::string fileName)
  // ################################
  // # MuMuVtxCL      = PreCuts[0]  #
  // # MuMuLsBS       = PreCuts[1]  #
  // # DCAMuMu        = PreCuts[2]  #
  // # DCAMuBS        = PreCuts[3]  #
  // # cosAlphaMuMuBS = PreCuts[4]  #
  // # MinMupT        = PreCuts[5]  #
  // # MuEta          = PreCuts[6]  #
  // # MuMupT         = PreCuts[7]  #
  // # MinMuMuMass    = PreCuts[8]  #
  // # MaxMuMuMass    = PreCuts[9]  #
  // # MinB0Mass      = PreCuts[10] #
  // # MaxB0Mass      = PreCuts[11] #
  // # B0VtxCL        = PreCuts[12] #
  // # KstMass        = PreCuts[13] #
  // # HadDCASBS      = PreCuts[14] #
  // # HadpT          = PreCuts[15] #
  // ################################
{
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ###########################
  // # Read pre-selection cuts #
  // ###########################
  std::cout << "\n[Utils::ReadPreselectionCut]\tPre-selection cuts from file : " << fileName.c_str() << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("precuts"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
  {
    PreCuts.push_back(atof(ParVector[i].c_str()));
    std::cout << "Pre-selection cut #" << i << " from config file : " << PreCuts[i] << std::endl;
  }


  ParVector.clear();
  delete ParameterFile;
}

bool Utils::SetPreCut (std::string cutName, double val)
{
  if      (cutName == "MuMuVtxCL")      PreCuts[0]  = val;
  else if (cutName == "MuMuLsBS")       PreCuts[1]  = val;
  else if (cutName == "DCAMuMu")        PreCuts[2]  = val;
  else if (cutName == "DCAMuBS")        PreCuts[3]  = val;
  else if (cutName == "cosAlphaMuMuBS") PreCuts[4]  = val;
  else if (cutName == "MinMupT")        PreCuts[5]  = val;
  else if (cutName == "MuEta")          PreCuts[6]  = val;
  else if (cutName == "MuMupT")         PreCuts[7]  = val;
  else if (cutName == "MinMuMuMass")    PreCuts[8]  = val;
  else if (cutName == "MaxMuMuMass")    PreCuts[9]  = val;
  else if (cutName == "MinB0Mass")      PreCuts[10] = val;
  else if (cutName == "MaxB0Mass")      PreCuts[11] = val;
  else if (cutName == "B0VtxCL")        PreCuts[12] = val;
  else if (cutName == "KstMass")        PreCuts[13] = val;
  else if (cutName == "HadDCASBS")      PreCuts[14] = val;
  else if (cutName == "HadpT")          PreCuts[15] = val;
  else return false;

  return true;
}

double Utils::GetPreCut (std::string cutName)
{
  if      (cutName == "MuMuVtxCL")      return PreCuts[0];
  else if (cutName == "MuMuLsBS")       return PreCuts[1];
  else if (cutName == "DCAMuMu")        return PreCuts[2];
  else if (cutName == "DCAMuBS")        return PreCuts[3];
  else if (cutName == "cosAlphaMuMuBS") return PreCuts[4];
  else if (cutName == "MinMupT")        return PreCuts[5];
  else if (cutName == "MuEta")          return PreCuts[6];
  else if (cutName == "MuMupT")         return PreCuts[7];
  else if (cutName == "MinMuMuMass")    return PreCuts[8];
  else if (cutName == "MaxMuMuMass")    return PreCuts[9];
  else if (cutName == "MinB0Mass")      return PreCuts[10];
  else if (cutName == "MaxB0Mass")      return PreCuts[11];
  else if (cutName == "B0VtxCL")        return PreCuts[12];
  else if (cutName == "KstMass")        return PreCuts[13];
  else if (cutName == "HadDCASBS")      return PreCuts[14];
  else if (cutName == "HadpT")          return PreCuts[15];
  else
  {
    std::cout << "[Utils::GetPreCut]\tPre-selection cut not valid : " << cutName << std::endl;
    exit (EXIT_FAILURE);
  }
}


void Utils::ReadGenericParam (std::string fileName)
  // #########################################
  // # NormJPSInotPSIP     = GenericPars[0]  #
  // # DegreeInterp        = GenericPars[1]  #
  // # TransfTolerance     = GenericPars[2]  #
  // # UseMINOS            = GenericPars[3]  #
  // # ApplyConstr         = GenericPars[4]  #
  // # CtrlFitWrkFlow      = GenericPars[5]  #
  // # CtrlMisTagWrkFlow   = GenericPars[6]  #
  // # SaveMisTagFrac      = GenericPars[7]  #
  // # UseSPwave           = GenericPars[8]  #
  // # doTransf            = GenericPars[9]  #
  // # B0MassIntervalLeft  = GenericPars[10] #
  // # B0MassIntervalRight = GenericPars[11] #
  // # NSigmaB0            = GenericPars[12] #
  // # NSigmaB0S           = GenericPars[13] #
  // # NSigmaB0B           = GenericPars[14] #
  // # NSigmaPsi           = GenericPars[15] #
  // # B&psiMassJpsiLo     = GenericPars[16] #
  // # B&psiMassJpsiHi     = GenericPars[17] #
  // # B&psiMassPsiPLo     = GenericPars[18] #
  // # B&psiMassPsiPHi     = GenericPars[19] #
  // # SIGMAS1             = GenericPars[20] #
  // # SIGMAS2             = GenericPars[21] #
  // # FRACMASSS           = GenericPars[22] #
  // #########################################
{
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ###########################
  // # Read generic parameters #
  // ###########################
  std::cout << "\n[Utils::ReadGenericParam]\tGeneric parameters from file : " << fileName.c_str() << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("genericpar"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
  {
    GenericPars.push_back(ParVector[i]);
    // std::cout << "Generic parameter #" << i << " from config file : " << GenericPars[i].c_str() << std::endl;
  }


  ParVector.clear();
  delete ParameterFile;
}

bool Utils::SetGenericParam (std::string parName, std::string val)
{
  if (parName == "UseMINOS")            GenericPars[0]  = val;
  else if (parName == "ApplyConstr")         GenericPars[1]  = val;
  else if (parName == "CtrlFitWrkFlow")      GenericPars[2]  = val;
  else if (parName == "CtrlMisTagWrkFlow")   GenericPars[3]  = val;
  else if (parName == "SaveMisTagFrac")      GenericPars[4]  = val;
  else if (parName == "B0MassIntervalLeft")  GenericPars[5] = val;
  else if (parName == "B0MassIntervalRight") GenericPars[6] = val;
  else if (parName == "NSigmaB0")            GenericPars[7] = val;
  else if (parName == "NSigmaB0S")           GenericPars[8] = val;
  else if (parName == "NSigmaB0B")           GenericPars[9] = val;
  else if (parName == "NSigmaPsi")           GenericPars[10] = val;
  else if (parName == "B&psiMassJpsiLo")     GenericPars[11] = val;
  else if (parName == "B&psiMassJpsiHi")     GenericPars[12] = val;
  else if (parName == "B&psiMassPsiPLo")     GenericPars[13] = val;
  else if (parName == "B&psiMassPsiPHi")     GenericPars[14] = val;
  else if (parName == "SIGMAS1")             GenericPars[15] = val;
  else if (parName == "SIGMAS2")             GenericPars[16] = val;
  else if (parName == "FRACMASSS")           GenericPars[17] = val;
  else return false;

  return true;
}

std::string Utils::GetGenericParam (std::string parName)
{
  if (parName == "UseMINOS")            return GenericPars[0];
  else if (parName == "ApplyConstr")         return GenericPars[1];
  else if (parName == "CtrlFitWrkFlow")      return GenericPars[2];
  else if (parName == "CtrlMisTagWrkFlow")   return GenericPars[3];
  else if (parName == "SaveMisTagFrac")      return GenericPars[4];
  else if (parName == "B0MassIntervalLeft")  return GenericPars[5];
  else if (parName == "B0MassIntervalRight") return GenericPars[6];
  else if (parName == "NSigmaB0")            return GenericPars[7];
  else if (parName == "NSigmaB0S")           return GenericPars[8];
  else if (parName == "NSigmaB0B")           return GenericPars[9];
  else if (parName == "NSigmaPsi")           return GenericPars[10];
  else if (parName == "B&psiMassJpsiLo")     return GenericPars[11];
  else if (parName == "B&psiMassJpsiHi")     return GenericPars[12];
  else if (parName == "B&psiMassPsiPLo")     return GenericPars[13];
  else if (parName == "B&psiMassPsiPHi")     return GenericPars[14];
  else if (parName == "SIGMAS1")             return GenericPars[15];
  else if (parName == "SIGMAS2")             return GenericPars[16];
  else if (parName == "FRACMASSS")           return GenericPars[17];
  else if (parName == "UseSPwave")           return GenericPars[18];
  else
  {
    std::cout << "[Utils::GetGenericParam]\tGeneric parameter not valid : " << parName << std::endl;
    exit (EXIT_FAILURE);
  }
}

double Utils::GetB0Width ()
{
  if (atof(GetGenericParam("FRACMASSS").c_str()) != 0.0)
    return sqrt(atof(GetGenericParam("FRACMASSS").c_str()) * pow(atof(GetGenericParam("SIGMAS1").c_str()),2.) + (1.-atof(GetGenericParam("FRACMASSS").c_str())) * pow(atof(GetGenericParam("SIGMAS2").c_str()),2.));
  else
    return atof(GetGenericParam("SIGMAS1").c_str());
}

double* Utils::MakeBinning (std::vector<double>* STLvec)
{
  double* vec = new double[STLvec->size()];

  for (unsigned int i = 0; i < STLvec->size(); i++) vec[i] = STLvec->operator[](i);

  return vec;
}

#if ROOFIT
std::string Utils::Transformer (std::string varName, bool doIt, double& varValOut, double& varValOutELo, double& varValOutEHi, RooFitResult* fitResult, RooRealVar* varValIn1, RooRealVar* varValIn2, RooRealVar* varValIn3, RooRealVar* varValIn4)
// #######################
// # varValIn1 = Fl, Fs  #
// # varValIn2 = Afb, Fs #
// # varValIn3 = As      #
// #######################
{
  double theorCoeff = 0.89; // @TMP@ : to be reviewed
  const TMatrixTSym<double>* CovM = (fitResult != NULL ? &fitResult->covarianceMatrix() : NULL);
  double val1,val2,val1ELo,val1EHi,val2ELo,val2EHi;
  std::string sVal1,sVal2;
  std::stringstream myString;
  myString.clear(); myString.str("");


  if (varValIn1 == NULL)
    {
      if (varName == "FlS")
	{
	  if (doIt == true) myString << "(1/2 + TMath::ATan(" << varName << ")/TMath::Pi())";
	  else              myString << "(" << varName << ")";
	  std::cout << "[Utils::Transformer]\tTransformer function: " << myString.str().c_str() << std::endl;
	  return myString.str();
	}
      else if (varName == "AfbS")
	{
	  sVal1 = Transformer("FlS",doIt,val1,val1ELo,val1EHi);

	  if (doIt == true) myString << "(3/4*(1 - " << sVal1 << ") * 2*TMath::ATan(" << varName << ")/TMath::Pi())";
	  else              myString << "(" << varName << ")";
	  std::cout << "[Utils::Transformer]\tTransformer function: " << myString.str().c_str() << std::endl;
	  return myString.str();
	}
      else if (varName == "FsS")
	{
	  myString << "(" << varName << ")";
	  std::cout << "[Utils::Transformer]\tTransformer function: " << myString.str().c_str() << std::endl;
	  return myString.str();
	}
      else if (varName == "AsS")
	{
	  sVal1 = Transformer("FlS",doIt,val1,val1ELo,val1EHi);
	  sVal2 = Transformer("FsS",doIt,val2,val2ELo,val2EHi);

	  myString << "(" << varName << "*2*" << theorCoeff << "*sqrt(3*" << sVal2 << "*(1-" << sVal2 << ")*" << sVal1 << "))";
	  std::cout << "[Utils::Transformer]\tTransformer function: " << myString.str().c_str() << std::endl;
	  return myString.str();
	}
      else if (varName == "As5S")
	{
	  sVal1 = Transformer("FlS",doIt,val1,val1ELo,val1EHi);
          sVal2 = Transformer("FsS",doIt,val2,val2ELo,val2EHi);

	  myString << "(" << varName << "*" << theorCoeff << "*sqrt(3*" << sVal2 << "*(1-" << sVal2 << ")*(1-" << sVal1 << ")*(1+P1S)))";
	  std::cout << "[Utils::Transformer]\tTransformer function: " << myString.str().c_str() << std::endl;
	  return myString.str();
	}
      else
	{
	  std::cout << "[Utils::Transformer]\tWrong parameter: " << varName << std::endl;
	  exit (EXIT_FAILURE);
	}
    }
  else if ((varName == "FlS") && (varValIn1 != NULL))
    {
      varValOut    = 1./2. + TMath::ATan(varValIn1->getVal()) / TMath::Pi();

      varValOutELo = 1./2. + TMath::ATan(varValIn1->getVal() + varValIn1->getErrorLo()) / TMath::Pi() - varValOut;
      varValOutEHi = 1./2. + TMath::ATan(varValIn1->getVal() + varValIn1->getErrorHi()) / TMath::Pi() - varValOut;

      if (doIt == false)
	{
	  varValOut    = varValIn1->getVal();
	  varValOutELo = varValIn1->getErrorLo();
	  varValOutEHi = varValIn1->getErrorHi();
	}
    }
  else if ((varName == "AfbS") && (varValIn1 != NULL) && (varValIn2 != NULL))
    {
      Transformer("FlS",doIt,val1,val1ELo,val1EHi,fitResult,varValIn1);
      val2 = 3./4. * (1. - val1);

      varValOut    = val2 * 2.*TMath::ATan(varValIn2->getVal()) / TMath::Pi();

      varValOutELo = val2 * 2.*TMath::ATan(varValIn2->getVal() + varValIn2->getErrorLo()) / TMath::Pi() - varValOut;
      varValOutEHi = val2 * 2.*TMath::ATan(varValIn2->getVal() + varValIn2->getErrorHi()) / TMath::Pi() - varValOut;

      varValOutELo = - sqrt( pow(3./4. * 2.*TMath::ATan(varValIn2->getVal()) / TMath::Pi() * val1ELo,2.) + pow(varValOutELo,2.) +

      			     (varValIn1->getErrorLo() != 0.0 && varValIn2->getErrorLo() != 0.0 ?
      			      2. *
      			      (3./4. * 2.*TMath::ATan(varValIn2->getVal()) / TMath::Pi() * val1ELo) / varValIn1->getErrorLo() *
      			      varValOutELo / varValIn2->getErrorLo() *
      			      (CovM != NULL && fitResult->floatParsFinal().index("FlS") != -1 && fitResult->floatParsFinal().index("AfbS") != -1 ?
      			       (*CovM)(fitResult->floatParsFinal().index("FlS"),fitResult->floatParsFinal().index("AfbS")) : 0.) : 0.));

      varValOutEHi = + sqrt( pow(3./4. * 2.*TMath::ATan(varValIn2->getVal()) / TMath::Pi() * val1EHi,2.) + pow(varValOutEHi,2.) +

      			     (varValIn1->getErrorHi() != 0.0 && varValIn2->getErrorHi() != 0.0 ?
      			      2. *
      			      (3./4. * 2.*TMath::ATan(varValIn2->getVal()) / TMath::Pi() * val1EHi)  / varValIn1->getErrorHi() *
      			      varValOutEHi / varValIn2->getErrorHi() *
      			      (CovM != NULL && fitResult->floatParsFinal().index("FlS") != -1 && fitResult->floatParsFinal().index("AfbS") != -1 ?
      			       (*CovM)(fitResult->floatParsFinal().index("FlS"),fitResult->floatParsFinal().index("AfbS")) : 0.) : 0.));

      if (doIt == false)
	{
	  varValOut    = varValIn2->getVal();
	  varValOutELo = varValIn2->getErrorLo();
	  varValOutEHi = varValIn2->getErrorHi();
	}
    }
  else if ((varName == "FsS") && (varValIn1 != NULL))
    {
      varValOut    = varValIn1->getVal();
      varValOutELo = varValIn1->getErrorLo();
      varValOutEHi = varValIn1->getErrorHi();
    }
  else if ((varName == "AsS") && (varValIn1 != NULL) && (varValIn2 != NULL) && (varValIn3 != NULL))
    {
      Transformer("FlS",doIt,val1,val1ELo,val1EHi,fitResult,varValIn1);
      Transformer("FsS",doIt,val2,val2ELo,val2EHi,fitResult,varValIn2);

      varValOut = varValIn3->getVal() * 2.*theorCoeff * sqrt(3. * val2 * (1. - val2) * val1);

      varValOutELo = - sqrt( pow(2.*theorCoeff * sqrt(3. * val2 * (1. - val2) * val1) * varValIn3->getErrorLo(),2.) +
			     pow(varValIn3->getVal() * 2.*theorCoeff * sqrt(3. * val2 * (1. - val2) / (4. * val1)) * varValIn1->getErrorLo(),2.) +
			     pow(varValIn3->getVal() * 2.*theorCoeff * sqrt(3. * val1 / (4. * val2 * (1. - val2))) * (1. - 2.*val2) * varValIn2->getErrorLo(),2.));

      varValOutEHi = + sqrt( pow(2.*theorCoeff * sqrt(3. * val2 * (1. - val2) * val1) * varValIn3->getErrorHi(),2.) +
			     pow(varValIn3->getVal() * 2.*theorCoeff * sqrt(3. * val2 * (1. - val2) / (4. * val1)) * varValIn1->getErrorHi(),2.) +
			     pow(varValIn3->getVal() * 2.*theorCoeff * sqrt(3. * val1 / (4. * val2 * (1. - val2))) * (1. - 2.*val2) * varValIn2->getErrorHi(),2.));
    }
  else if ((varName == "As5S") && (varValIn1 != NULL) && (varValIn2 != NULL) && (varValIn3 != NULL) && (varValIn4 != NULL))
    {
      Transformer("FlS",doIt,val1,val1ELo,val1EHi,fitResult,varValIn1);
      Transformer("FsS",doIt,val2,val2ELo,val2EHi,fitResult,varValIn2);

      varValOut = varValIn3->getVal() * theorCoeff * sqrt(3. * val2 * (1. - val2) * (1. - val1) * (1. + varValIn4->getVal()));

      varValOutELo = - sqrt( pow(varValOut / varValIn3->getVal() * varValIn3->getErrorLo(),2.) +
			     pow(varValOut / 2 * varValIn4->getErrorLo(),2.)/(1. + varValIn4->getVal()));
			     
      varValOutEHi = + sqrt( pow(varValOut / varValIn3->getVal() * varValIn3->getErrorHi(),2.) +
                             pow(varValOut / 2 * varValIn4->getErrorHi(),2.)/(1. + varValIn4->getVal()));
    }
  else
    {
      std::cout << "[Utils::Transformer]\tWrong parameter: " << varName << std::endl;
      exit (EXIT_FAILURE);
    }


  return "";
}


void Utils::AntiTransformer (std::string varName, bool doIt, double& varValOut, double& varValOutELo, double& varValOutEHi, RooRealVar* varValIn1, RooRealVar* varValIn2)
// ##########################
// # varValIn1 = Fl, Fs, As #
// # varValIn2 = Afb        #
// ##########################
{
  double tmpVal,limit;
  double tolerance = atof(GetGenericParam("TransfTolerance").c_str());


  if ((varName == "FlS") && (varValIn1 != NULL))
    {
      varValOut = TMath::Tan((varValIn1->getVal() - 1./2.) * TMath::Pi());

      if (varValIn1->getErrorLo() == -1.0) varValOutELo = varValIn1->getErrorLo();
      else
	{
	  if ((varValIn1->getVal() + varValIn1->getErrorLo()) <= 0.) tmpVal = tolerance;
	  else                                                       tmpVal = varValIn1->getVal() + varValIn1->getErrorLo();
	  varValOutELo = TMath::Tan((tmpVal - 1./2.) * TMath::Pi()) - varValOut;
	}

      if (varValIn1->getErrorHi() == 1.0) varValOutEHi = varValIn1->getErrorHi();
      else
	{
	  if ((varValIn1->getVal() + varValIn1->getErrorHi()) >= 1.) tmpVal = 1. - tolerance;
	  else                                                       tmpVal = varValIn1->getVal() + varValIn1->getErrorHi();
	  varValOutEHi = TMath::Tan((tmpVal - 1./2.) * TMath::Pi()) - varValOut;
	}

      if (doIt == false)
	{
	  varValOut    = varValIn1->getVal();
	  varValOutELo = varValIn1->getErrorLo();
	  varValOutEHi = varValIn1->getErrorHi();
	}
    }
  else if ((varName == "AfbS") && (varValIn1 != NULL) && (varValIn2 != NULL))
    {
      limit = 3./4. * (1. - varValIn1->getVal());
      
      varValOut = TMath::Tan(varValIn2->getVal() / limit / 2. * TMath::Pi());

      if (varValIn2->getErrorLo() == -1.0) varValOutELo = varValIn2->getErrorLo();
      else
	{
	  if ((varValIn2->getVal() + varValIn2->getErrorLo()) <= -limit) tmpVal = (limit >= 2.*tolerance ? -limit + tolerance : -limit + tolerance*limit);
	  else                                                           tmpVal = varValIn2->getVal() + varValIn2->getErrorLo();
	  varValOutELo = TMath::Tan(tmpVal / limit / 2. * TMath::Pi()) - varValOut;
	}

      if (varValIn2->getErrorHi() == 1.0) varValOutEHi = varValIn2->getErrorHi();
      else
	{
	  if ((varValIn2->getVal() + varValIn2->getErrorHi()) >= limit) tmpVal = (limit >= 2.*tolerance ? limit - tolerance : limit - tolerance*limit);
	  else                                                          tmpVal = varValIn2->getVal() + varValIn2->getErrorHi();
	  varValOutEHi = TMath::Tan(tmpVal / limit / 2. * TMath::Pi()) - varValOut;
	}

      if (doIt == false)
	{
	  varValOut    = varValIn2->getVal();
	  varValOutELo = varValIn2->getErrorLo();
	  varValOutEHi = varValIn2->getErrorHi();
	}
    }
  else if ((varName == "FsS") && (varValIn1 != NULL))
    {
      varValOut    = varValIn1->getVal();
      varValOutELo = varValIn1->getErrorLo();
      varValOutEHi = varValIn1->getErrorHi();
    }
  else if ((varName == "AsS") && (varValIn1 != NULL))
    {
      varValOut    = varValIn1->getVal();
      varValOutELo = varValIn1->getErrorLo();
      varValOutEHi = varValIn1->getErrorHi();
    }
  else
    {
      std::cout << "[Utils::AntiTransformer]\tWrong parameter: " << varName << std::endl;
      exit (EXIT_FAILURE);
    }
}
#endif
