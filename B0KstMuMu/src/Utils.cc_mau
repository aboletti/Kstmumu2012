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
  nFitObserv   = 5; // FL --- AFB --- P1 --- P2 --- BF

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
  DirEfficiency  = "/efficiency/";

  Histo2DEffNameOkTagSig   = "H2Deff_OkTag_q2Bin_interp";
  Histo2DEffNameOkTagJPsi  = "H2Deff_OkTagJPsi_q2Bin_interp";
  Histo2DEffNameOkTagPsi2S = "H2Deff_OkTagPsiP_q2Bin_interp";

  Histo3DEffNameOkTagSig   = "H3Deff_OkTag_q2Bin_interp";
  Histo3DEffNameOkTagJPsi  = "H3Deff_OkTagJPsi_q2Bin_interp";
  Histo3DEffNameOkTagPsi2S = "H3Deff_OkTagPsiP_q2Bin_interp";

  Histo2DEffNameMisTagSig   = "H2Deff_MisTag_q2Bin_interp";
  Histo2DEffNameMisTagJPsi  = "H2Deff_MisTagJPsi_q2Bin_interp";
  Histo2DEffNameMisTagPsi2S = "H2Deff_MisTagPsiP_q2Bin_interp";

  Histo3DEffNameMisTagSig   = "H3Deff_MisTag_q2Bin_interp";
  Histo3DEffNameMisTagJPsi  = "H3Deff_MisTagJPsi_q2Bin_interp";
  Histo3DEffNameMisTagPsi2S = "H3Deff_MisTagPsiP_q2Bin_interp";

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
  std::cout << "\n@@@@@@ Utils class settings : private @@@@@@" << std::endl;
  std::cout << "Analysis environment variable: " << ANALYPATH << std::endl;
  std::cout << "nFitObserv: "                << nFitObserv << std::endl;
  std::cout << "ProbThreshold: "             << ProbThreshold << std::endl;
  std::cout << "scrambleFraction: "          << scrambleFraction << std::endl;
  std::cout << "DirEfficiency: "             << DirEfficiency << std::endl;

  std::cout << "\nHisto2DEffNameOkTagSig: "  << Histo2DEffNameOkTagSig << std::endl;
  std::cout << "Histo2DEffNameOkTagJPsi: "   << Histo2DEffNameOkTagJPsi << std::endl;
  std::cout << "Histo2DEffNameOkTagPsi2S: "  << Histo2DEffNameOkTagPsi2S << std::endl;

  std::cout << "\nHisto3DEffNameOkTagSig: "  << Histo3DEffNameOkTagSig << std::endl;
  std::cout << "Histo3DEffNameOkTagJPsi: "   << Histo3DEffNameOkTagJPsi << std::endl;
  std::cout << "Histo3DEffNameOkTagPsi2S: "  << Histo3DEffNameOkTagPsi2S << std::endl;

  std::cout << "\nHisto2DEffNameMisTagSig: " << Histo2DEffNameMisTagSig << std::endl;
  std::cout << "Histo2DEffNameMisTagJPsi: "  << Histo2DEffNameMisTagJPsi << std::endl;
  std::cout << "Histo2DEffNameMisTagPsi2S: " << Histo2DEffNameMisTagPsi2S << std::endl;

  std::cout << "\nHisto3DEffNameMisTagSig: " << Histo3DEffNameMisTagSig << std::endl;
  std::cout << "Histo3DEffNameMisTagJPsi: "  << Histo3DEffNameMisTagJPsi << std::endl;
  std::cout << "Histo3DEffNameMisTagPsi2S: " << Histo3DEffNameMisTagPsi2S << std::endl;

  std::cout << "\n@@@@@@ Utils class settings : public  @@@@@@" << std::endl;
  std::cout << "NcoeffThetaL: "              << NcoeffThetaL << std::endl;
  std::cout << "NcoeffThetaK: "              << NcoeffThetaK << std::endl;
  std::cout << "NcoeffPhi: "                 << NcoeffPhi << std::endl;
  std::cout << "RIGHTflavorTAG: "            << RIGHTflavorTAG << std::endl;
  std::cout << "B0ToKstMuMu: "               << B0ToKstMuMu << std::endl;
  std::cout << "B0ToJPsiKst: "               << B0ToJPsiKst << std::endl;
  std::cout << "B0ToPsi2SKst: "              << B0ToPsi2SKst << std::endl;
  std::cout << "nFitParam: "                 << nFitParam << std::endl;
  std::cout << "nConfigParam: "              << nConfigParam << std::endl;
  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
  std::cout << "@@@ Consider to double-check values for: @@@" << std::endl;
  std::cout << "- Utils::AddConstraintThetaL" << std::endl;
  std::cout << "- Utils::AddConstraintThetaKThetaL" << std::endl;
  std::cout << "- Utils::AddConstraintThetaKThetaLPhi" << std::endl;
  std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
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
  std::cout << "[Utils::ReadAllBins]\tAll bins from file : " << fileName.c_str() << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("q2"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
    {
      q2Bins->push_back(atof(ParVector[i].c_str()));
      std::cout << "Bin " << i << "\t" << q2Bins->back() << std::endl;
    }


  if (signalType == "goodTag")
    {
      // #######################
      // # Read cosThetaK bins #
      // #######################
      std::cout << "@@@ cos(theta_K) bins from file @@@" << std::endl;
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
      std::cout << "@@@ cos(theta_l) bins from file @@@" << std::endl;
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
      std::cout << "@@@ phi bins from file @@@" << std::endl;
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
      std::cout << "@@@ cos(theta_K) bins from file @@@" << std::endl;
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
      std::cout << "@@@ cos(theta_l) bins from file @@@" << std::endl;
      ParVector.clear();
      ParameterFile->ReadFromFile(ParFileBlockN("thetaLmisTag"),&ParVector);
      for (unsigned int i = 0; i < ParVector.size(); i++)
	{
	  cosThetaLBins->push_back(atof(ParVector[i].c_str()));
	  std::cout << "Bin " << i << "\t" << cosThetaLBins->back() << std::endl;
	}


      // #################
      // # Read phi bins #
      // #################
      std::cout << "@@@ phi bins from file @@@" << std::endl;
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


  ParVector.clear();
  delete ParameterFile;
}

void Utils::Readq2Bins (std::string fileName, std::vector<double>* q2Bins)
{
  std::vector<std::string>* ParVector = NULL;
  ReadParVsq2Bins(fileName,"q2",&ParVector);


  std::cout << "[Utils::Readq2Bins]\tq^2 bins from file : " << fileName.c_str() << std::endl;
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
  std::cout << "[Utils::ReadHLTpaths]\tHLT-trigger table from file : " << fileName.c_str() << std::endl;
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

void Utils::GenEfficiency (effStruct* myEff, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins)
{
  myEff->Num1 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
  myEff->Num2 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
  myEff->Den1 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
  myEff->Den2 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];

  myEff->Err2PoisNum1 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
  myEff->Err2PoisNum2 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
  myEff->Err2PoisDen1 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
  myEff->Err2PoisDen2 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];

  myEff->Err2WeigNum1 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
  myEff->Err2WeigNum2 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
  myEff->Err2WeigDen1 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
  myEff->Err2WeigDen2 = new double[(q2Bins->size()-1) * (cosThetaKBins->size()-1) * (cosThetaLBins->size()-1) * (phiBins->size()-1)];
}

void Utils::InitEfficiency (effValue myEffVal, effStruct myEff, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins)
{
  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	for (unsigned int l = 0; l < phiBins->size()-1; l++)
	  {
	    myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Num1;
	    myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Num2;
	    myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Den1;
	    myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Den2;

	    myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2PoisNum1;
	    myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2PoisNum2;
	    myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2PoisDen1;
	    myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2PoisDen2;

	    myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2WeigNum1;
	    myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2WeigNum2;
	    myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2WeigDen1;
	    myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2WeigDen2;
	  }
}

void Utils::SaveEfficiency (std::string fileName, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, effStruct myEff)
{
  effValue myEffVal;

  std::ofstream out;
  out.open(fileName.c_str(), std::ofstream::out);
  if (out.good() == false)
    {
      std::cout << "[Utils::SaveMatrices]\tError opening file : " << fileName << std::endl;
      exit (EXIT_FAILURE);
    }

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	for (unsigned int l = 0; l < phiBins->size()-1; l++)
	  {
	    myEffVal.Num1 = myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Num2 = myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Den1 = myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Den2 = myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	    myEffVal.Err2PoisNum1 = myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2PoisNum2 = myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2PoisDen1 = myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2PoisDen2 = myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	    myEffVal.Err2WeigNum1 = myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2WeigNum2 = myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2WeigDen1 = myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2WeigDen2 = myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	    out << myEffVal.Num1 << "   " << myEffVal.Num2 << "   " << myEffVal.Den1 << "   " << myEffVal.Den2 << "   ";
	    out << myEffVal.Err2PoisNum1 << "   " << myEffVal.Err2PoisNum2 << "   " << myEffVal.Err2PoisDen1 << "   " << myEffVal.Err2PoisDen2 << "   ";
	    out << myEffVal.Err2WeigNum1 << "   " << myEffVal.Err2WeigNum2 << "   " << myEffVal.Err2WeigDen1 << "   " << myEffVal.Err2WeigDen2 << std::endl;
	  }

  out.close();
}

void Utils::ReadEfficiency (std::string fileName, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, effStruct* myEff)
{
  std::string tmpString;
  effValue myEffVal;

  std::ifstream in;
  in.open(fileName.c_str(), std::ofstream::in);
  if (in.good() == false)
    {
      std::cout << "[Utils::ReadEfficiency]\tError opening file : " << fileName << std::endl;
      exit (EXIT_FAILURE);
    }
  
  GenEfficiency(myEff,q2Bins,cosThetaKBins,cosThetaLBins,phiBins);
  ResetEffValue(&myEffVal,0.0);
  InitEfficiency(myEffVal,*myEff,q2Bins,cosThetaKBins,cosThetaLBins,phiBins);

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	for (unsigned int l = 0; l < phiBins->size()-1; l++)
	  {
	    getline(in,tmpString);
	    std::stringstream rawString(tmpString);
	    rawString >> myEffVal.Num1 >> myEffVal.Num2 >> myEffVal.Den1 >> myEffVal.Den2 >> myEffVal.Err2PoisNum1 >> myEffVal.Err2PoisNum2 >> myEffVal.Err2PoisDen1 >> myEffVal.Err2PoisDen2 >> myEffVal.Err2WeigNum1 >> myEffVal.Err2WeigNum2 >> myEffVal.Err2WeigDen1 >> myEffVal.Err2WeigDen2;

	    myEff->Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Num1;
	    myEff->Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Num2;
	    myEff->Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Den1;
	    myEff->Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Den2;

	    myEff->Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2PoisNum1;
	    myEff->Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2PoisNum2;
	    myEff->Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2PoisDen1;
	    myEff->Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2PoisDen2;

	    myEff->Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2WeigNum1;
	    myEff->Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2WeigNum2;
	    myEff->Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2WeigDen1;
	    myEff->Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] = myEffVal.Err2WeigDen2;
	  }

  in.close();
}

void Utils::GetEffq2Bin (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, unsigned int q2Indx, unsigned int cosThetaKIndx, unsigned int cosThetaMuIndx, unsigned int phiIndx, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;
  
  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  myEffVal.Num1 = myEff.Num1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx];
  myEffVal.Num2 = myEff.Num2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx];
  myEffVal.Den1 = myEff.Den1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx];
  myEffVal.Den2 = myEff.Den2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx];

  myEffVal.Err2PoisNum1 = myEff.Err2PoisNum1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx];
  myEffVal.Err2PoisNum2 = myEff.Err2PoisNum2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx];
  myEffVal.Err2PoisDen1 = myEff.Err2PoisDen1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx];
  myEffVal.Err2PoisDen2 = myEff.Err2PoisDen2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx];

  if (myEff.Err2PoisNum1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx] > 0.0)
    myEffValOrg.Num1 = pow(myEff.Num1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx],2.0) /
      myEff.Err2PoisNum1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx];
  
  if (myEff.Err2PoisNum2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx] > 0.0)
    myEffValOrg.Num2 = pow(myEff.Num2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx],2.0) /
      myEff.Err2PoisNum2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx];
  
  if (myEff.Err2PoisDen1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx] > 0.0)
    myEffValOrg.Den1 = pow(myEff.Den1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx],2.0) /
      myEff.Err2PoisDen1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx];
  
  if (myEff.Err2PoisDen2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx] > 0.0)
    myEffValOrg.Den2 = pow(myEff.Den2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx],2.0) /
      myEff.Err2PoisDen2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx];
  
  myEffVal.Err2WeigNum1 = myEff.Err2WeigNum1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx];
  myEffVal.Err2WeigNum2 = myEff.Err2WeigNum2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx];
  myEffVal.Err2WeigDen1 = myEff.Err2WeigDen1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx];
  myEffVal.Err2WeigDen2 = myEff.Err2WeigDen2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaMuIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKIndx*(q2Bins->size()-1) + q2Indx];

  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && (myEffValOrg.Den1-1) > 0.0 && myEffValOrg.Num2 > 0.0 && (myEffValOrg.Den2-1) > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / pow(myEffVal.Err2PoisNum1 / myEffVal.Num1,2.);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / pow(myEffVal.Err2PoisNum2 / myEffVal.Num2,2.);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / pow(myEffVal.Err2PoisDen1 / myEffVal.Den1,2.);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / pow(myEffVal.Err2PoisDen2 / myEffVal.Den2,2.);

      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((myEffValOrg.Den1 / (myEffValOrg.Den1 - 1) * (1./myEffValOrg.Num1 - 1./myEffValOrg.Den1) + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) +
			    (myEffValOrg.Den2 / (myEffValOrg.Den2 - 1) * (1./myEffValOrg.Num2 - 1./myEffValOrg.Den2) + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

TH2D* Utils::Get2DEffHistoq2Bin (std::string histoName, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int q2Indx, effStruct myEff)
{
  double Eff    = 0.0;
  double EffErr = 0.0;
  double* cosThetaKBins_ = MakeBinning(cosThetaKBins);
  double* cosThetaLBins_ = MakeBinning(cosThetaLBins);

  TH2D* Histo = new TH2D(histoName.c_str(), histoName.c_str(), cosThetaKBins->size()-1, cosThetaKBins_, cosThetaLBins->size()-1, cosThetaLBins_);
  Histo->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
  Histo->GetXaxis()->SetTitleOffset(1.8);
  Histo->SetYTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
  Histo->GetYaxis()->SetTitleOffset(1.8);
  Histo->SetZTitle("Efficiency");


  // ##########################
  // # Read binned efficiency #
  // ##########################
  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
    for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
      {
	IntegrateEffPhi(q2Bins,cosThetaKBins,cosThetaLBins,phiBins,q2Bins->operator[](q2Indx),cosThetaKBins->operator[](j),cosThetaLBins->operator[](k),myEff,&Eff,&EffErr);
	
	if ((Eff != 0.0 ) && (EffErr != 0.0))
	  {
	    Histo->SetBinContent(j+1,k+1,Eff);
	    Histo->SetBinError(j+1,k+1,EffErr);
	  }
      }
  

  return Histo;
}

TH3D* Utils::Get3DEffHistoq2Bin (std::string histoName, std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int q2Indx, effStruct myEff)
{
  double Eff    = 0.0;
  double EffErr = 0.0;
  double* cosThetaKBins_ = MakeBinning(cosThetaKBins);
  double* cosThetaLBins_ = MakeBinning(cosThetaLBins);
  double* phiBins_       = MakeBinning(phiBins);

  TH3D* Histo = new TH3D(histoName.c_str(), histoName.c_str(), cosThetaKBins->size()-1, cosThetaKBins_, cosThetaLBins->size()-1, cosThetaLBins_, phiBins->size()-1, phiBins_);
  Histo->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
  Histo->GetXaxis()->SetTitleOffset(1.8);
  Histo->SetYTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
  Histo->GetYaxis()->SetTitleOffset(1.8);
  Histo->SetZTitle("#phi");
  Histo->GetZaxis()->SetTitleOffset(1.8);


  // ##########################
  // # Read binned efficiency #
  // ##########################
  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
    for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  GetEffq2Bin(q2Bins,cosThetaKBins,cosThetaLBins,q2Indx,j,k,l,myEff,&Eff,&EffErr);

	  if ((Eff != 0.0 ) && (EffErr != 0.0))
	    {
	      Histo->SetBinContent(j+1,k+1,l+1,Eff);
	      Histo->SetBinError(j+1,k+1,l+1,EffErr);
	    }
	}

  
  return Histo;
}

TH2D* Utils::Get2DEffHistoq2Bin (std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, unsigned int q2Indx, int SignalType, bool giveMeOriginal,
				 std::pair <double,double> cosThetaKRange, std::pair <double,double> cosThetaLRange)
{
  std::ifstream inputFile;
  std::stringstream myString;
  double xx, xw, yy, yw, cont, err, tmp;
  std::vector<double> Xbins;
  std::vector<double> Ybins;


  myString.clear(); myString.str("");
  myString << (DirEfficiency != "" ? MakeAnalysisPATH(DirEfficiency).c_str() : "") << GetHisto2DEffName(SignalType) << "_" << q2Indx << ".txt";
  std::cout << "[Utils::Get2DEffHistoq2Bin]\tReading 2D binned efficiency file : " << myString.str().c_str() << std::endl;
  inputFile.open(myString.str().c_str(), std::ifstream::in);
  if (inputFile.good() == false)
    {
      std::cout << "[Utils::Get2DEffHistoq2Bin]\tError opening file : " << myString.str().c_str() << std::endl;
      exit (EXIT_FAILURE);
    }

  // ##################
  // # Reading Y bins #
  // ##################
  inputFile >> xx >> xw >> yy >> yw >> cont >> err;
  tmp = xx;
  while (xx == tmp)
    {
      if ((yy >= cosThetaLRange.first) && (yy < cosThetaLRange.second)) Ybins.push_back(yy);
      inputFile >> xx >> xw >> yy >> yw >> cont >> err;
    }
  if (cosThetaLBins->operator[](cosThetaLBins->size()-1) < cosThetaLRange.second) Ybins.push_back(cosThetaLBins->operator[](cosThetaLBins->size()-1));
  else if (Ybins.back() < cosThetaLRange.second)                                  Ybins.push_back(cosThetaLRange.second);
  inputFile.clear();
  inputFile.seekg(0, std::ios::beg);

  // ##################
  // # Reading X bins #
  // ##################
  inputFile >> xx >> xw >> yy >> yw >> cont >> err;
  tmp = xx;
  while (inputFile.eof() == false)
    {
      if ((xx >= cosThetaKRange.first) && (xx < cosThetaKRange.second)) Xbins.push_back(xx);
      while ((xx == tmp) && (inputFile.eof() == false)) inputFile >> xx >> xw >> yy >> yw >> cont >> err;
      tmp = xx;
    }
  if (cosThetaKBins->operator[](cosThetaKBins->size()-1) < cosThetaKRange.second) Xbins.push_back(cosThetaKBins->operator[](cosThetaKBins->size()-1));
  else if (Xbins.back() < cosThetaKRange.second)                                  Xbins.push_back(cosThetaKRange.second);
  inputFile.clear();
  inputFile.seekg(0, std::ios::beg);


  std::cout << "[Utils::Get2DEffHistoq2Bin]\tNew X-axis binning" << std::endl;
  for (unsigned int j = 0; j < Xbins.size(); j++) std::cout << "bin #" << j << " --> " << Xbins[j] << std::endl;

  std::cout << "[Utils::Get2DEffHistoq2Bin]\tNew Y-axis binning" << std::endl;
  for (unsigned int k = 0; k < Ybins.size(); k++) std::cout << "bin #" << k << " --> " << Ybins[k] << std::endl;


  double* Xbins_ = MakeBinning(&Xbins);
  double* Ybins_ = MakeBinning(&Ybins);

  TH2D* Histo = new TH2D(GetHisto2DEffName(SignalType).c_str(), GetHisto2DEffName(SignalType).c_str(), Xbins.size()-1, Xbins_, Ybins.size()-1, Ybins_);
  Histo->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
  Histo->GetXaxis()->SetTitleOffset(1.8);
  Histo->SetYTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
  Histo->GetYaxis()->SetTitleOffset(1.8);
  Histo->SetZTitle("Efficiency");
  TH2D* Histo_clone = (TH2D*)Histo->Clone();


  // ##########################
  // # Read binned efficiency #
  // ##########################
  unsigned int j = 1;
  while (j <= Xbins.size()-1)
    {
      unsigned int k = 1;
      while (k <= Ybins.size()-1)
	{
	  inputFile >> xx >> xw >> yy >> yw >> cont >> err;

	  if ((xx >= Xbins[0]) && (xx <= Xbins[Xbins.size()-1]) &&
	      (yy >= Ybins[0]) && (yy <= Ybins[Ybins.size()-1]))
	    {
	      Histo->SetBinContent(j,k,cont);
	      Histo->SetBinError(j,k,err);

	      if (RIGHTflavorTAG == true) Histo_clone->SetBinContent(j,k,cont);
	      else                        Histo_clone->SetBinContent(Histo->GetNbinsX()-j+1,Histo->GetNbinsY()-k+1,cont);
	      k++;
	    }
	}
      if (k != 1 ) j++;
    }


  Xbins.clear();
  Ybins.clear();
  inputFile.close();
  if (giveMeOriginal == true)
    {
      delete Histo_clone;
      return Histo;
    }
  else
    {
      delete Histo;
      return Histo_clone;
    }
}

TH3D* Utils::Get3DEffHistoq2Bin (std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int q2Indx, int SignalType, bool giveMeOriginal,
				 std::pair <double,double> cosThetaKRange, std::pair <double,double> cosThetaLRange, std::pair <double,double> phiRange)
{
  std::ifstream inputFile;
  std::stringstream myString;
  double xx, xw, yy, yw, zz, zw, cont, err, tmp;
  std::vector<double> Xbins;
  std::vector<double> Ybins;
  std::vector<double> Zbins;


  myString.clear(); myString.str("");
  myString << MakeAnalysisPATH(DirEfficiency).c_str() << GetHisto3DEffName(SignalType) << "_" << q2Indx << ".txt";
  std::cout << "[Utils::Get3DEffHistoq2Bin]\tReading 3D binned efficiency file : " << myString.str().c_str() << std::endl;
  inputFile.open(myString.str().c_str(), std::ifstream::in);
  if (inputFile.good() == false)
    {
      std::cout << "[Utils::Get3DEffHistoq2Bin]\tError opening file : " << myString.str().c_str() << std::endl;
      exit (EXIT_FAILURE);
    }

  // ##################
  // # Reading Z bins #
  // ##################
  inputFile >> xx >> xw >> yy >> yw >> zz >> zw >> cont >> err;
  tmp = yy;
  while (yy == tmp)
    {
      if ((zz >= phiRange.first) && (zz < phiRange.second)) Zbins.push_back(zz);
      inputFile >> xx >> xw >> yy >> yw >> zz >> zw >> cont >> err;
    }
  if (phiBins->operator[](phiBins->size()-1) < phiRange.second) Zbins.push_back(phiBins->operator[](phiBins->size()-1));
  else if (Zbins.back() < phiRange.second)                      Zbins.push_back(phiRange.second);
  inputFile.clear();
  inputFile.seekg(0, std::ios::beg);

  // ##################
  // # Reading Y bins #
  // ##################
  inputFile >> xx >> xw >> yy >> yw >> zz >> zw >> cont >> err;
  tmp = xx;
  while (xx == tmp)
    {
      if ((yy >= cosThetaLRange.first) && (yy < cosThetaLRange.second)) Ybins.push_back(yy);
      Ybins.push_back(yy);
      inputFile >> xx >> xw >> yy >> yw >> zz >> zw >> cont >> err;
    }
  if (cosThetaLBins->operator[](cosThetaLBins->size()-1) < cosThetaLRange.second) Ybins.push_back(cosThetaLBins->operator[](cosThetaLBins->size()-1));
  else if (Ybins.back() < cosThetaLRange.second)                                  Ybins.push_back(cosThetaLRange.second);
  inputFile.clear();
  inputFile.seekg(0, std::ios::beg);

  // ##################
  // # Reading X bins #
  // ##################
  inputFile >> xx >> xw >> yy >> yw >> zz >> zw >> cont >> err;
  tmp = xx;
  while (inputFile.eof() == false)
    {
      if ((xx >= cosThetaKRange.first) && (xx < cosThetaKRange.second)) Xbins.push_back(xx);
      while ((xx == tmp) && (inputFile.eof() == false)) inputFile >> xx >> xw >> yy >> yw >> zz >> zw >> cont >> err;
      tmp = xx;
    }
  if (cosThetaKBins->operator[](cosThetaKBins->size()-1) < cosThetaKRange.second) Xbins.push_back(cosThetaKBins->operator[](cosThetaKBins->size()-1));
  else if (Xbins.back() < cosThetaKRange.second)                                  Xbins.push_back(cosThetaKRange.second);
  inputFile.clear();
  inputFile.seekg(0, std::ios::beg);


  std::cout << "[Utils::Get3DEffHistoq2Bin]\tNew X-axis binning" << std::endl;
  for (unsigned int j = 0; j < Xbins.size(); j++) std::cout << "bin #" << j << " --> " << Xbins[j] << std::endl;

  std::cout << "[Utils::Get3DEffHistoq2Bin]\tNew Y-axis binning" << std::endl;
  for (unsigned int k = 0; k < Ybins.size(); k++) std::cout << "bin #" << k << " --> " << Ybins[k] << std::endl;

  std::cout << "[Utils::Get3DEffHistoq2Bin]\tNew Z-axis binning" << std::endl;
  for (unsigned int l = 0; l < Zbins.size(); l++) std::cout << "bin #" << l << " --> " << Zbins[l] << std::endl;


  double* Xbins_ = MakeBinning(&Xbins);
  double* Ybins_ = MakeBinning(&Ybins);
  double* Zbins_ = MakeBinning(&Zbins);

  TH3D* Histo = new TH3D(GetHisto3DEffName(SignalType).c_str(), GetHisto3DEffName(SignalType).c_str(), Xbins.size()-1, Xbins_, Ybins.size()-1, Ybins_, Zbins.size()-1, Zbins_);
  Histo->SetXTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
  Histo->GetXaxis()->SetTitleOffset(1.8);
  Histo->SetYTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
  Histo->GetYaxis()->SetTitleOffset(1.8);
  Histo->SetZTitle("#phi");
  Histo->GetZaxis()->SetTitleOffset(1.8);
  TH3D* Histo_clone = (TH3D*)Histo->Clone();


  // ##########################
  // # Read binned efficiency #
  // ##########################
  unsigned int j = 1;
  while (j <= Xbins.size()-1)
    {
      unsigned int k = 1;
      while (k <= Ybins.size()-1)
	{
	  unsigned int l = 1;
	  while (l <= Zbins.size()-1)
	    {
	      inputFile >> xx >> xw >> yy >> yw >> zz >> zw >> cont >> err;

	      if ((xx >= Xbins[0]) && (xx <= Xbins[Xbins.size()-1]) &&
		  (yy >= Ybins[0]) && (yy <= Ybins[Ybins.size()-1]) &&
		  (zz >= Zbins[0]) && (zz <= Zbins[Zbins.size()-1]))
		{
		  Histo->SetBinContent(j,k,l,cont);
		  Histo->SetBinError(j,k,l,err);

		  if (RIGHTflavorTAG == true) Histo_clone->SetBinContent(j,k,l,cont);
		  else                        Histo_clone->SetBinContent(Histo->GetNbinsX()-j+1,Histo->GetNbinsY()-k+1,Histo->GetNbinsZ()-l+1,cont);
		  l++;
		}
	    }
	  if (l != 1 ) k++;
	}
      if (k != 1 ) j++;
    }


  Xbins.clear();
  Ybins.clear();
  Zbins.clear();
  inputFile.close();
  if (giveMeOriginal == true)
    {
      delete Histo_clone;
      return Histo;
    }
  else
    {
      delete Histo;
      return Histo_clone;
    }
}

void Utils::Put2DEffHistoq2Bin (std::string fileName, TH2D* histo)
{
  std::ofstream outputFile;
  TAxis* XAxis = histo->GetXaxis();
  TAxis* YAxis = histo->GetYaxis();


  outputFile.open(fileName.c_str(), std::ofstream::out);
  if (outputFile.good() == false)
    {
      std::cout << "[Utils::Put2DEffHistoq2Bin]\tError opening file : " << fileName.c_str() << std::endl;
      exit (EXIT_FAILURE);
    }
  
  
  // ##########################
  // # Save binned efficiency #
  // ##########################
  for (int i = 1; i <= histo->GetNbinsX(); i++)
    for (int j = 1; j <= histo->GetNbinsY(); j++)
      outputFile << XAxis->GetBinLowEdge(i) << "   " << XAxis->GetBinWidth(i) << "   " << YAxis->GetBinLowEdge(j) << "   " << YAxis->GetBinWidth(j) << "   " << histo->GetBinContent(i,j) << "   " << histo->GetBinError(i,j) << std::endl;
  
  
  outputFile.close();
}

void Utils::Put3DEffHistoq2Bin (std::string fileName, TH3D* histo)
{
  std::ofstream outputFile;
  TAxis* XAxis = histo->GetXaxis();
  TAxis* YAxis = histo->GetYaxis();
  TAxis* ZAxis = histo->GetZaxis();


  outputFile.open(fileName.c_str(), std::ofstream::out);
  if (outputFile.good() == false)
    {
      std::cout << "[Utils::Put3DEffHistoq2Bin]\tError opening file : " << fileName.c_str() << std::endl;
      exit (EXIT_FAILURE);
    }
  
  
  // ##########################
  // # Save binned efficiency #
  // ##########################
  for (int i = 1; i <= histo->GetNbinsX(); i++)
    for (int j = 1; j <= histo->GetNbinsY(); j++)
      for (int k = 1; k <= histo->GetNbinsZ(); k++)
	outputFile << XAxis->GetBinLowEdge(i) << "   " << XAxis->GetBinWidth(i) << "   " << YAxis->GetBinLowEdge(j) << "   " << YAxis->GetBinWidth(j) << "   " << ZAxis->GetBinLowEdge(k) << "   " << ZAxis->GetBinWidth(k) << "   " << histo->GetBinContent(i,j,k) << "   " << histo->GetBinError(i,j,k) << std::endl;
  

  outputFile.close();
}

void Utils::DeleteEfficiency (effStruct myEff)
{
  delete myEff.Num1;
  delete myEff.Num2;
  delete myEff.Den1;
  delete myEff.Den2;
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

void Utils::IntegrateEffAll (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;

  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	for (unsigned int l = 0; l < phiBins->size()-1; l++)
	  {
	    myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	    myEffVal.Err2PoisNum1 = myEffVal.Err2PoisNum1 + myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2PoisNum2 = myEffVal.Err2PoisNum2 + myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2PoisDen1 = myEffVal.Err2PoisDen1 + myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2PoisDen2 = myEffVal.Err2PoisDen2 + myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	    if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	      myEffValOrg.Num1 = myEffValOrg.Num1 + pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	    if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	      myEffValOrg.Num2 = myEffValOrg.Num2 + pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	    if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	      myEffValOrg.Den1 = myEffValOrg.Den1 + pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	    if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	      myEffValOrg.Den2 = myEffValOrg.Den2 + pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	    myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 + myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 + myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 + myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 + myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  }
  
  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && (myEffValOrg.Den1-1) > 0.0 && myEffValOrg.Num2 > 0.0 && (myEffValOrg.Den2-1) > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / pow(myEffVal.Err2PoisNum1 / myEffVal.Num1,2.);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / pow(myEffVal.Err2PoisNum2 / myEffVal.Num2,2.);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / pow(myEffVal.Err2PoisDen1 / myEffVal.Den1,2.);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / pow(myEffVal.Err2PoisDen2 / myEffVal.Den2,2.);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((myEffValOrg.Den1 / (myEffValOrg.Den1 - 1) * (1./myEffValOrg.Num1 - 1./myEffValOrg.Den1) + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) +
			    (myEffValOrg.Den2 / (myEffValOrg.Den2 - 1) * (1./myEffValOrg.Num2 - 1./myEffValOrg.Den2) + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffButPsi (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;

  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	for (unsigned int l = 0; l < phiBins->size()-1; l++)
	  if (!(((q2Bins->operator[](i) < (JPsiMass*JPsiMass)) && (q2Bins->operator[](i+1) > (JPsiMass*JPsiMass))) ||
		((q2Bins->operator[](i) < (PsiPMass*PsiPMass)) && (q2Bins->operator[](i+1) > (PsiPMass*PsiPMass)))))
	    {
	      myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      myEffVal.Err2PoisNum1 = myEffVal.Err2PoisNum1 + myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2PoisNum2 = myEffVal.Err2PoisNum2 + myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2PoisDen1 = myEffVal.Err2PoisDen1 + myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2PoisDen2 = myEffVal.Err2PoisDen2 + myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Num1 = myEffValOrg.Num1 + pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	      if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Num2 = myEffValOrg.Num2 + pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Den1 = myEffValOrg.Den1 + pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	      if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Den2 = myEffValOrg.Den2 + pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 + myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 + myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 + myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 + myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    }
  
  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && (myEffValOrg.Den1-1) > 0.0 && myEffValOrg.Num2 > 0.0 && (myEffValOrg.Den2-1) > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / pow(myEffVal.Err2PoisNum1 / myEffVal.Num1,2.);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / pow(myEffVal.Err2PoisNum2 / myEffVal.Num2,2.);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / pow(myEffVal.Err2PoisDen1 / myEffVal.Den1,2.);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / pow(myEffVal.Err2PoisDen2 / myEffVal.Den2,2.);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((myEffValOrg.Den1 / (myEffValOrg.Den1 - 1) * (1./myEffValOrg.Num1 - 1./myEffValOrg.Den1) + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) +
			    (myEffValOrg.Den2 / (myEffValOrg.Den2 - 1) * (1./myEffValOrg.Num2 - 1./myEffValOrg.Den2) + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffInJPsi (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;

  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	for (unsigned int l = 0; l < phiBins->size()-1; l++)
	  if ((q2Bins->operator[](i) < (JPsiMass*JPsiMass)) && (q2Bins->operator[](i+1) > (JPsiMass*JPsiMass)))
	    {
	      myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      myEffVal.Err2PoisNum1 = myEffVal.Err2PoisNum1 + myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2PoisNum2 = myEffVal.Err2PoisNum2 + myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2PoisDen1 = myEffVal.Err2PoisDen1 + myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2PoisDen2 = myEffVal.Err2PoisDen2 + myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Num1 = myEffValOrg.Num1 + pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	      if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Num2 = myEffValOrg.Num2 + pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Den1 = myEffValOrg.Den1 + pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	      if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Den2 = myEffValOrg.Den2 + pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 + myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 + myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 + myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 + myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    }
  
  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && (myEffValOrg.Den1-1) > 0.0 && myEffValOrg.Num2 > 0.0 && (myEffValOrg.Den2-1) > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / pow(myEffVal.Err2PoisNum1 / myEffVal.Num1,2.);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / pow(myEffVal.Err2PoisNum2 / myEffVal.Num2,2.);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / pow(myEffVal.Err2PoisDen1 / myEffVal.Den1,2.);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / pow(myEffVal.Err2PoisDen2 / myEffVal.Den2,2.);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((myEffValOrg.Den1 / (myEffValOrg.Den1 - 1) * (1./myEffValOrg.Num1 - 1./myEffValOrg.Den1) + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) +
			    (myEffValOrg.Den2 / (myEffValOrg.Den2 - 1) * (1./myEffValOrg.Num2 - 1./myEffValOrg.Den2) + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffInPsiP (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;

  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	for (unsigned int l = 0; l < phiBins->size()-1; l++)
	  if ((q2Bins->operator[](i) < (PsiPMass*PsiPMass)) && (q2Bins->operator[](i+1) > (PsiPMass*PsiPMass)))
	    {
	      myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      myEffVal.Err2PoisNum1 = myEffVal.Err2PoisNum1 + myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2PoisNum2 = myEffVal.Err2PoisNum2 + myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2PoisDen1 = myEffVal.Err2PoisDen1 + myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2PoisDen2 = myEffVal.Err2PoisDen2 + myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Num1 = myEffValOrg.Num1 + pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	      if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Num2 = myEffValOrg.Num2 + pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Den1 = myEffValOrg.Den1 + pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	      if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
		myEffValOrg.Den2 = myEffValOrg.Den2 + pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
		  myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 + myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 + myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 + myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 + myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    }

  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && (myEffValOrg.Den1-1) > 0.0 && myEffValOrg.Num2 > 0.0 && (myEffValOrg.Den2-1) > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / pow(myEffVal.Err2PoisNum1 / myEffVal.Num1,2.);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / pow(myEffVal.Err2PoisNum2 / myEffVal.Num2,2.);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / pow(myEffVal.Err2PoisDen1 / myEffVal.Den1,2.);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / pow(myEffVal.Err2PoisDen2 / myEffVal.Den2,2.);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((myEffValOrg.Den1 / (myEffValOrg.Den1 - 1) * (1./myEffValOrg.Num1 - 1./myEffValOrg.Den1) + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) +
			    (myEffValOrg.Den2 / (myEffValOrg.Den2 - 1) * (1./myEffValOrg.Num2 - 1./myEffValOrg.Den2) + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffButq2 (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int q2Indx, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;

  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
    for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	  myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	  myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	  myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	  myEffVal.Err2PoisNum1 = myEffVal.Err2PoisNum1 + myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	  myEffVal.Err2PoisNum2 = myEffVal.Err2PoisNum2 + myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	  myEffVal.Err2PoisDen1 = myEffVal.Err2PoisDen1 + myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	  myEffVal.Err2PoisDen2 = myEffVal.Err2PoisDen2 + myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	  if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	    myEffValOrg.Num1 = myEffValOrg.Num1 + pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	      myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	    
	  if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	    myEffValOrg.Num2 = myEffValOrg.Num2 + pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	      myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	  if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	    myEffValOrg.Den1 = myEffValOrg.Den1 + pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	      myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	    
	  if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	    myEffValOrg.Den2 = myEffValOrg.Den2 + pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	      myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	  myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 + myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	  myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 + myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	  myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 + myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	  myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 + myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	}
  
  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && (myEffValOrg.Den1-1) > 0.0 && myEffValOrg.Num2 > 0.0 && (myEffValOrg.Den2-1) > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / pow(myEffVal.Err2PoisNum1 / myEffVal.Num1,2.);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / pow(myEffVal.Err2PoisNum2 / myEffVal.Num2,2.);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / pow(myEffVal.Err2PoisDen1 / myEffVal.Den1,2.);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / pow(myEffVal.Err2PoisDen2 / myEffVal.Den2,2.);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((myEffValOrg.Den1 / (myEffValOrg.Den1 - 1) * (1./myEffValOrg.Num1 - 1./myEffValOrg.Den1) + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) +
			    (myEffValOrg.Den2 / (myEffValOrg.Den2 - 1) * (1./myEffValOrg.Num2 - 1./myEffValOrg.Den2) + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffButCosThetaK (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int cosThetaKBinIndx, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;

  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	  myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	  myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	  myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];

	  myEffVal.Err2PoisNum1 = myEffVal.Err2PoisNum1 + myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	  myEffVal.Err2PoisNum2 = myEffVal.Err2PoisNum2 + myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	  myEffVal.Err2PoisDen1 = myEffVal.Err2PoisDen1 + myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	  myEffVal.Err2PoisDen2 = myEffVal.Err2PoisDen2 + myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];

	  if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Num1 = myEffValOrg.Num1 + pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	    
	  if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Num2 = myEffValOrg.Num2 + pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];

	  if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Den1 = myEffValOrg.Den1 + pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	    
	  if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Den2 = myEffValOrg.Den2 + pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];

	  myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 + myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	  myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 + myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	  myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 + myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	  myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 + myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + i];
	}
  
  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && (myEffValOrg.Den1-1) > 0.0 && myEffValOrg.Num2 > 0.0 && (myEffValOrg.Den2-1) > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / pow(myEffVal.Err2PoisNum1 / myEffVal.Num1,2.);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / pow(myEffVal.Err2PoisNum2 / myEffVal.Num2,2.);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / pow(myEffVal.Err2PoisDen1 / myEffVal.Den1,2.);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / pow(myEffVal.Err2PoisDen2 / myEffVal.Den2,2.);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((myEffValOrg.Den1 / (myEffValOrg.Den1 - 1) * (1./myEffValOrg.Num1 - 1./myEffValOrg.Den1) + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) +
			    (myEffValOrg.Den2 / (myEffValOrg.Den2 - 1) * (1./myEffValOrg.Num2 - 1./myEffValOrg.Den2) + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffButCosThetaL (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int cosThetaLBinIndx, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;

  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	  myEffVal.Err2PoisNum1 = myEffVal.Err2PoisNum1 + myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Err2PoisNum2 = myEffVal.Err2PoisNum2 + myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Err2PoisDen1 = myEffVal.Err2PoisDen1 + myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Err2PoisDen2 = myEffVal.Err2PoisDen2 + myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	  if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Num1 = myEffValOrg.Num1 + pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	  if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Num2 = myEffValOrg.Num2 + pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	  if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Den1 = myEffValOrg.Den1 + pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	  if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Den2 = myEffValOrg.Den2 + pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	  myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 + myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 + myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 + myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 + myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	}
  
  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && (myEffValOrg.Den1-1) > 0.0 && myEffValOrg.Num2 > 0.0 && (myEffValOrg.Den2-1) > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / pow(myEffVal.Err2PoisNum1 / myEffVal.Num1,2.);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / pow(myEffVal.Err2PoisNum2 / myEffVal.Num2,2.);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / pow(myEffVal.Err2PoisDen1 / myEffVal.Den1,2.);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / pow(myEffVal.Err2PoisDen2 / myEffVal.Den2,2.);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((myEffValOrg.Den1 / (myEffValOrg.Den1 - 1) * (1./myEffValOrg.Num1 - 1./myEffValOrg.Den1) + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) +
			    (myEffValOrg.Den2 / (myEffValOrg.Den2 - 1) * (1./myEffValOrg.Num2 - 1./myEffValOrg.Den2) + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffButPhi (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, unsigned int phiIndx, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;

  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  for (unsigned int i = 0; i < q2Bins->size()-1; i++)
    for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
      for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
	{
	  myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	  myEffVal.Err2PoisNum1 = myEffVal.Err2PoisNum1 + myEff.Err2PoisNum1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Err2PoisNum2 = myEffVal.Err2PoisNum2 + myEff.Err2PoisNum2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Err2PoisDen1 = myEffVal.Err2PoisDen1 + myEff.Err2PoisDen1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Err2PoisDen2 = myEffVal.Err2PoisDen2 + myEff.Err2PoisDen2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	  if (myEff.Err2PoisNum1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Num1 = myEffValOrg.Num1 + pow(myEff.Num1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisNum1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	  if (myEff.Err2PoisNum2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Num2 = myEffValOrg.Num2 + pow(myEff.Num2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisNum2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	  if (myEff.Err2PoisDen1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Den1 = myEffValOrg.Den1 + pow(myEff.Den1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisDen1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	    
	  if (myEff.Err2PoisDen2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i] > 0.0)
	    myEffValOrg.Den2 = myEffValOrg.Den2 + pow(myEff.Den2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i],2.0) /
	      myEff.Err2PoisDen2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];

	  myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 + myEff.Err2WeigNum1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 + myEff.Err2WeigNum2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 + myEff.Err2WeigDen1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	  myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 + myEff.Err2WeigDen2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + i];
	}

  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && (myEffValOrg.Den1-1) > 0.0 && myEffValOrg.Num2 > 0.0 && (myEffValOrg.Den2-1) > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / pow(myEffVal.Err2PoisNum1 / myEffVal.Num1,2.);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / pow(myEffVal.Err2PoisNum2 / myEffVal.Num2,2.);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / pow(myEffVal.Err2PoisDen1 / myEffVal.Den1,2.);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / pow(myEffVal.Err2PoisDen2 / myEffVal.Den2,2.);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((myEffValOrg.Den1 / (myEffValOrg.Den1 - 1) * (1./myEffValOrg.Num1 - 1./myEffValOrg.Den1) + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) +
			    (myEffValOrg.Den2 / (myEffValOrg.Den2 - 1) * (1./myEffValOrg.Num2 - 1./myEffValOrg.Den2) + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffPhiCosThetaL (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int q2Indx, unsigned int cosThetaKBinIndx, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;
  
  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
    for (unsigned int l = 0; l < phiBins->size()-1; l++)
      {
	myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];

	myEffVal.Err2PoisNum1 = myEffVal.Err2PoisNum1 + myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2PoisNum2 = myEffVal.Err2PoisNum2 + myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2PoisDen1 = myEffVal.Err2PoisDen1 + myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2PoisDen2 = myEffVal.Err2PoisDen2 + myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];

	if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Num1 = myEffValOrg.Num1 + pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	    
	if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Num2 = myEffValOrg.Num2 + pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];

	if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Den1 = myEffValOrg.Den1 + pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	    
	if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Den2 = myEffValOrg.Den2 + pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];

	myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 + myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 + myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 + myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 + myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + q2Indx];
      }

  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && (myEffValOrg.Den1-1) > 0.0 && myEffValOrg.Num2 > 0.0 && (myEffValOrg.Den2-1) > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / pow(myEffVal.Err2PoisNum1 / myEffVal.Num1,2.);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / pow(myEffVal.Err2PoisNum2 / myEffVal.Num2,2.);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / pow(myEffVal.Err2PoisDen1 / myEffVal.Den1,2.);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / pow(myEffVal.Err2PoisDen2 / myEffVal.Den2,2.);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((myEffValOrg.Den1 / (myEffValOrg.Den1 - 1) * (1./myEffValOrg.Num1 - 1./myEffValOrg.Den1) + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) +
			    (myEffValOrg.Den2 / (myEffValOrg.Den2 - 1) * (1./myEffValOrg.Num2 - 1./myEffValOrg.Den2) + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffPhiCosThetaK (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, unsigned int q2Indx, unsigned int cosThetaLBinIndx, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;
  
  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);
  
  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
    for (unsigned int l = 0; l < phiBins->size()-1; l++)
      {
	myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	myEffVal.Err2PoisNum1 = myEffVal.Err2PoisNum1 + myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2PoisNum2 = myEffVal.Err2PoisNum2 + myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2PoisDen1 = myEffVal.Err2PoisDen1 + myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2PoisDen2 = myEffVal.Err2PoisDen2 + myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Num1 = myEffValOrg.Num1 + pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	    
	if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Num2 = myEffValOrg.Num2 + pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Den1 = myEffValOrg.Den1 + pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	    
	if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Den2 = myEffValOrg.Den2 + pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 + myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 + myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 + myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 + myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
      }
  
  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && (myEffValOrg.Den1-1) > 0.0 && myEffValOrg.Num2 > 0.0 && (myEffValOrg.Den2-1) > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / pow(myEffVal.Err2PoisNum1 / myEffVal.Num1,2.);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / pow(myEffVal.Err2PoisNum2 / myEffVal.Num2,2.);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / pow(myEffVal.Err2PoisDen1 / myEffVal.Den1,2.);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / pow(myEffVal.Err2PoisDen2 / myEffVal.Den2,2.);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((myEffValOrg.Den1 / (myEffValOrg.Den1 - 1) * (1./myEffValOrg.Num1 - 1./myEffValOrg.Den1) + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) +
			    (myEffValOrg.Den2 / (myEffValOrg.Den2 - 1) * (1./myEffValOrg.Num2 - 1./myEffValOrg.Den2) + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

void Utils::IntegrateEffCosThetaKCosThetaL (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, unsigned int q2Indx, unsigned int phiIndx, effStruct myEff, double* Eff, double* EffErr)
{
  effValue myEffVal;
  effValue myEffValOrg;
  
  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);
  
  for (unsigned int j = 0; j < cosThetaKBins->size()-1; j++)
    for (unsigned int k = 0; k < cosThetaLBins->size()-1; k++)
      {
	myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	myEffVal.Err2PoisNum1 = myEffVal.Err2PoisNum1 + myEff.Err2PoisNum1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2PoisNum2 = myEffVal.Err2PoisNum2 + myEff.Err2PoisNum2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2PoisDen1 = myEffVal.Err2PoisDen1 + myEff.Err2PoisDen1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2PoisDen2 = myEffVal.Err2PoisDen2 + myEff.Err2PoisDen2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	if (myEff.Err2PoisNum1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Num1 = myEffValOrg.Num1 + pow(myEff.Num1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisNum1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	    
	if (myEff.Err2PoisNum2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Num2 = myEffValOrg.Num2 + pow(myEff.Num2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisNum2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	if (myEff.Err2PoisDen1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Den1 = myEffValOrg.Den1 + pow(myEff.Den1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisDen1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	    
	if (myEff.Err2PoisDen2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx] > 0.0)
	  myEffValOrg.Den2 = myEffValOrg.Den2 + pow(myEff.Den2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx],2.0) /
	    myEff.Err2PoisDen2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];

	myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 + myEff.Err2WeigNum1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 + myEff.Err2WeigNum2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 + myEff.Err2WeigDen1[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
	myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 + myEff.Err2WeigDen2[phiIndx*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + k*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + j*(q2Bins->size()-1) + q2Indx];
      }
  
  if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && (myEffValOrg.Den1-1) > 0.0 && myEffValOrg.Num2 > 0.0 && (myEffValOrg.Den2-1) > 0.0)
    {
      myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / pow(myEffVal.Err2PoisNum1 / myEffVal.Num1,2.);
      myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / pow(myEffVal.Err2PoisNum2 / myEffVal.Num2,2.);
      myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / pow(myEffVal.Err2PoisDen1 / myEffVal.Den1,2.);
      myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / pow(myEffVal.Err2PoisDen2 / myEffVal.Den2,2.);
      
      *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
      *EffErr = *Eff * sqrt((myEffValOrg.Den1 / (myEffValOrg.Den1 - 1) * (1./myEffValOrg.Num1 - 1./myEffValOrg.Den1) + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) +
			    (myEffValOrg.Den2 / (myEffValOrg.Den2 - 1) * (1./myEffValOrg.Num2 - 1./myEffValOrg.Den2) + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
    }
  else
    {
      *Eff    = 0.0;
      *EffErr = 0.0;
    }
}

bool Utils::IntegrateEffPhi (std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, double mumuMass2, double cosThetaK, double cosThetaMu, effStruct myEff, double* Eff, double* EffErr)
{
  int mumuq2Indx;
  int cosThetaKBinIndx;
  int cosThetaLBinIndx;
  
  effValue myEffVal;
  effValue myEffValOrg;
  
  ResetEffValue(&myEffVal,0.0);
  ResetEffValue(&myEffValOrg,0.0);

  mumuq2Indx       = SearchBin(mumuMass2,q2Bins);
  cosThetaKBinIndx = SearchBin(cosThetaK,cosThetaKBins);
  cosThetaLBinIndx = SearchBin(cosThetaMu,cosThetaLBins);

  if ((mumuq2Indx != -1) && (cosThetaKBinIndx != -1) && (cosThetaLBinIndx != -1))
    {
      mumuq2Indx       = SearchBin(mumuMass2,q2Bins);
      cosThetaKBinIndx = SearchBin(cosThetaK,cosThetaKBins);
      cosThetaLBinIndx = SearchBin(cosThetaMu,cosThetaLBins);

      for (unsigned int l = 0; l < phiBins->size()-1; l++)
	{
	  myEffVal.Num1 = myEffVal.Num1 + myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
						     cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
						     cosThetaKBinIndx*(q2Bins->size()-1) +
						     mumuq2Indx];
	  myEffVal.Num2 = myEffVal.Num2 + myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
						     cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
						     cosThetaKBinIndx*(q2Bins->size()-1) +
						     mumuq2Indx];
	  myEffVal.Den1 = myEffVal.Den1 + myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
						     cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
						     cosThetaKBinIndx*(q2Bins->size()-1) +
						     mumuq2Indx];
	  myEffVal.Den2 = myEffVal.Den2 + myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
						     cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
						     cosThetaKBinIndx*(q2Bins->size()-1) +
						     mumuq2Indx];

	  myEffVal.Err2PoisNum1 = myEffVal.Err2PoisNum1 + myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
									     cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
									     cosThetaKBinIndx*(q2Bins->size()-1) +
									     mumuq2Indx];
	  myEffVal.Err2PoisNum2 = myEffVal.Err2PoisNum2 + myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
									     cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
									     cosThetaKBinIndx*(q2Bins->size()-1) +
									     mumuq2Indx];
	  myEffVal.Err2PoisDen1 = myEffVal.Err2PoisDen1 + myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
									     cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
									     cosThetaKBinIndx*(q2Bins->size()-1) +
									     mumuq2Indx];
	  myEffVal.Err2PoisDen2 = myEffVal.Err2PoisDen2 + myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
									     cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
									     cosThetaKBinIndx*(q2Bins->size()-1) +
									     mumuq2Indx];

	  if (myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaKBinIndx*(q2Bins->size()-1) +
				 mumuq2Indx] > 0.0)
	    myEffValOrg.Num1 = myEffValOrg.Num1 + pow(myEff.Num1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
								 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
								 cosThetaKBinIndx*(q2Bins->size()-1) +
								 mumuq2Indx],2.0) /
	      myEff.Err2PoisNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaKBinIndx*(q2Bins->size()-1) +
				 mumuq2Indx];
	    
	  if (myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaKBinIndx*(q2Bins->size()-1) +
				 mumuq2Indx] > 0.0)
	    myEffValOrg.Num2 = myEffValOrg.Num2 + pow(myEff.Num2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
								 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
								 cosThetaKBinIndx*(q2Bins->size()-1) +
								 mumuq2Indx],2.0) /
	      myEff.Err2PoisNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaKBinIndx*(q2Bins->size()-1) +
				 mumuq2Indx];

	  if (myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaKBinIndx*(q2Bins->size()-1) +
				 mumuq2Indx] > 0.0)
	    myEffValOrg.Den1 = myEffValOrg.Den1 + pow(myEff.Den1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
								 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
								 cosThetaKBinIndx*(q2Bins->size()-1) +
								 mumuq2Indx],2.0) /
	      myEff.Err2PoisDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaKBinIndx*(q2Bins->size()-1) +
				 mumuq2Indx];
	    
	  if (myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaKBinIndx*(q2Bins->size()-1) +
				 mumuq2Indx] > 0.0)
	    myEffValOrg.Den2 = myEffValOrg.Den2 + pow(myEff.Den2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
								 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
								 cosThetaKBinIndx*(q2Bins->size()-1) +
								 mumuq2Indx],2.0) /
	      myEff.Err2PoisDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) +
				 cosThetaKBinIndx*(q2Bins->size()-1) +
				 mumuq2Indx];

	  myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 + myEff.Err2WeigNum1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + mumuq2Indx];
	  myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 + myEff.Err2WeigNum2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + mumuq2Indx];
	  myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 + myEff.Err2WeigDen1[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + mumuq2Indx];
	  myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 + myEff.Err2WeigDen2[l*(cosThetaLBins->size()-1)*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaLBinIndx*(cosThetaKBins->size()-1)*(q2Bins->size()-1) + cosThetaKBinIndx*(q2Bins->size()-1) + mumuq2Indx];
	}
      
      if (myEffVal.Num1 > 0.0 && myEffVal.Den1 > 0.0 && myEffVal.Num2 > 0.0 && myEffVal.Den2 > 0.0 && myEffValOrg.Num1 > 0.0 && (myEffValOrg.Den1-1) > 0.0 && myEffValOrg.Num2 > 0.0 && (myEffValOrg.Den2-1) > 0.0)
	{
	  myEffVal.Err2WeigNum1 = myEffVal.Err2WeigNum1 / pow(myEffVal.Err2PoisNum1 / myEffVal.Num1,2.);
	  myEffVal.Err2WeigNum2 = myEffVal.Err2WeigNum2 / pow(myEffVal.Err2PoisNum2 / myEffVal.Num2,2.);
	  myEffVal.Err2WeigDen1 = myEffVal.Err2WeigDen1 / pow(myEffVal.Err2PoisDen1 / myEffVal.Den1,2.);
	  myEffVal.Err2WeigDen2 = myEffVal.Err2WeigDen2 / pow(myEffVal.Err2PoisDen2 / myEffVal.Den2,2.);
	  
	  *Eff    = myEffVal.Num1*myEffVal.Num2 / (myEffVal.Den1*myEffVal.Den2);
	  *EffErr = *Eff * sqrt((myEffValOrg.Den1 / (myEffValOrg.Den1 - 1) * (1./myEffValOrg.Num1 - 1./myEffValOrg.Den1) + myEffVal.Err2WeigNum1 + myEffVal.Err2WeigDen1) +
				(myEffValOrg.Den2 / (myEffValOrg.Den2 - 1) * (1./myEffValOrg.Num2 - 1./myEffValOrg.Den2) + myEffVal.Err2WeigNum2 + myEffVal.Err2WeigDen2));
	}
      else
	{
	  *Eff    = 0.0;
	  *EffErr = 0.0;
	}
      
      return true;
    }
  else return false;
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
  std::cout << "[Utils::ReadLumi]\t@@@ Recorded luminosity: " << val << " fb-1 @@@" << std::endl;


  ParVector.clear();
  delete ParameterFile;
  return val;
}

void Utils::ReadNLLval (std::string fileName, std::vector<std::vector<double>*>* vecParam)
// ######################################
// # vecParam[0] --> Fl                 #
// # vecParam[1] --> Afb                #
// # vecParam[2] --> P1                 #
// # vecParam[3] --> P2                 #
// # vecParam[4] --> Branching-Fraction #
// ######################################
{
  double val;
  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileName.c_str());


  // ###################
  // # Clear vector(s) #
  // ###################
  vecParam->clear();
  

  // ###################
  // # Read NLL values #
  // ###################
  ParameterFile->ReadFromFile(ParFileBlockN("fitNLL"),&ParVector);

  for (unsigned int j = 0; j < nFitObserv*2; j++) vecParam->push_back(new std::vector<double>);
  
  for (unsigned int i = 0; i < ParVector.size(); i=i+nFitObserv)
    {

      std::cout << "Read set-" << static_cast<int>(static_cast<double>(i)/static_cast<double>(nFitObserv)) << " of fit-observable NLL" << std::endl;

      for (unsigned int j = 0; j < nFitObserv; j++)
	{
	  std::stringstream rawString(ParVector[i+j]);
	  rawString >> val;
	  vecParam->operator[](j)->push_back(val);
	  std::cout << "Fit observable-" << j << " NLL: " << vecParam->operator[](j)->back() << std::endl;
	}
    }

  
  ParVector.clear();
  delete ParameterFile;
}

double Utils::GetNLLval (std::vector<std::vector<double>*>* NLLvals, std::string varName, unsigned int q2Indx)
// ###################
// # varName = "Fl"  #
// # varName = "Afb" #
// # varName = "P1"  #
// # varName = "P2"  #
// # varName = "BF"  #
// ###################
{
  if      (varName == "Fl")  return (*NLLvals)[0]->operator[](q2Indx);
  else if (varName == "Afb") return (*NLLvals)[1]->operator[](q2Indx);
  else if (varName == "P1")  return (*NLLvals)[2]->operator[](q2Indx);
  else if (varName == "P2")  return (*NLLvals)[3]->operator[](q2Indx);
  else if (varName == "BF")  return (*NLLvals)[4]->operator[](q2Indx);
  else
    {
      std::cout << "[Utils::GetNLLval]\tNLL parameter not valid : " << varName << std::endl;
      exit (EXIT_FAILURE);
    }
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

      std::cout << "Read trigger path from config file : " << HLTpath.back() << std::endl;
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

      std::cout << "Read set-" << static_cast<int>(static_cast<double>(i)/static_cast<double>(nFitParam+nConfigParam)) << " of fit-parameter starting values" << std::endl;

      for (unsigned int j = 0; j < nFitParam; j++)
	{
	  std::stringstream rawString(ParVector[i+j]);
	  vecParam->operator[](j)->push_back(rawString.str());
	  std::cout << "Fit parameter-" << j << " value: " << vecParam->operator[](j)->back() << std::endl;
	}
      for (unsigned int j = 0; j < nConfigParam; j++)
	{
	  configParam->operator[](j)->push_back(atoi(ParVector[i+nFitParam+j].c_str()));
	  std::cout << "Config. parameter-" << j << " value: " << configParam->operator[](j)->back() << std::endl;
	}
    }

  
  ParVector.clear();
  delete ParameterFile;
}

void Utils::ReadFitSystematics (std::string fileName, std::vector<std::vector<double>*>* vecParam)
// ###########################################
// # vecParam[0] --> Fl +err                 #
// # vecParam[1] --> Fl -err                 #
// # vecParam[2] --> Afb +err                #
// # vecParam[3] --> Afb -err                #
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

      std::cout << "Read set-" << static_cast<int>(static_cast<double>(i)/static_cast<double>(nFitObserv)) << " of fit-observable systematic errors" << std::endl;

      for (unsigned int j = 0; j < nFitObserv; j++)
	{
	  std::stringstream rawString(ParVector[i+j]);
	  vecParam->operator[](j*2)->push_back(0.0);
	  vecParam->operator[](j*2+1)->push_back(0.0);

	  while (true)
	    {
	      val = -1.0;
	      rawString >> val;
	      if (val == -1.0) break;
	      vecParam->operator[](j*2)->back() += val*val;

	      val = -1.0;
	      rawString >> val;
	      if (val == -1.0)
		{
		  std::cout << "[Utils::ReadFitSystematics]\tData format not valid in systematic error" << std::endl;
		  exit (EXIT_FAILURE);
		}
	      vecParam->operator[](j*2+1)->back() += val*val;
	    }

	  vecParam->operator[](j*2)->back()   = sqrt(vecParam->operator[](j*2)->back());
	  vecParam->operator[](j*2+1)->back() = sqrt(vecParam->operator[](j*2+1)->back());
	  std::cout << "[Utils::ReadFitSystematics]\tFit observable-" << j << " systematic error: +" << vecParam->operator[](j*2)->back() << "/-" << vecParam->operator[](j*2+1)->back() << std::endl;
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

  std::cout << "[Utils::ReadParVsq2Bins]\tReading parameters vs q^2 from file : " << fileName << std::endl;
  for (unsigned int i = 0; i < (*vecParam)->size(); i++)
    std::cout << "Parameter value and errors for q2 bin " << i << ": " << (*vecParam)->operator[](i) << std::endl;
  
  delete ParameterFile;
}

void Utils::SaveAnalyticalEff (std::string fileName, TF2* effFunc, double q2Val)
{
  std::ofstream fileOutput;
  fileOutput.open(fileName.c_str(),std::ofstream::app);
  if (fileOutput.good() == false)
    {
      std::cout << "[Utils::SaveAnalyticalEff]\tError opening file : " << fileName << std::endl;
      exit (EXIT_FAILURE);
    }
  
  for (unsigned int k = 0; k < NcoeffThetaL; k++)
    {
      fileOutput << q2Val;
      for (unsigned int j = 0; j < NcoeffThetaK; j++) fileOutput << "   " << effFunc->GetParameter(j+k*NcoeffThetaK) << "   " << effFunc->GetParError(j+k*NcoeffThetaK);
      fileOutput << std::endl;
    }

  fileOutput.close();
}

void Utils::SaveAnalyticalEff (std::string fileName, TF3* effFunc, double q2Val)
{
  std::ofstream fileOutput;
  fileOutput.open(fileName.c_str(),std::ofstream::app);
  if (fileOutput.good() == false)
    {
      std::cout << "[Utils::SaveAnalyticalEff]\tError opening file : " << fileName << std::endl;
      exit (EXIT_FAILURE);
    }
  
  for (unsigned int k = 0; k < NcoeffThetaL; k++)
    {
      fileOutput << q2Val;
      for (unsigned int j = 0; j < NcoeffThetaK; j++) fileOutput << "   " << effFunc->GetParameter(j+k*NcoeffThetaK) << "   " << effFunc->GetParError(j+k*NcoeffThetaK);
      fileOutput << std::endl;
    }

  fileOutput << q2Val;
  for (unsigned int l = 0; l < NcoeffPhi; l++) fileOutput << "   " << effFunc->GetParameter(l+NcoeffThetaK*NcoeffThetaL) << "   " << effFunc->GetParError(l+NcoeffThetaK*NcoeffThetaL);
  fileOutput << std::endl;

  fileOutput.close();
}

void Utils::SaveAnalyticalEffFullCovariance (std::string fileName, TMatrixTSym<double>* covMatrix, double q2Val)
{
  std::ofstream fileOutput;
  fileOutput.open(fileName.c_str(),std::ofstream::app);
  if (fileOutput.good() == false)
    {
      std::cout << "[Utils::SaveAnalyticalEffFullCovariance]\tError opening file : " << fileName << std::endl;
      exit (EXIT_FAILURE);
    }
  
  for (int i = 0; i < covMatrix->GetNrows(); i++)
    {
      fileOutput << q2Val;
      for (int j = 0; j < covMatrix->GetNrows(); j++) fileOutput << "   " << (*covMatrix)(i,j);
      fileOutput << std::endl;
    }

  fileOutput.close();
}

std::string Utils::TellMeEffFuncThetaKThetaLPhi ()
{
  std::string myString = "([0]+[1]*x+[2]*x*x+[3]*x*x*x) + ([4]+[5]*x+[6]*x*x+[7]*x*x*x)*y + ([8]+[9]*x+[10]*x*x+[11]*x*x*x)*y*y + ([12]+[13]*x+[14]*x*x+[15]*x*x*x)*y*y*y + ([16]+[17]*x+[18]*x*x+[19]*x*x*x)*y*y*y*y + ([20]+[21]*x+[22]*x*x+[23]*x*x*x)*y*y*y*y*y +  ([24] + [25]*x*x + [26]*y+[27]*y*y)*z*z";
  std::cout << "[Utils::TellMeEffFuncThetaKThetaLPhi]\tEfficiency shape: " << myString << std::endl;
  return myString;
}

std::string Utils::TellMeEffFuncThetaKThetaL ()
{
  std::string myString = "([0]+[1]*x+[2]*x*x+[3]*x*x*x) + ([4]+[5]*x+[6]*x*x+[7]*x*x*x)*y + ([8]+[9]*x+[10]*x*x+[11]*x*x*x)*y*y + ([12]+[13]*x+[14]*x*x+[15]*x*x*x)*y*y*y + ([16]+[17]*x+[18]*x*x+[19]*x*x*x)*y*y*y*y + ([20]+[21]*x+[22]*x*x+[23]*x*x*x)*y*y*y*y*y";
  std::cout << "[Utils::TellMeEffFuncThetaKThetaL]\tEfficiency shape: " << myString << std::endl;
  return myString;
}

std::string Utils::TellMeEffFuncThetaK ()
{
  std::string myString = "[0] + [1]*x + [2]*x*x + [3]*x*x*x";
  std::cout << "[Utils::TellMeEffFuncThetaK]\tEfficiency shape: " << myString << std::endl;
  return myString;
}

std::string Utils::TellMeEffFuncThetaL ()
{
  std::string myString = "[0]  + [1]*x + [2]*x*x + [3]*x*x*x + [4]*x*x*x*x + [5]*x*x*x*x*x";
  std::cout << "[Utils::TellMeEffFuncThetaL]\tEfficiency shape: " << myString << std::endl;
  return myString;
}

std::string Utils::TellMeEffFuncPhi ()
{
  std::string myString = "[0]  + [1]*x + [2]*x*x + [3]*x*x*x";
  std::cout << "[Utils::TellMeEffFuncPhi]\tEfficiency shape: " << myString << std::endl;
  return myString;
}

void Utils::ReadAnalyticalEff (std::string fileNameEffParams,
			       std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins,
			       std::vector<TF2*>* effFuncs, std::string effID, const unsigned int dataBlockN)
{
  unsigned int indx;
  std::stringstream myString;
  double* coeffVec = new double[2*NcoeffThetaK*NcoeffThetaL];

  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileNameEffParams.c_str());


  // ######################################
  // # Read block of parameters from file #
  // ######################################
  ParameterFile->ReadFromFile(dataBlockN,&ParVector);

  for (unsigned int q2Indx = 0; q2Indx < q2Bins->size()-1; q2Indx++)
    {
      myString.clear(); myString.str("");
      myString << effID << "_" << q2Indx;
      effFuncs->push_back(new TF2(myString.str().c_str(),TellMeEffFuncThetaKThetaL().c_str(),
				  cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),
				  cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1)));
      
      std::cout << "@@@ Reading coefficients for analytical efficiency for set-" << q2Indx << " from file : " << fileNameEffParams.c_str() << " @@@" << std::endl;
      
      for (unsigned int k = 0; k < NcoeffThetaL; k++)
	{
	  std::stringstream rawStringK(ParVector[k+q2Indx*NcoeffThetaL]);
	  rawStringK >> coeffVec[0]; // Discard q2 bin value
	  indx = 0;
	  for (unsigned int j = 0; j < NcoeffThetaK*2; j = j+2)
	    {
	      rawStringK >> coeffVec[indx];
	      effFuncs->operator[](q2Indx)->SetParameter(indx+NcoeffThetaK*k,coeffVec[indx]);
	      rawStringK >> coeffVec[indx];
	      effFuncs->operator[](q2Indx)->SetParError(indx+NcoeffThetaK*k,coeffVec[indx]);

	      std::cout << "Theta_L coef. " << k << " --> reading coef. " << indx << " for var. theta_K: ";
	      std::cout << effFuncs->operator[](q2Indx)->GetParameter(indx+NcoeffThetaK*k) << " +/- ";
	      std::cout << effFuncs->operator[](q2Indx)->GetParError(indx+NcoeffThetaK*k) << std::endl;

	      indx++;
	    }
	}
        
      effFuncs->operator[](q2Indx)->GetXaxis()->SetTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
      effFuncs->operator[](q2Indx)->GetXaxis()->SetTitleOffset(1.8);
      effFuncs->operator[](q2Indx)->GetYaxis()->SetTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
      effFuncs->operator[](q2Indx)->GetYaxis()->SetTitleOffset(1.8);
      effFuncs->operator[](q2Indx)->GetZaxis()->SetTitle("Efficiency");
    }


  delete[] coeffVec;
  ParVector.clear();
  delete ParameterFile;
}

void Utils::ReadAnalyticalEff (std::string fileNameEffParams,
			       std::vector<double>* q2Bins, std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins,
			       std::vector<TF3*>* effFuncs, std::string effID, const unsigned int dataBlockN)
{
  unsigned int indx;
  std::stringstream myString;
  double* coeffVec = new double[2*(NcoeffThetaK*NcoeffThetaL+NcoeffPhi)];

  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileNameEffParams.c_str());


  // ######################################
  // # Read block of parameters from file #
  // ######################################
  ParameterFile->ReadFromFile(dataBlockN,&ParVector);

  for (unsigned int q2Indx = 0; q2Indx < q2Bins->size()-1; q2Indx++)
    {
      myString.clear(); myString.str("");
      myString << effID << "_" << q2Indx;
      effFuncs->push_back(new TF3(myString.str().c_str(),TellMeEffFuncThetaKThetaLPhi().c_str(),
				  cosThetaKBins->operator[](0),cosThetaKBins->operator[](cosThetaKBins->size()-1),
				  cosThetaLBins->operator[](0),cosThetaLBins->operator[](cosThetaLBins->size()-1),
				  phiBins->operator[](0),      phiBins->operator[](cosThetaLBins->size()-1)));
      
      std::cout << "@@@ Reading coefficients for analytical efficiency for set-" << q2Indx << " from file : " << fileNameEffParams.c_str() << " @@@" << std::endl;
      
      for (unsigned int k = 0; k < NcoeffThetaL; k++)
	{
	  std::stringstream rawStringK(ParVector[k+q2Indx*(NcoeffThetaL+1)]);
	  rawStringK >> coeffVec[0]; // Discard q2 bin value
	  indx = 0;
	  for (unsigned int j = 0; j < NcoeffThetaK*2; j = j+2)
	    {
	      rawStringK >> coeffVec[indx];
	      effFuncs->operator[](q2Indx)->SetParameter(indx+NcoeffThetaK*k,coeffVec[indx]);
	      rawStringK >> coeffVec[indx];
	      effFuncs->operator[](q2Indx)->SetParError(indx+NcoeffThetaK*k,coeffVec[indx]);

	      std::cout << "Theta_L coef. " << k << " --> reading coef. " << indx << " for var. theta_K: ";
	      std::cout << effFuncs->operator[](q2Indx)->GetParameter(indx+NcoeffThetaK*k) << " +/- ";
	      std::cout << effFuncs->operator[](q2Indx)->GetParError(indx+NcoeffThetaK*k) << std::endl;

	      indx++;
	    }
	}

      std::stringstream rawStringK(ParVector[NcoeffThetaL+q2Indx*(NcoeffThetaL+1)]);
      rawStringK >> coeffVec[0]; // Discard q2 bin value
      indx = 0;
      for (unsigned int l = 0; l < NcoeffPhi*2; l = l+2)
	{
	  rawStringK >> coeffVec[indx];
	  effFuncs->operator[](q2Indx)->SetParameter(indx+NcoeffThetaK*NcoeffThetaL,coeffVec[indx]);
	  rawStringK >> coeffVec[indx];
	  effFuncs->operator[](q2Indx)->SetParError(indx+NcoeffThetaK*NcoeffThetaL,coeffVec[indx]);
	  
	  std::cout << "Reading coef. " << indx << " for var. phi: ";
	  std::cout << effFuncs->operator[](q2Indx)->GetParameter(indx+NcoeffThetaK*NcoeffThetaL) << " +/- ";
	  std::cout << effFuncs->operator[](q2Indx)->GetParError(indx+NcoeffThetaK*NcoeffThetaL) << std::endl;

	  indx++;
	}
      
      effFuncs->operator[](q2Indx)->GetXaxis()->SetTitle("cos(#theta#lower[-0.4]{_{#font[122]{K}}})");
      effFuncs->operator[](q2Indx)->GetXaxis()->SetTitleOffset(1.8);
      effFuncs->operator[](q2Indx)->GetYaxis()->SetTitle("cos(#theta#lower[-0.4]{_{#font[12]{l}}})");
      effFuncs->operator[](q2Indx)->GetYaxis()->SetTitleOffset(1.8);
      effFuncs->operator[](q2Indx)->GetZaxis()->SetTitle("#phi");
      effFuncs->operator[](q2Indx)->GetZaxis()->SetTitleOffset(1.8);
    }


  delete[] coeffVec;
  ParVector.clear();
  delete ParameterFile;
}

void Utils::ReadAnalyticalEffFullCovariance (std::string fileNameEffParams, std::vector<TMatrixTSym<double>*>* covMatrices, std::string dimensions, const unsigned int dataBlockN)
// ##################################################################
// # dimensions : "2D" (cos(theta_K), cos(theta_l)) efficiency      #
// # dimensions : "3D" (cos(theta_K), cos(theta_l), phi) efficiency #
// ##################################################################
{
  unsigned int Ncoeff;
  if      (dimensions == "2D") Ncoeff = NcoeffThetaK*NcoeffThetaL;
  else if (dimensions == "3D") Ncoeff = NcoeffThetaK*NcoeffThetaL+NcoeffPhi;
  else
    {
      std::cout << "[Utils::ReadAnalyticalEffFullCovariance]\tError wrong parameter name : " << dimensions << std::endl;
      exit (EXIT_FAILURE);
    }

  double* coeffVec = new double[Ncoeff];

  std::vector<std::string> ParVector;
  ReadParameters* ParameterFile = new ReadParameters(fileNameEffParams.c_str());


  // ######################################
  // # Read block of parameters from file #
  // ######################################
  ParameterFile->ReadFromFile(dataBlockN,&ParVector);

  for (int q2Indx = 0; q2Indx < static_cast<int>(rint(ParVector.size()/Ncoeff)); q2Indx++)
    {
      std::cout << "@@@ Reading covariance matrix for analytical efficiency for set-" << q2Indx << " from file : " << fileNameEffParams.c_str() << " @@@" << std::endl;
      
      covMatrices->push_back(new TMatrixTSym<double>(Ncoeff));
      
      for (unsigned int j = 0; j < Ncoeff; j++)
	{
	  std::stringstream rawString(ParVector[j+q2Indx*Ncoeff]);
	  std::cout << "Row #" << j << ": " << rawString.str().c_str() << std::endl;
	  rawString >> coeffVec[0]; // Discard q2 bin value

	  for (unsigned int k = 0; k < Ncoeff; k++)
	    {
	      rawString >> coeffVec[k];
	      (*covMatrices->operator[](q2Indx))[j][k] = coeffVec[k];
	      std::cout << "Covariance [" << j << "][" << k << "] --> " << coeffVec[k] << std::endl;
	    }
	}
    }


  delete[] coeffVec;
  ParVector.clear();
  delete ParameterFile;
}

double Utils::EffMinValue1D (double minX, double maxX, TF1* effFunc)
{
  const unsigned int nsteps = 2000;
  unsigned int iMem = 0;
  double minVal = 0.0;

  // ##############################################################################################
  // # Search for the minimal value of the 1D-function scanning the domain divided in nsteps grid #
  // ##############################################################################################
  for (unsigned int i = 0; i <= nsteps; i++)
    if (effFunc->Eval(maxX - (maxX - minX) / static_cast<double>(nsteps) * static_cast<double>(i)) < minVal)
      {
	minVal = effFunc->Eval(maxX - (maxX - minX) / static_cast<double>(nsteps) * static_cast<double>(i));
	iMem = i;
      }
  
  std::cout << "@@@ Efficiency minimal value (if less than zero): " << minVal << " at: X = " << maxX - (maxX - minX) / static_cast<double>(nsteps) * static_cast<double>(iMem);
  std::cout << " (step along X: " << (maxX - minX) / static_cast<double>(nsteps) << ") @@@" << std::endl;
  
  return minVal;
}

double Utils::EffMinValue2D (std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, TF2* effFunc)
{
  const unsigned int nsteps = 2000;
  unsigned int kMem = 0;
  unsigned int jMem = 0;
  double minVal = 0.0;
  double tmpVal;

  double valj;
  double valk;

  // ############################################################################################################
  // # Search for the minimal value of the 2D-function scanning the domain divided in nsteps*nsteps grid matrix #
  // ############################################################################################################
  for (unsigned int k = 0; k <= nsteps; k++)
    for (unsigned int j = 0; j <= nsteps; j++)
      {
	tmpVal = effFunc->Eval(cosThetaKBins->operator[](cosThetaKBins->size()-1) - (cosThetaKBins->operator[](cosThetaKBins->size()-1) - cosThetaKBins->operator[](0)) / static_cast<double>(nsteps) * static_cast<double>(j),
			       cosThetaLBins->operator[](cosThetaLBins->size()-1) - (cosThetaLBins->operator[](cosThetaLBins->size()-1) - cosThetaLBins->operator[](0)) / static_cast<double>(nsteps) * static_cast<double>(k));
	if (tmpVal < minVal)
	  {
	    minVal = tmpVal;
	    
	    kMem = k;
	    jMem = j;
	  }
      }
  
  valj = cosThetaKBins->operator[](cosThetaKBins->size()-1) - (cosThetaKBins->operator[](cosThetaKBins->size()-1) - cosThetaKBins->operator[](0)) / static_cast<double>(nsteps) * static_cast<double>(jMem);
  valk = cosThetaLBins->operator[](cosThetaLBins->size()-1) - (cosThetaLBins->operator[](cosThetaLBins->size()-1) - cosThetaLBins->operator[](0)) / static_cast<double>(nsteps) * static_cast<double>(kMem);

  std::cout << "[Utils::EffMinValue2D]\t@@@ Efficiency minimal value (if less than zero): " << minVal << " at: (theta_K,theta_l = " << valj << "," << valk << ")" << std::endl;

  std::cout << "Corresponding to bin [#bin convention 0...N-1 (-1 = upper bin)] (theta_K,theta_l): ";
  std::cout << SearchBin(valj,cosThetaKBins) << ",";
  std::cout << SearchBin(valk,cosThetaLBins) << std::endl;

  std::cout << "Grid step along cos(theta_K): ";
  std::cout << (cosThetaKBins->operator[](cosThetaKBins->size()-1) - cosThetaKBins->operator[](0)) / static_cast<double>(nsteps);
  std::cout << "; grid step along cos(theta_l): ";
  std::cout << (cosThetaLBins->operator[](cosThetaLBins->size()-1) - cosThetaLBins->operator[](0)) / static_cast<double>(nsteps) << std::endl;

  return minVal;
}

double Utils::EffMinValue3D (std::vector<double>* cosThetaKBins, std::vector<double>* cosThetaLBins, std::vector<double>* phiBins, TF3* effFunc)
{
  const unsigned int nsteps = 200;
  unsigned int kMem = 0;
  unsigned int jMem = 0;
  unsigned int lMem = 0;
  double minVal = 0.0;
  double tmpVal;

  double valj;
  double valk;
  double vall;

  // ###################################################################################################################
  // # Search for the minimal value of the 3D-function scanning the domain divided in nsteps*nsteps*nsteps grid matrix #
  // ###################################################################################################################
  for (unsigned int k = 0; k <= nsteps; k++)
    for (unsigned int j = 0; j <= nsteps; j++)
      for (unsigned int l = 0; l <= nsteps; l++)
	{
	  tmpVal = effFunc->Eval(cosThetaKBins->operator[](cosThetaKBins->size()-1) - (cosThetaKBins->operator[](cosThetaKBins->size()-1) - cosThetaKBins->operator[](0)) / static_cast<double>(nsteps) * static_cast<double>(j),
				 cosThetaLBins->operator[](cosThetaLBins->size()-1) - (cosThetaLBins->operator[](cosThetaLBins->size()-1) - cosThetaLBins->operator[](0)) / static_cast<double>(nsteps) * static_cast<double>(k),
				 phiBins->operator[](phiBins->size()-1)             - (phiBins->operator[](phiBins->size()-1)             - phiBins->operator[](0))       / static_cast<double>(nsteps) * static_cast<double>(l));
	  if (tmpVal < minVal)
	    {
	      minVal = tmpVal;

	      kMem = k;
	      jMem = j;
	      lMem = l;
	    }
	}

  valj = cosThetaKBins->operator[](cosThetaKBins->size()-1) - (cosThetaKBins->operator[](cosThetaKBins->size()-1) - cosThetaKBins->operator[](0)) / static_cast<double>(nsteps) * static_cast<double>(jMem);
  valk = cosThetaLBins->operator[](cosThetaLBins->size()-1) - (cosThetaLBins->operator[](cosThetaLBins->size()-1) - cosThetaLBins->operator[](0)) / static_cast<double>(nsteps) * static_cast<double>(kMem);
  vall = phiBins->operator[](phiBins->size()-1)             - (phiBins->operator[](phiBins->size()-1)             - phiBins->operator[](0))       / static_cast<double>(nsteps) * static_cast<double>(lMem);

  std::cout << "[Utils::EffMinValue3D]\t@@@ Efficiency minimal value (if less than zero): " << minVal << " at: (theta_K,theta_l,phi = " << valj << "," << valk << "," << vall << ")" << std::endl;

  std::cout << "Corresponding to bin [#bin convention 0...N-1 (-1 = upper bin)] (theta_K,theta_l,phi): ";
  std::cout << SearchBin(valj,cosThetaKBins) << ",";
  std::cout << SearchBin(valk,cosThetaLBins) << ",";
  std::cout << SearchBin(vall,phiBins) << std::endl;

  std::cout << "Grid step along cos(theta_K): ";
  std::cout << (cosThetaKBins->operator[](cosThetaKBins->size()-1) - cosThetaKBins->operator[](0)) / static_cast<double>(nsteps);
  std::cout << "; grid step along cos(theta_l): ";
  std::cout << (cosThetaLBins->operator[](cosThetaLBins->size()-1) - cosThetaLBins->operator[](0)) / static_cast<double>(nsteps);
  std::cout << "; grid step along phi: ";
  std::cout << (phiBins->operator[](phiBins->size()-1)             - phiBins->operator[](0))       / static_cast<double>(nsteps) << std::endl;

  return minVal;
}

void Utils::MakeGraphVar (std::string fileName, TGraphAsymmErrors** graph, std::string varName, double offset)
// ###################
// # varName = "Fl"  #
// # varName = "Afb" #
// # varName = "P1"  #
// # varName = "P2"  #
// # varName = "FS"  #
// # varName = "AS"  #
// # varName = "BF"  #
// ###################
{
  double tmpVar;

  std::vector<std::vector<std::string>*> vecParam;
  std::vector<std::vector<unsigned int>*> configParam;
  std::vector<std::string>* ParVector = NULL;
  std::vector<double> vxs;
  std::vector<double> vys;
  std::vector<double> vxel;
  std::vector<double> vxeh;
  std::vector<double> vyel;
  std::vector<double> vyeh;
  

  // ################################
  // # Read values from config file #
  // ################################
  std::vector<double> q2Bins;
  Readq2Bins(fileName,&q2Bins);
  ReadFitStartingValues(fileName,&vecParam,&configParam,ParFileBlockN("fitValBins"));
  ReadParVsq2Bins(fileName,"BF",&ParVector);


  for (unsigned int i = 0; i < q2Bins.size()-1; i++)
    {
      std::stringstream rawString;
      if      (varName == "BF")  rawString << ParVector->operator[](i);
      else if (varName == "Fl")  rawString << vecParam[GetFitParamIndx("FlS")]->operator[](i);
      else if (varName == "Afb") rawString << vecParam[GetFitParamIndx("AfbS")]->operator[](i);
      else if (varName == "P1")  rawString << vecParam[GetFitParamIndx("P1S")]->operator[](i);
      else if (varName == "P2")  rawString << vecParam[GetFitParamIndx("P2S")]->operator[](i);
      else if (varName == "FS")  rawString << vecParam[GetFitParamIndx("FsS")]->operator[](i);
      else if (varName == "AS")  rawString << vecParam[GetFitParamIndx("AsS")]->operator[](i);
      else { std::cout << "[Utils::MakeGraphVar]\tVariable name unknown: " << varName << std::endl; exit (EXIT_FAILURE); }
      std::cout << "[Utils::MakeGraphVar]\tI've read the raw value: " << rawString.str().c_str() << std::endl;

      if (ValIsInPsi(&q2Bins,(q2Bins[i+1]+q2Bins[i])/2.) == false)
	{
	  vxs.push_back((q2Bins[i+1] + q2Bins[i]) / 2. + offset);
	  vxeh.push_back(q2Bins[i+1] - vxs.back() - offset);
	  vxel.push_back(vxs.back() - q2Bins[i] + offset);
 
	  rawString >> tmpVar;
	  vys.push_back(tmpVar);
	  rawString >> tmpVar;
	  vyel.push_back(fabs(tmpVar));
	  rawString >> tmpVar;
	  vyeh.push_back(fabs(tmpVar));
	}
      else
	{
	  vxs.push_back((q2Bins[i+1] + q2Bins[i]) / 2.);
	  vxeh.push_back(0.0);
	  vxel.push_back(0.0);

	  vys.push_back(YvalueOutsideLimits);
	  vyeh.push_back(0.0);
	  vyel.push_back(0.0);
	}
    }
  *graph = new TGraphAsymmErrors(vxs.size(), &vxs[0], &vys[0], &vxel[0], &vxeh[0], &vyel[0], &vyeh[0]);


  // #################
  // # Clear vectors #
  // #################
  q2Bins.clear();
  ParVector->clear();
  delete ParVector;
  for (unsigned int i = 0; i < vecParam.size(); i++) vecParam[i]->clear();
  vecParam.clear();
  vxs.clear();
  vys.clear();
  vxel.clear();
  vxeh.clear();
  vyel.clear();
  vyeh.clear();
}

void Utils::InitEffFuncThetaL (TF1* fitFun)
{
  fitFun->ReleaseParameter(0);
  fitFun->ReleaseParameter(1);
  fitFun->ReleaseParameter(2);
  fitFun->ReleaseParameter(3);
  fitFun->ReleaseParameter(4);
  fitFun->ReleaseParameter(5);

  fitFun->SetParameter(0,1e-3);
  fitFun->SetParameter(1,0.0);
  fitFun->SetParameter(2,-1.0);
  fitFun->SetParameter(3,0.0);
  fitFun->SetParameter(4,0.0);
  fitFun->SetParameter(5,0.0);
}

void Utils::InitEffFuncThetaK (TF1* fitFun)
{
  fitFun->ReleaseParameter(0);
  fitFun->ReleaseParameter(1);
  fitFun->ReleaseParameter(2);
  fitFun->ReleaseParameter(3);

  fitFun->SetParameter(0,0.0);
  fitFun->SetParameter(1,0.0);
  fitFun->SetParameter(2,0.0);
  fitFun->SetParameter(3,0.0);
}

void Utils::InitEffFuncPhi (TF1* fitFun)
{
  fitFun->ReleaseParameter(0);
  fitFun->ReleaseParameter(1);
  fitFun->ReleaseParameter(2);

  fitFun->SetParameter(0,0.0);
  fitFun->SetParameter(1,0.0);
  fitFun->SetParameter(2,0.0);
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


  std::cout << "[Utils::AddConstraint3D]" << std::endl;
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


#if ROOFIT
std::string Utils::Transformer (std::string varName, bool doIt, double& varValOut, double& varValOutELo, double& varValOutEHi, RooFitResult* fitResult, RooRealVar* varValIn1, RooRealVar* varValIn2, RooRealVar* varValIn3)
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

  else if (blockName == "fitValGlob")     return 12;
  else if (blockName == "fitValBins")     return 13;

  else if (blockName == "I[S*E]okTag")    return 14;
  else if (blockName == "I[S*E]misTag")   return 15;

  else if (blockName == "BF")             return 16;
  else if (blockName == "fitSyst")        return 17;
  else if (blockName == "fitNLL")         return 18;

  else if (blockName == "analyEffokTag")  return 19;
  else if (blockName == "analyEffmisTag") return 20;

  else if (blockName == "genericpar")     return 21;
  else if (blockName == "lumi")           return 22;
  else if (blockName == "dtype")          return 23;

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

  else if (varName == "nBkgComb")       return 21;
  else if (varName == "nMisTagFrac")    return 22;
  else if (varName == "nBkgPeak")       return 23;
  else if (varName == "nSig")           return 24;

  else if (varName == "nPolyP1")        return 25;
  else if (varName == "p1Poly0")        return 26;
  else if (varName == "p1Poly1")        return 27;
  else if (varName == "p1Poly2")        return 28;
  else if (varName == "p1Poly3")        return 29;
  else if (varName == "p1Poly4")        return 30;

  else if (varName == "nPolyC1")        return 31;
  else if (varName == "c1Poly0")        return 32;
  else if (varName == "c1Poly1")        return 33;
  else if (varName == "c1Poly2")        return 34;
  else if (varName == "c1Poly3")        return 35;
  else if (varName == "c1Poly4")        return 36;

  else if (varName == "nPolyP2")        return 37;
  else if (varName == "p2Poly0")        return 38;
  else if (varName == "p2Poly1")        return 39;
  else if (varName == "p2Poly2")        return 40;
  else if (varName == "p2Poly3")        return 41;
  else if (varName == "p2Poly4")        return 42;

  else if (varName == "nPolyC2")        return 43;
  else if (varName == "c2Poly0")        return 44;
  else if (varName == "c2Poly1")        return 45;
  else if (varName == "c2Poly2")        return 46;
  else if (varName == "c2Poly3")        return 47;
  else if (varName == "c2Poly4")        return 48;

  else if (varName == "nPolyP3")        return 49;
  else if (varName == "p3Poly0")        return 50;
  else if (varName == "p3Poly1")        return 51;
  else if (varName == "p3Poly2")        return 52;
  else if (varName == "p3Poly3")        return 53;
  else if (varName == "p3Poly4")        return 54;

  else if (varName == "nPolyC3")        return 55;
  else if (varName == "c3Poly0")        return 56;
  else if (varName == "c3Poly1")        return 57;
  else if (varName == "c3Poly2")        return 58;
  else if (varName == "c3Poly3")        return 59;
  else if (varName == "c3Poly4")        return 60;

  else if (varName == "FlS")            return 61;
  else if (varName == "AfbS")           return 62;
  else if (varName == "P1S")            return 63;
  else if (varName == "P2S")            return 64;
  else if (varName == "FsS")            return 65;
  else if (varName == "AsS")            return 66;

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
  std::cout << "[Utils::ReadSelectionCuts]\tSelection cuts from file : " << fileName.c_str() << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("selecuts"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
    {
      SeleCuts.push_back(atof(ParVector[i].c_str()));
      std::cout << "Selection cut #" << i << " from config file : " << SeleCuts[i] << std::endl;
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
  std::cout << "[Utils::ReadPreselectionCut]\tPre-selection cuts from file : " << fileName.c_str() << std::endl;
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
  std::cout << "[Utils::ReadGenericParam]\tGeneric parameters from file : " << fileName.c_str() << std::endl;
  ParVector.clear();
  ParameterFile->ReadFromFile(ParFileBlockN("genericpar"),&ParVector);
  for (unsigned int i = 0; i < ParVector.size(); i++)
    {
      GenericPars.push_back(ParVector[i]);
      std::cout << "Generic parameter #" << i << " from config file : " << GenericPars[i].c_str() << std::endl;
    }


  ParVector.clear();
  delete ParameterFile;
}

bool Utils::SetGenericParam (std::string parName, std::string val)
{
       if (parName == "NormJPSInotPSIP")     GenericPars[0]  = val;
  else if (parName == "DegreeInterp")        GenericPars[1]  = val;
  else if (parName == "TransfTolerance")     GenericPars[2]  = val;
  else if (parName == "UseMINOS")            GenericPars[3]  = val;
  else if (parName == "ApplyConstr")         GenericPars[4]  = val;
  else if (parName == "CtrlFitWrkFlow")      GenericPars[5]  = val;
  else if (parName == "CtrlMisTagWrkFlow")   GenericPars[6]  = val;
  else if (parName == "SaveMisTagFrac")      GenericPars[7]  = val;
  else if (parName == "UseSPwave")           GenericPars[8]  = val;
  else if (parName == "doTransf")            GenericPars[9]  = val;
  else if (parName == "B0MassIntervalLeft")  GenericPars[10] = val;
  else if (parName == "B0MassIntervalRight") GenericPars[11] = val;
  else if (parName == "NSigmaB0")            GenericPars[12] = val;
  else if (parName == "NSigmaB0S")           GenericPars[13] = val;
  else if (parName == "NSigmaB0B")           GenericPars[14] = val;
  else if (parName == "NSigmaPsi")           GenericPars[15] = val;
  else if (parName == "B&psiMassJpsiLo")     GenericPars[16] = val;
  else if (parName == "B&psiMassJpsiHi")     GenericPars[17] = val;
  else if (parName == "B&psiMassPsiPLo")     GenericPars[18] = val;
  else if (parName == "B&psiMassPsiPHi")     GenericPars[19] = val;
  else if (parName == "SIGMAS1")             GenericPars[20] = val;
  else if (parName == "SIGMAS2")             GenericPars[21] = val;
  else if (parName == "FRACMASSS")           GenericPars[22] = val;
  else return false;

  return true;
}

std::string Utils::GetGenericParam (std::string parName)
{
       if (parName == "NormJPSInotPSIP")     return GenericPars[0];
  else if (parName == "DegreeInterp")        return GenericPars[1];
  else if (parName == "TransfTolerance")     return GenericPars[2];
  else if (parName == "UseMINOS")            return GenericPars[3];
  else if (parName == "ApplyConstr")         return GenericPars[4];
  else if (parName == "CtrlFitWrkFlow")      return GenericPars[5];
  else if (parName == "CtrlMisTagWrkFlow")   return GenericPars[6];
  else if (parName == "SaveMisTagFrac")      return GenericPars[7];
  else if (parName == "UseSPwave")           return GenericPars[8];
  else if (parName == "doTransf")            return GenericPars[9];
  else if (parName == "B0MassIntervalLeft")  return GenericPars[10];
  else if (parName == "B0MassIntervalRight") return GenericPars[11];
  else if (parName == "NSigmaB0")            return GenericPars[12];
  else if (parName == "NSigmaB0S")           return GenericPars[13];
  else if (parName == "NSigmaB0B")           return GenericPars[14];
  else if (parName == "NSigmaPsi")           return GenericPars[15];
  else if (parName == "B&psiMassJpsiLo")     return GenericPars[16];
  else if (parName == "B&psiMassJpsiHi")     return GenericPars[17];
  else if (parName == "B&psiMassPsiPLo")     return GenericPars[18];
  else if (parName == "B&psiMassPsiPHi")     return GenericPars[19];
  else if (parName == "SIGMAS1")             return GenericPars[20];
  else if (parName == "SIGMAS2")             return GenericPars[21];
  else if (parName == "FRACMASSS")           return GenericPars[22];
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

void Utils::ResetEffValue (effValue* myEffVal, double value)
{
  myEffVal->Num1 = value;
  myEffVal->Num2 = value;
  myEffVal->Den1 = value;
  myEffVal->Den2 = value;
  
  // Poissonian errors
  myEffVal->Err2PoisNum1 = value;
  myEffVal->Err2PoisNum2 = value;
  myEffVal->Err2PoisDen1 = value;
  myEffVal->Err2PoisDen2 = value;
  
  // Weight errors
  myEffVal->Err2WeigNum1 = value;
  myEffVal->Err2WeigNum2 = value;
  myEffVal->Err2WeigDen1 = value;
  myEffVal->Err2WeigDen2 = value;
}

std::string Utils::GetHisto2DEffName (int SignalType)
{
  std::string name;

  if (RIGHTflavorTAG == true)
    {
      if      (SignalType == B0ToKstMuMu)  name = Histo2DEffNameOkTagSig;
      else if (SignalType == B0ToJPsiKst)  name = Histo2DEffNameOkTagJPsi;
      else if (SignalType == B0ToPsi2SKst) name = Histo2DEffNameOkTagPsi2S;
      else
	{
	  std::cout << "[Utils::GetHisto2DEffName]\tSignalType not valid : " << SignalType << std::endl;
	  exit (EXIT_FAILURE);
	}
    }
  else
    {
      if      (SignalType == B0ToKstMuMu)  name = Histo2DEffNameMisTagSig;
      else if (SignalType == B0ToJPsiKst)  name = Histo2DEffNameMisTagJPsi;
      else if (SignalType == B0ToPsi2SKst) name = Histo2DEffNameMisTagPsi2S;
      else
	{
	  std::cout << "[Utils::GetHisto2DEffName]\tSignalType not valid : " << SignalType << std::endl;
	  exit (EXIT_FAILURE);
	}
    }

  std::cout << "[Utils::GetHisto2DEffName]\tChosen name : " << name.c_str() << std::endl;
  return name;
}

void Utils::SetHisto2DEffName (int SignalType, std::string newName)
{
  if (RIGHTflavorTAG == true)
    {
      if      (SignalType == B0ToKstMuMu)  Histo2DEffNameOkTagSig   = newName;
      else if (SignalType == B0ToJPsiKst)  Histo2DEffNameOkTagJPsi  = newName;
      else if (SignalType == B0ToPsi2SKst) Histo2DEffNameOkTagPsi2S = newName;
      else
	{
	  std::cout << "[Utils::SetHisto2DEffName]\tSignalType not valid : " << SignalType << std::endl;
	  exit (EXIT_FAILURE);
	}
    }
  else
    {
      if      (SignalType == B0ToKstMuMu)  Histo2DEffNameMisTagSig   = newName;
      else if (SignalType == B0ToJPsiKst)  Histo2DEffNameMisTagJPsi  = newName;
      else if (SignalType == B0ToPsi2SKst) Histo2DEffNameMisTagPsi2S = newName;
      else
	{
	  std::cout << "[Utils::GetHisto2DEffName]\tSignalType not valid : " << SignalType << std::endl;
	  exit (EXIT_FAILURE);
	}
    }
}

std::string Utils::GetHisto3DEffName (int SignalType)
{
  std::string name;

  if (RIGHTflavorTAG == true)
    {
      if      (SignalType == B0ToKstMuMu)  name = Histo3DEffNameOkTagSig;
      else if (SignalType == B0ToJPsiKst)  name = Histo3DEffNameOkTagJPsi;
      else if (SignalType == B0ToPsi2SKst) name = Histo3DEffNameOkTagPsi2S;
      else
	{
	  std::cout << "[Utils::GetHisto3DEffName]\tSignalType not valid : " << SignalType << std::endl;
	  exit (EXIT_FAILURE);
	}
    }
  else
    {
      if      (SignalType == B0ToKstMuMu)  name = Histo3DEffNameMisTagSig;
      else if (SignalType == B0ToJPsiKst)  name = Histo3DEffNameMisTagJPsi;
      else if (SignalType == B0ToPsi2SKst) name = Histo3DEffNameMisTagPsi2S;
      else
	{
	  std::cout << "[Utils::GetHisto3DEffName]\tSignalType not valid : " << SignalType << std::endl;
	  exit (EXIT_FAILURE);
	}
    }

  std::cout << "[Utils::GetHisto3DEffName]\tChosen name : " << name.c_str() << std::endl;
  return name;
}

void Utils::SetHisto3DEffName (int SignalType, std::string newName)
{
  if (RIGHTflavorTAG == true)
    {
      if      (SignalType == B0ToKstMuMu)  Histo3DEffNameOkTagSig   = newName;
      else if (SignalType == B0ToJPsiKst)  Histo3DEffNameOkTagJPsi  = newName;
      else if (SignalType == B0ToPsi2SKst) Histo3DEffNameOkTagPsi2S = newName;
      else
	{
	  std::cout << "[Utils::SetHisto3DEffName]\tSignalType not valid : " << SignalType << std::endl;
	  exit (EXIT_FAILURE);
	}
    }
  else
    {
      if      (SignalType == B0ToKstMuMu)  Histo3DEffNameMisTagSig   = newName;
      else if (SignalType == B0ToJPsiKst)  Histo3DEffNameMisTagJPsi  = newName;
      else if (SignalType == B0ToPsi2SKst) Histo3DEffNameMisTagPsi2S = newName;
      else
	{
	  std::cout << "[Utils::GetHisto3DEffName]\tSignalType not valid : " << SignalType << std::endl;
	  exit (EXIT_FAILURE);
	}
    }
}

void Utils::SetDirEfficiency (std::string newName)
// ###############################################
// # If DirEfficiency = "" then the program will #
// # use a relative path to the efficiency files #
// ###############################################
{
  DirEfficiency = newName;
}
