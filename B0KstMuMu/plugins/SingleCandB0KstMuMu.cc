// ##########################################################################
// # Program to select one single candidate the B0 --> K*0 mu+ mu- analysis #
// ##########################################################################
// # Author: Mauro Dinardo                                                  #
// ##########################################################################

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGaxis.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLorentzVector.h>

#include <cmath>
#include <iostream>
#include <sstream>

#include "Utils.h"
#include "B0KstMuMuSingleCandTreeContent.h"

using std::cout;
using std::endl;
using std::string;


// ####################
// # Global constants #
// ####################
#define DoTrigCheck     1
#define DoMCTruth       false // Compute the single candidate variables from MC-GEN values
#define TagFromTruth    false // Assign the CP-eigenstate from MC-GEN
#define PARAMETERFILEIN "/python/ParameterFile.txt"


// ####################
// # Global variables #
// ####################
Utils* Utility;
TTree* theTreeIn;
TTree* theTreeOut;
B0KstMuMuSingleCandTreeContent* NTupleIn;
B0KstMuMuSingleCandTreeContent* NTupleOut;


// #######################
// # Function Definition #
// #######################
void SelectBestCand      (int SignalType);
void BestCandPerformance (int SignalType);


// ###########################
// # Function Implementation #
// ###########################
void SelectBestCand (int SignalType)
{
  unsigned int countCands;
  bool B0notB0bar;
  int BestCandIndx;
  int TrigCat;

  NTupleOut->ClearNTuple();
  NTupleOut->MakeTreeBranches(theTreeOut);

  NTupleIn->ClearNTuple();
  NTupleIn->SetBranchAddresses(theTreeIn);
  int nEntries = theTreeIn->GetEntries();
  cout << "\n[SingleCandB0KstMuMu::SelectBestCand]\t@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;


  // ###########################
  // # Loop over B0 Candidates #
  // ###########################
  for (int entry = 0; entry < nEntries; entry++)
    {
      theTreeIn->GetEntry(entry);

      if (Utility->ChooseBestCand(NTupleIn, DoTrigCheck, static_cast<double>(entry)/static_cast<double>(nEntries), &BestCandIndx, &B0notB0bar, &TrigCat, &countCands) == true)
	{
	  NTupleOut->CopyCandidate(NTupleIn, BestCandIndx);
	  
	  // #############################################
	  // # Adding new variables for angular analysis #
	  // #############################################
	  if (((TagFromTruth == false) && (B0notB0bar == true)) || ((TagFromTruth == true) && (NTupleIn->genSignal == SignalType)))
	    {
	      double cosThetaMup, cosThetaMupErr;
	      double cosThetaK, cosThetaKErr;
	      double phiKstMuMuPlane;
	      if (DoMCTruth == false)
		{
		  // ############
		  // # B0 boost #
		  // ############
		  TLorentzVector LoreVecB0;


		  TLorentzVector LoreVecMuMu;
		  LoreVecMuMu.SetXYZM(NTupleIn->mumuPx->at(BestCandIndx),NTupleIn->mumuPy->at(BestCandIndx),NTupleIn->mumuPz->at(BestCandIndx),NTupleIn->mumuMass->at(BestCandIndx));
		  TVector3 boostMuMu = LoreVecMuMu.BoostVector();
		  TLorentzVector LoreVecMup;
		  LoreVecMup.SetXYZM(NTupleIn->mupPx->at(BestCandIndx),NTupleIn->mupPy->at(BestCandIndx),NTupleIn->mupPz->at(BestCandIndx),Utility->muonMass);
		  // ################################
		  // # Boost mu+ to mumu ref. frame #
		  // ################################
		  LoreVecMup.Boost(-boostMuMu);
		  // ###############################
		  // # Boost B0 to mumu ref. frame #
		  // ###############################
		  LoreVecB0.SetXYZM(NTupleIn->bPx->at(BestCandIndx),NTupleIn->bPy->at(BestCandIndx),NTupleIn->bPz->at(BestCandIndx),NTupleIn->bMass->at(BestCandIndx));
		  LoreVecB0.Boost(-boostMuMu);
		  // ##################
		  // # Compute angles #
		  // ##################
		  Utility->computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
					   LoreVecMup.Px(),LoreVecMup.Py(),LoreVecMup.Pz(),
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   &cosThetaMup,&cosThetaMupErr);


		  TLorentzVector LoreVecKst;
		  LoreVecKst.SetXYZM(NTupleIn->kstPx->at(BestCandIndx),NTupleIn->kstPy->at(BestCandIndx),NTupleIn->kstPz->at(BestCandIndx),NTupleIn->kstMass->at(BestCandIndx));
		  TVector3 boostKst = LoreVecKst.BoostVector();
		  TLorentzVector LoreVecK;
		  LoreVecK.SetXYZM(NTupleIn->kstTrkpPx->at(BestCandIndx),NTupleIn->kstTrkpPy->at(BestCandIndx),NTupleIn->kstTrkpPz->at(BestCandIndx),Utility->kaonMass);
		  // ##############################
		  // # Boost K+ to K*0 ref. frame #
		  // ##############################
		  LoreVecK.Boost(-boostKst);
		  // ##############################
		  // # Boost B0 to K*0 ref. frame #
		  // ##############################
		  LoreVecB0.SetXYZM(NTupleIn->bPx->at(BestCandIndx),NTupleIn->bPy->at(BestCandIndx),NTupleIn->bPz->at(BestCandIndx),NTupleIn->bMass->at(BestCandIndx));
		  LoreVecB0.Boost(-boostKst);
		  // ##################
		  // # Compute angles #
		  // ##################
		  Utility->computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
					   LoreVecK.Px(),LoreVecK.Py(),LoreVecK.Pz(),
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   &cosThetaK,&cosThetaKErr);


		  // ######################################################################
		  // # Angle between [mu+ - mu-] and [K - pi] planes in the B0 ref. frame #
		  // ######################################################################
		  LoreVecB0.SetXYZM(NTupleIn->bPx->at(BestCandIndx),NTupleIn->bPy->at(BestCandIndx),NTupleIn->bPz->at(BestCandIndx),NTupleIn->bBarMass->at(BestCandIndx));
		  TVector3 boostB0 = LoreVecB0.BoostVector();
		  LoreVecMup.SetXYZM(NTupleIn->mupPx->at(BestCandIndx),NTupleIn->mupPy->at(BestCandIndx),NTupleIn->mupPz->at(BestCandIndx),Utility->muonMass);
		  TLorentzVector LoreVecMum;
		  LoreVecMum.SetXYZM(NTupleIn->mumPx->at(BestCandIndx),NTupleIn->mumPy->at(BestCandIndx),NTupleIn->mumPz->at(BestCandIndx),Utility->muonMass);
		  LoreVecK.SetXYZM(NTupleIn->kstTrkpPx->at(BestCandIndx),NTupleIn->kstTrkpPy->at(BestCandIndx),NTupleIn->kstTrkpPz->at(BestCandIndx),Utility->kaonMass);
		  TLorentzVector LoreVecPi;
		  LoreVecPi.SetXYZM(NTupleIn->kstTrkmPx->at(BestCandIndx),NTupleIn->kstTrkmPy->at(BestCandIndx),NTupleIn->kstTrkmPz->at(BestCandIndx),Utility->pionMass);

		  LoreVecMum.Boost(-boostB0);
		  LoreVecMup.Boost(-boostB0);
		  LoreVecK.Boost(-boostB0);
		  LoreVecPi.Boost(-boostB0);
		  TVector3 MuMuPlane = LoreVecMup.Vect().Cross(LoreVecMum.Vect());
		  TVector3 KstPlane  = LoreVecK.Vect().Cross(LoreVecPi.Vect());
		  if (MuMuPlane.Cross(KstPlane).Dot(-LoreVecB0.Vect()) > 0.0) phiKstMuMuPlane = MuMuPlane.Angle(KstPlane);
		  else                                                        phiKstMuMuPlane = -MuMuPlane.Angle(KstPlane);
 		}
	      else if (DoMCTruth == true)
		{
		  // ############
		  // # B0 boost #
		  // ############
		  TLorentzVector LoreVecB0;
	

		  TLorentzVector LoreVecMuMu;
		  LoreVecMuMu.SetXYZM(NTupleIn->genMumPx+NTupleIn->genMupPx,NTupleIn->genMumPy+NTupleIn->genMupPy,NTupleIn->genMumPz+NTupleIn->genMupPz,
				      Utility->computeInvMass(NTupleIn->genMumPx,NTupleIn->genMumPy,NTupleIn->genMumPz,Utility->muonMass,
							      NTupleIn->genMupPx,NTupleIn->genMupPy,NTupleIn->genMupPz,Utility->muonMass));
		  TVector3 boostMuMu = LoreVecMuMu.BoostVector();
		  TLorentzVector LoreVecMup;
		  LoreVecMup.SetXYZM(NTupleIn->genMupPx,NTupleIn->genMupPy,NTupleIn->genMupPz,Utility->muonMass);
		  // ################################
		  // # Boost mu+ to mumu ref. frame #
		  // ################################
		  LoreVecMup.Boost(-boostMuMu);
		  // ###############################
		  // # Boost B0 to mumu ref. frame #
		  // ###############################
		  LoreVecB0.SetXYZM(NTupleIn->genB0Px,NTupleIn->genB0Py,NTupleIn->genB0Pz,NTupleIn->genB0Mass);
		  LoreVecB0.Boost(-boostMuMu);
		  // ##################
		  // # Compute angles #
		  // ##################
		  Utility->computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
					   LoreVecMup.Px(),LoreVecMup.Py(),LoreVecMup.Pz(),
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   &cosThetaMup,&cosThetaMupErr);


		  TLorentzVector LoreVecKst;
		  LoreVecKst.SetXYZM(NTupleIn->genKstPx,NTupleIn->genKstPy,NTupleIn->genKstPz,NTupleIn->genKstMass);
		  TVector3 boostKst = LoreVecKst.BoostVector();
		  TLorentzVector LoreVecK;
		  LoreVecK.SetXYZM(NTupleIn->genKstTrkpPx,NTupleIn->genKstTrkpPy,NTupleIn->genKstTrkpPz,Utility->kaonMass);
		  // ##############################
		  // # Boost K+ to K*0 ref. frame #
		  // ##############################
		  LoreVecK.Boost(-boostKst);
		  // ##############################
		  // # Boost B0 to K*0 ref. frame #
		  // ##############################
		  LoreVecB0.SetXYZM(NTupleIn->genB0Px,NTupleIn->genB0Py,NTupleIn->genB0Pz,NTupleIn->genB0Mass);
		  LoreVecB0.Boost(-boostKst);
		  // ##################
		  // # Compute angles #
		  // ##################
		  Utility->computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
					   LoreVecK.Px(),LoreVecK.Py(),LoreVecK.Pz(),
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   &cosThetaK,&cosThetaKErr);


		  // ######################################################################
		  // # Angle between [mu+ - mu-] and [K - pi] planes in the B0 ref. frame #
		  // ######################################################################
		  LoreVecB0.SetXYZM(NTupleIn->genB0Px,NTupleIn->genB0Py,NTupleIn->genB0Pz,NTupleIn->genB0Mass);
		  TVector3 boostB0 = LoreVecB0.BoostVector();
		  LoreVecMup.SetXYZM(NTupleIn->genMupPx,NTupleIn->genMupPy,NTupleIn->genMupPz,Utility->muonMass);
		  TLorentzVector LoreVecMum;
		  LoreVecMum.SetXYZM(NTupleIn->genMumPx,NTupleIn->genMumPy,NTupleIn->genMumPz,Utility->muonMass);
		  LoreVecK.SetXYZM(NTupleIn->genKstTrkpPx,NTupleIn->genKstTrkpPy,NTupleIn->genKstTrkpPz,Utility->kaonMass);
		  TLorentzVector LoreVecPi;
		  LoreVecPi.SetXYZM(NTupleIn->genKstTrkmPx,NTupleIn->genKstTrkmPy,NTupleIn->genKstTrkmPz,Utility->pionMass);

		  LoreVecMum.Boost(-boostB0);
		  LoreVecMup.Boost(-boostB0);
		  LoreVecK.Boost(-boostB0);
		  LoreVecPi.Boost(-boostB0);
		  TVector3 MuMuPlane = LoreVecMup.Vect().Cross(LoreVecMum.Vect());
		  TVector3 KstPlane  = LoreVecK.Vect().Cross(LoreVecPi.Vect());
		  if (MuMuPlane.Cross(KstPlane).Dot(-LoreVecB0.Vect()) > 0.0) phiKstMuMuPlane = MuMuPlane.Angle(KstPlane);
		  else                                                        phiKstMuMuPlane = -MuMuPlane.Angle(KstPlane);
		}

	      NTupleOut->B0MassArb          = NTupleIn->bMass->at(BestCandIndx);
	      NTupleOut->kstMassArb         = NTupleIn->kstMass->at(BestCandIndx);
	      NTupleOut->CosThetaMuArb      = cosThetaMup;
	      NTupleOut->CosThetaKArb       = cosThetaK;
	      NTupleOut->PhiKstMuMuPlaneArb = phiKstMuMuPlane;
	    }
	  else if (((TagFromTruth == false) && (B0notB0bar == false)) || ((TagFromTruth == true) && (NTupleIn->genSignal == SignalType+1)))
	    {
	      double cosThetaMum, cosThetaMumErr;
	      double cosThetaK, cosThetaKErr;
	      double phiKstMuMuPlane;
	      if (DoMCTruth == false)
		{
		  // ############
		  // # B0 boost #
		  // ############
		  TLorentzVector LoreVecB0;
	

		  TLorentzVector LoreVecMuMu;
		  LoreVecMuMu.SetXYZM(NTupleIn->mumuPx->at(BestCandIndx),NTupleIn->mumuPy->at(BestCandIndx),NTupleIn->mumuPz->at(BestCandIndx),NTupleIn->mumuMass->at(BestCandIndx));
		  TVector3 boostMuMu = LoreVecMuMu.BoostVector();
		  TLorentzVector LoreVecMum;
		  LoreVecMum.SetXYZM(NTupleIn->mumPx->at(BestCandIndx),NTupleIn->mumPy->at(BestCandIndx),NTupleIn->mumPz->at(BestCandIndx),Utility->muonMass);
		  // ################################
		  // # Boost mu- to mumu ref. frame #
		  // ################################
		  LoreVecMum.Boost(-boostMuMu);
		  // ###############################
		  // # Boost B0 to mumu ref. frame #
		  // ###############################
		  LoreVecB0.SetXYZM(NTupleIn->bPx->at(BestCandIndx),NTupleIn->bPy->at(BestCandIndx),NTupleIn->bPz->at(BestCandIndx),NTupleIn->bBarMass->at(BestCandIndx));
		  LoreVecB0.Boost(-boostMuMu);
		  // ##################
		  // # Compute angles #
		  // ##################
		  Utility->computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
					   LoreVecMum.Px(),LoreVecMum.Py(),LoreVecMum.Pz(),
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   &cosThetaMum,&cosThetaMumErr);


		  TLorentzVector LoreVecKst;
		  LoreVecKst.SetXYZM(NTupleIn->kstPx->at(BestCandIndx),NTupleIn->kstPy->at(BestCandIndx),NTupleIn->kstPz->at(BestCandIndx),NTupleIn->kstBarMass->at(BestCandIndx));
		  TVector3 boostKst = LoreVecKst.BoostVector();
		  TLorentzVector LoreVecK;
		  LoreVecK.SetXYZM(NTupleIn->kstTrkmPx->at(BestCandIndx),NTupleIn->kstTrkmPy->at(BestCandIndx),NTupleIn->kstTrkmPz->at(BestCandIndx),Utility->kaonMass);
		  // ##############################
		  // # Boost K- to K*0 ref. frame #
		  // ##############################
		  LoreVecK.Boost(-boostKst);
		  // ##############################
		  // # Boost B0 to K*0 ref. frame #
		  // ##############################
		  LoreVecB0.SetXYZM(NTupleIn->bPx->at(BestCandIndx),NTupleIn->bPy->at(BestCandIndx),NTupleIn->bPz->at(BestCandIndx),NTupleIn->bBarMass->at(BestCandIndx));
		  LoreVecB0.Boost(-boostKst);
		  // ##################
		  // # Compute angles #
		  // ##################
		  Utility->computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
					   LoreVecK.Px(),LoreVecK.Py(),LoreVecK.Pz(),
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   &cosThetaK,&cosThetaKErr);


		  // ######################################################################
		  // # Angle between [mu+ - mu-] and [K - pi] planes in the B0 ref. frame #
		  // ######################################################################
		  LoreVecB0.SetXYZM(NTupleIn->bPx->at(BestCandIndx),NTupleIn->bPy->at(BestCandIndx),NTupleIn->bPz->at(BestCandIndx),NTupleIn->bBarMass->at(BestCandIndx));
		  TVector3 boostB0 = LoreVecB0.BoostVector();
		  LoreVecMum.SetXYZM(NTupleIn->mumPx->at(BestCandIndx),NTupleIn->mumPy->at(BestCandIndx),NTupleIn->mumPz->at(BestCandIndx),Utility->muonMass);
		  TLorentzVector LoreVecMup;
		  LoreVecMup.SetXYZM(NTupleIn->mupPx->at(BestCandIndx),NTupleIn->mupPy->at(BestCandIndx),NTupleIn->mupPz->at(BestCandIndx),Utility->muonMass);
		  LoreVecK.SetXYZM(NTupleIn->kstTrkmPx->at(BestCandIndx),NTupleIn->kstTrkmPy->at(BestCandIndx),NTupleIn->kstTrkmPz->at(BestCandIndx),Utility->kaonMass);
		  TLorentzVector LoreVecPi;
		  LoreVecPi.SetXYZM(NTupleIn->kstTrkpPx->at(BestCandIndx),NTupleIn->kstTrkpPy->at(BestCandIndx),NTupleIn->kstTrkpPz->at(BestCandIndx),Utility->pionMass);

		  LoreVecMum.Boost(-boostB0);
		  LoreVecMup.Boost(-boostB0);
		  LoreVecK.Boost(-boostB0);
		  LoreVecPi.Boost(-boostB0);
		  TVector3 MuMuPlane = LoreVecMum.Vect().Cross(LoreVecMup.Vect());
		  TVector3 KstPlane  = LoreVecK.Vect().Cross(LoreVecPi.Vect());
		  if (MuMuPlane.Cross(KstPlane).Dot(-LoreVecB0.Vect()) > 0.0) phiKstMuMuPlane = MuMuPlane.Angle(KstPlane);
		  else                                                        phiKstMuMuPlane = -MuMuPlane.Angle(KstPlane);
 		}
	      else if (DoMCTruth == true)
		{
		  // ############
		  // # B0 boost #
		  // ############
		  TLorentzVector LoreVecB0;
	

		  TLorentzVector LoreVecMuMu;
		  LoreVecMuMu.SetXYZM(NTupleIn->genMumPx+NTupleIn->genMupPx,NTupleIn->genMumPy+NTupleIn->genMupPy,NTupleIn->genMumPz+NTupleIn->genMupPz,
				      Utility->computeInvMass(NTupleIn->genMumPx,NTupleIn->genMumPy,NTupleIn->genMumPz,Utility->muonMass,
							      NTupleIn->genMupPx,NTupleIn->genMupPy,NTupleIn->genMupPz,Utility->muonMass));
		  TVector3 boostMuMu = LoreVecMuMu.BoostVector();
		  TLorentzVector LoreVecMum;
		  LoreVecMum.SetXYZM(NTupleIn->genMumPx,NTupleIn->genMumPy,NTupleIn->genMumPz,Utility->muonMass);
		  // ################################
		  // # Boost mu- to mumu ref. frame #
		  // ################################
		  LoreVecMum.Boost(-boostMuMu);
		  // ###############################
		  // # Boost B0 to mumu ref. frame #
		  // ###############################
		  LoreVecB0.SetXYZM(NTupleIn->genB0Px,NTupleIn->genB0Py,NTupleIn->genB0Pz,NTupleIn->genB0Mass);
		  LoreVecB0.Boost(-boostMuMu);
		  // ##################
		  // # Compute angles #
		  // ##################
		  Utility->computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
					   LoreVecMum.Px(),LoreVecMum.Py(),LoreVecMum.Pz(),
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   &cosThetaMum,&cosThetaMumErr);


		  TLorentzVector LoreVecKst;
		  LoreVecKst.SetXYZM(NTupleIn->genKstPx,NTupleIn->genKstPy,NTupleIn->genKstPz,NTupleIn->genKstMass);
		  TVector3 boostKst = LoreVecKst.BoostVector();
		  TLorentzVector LoreVecK;
		  LoreVecK.SetXYZM(NTupleIn->genKstTrkmPx,NTupleIn->genKstTrkmPy,NTupleIn->genKstTrkmPz,Utility->kaonMass);
		  // ##############################
		  // # Boost K- to K*0 ref. frame #
		  // ##############################
		  LoreVecK.Boost(-boostKst);
		  // ##############################
		  // # Boost B0 to K*0 ref. frame #
		  // ##############################
		  LoreVecB0.SetXYZM(NTupleIn->genB0Px,NTupleIn->genB0Py,NTupleIn->genB0Pz,NTupleIn->genB0Mass);
		  LoreVecB0.Boost(-boostKst);
		  // ##################
		  // # Compute angles #
		  // ##################
		  Utility->computeCosAlpha(-LoreVecB0.Px(),-LoreVecB0.Py(),-LoreVecB0.Pz(),
					   LoreVecK.Px(),LoreVecK.Py(),LoreVecK.Pz(),
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   0.0,0.0,0.0,
					   &cosThetaK,&cosThetaKErr);


		  // ######################################################################
		  // # Angle between [mu+ - mu-] and [K - pi] planes in the B0 ref. frame #
		  // ######################################################################
		  LoreVecB0.SetXYZM(NTupleIn->genB0Px,NTupleIn->genB0Py,NTupleIn->genB0Pz,NTupleIn->genB0Mass);
		  TVector3 boostB0 = LoreVecB0.BoostVector();
		  LoreVecMum.SetXYZM(NTupleIn->genMumPx,NTupleIn->genMumPy,NTupleIn->genMumPz,Utility->muonMass);
		  TLorentzVector LoreVecMup;
		  LoreVecMup.SetXYZM(NTupleIn->genMupPx,NTupleIn->genMupPy,NTupleIn->genMupPz,Utility->muonMass);
		  LoreVecK.SetXYZM(NTupleIn->genKstTrkmPx,NTupleIn->genKstTrkmPy,NTupleIn->genKstTrkmPz,Utility->kaonMass);
		  TLorentzVector LoreVecPi;
		  LoreVecPi.SetXYZM(NTupleIn->genKstTrkpPx,NTupleIn->genKstTrkpPy,NTupleIn->genKstTrkpPz,Utility->pionMass);

		  LoreVecMum.Boost(-boostB0);
		  LoreVecMup.Boost(-boostB0);
		  LoreVecK.Boost(-boostB0);
		  LoreVecPi.Boost(-boostB0);
		  TVector3 MuMuPlane = LoreVecMum.Vect().Cross(LoreVecMup.Vect());
		  TVector3 KstPlane  = LoreVecK.Vect().Cross(LoreVecPi.Vect());
		  if (MuMuPlane.Cross(KstPlane).Dot(-LoreVecB0.Vect()) > 0.0) phiKstMuMuPlane = MuMuPlane.Angle(KstPlane);
		  else                                                        phiKstMuMuPlane = -MuMuPlane.Angle(KstPlane);
 		}

	      NTupleOut->B0MassArb          = NTupleIn->bBarMass->at(BestCandIndx);
	      NTupleOut->kstMassArb         = NTupleIn->kstBarMass->at(BestCandIndx);
	      NTupleOut->CosThetaMuArb      = cosThetaMum;
	      NTupleOut->CosThetaKArb       = cosThetaK;
	      NTupleOut->PhiKstMuMuPlaneArb = phiKstMuMuPlane;
	    }

	  NTupleOut->TrigCat = TrigCat;

	  NTupleOut->B0notB0bar = B0notB0bar;
	  if (((NTupleIn->genSignal == SignalType)   && (B0notB0bar == true)) ||
	      ((NTupleIn->genSignal == SignalType+1) && (B0notB0bar == false))) NTupleOut->rightFlavorTag = true;
	  else                                                                  NTupleOut->rightFlavorTag = false;
	  NTupleOut->B0pT  = sqrt(NTupleIn->bPx->at(BestCandIndx)*NTupleIn->bPx->at(BestCandIndx) + NTupleIn->bPy->at(BestCandIndx)*NTupleIn->bPy->at(BestCandIndx));
	  NTupleOut->B0Eta = Utility->computeEta (NTupleIn->bPx->at(BestCandIndx),NTupleIn->bPy->at(BestCandIndx),NTupleIn->bPz->at(BestCandIndx));
	  NTupleOut->B0Phi = Utility->computePhi (NTupleIn->bPx->at(BestCandIndx),NTupleIn->bPy->at(BestCandIndx));

	  theTreeOut->Fill();
	  NTupleOut->ClearNTuple();
	}
    }
}


void BestCandPerformance (int SignalType)
// ##############################################################
// # Compute efficiency, purity and B-flavor tagging efficiency #
// ##############################################################
{
  unsigned int coutGoodEv          = 0;
  unsigned int countSeleEv         = 0;
  unsigned int countGoodBTag       = 0;
  unsigned int countGoodBbarTag    = 0;
  unsigned int countGoodBSeleEv    = 0;
  unsigned int countGoodBbarSeleEv = 0;
  unsigned int countGoodSeleEv     = 0;

  unsigned int countCands;
  bool B0notB0bar;
  int BestCandIndx;
  int TrigCat;

  unsigned int counter = 0;
  double countCandsAVG = 0.0;

  NTupleIn->ClearNTuple();
  NTupleIn->SetBranchAddresses(theTreeIn);
  int nEntries = theTreeIn->GetEntries();
  cout << "\n[SingleCandB0KstMuMu::BestCandPerformance]\t@@@ Total number of events in the tree: " << nEntries << " @@@" << endl;


  // ###########################
  // # Loop over B0 Candidates #
  // ###########################
  for (int entry = 0; entry < nEntries; entry++)
    {
      theTreeIn->GetEntry(entry);


      // #####################
      // # Count good events #
      // #####################
      for (unsigned int i = 0; i < NTupleIn->bMass->size(); i++)
	if (((NTupleIn->genSignal == SignalType) || (NTupleIn->genSignal == SignalType+1)) && (NTupleIn->truthMatchSignal->at(i) == true))
	  {
	    coutGoodEv++;
	    break;
	  }


      if (Utility->ChooseBestCand(NTupleIn, DoTrigCheck, static_cast<double>(entry)/static_cast<double>(nEntries), &BestCandIndx, &B0notB0bar, &TrigCat, &countCands) == true)
	{
	  countSeleEv++;

	  if ((NTupleIn->genSignal == SignalType) && (NTupleIn->truthMatchSignal->at(BestCandIndx) == true) && (B0notB0bar == true)) countGoodBTag++;
	  if ((NTupleIn->genSignal == SignalType) && (NTupleIn->truthMatchSignal->at(BestCandIndx) == true))                         countGoodBSeleEv++;

	  if ((NTupleIn->genSignal == SignalType+1) && (NTupleIn->truthMatchSignal->at(BestCandIndx) == true) && (B0notB0bar == false)) countGoodBbarTag++;
	  if ((NTupleIn->genSignal == SignalType+1) && (NTupleIn->truthMatchSignal->at(BestCandIndx) == true))                          countGoodBbarSeleEv++;

	  if (((NTupleIn->genSignal == SignalType) || (NTupleIn->genSignal == SignalType+1)) && (NTupleIn->truthMatchSignal->at(BestCandIndx) == true)) countGoodSeleEv++;
	}


      if (countCands != 0)
	{
	  counter++;
	  countCandsAVG += countCands;
	}
    }
  countCandsAVG = countCandsAVG / static_cast<double>(counter);

  cout << "\n@@@ Efficiency: " << static_cast<double>(countGoodSeleEv) / static_cast<double>(coutGoodEv) * 100. << " @@@" << endl;
  cout << "@@@ Purity: " << static_cast<double>(countGoodSeleEv) / static_cast<double>(countSeleEv) * 100. << " @@@" << endl;
  cout << "@@@ Total B-flavor tagging efficiency: " << static_cast<double>(countGoodBTag+countGoodBbarTag) / static_cast<double>(countGoodSeleEv) * 100. << " @@@" << endl;
  cout << "@@@ B-flavor tagging efficiency: " << static_cast<double>(countGoodBbarTag) / static_cast<double>(countGoodBSeleEv) * 100. << " @@@" << endl;
  cout << "@@@ Bbar-flavor tagging efficiency: " << static_cast<double>(countGoodBbarTag) / static_cast<double>(countGoodBbarSeleEv) * 100. << " @@@" << endl;
  cout << "@@@ Average number of candidates after all cuts applied but B0-vtx CL: " << countCandsAVG << " @@@" << endl;
}


int main (int argc, char** argv)
{
  if (argc >= 4)
    {
      string option = argv[1];
      if ((option == "singlecand") || (option == "eff"))
	{
	  string fileNameIn  = argv[2];
	  string fileNameOut = "";
	  string localVar    = "1";

	  Utility = new Utils();
	  Utility->ReadTriggerPathsANDCutsANDEntries(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str());
	  Utility->ReadSelectionCuts(Utility->MakeAnalysisPATH(PARAMETERFILEIN).c_str());


	  cout << "\n[SingleCandB0KstMuMu::main]\t@@@ Settings @@@" << endl;
	  cout << "Do Trigger Check: "   << DoTrigCheck << endl;
	  cout << "DoMCTruth: "          << DoMCTruth << endl;
	  cout << "TagFromTruth: "       << TagFromTruth << endl;
	  cout << "PARAMETERFILEIN: "    << PARAMETERFILEIN << endl;
	  cout << "Default SignalType: " << localVar << endl;


	  TFile* NtplFileIn  = new TFile(fileNameIn.c_str(), "READ");
	  TFile* NtplFileOut = NULL;
	  theTreeIn          = (TTree*) NtplFileIn->Get("B0KstMuMu/B0KstMuMuNTuple");
	  NTupleIn           = new B0KstMuMuSingleCandTreeContent();
	  NTupleIn->Init();

	  if (option == "singlecand")
	    {
	      fileNameOut = argv[3];
	      if (argc == 5) localVar = argv[4];

	      NtplFileOut = new TFile(fileNameOut.c_str(), "RECREATE");
	      NtplFileOut->mkdir("B0KstMuMu");
	      NtplFileOut->cd("B0KstMuMu");
	      theTreeOut  = new TTree("B0KstMuMuNTuple","B0KstMuMuNTuple");
	      NTupleOut   = new B0KstMuMuSingleCandTreeContent();
	      NTupleOut->Init();

	      SelectBestCand(atoi(localVar.c_str()));
	    }
	  else if (option == "eff")
	    {
	      if (argc == 4) localVar = argv[3];
	      BestCandPerformance(atoi(localVar.c_str()));
	    }
	  else
	    {
	      cout << "[SingleCandB0KstMuMu::main]\tWrong option: " << option << endl;
	      exit (EXIT_FAILURE);
	    }
	  cout << "\n@@@ Single candidate program done @@@" << endl;


	  // #############
	  // # Save Tree #
	  // #############
	  if (option == "singlecand")
	    {
	      NtplFileOut->cd("B0KstMuMu");
	      theTreeOut->Write();
	      NtplFileOut->Close();
	      NtplFileIn->Close();

	      delete NTupleOut;
	    }
	  return EXIT_SUCCESS;
	}
      else
	{
	  cout << "Wrong option: " << endl;
	  cout << "./SingleCandB0KstMuMu [singlecand eff] inputFile.root [[if singlecand]outputFile.root] [SignalType]" << endl;
	  cout << "- SignalType : if B0 --> K*0 mumu : 1; if B0 --> J/psi K*0 : 3; if B0 --> psi(2S) K*0 : 5" << endl;

	  return EXIT_FAILURE;
	}
    }
  else
    {
      cout << "Parameter missing: " << endl;
      cout << "./SingleCandB0KstMuMu [singlecand eff] inputFile.root [[if singlecand]outputFile.root] [SignalType]" << endl;
      cout << "- SignalType : if B0 --> K*0 mumu : 1; if B0 --> J/psi K*0 : 3; if B0 --> psi(2S) K*0 : 5" << endl;

      return EXIT_FAILURE;
    }

  return EXIT_SUCCESS;
}
