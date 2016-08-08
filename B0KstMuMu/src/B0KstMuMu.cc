// ##################################################
// # Description:                                   #
// # Make rootTuple for b --> s mu+ mu- analysis    #
// # Original Author:  Mauro Dinardo                #
// #         Created:  Mon Apr 27 09:53:19 MDT 2011 #
// ##################################################

// System include files
#include <memory>

// User include files
#include "../interface/B0KstMuMu.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"

#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "TMath.h"

#include <sstream>
#include <utility>


// ####################
// # Global constants #
// ####################
#define TRKMAXR 110.0 // [cm]
#define TRKMAXZ 280.0 // [cm]

#define PRIVTXNDOF  4.0
#define PRIVTXMAXZ 50.0 // [cm]
#define PRIVTXMAXR  2.0 // [cm]

#define MUVARTOLE 0.01  // [From 0 to 1]
#define HADVARTOLE 0.10 // [From 0 to 1]


// #######################
// # Truth matching cuts #
// #######################
#define RCUTMU 0.004 // [eta-phi]
#define RCUTTRK 0.3  // [eta-phi]


B0KstMuMu::B0KstMuMu (const edm::ParameterSet& iConfig) :
  hltTriggerResults_ ( iConfig.getUntrackedParameter<std::string>("HLTriggerResults", std::string("HLT")) ),
  vtxSample_         ( iConfig.getUntrackedParameter<std::string>("VtxSample",        std::string("offlinePrimaryVertices")) ),
  beamSpot_          ( iConfig.getUntrackedParameter<std::string>("BeamSpot",         std::string("offlineBeamSpot")) ),
  genParticles_      ( iConfig.getUntrackedParameter<std::string>("GenParticles",     std::string("genParticles")) ),
  muonType_          ( iConfig.getUntrackedParameter<std::string>("MuonType",         std::string("cleanPatMuonsTriggerMatch")) ),
  trackType_         ( iConfig.getUntrackedParameter<std::string>("TrackType",        std::string("cleanPatTrackCands")) ),
  parameterFile_     ( iConfig.getUntrackedParameter<std::string>("ParameterFile",    std::string("ParameterFile.txt")) ),
  doGenReco_         ( iConfig.getUntrackedParameter<unsigned int>("doGenReco",       1) ),
  printMsg           ( iConfig.getUntrackedParameter<bool>("printMsg",                false) ),
  
  theTree(0)
{
  std::cout << "\n@@@ Ntuplizer configuration parameters @@@" << std::endl;
  std::cout << __LINE__ << " : hltTriggerResults    = " << hltTriggerResults_ << std::endl;
  std::cout << __LINE__ << " : vtxSample     = " << vtxSample_ << std::endl;
  std::cout << __LINE__ << " : beamSpot      = " << beamSpot_ << std::endl;
  std::cout << __LINE__ << " : genParticles  = " << genParticles_ << std::endl;
  std::cout << __LINE__ << " : trackType     = " << trackType_ << std::endl;
  std::cout << __LINE__ << " : parameterFile = " << parameterFile_ << std::endl;
  std::cout << __LINE__ << " : doGenReco     = " << doGenReco_ << std::endl;
  std::cout << __LINE__ << " : printMsg      = " << printMsg << std::endl;

  NTuple = new B0KstMuMuTreeContent();
  NTuple->Init();
  NTuple->ClearNTuple();

  Utility = new Utils();

  ParameterFile = new ReadParameters(parameterFile_.c_str());
}


B0KstMuMu::~B0KstMuMu ()
{
  delete NTuple;
  delete Utility;
  delete ParameterFile;
}


void B0KstMuMu::analyze (const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // ######################
  // # Internal Variables #
  // ######################
  std::stringstream myString;
  std::string tmpString;

  std::string MuMCat;
  std::string MuPCat;

  const ParticleMass muonMass = Utility->muonMass;
  const ParticleMass pionMass = Utility->pionMass;
  const ParticleMass kaonMass = Utility->kaonMass;
  float muonMassErr           = Utility->muonMassErr;
  float pionMassErr           = Utility->pionMassErr;
  float kaonMassErr           = Utility->kaonMassErr;
  
  double chi;
  double ndf;

  double LSVtx;
  double LSVtxErr;
  double LSBS;
  double LSBSErr;

  double cosAlphaVtx;
  double cosAlphaVtxErr;
  double cosAlphaBS;
  double cosAlphaBSErr;

  double deltaEtaPhi;
  double pT;
  double eta;

  KinematicParticleFactoryFromTransientTrack partFactory;

  AdaptiveVertexFitter theVtxFitter;                              // Vertex fitter in nominal reconstruction
  KinematicParticleVertexFitter PartVtxFitter;                    // Vertex fit with vtx constraint
  // KinematicParticleFitter PartVtxFitterMassVtxConst;              // Vertex fit with constraint on the vtx and on the vtx mass
  // KinematicConstrainedVertexFitter PartVtxFitterMassPartConst;    // Vertex fit with more than 3 tracks and with a constraint on the vtx and on the mass of 2 of the tracks

  KinematicConstraint* B0MassConstraint = new MassKinematicConstraint(Utility->B0Mass, Utility->B0MassErr);

  ClosestApproachInRPhi ClosestApp;
  GlobalPoint XingPoint;

  reco::TrackRef muTrackTmp;
  reco::TrackRef muTrackm;
  reco::TrackRef muTrackp;
  reco::TrackRef Trackm;
  reco::TrackRef Trackp;

  std::pair <bool,Measurement1D> theDCAXVtx;
  TrajectoryStateClosestToPoint theDCAXBS;


  if (printMsg == true) std::cout << "\n\n" << __LINE__ << " : @@@@@@ Start Analyzer @@@@@@" << std::endl;


  // Get Gen-Particles
  edm::Handle<reco::GenParticleCollection> genParticles;
  if ((doGenReco_ == 2) || (doGenReco_ == 3)) iEvent.getByLabel(genParticles_, genParticles);

  // Get magnetic field
  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);      
  
  // Get HLT results
  edm::Handle<edm::TriggerResults> hltTriggerResults;
  try
    {
      iEvent.getByLabel(edm::InputTag(std::string("TriggerResults::")+hltTriggerResults_), hltTriggerResults);
    }
  catch ( ... )
    {
      if (printMsg == true) std::cout << __LINE__ << " : couldn't get handle on HLT Trigger" << std::endl;
    }
  

  // ############################
  // # Save trigger information #
  // ############################
  HLTConfigProvider hltConfig_;
  bool changed = true;
  if (((hltTriggerResults.isValid() == false) || (hltConfig_.init(iEvent.getRun(), iSetup, hltTriggerResults_, changed) == false)) &&
      (printMsg == true)) std::cout << __LINE__ << " : no trigger results" << std::endl;
  else
    {
      // Get hold of trigger names - based on TriggerResults object
      const edm::TriggerNames& triggerNames_ = iEvent.triggerNames(*hltTriggerResults);
      
      for (unsigned int itrig = 0; itrig < hltTriggerResults->size(); itrig++)
	{
	  if ((*hltTriggerResults)[itrig].accept() == 1)
	    {
	      std::string trigName = triggerNames_.triggerName(itrig);
	      int trigPrescale = hltConfig_.prescaleValue(itrig, trigName);
	      
	      if (printMsg == true) std::cout << __LINE__ << " : Trigger name in the event: "  << trigName << "\twith prescale: " << trigPrescale << std::endl;
	      
	      
	      // ################################
	      // # Save HLT trigger information #
	      // ################################
	      for (unsigned int it = 0; it < TrigTable.size(); it++)
		if (trigName.find(TrigTable[it]) != std::string::npos)
		  {
		    NTuple->TrigTable->push_back(trigName);
		    NTuple->TrigPrescales->push_back(trigPrescale);
		    
		    break;
		  }
	    }
	}
      if (NTuple->TrigTable->size() == 0)
	{
	  NTuple->TrigTable->push_back("NotInTable");
	  NTuple->TrigPrescales->push_back(-1);
	  
	  if (printMsg == true) std::cout << __LINE__ << " : no trigger path in " << parameterFile_ << " have been found in the event" << std::endl;
	}
      else if (printMsg == true) std::cout << __LINE__ << " : some trigger paths in " << parameterFile_ << " have been found in the event" << std::endl;
    }


  if ((doGenReco_ == 1) || (doGenReco_ == 2))
    {
      // Get BeamSpot
      edm::Handle<reco::BeamSpot> beamSpotH;
      iEvent.getByLabel(beamSpot_, beamSpotH);
      reco::BeamSpot beamSpot = *beamSpotH;
      
      // Get primary vertex with BeamSpot constraint: ordered by the scalar sum of the |pT|^2 of the tracks that compose the vertex
      edm::Handle<reco::VertexCollection> recVtx;
      iEvent.getByLabel(vtxSample_, recVtx);
      reco::Vertex bestVtx = *(recVtx->begin());
      reco::Vertex bestVtxReFit;
      
      // Get PAT Muons
      edm::Handle< std::vector<pat::Muon> > thePATMuonHandle;
      iEvent.getByLabel(muonType_, thePATMuonHandle);
      
      // Get PAT Tracks
      edm::Handle< std::vector<pat::GenericParticle> > thePATTrackHandle;
      iEvent.getByLabel(trackType_, thePATTrackHandle);


      if (printMsg == true) std::cout << __LINE__ << " : the event has: " << thePATMuonHandle->size() << " muons AND " << thePATTrackHandle->size() << " tracks" << std::endl;


      // #################################
      // # Search for first valid vertex #
      // #################################
      for (std::vector<reco::Vertex>::const_iterator iVertex = recVtx->begin(); iVertex != recVtx->end(); iVertex++) { bestVtx = *(iVertex); if (bestVtx.isValid() == true) break; }


      if (bestVtx.isValid() == true)
	{
	  // ###########
	  // # Get mu- #
	  // ###########
	  for (std::vector<pat::Muon>::const_iterator iMuonM = thePATMuonHandle->begin(); iMuonM != thePATMuonHandle->end(); iMuonM++)
	    {
	      bool skip = false;


	      // ########################
	      // # Check mu- kinematics #
	      // ########################
	      muTrackm = iMuonM->innerTrack();
	      if ((muTrackm.isNull() == true) || (muTrackm->charge() != -1)) continue;
 
	      const reco::TransientTrack muTrackmTT(muTrackm, &(*bFieldHandle));

		
	      // #########################
	      // # Muon- pT and eta cuts #
	      // #########################
	      pT  = muTrackm->pt();
	      eta = Utility->computeEta (muTrackmTT.track().momentum().x(),muTrackmTT.track().momentum().y(),muTrackmTT.track().momentum().z());
	      if ((pT < (MUMINPT*(1.0-MUVARTOLE))) || (fabs(eta) > (MUMAXETA*(1.0+MUVARTOLE))))
		{
		  if (printMsg == true) std::cout << __LINE__ << " : break --> too low pT of mu- : " << pT << " or too high eta : " << eta << std::endl;
		  break;
		}


	      // ###############################
	      // # Compute mu- DCA to BeamSpot #
	      // ###############################
	      theDCAXBS = muTrackmTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
	      if (theDCAXBS.isValid() == false)
		{
		  if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for mu-" << std::endl;
		  continue;
		}
	      double DCAmumBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
	      double DCAmumBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
	      if (fabs(DCAmumBS) > DCAMUBS)
		{
		  if (printMsg == true) std::cout << __LINE__ << " : continue --> bad absolute impact parameter 2D for mu- : " << DCAmumBS << std::endl;
		  continue;
		}


	      // ###########
	      // # Get mu+ #
	      // ###########
	      for (std::vector<pat::Muon>::const_iterator iMuonP = thePATMuonHandle->begin(); iMuonP != thePATMuonHandle->end(); iMuonP++)
		{
		  // ########################
		  // # Check mu+ kinematics #
		  // ########################
		  muTrackp = iMuonP->innerTrack();
		  if ((muTrackp.isNull() == true) || (muTrackp->charge() != 1)) continue;

		  const reco::TransientTrack muTrackpTT(muTrackp, &(*bFieldHandle));


		  // #########################
		  // # Muon+ pT and eta cuts #
		  // #########################
		  pT  = muTrackp->pt();
		  eta = Utility->computeEta (muTrackpTT.track().momentum().x(),muTrackpTT.track().momentum().y(),muTrackpTT.track().momentum().z());
		  if ((pT < (MUMINPT*(1.0-MUVARTOLE))) || (fabs(eta) > (MUMAXETA*(1.0+MUVARTOLE))))
		    {
		      if (printMsg == true) std::cout << __LINE__ << " : break --> too low pT of mu+ : " << pT << " or too high eta : " << eta << std::endl;
		      break;
		    }


		  // ###############################
		  // # Compute mu+ DCA to BeamSpot #
		  // ###############################
		  theDCAXBS = muTrackpTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
		  if (theDCAXBS.isValid() == false)
		    {
		      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for mu+" << std::endl;
		      continue;
		    }
		  double DCAmupBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
		  double DCAmupBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
		  if (fabs(DCAmupBS) > DCAMUBS)
		    {
		      if (printMsg == true) std::cout << __LINE__ << " : continue --> bad absolute impact parameter 2D for mu+: " << DCAmupBS << std::endl;
		      continue;
		    }


		  // ############################################
		  // # Check goodness of muons closest approach #
		  // ############################################
		  ClosestApp.calculate(muTrackpTT.initialFreeState(),muTrackmTT.initialFreeState());
		  if (ClosestApp.status() == false)
		    {
		      if (printMsg == true) std::cout << __LINE__ << " : continue --> bad status of closest approach" << std::endl;
		      continue;
		    }
		  XingPoint = ClosestApp.crossingPoint();
		  if ((sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y()) > TRKMAXR) || (fabs(XingPoint.z()) > TRKMAXZ))
		    {
		      if (printMsg == true) std::cout << __LINE__ << " : continue --> closest approach crossing point outside the tracker volume" << std::endl;
		      continue;
		    }


		  // #####################################################
		  // # Cut on the mumu 3D-DCA with respect to each other #
		  // #####################################################
		  double mumuDCA = ClosestApp.distance();
		  if (mumuDCA > DCAMUMU)
		    {
		      if (printMsg == true) std::cout << __LINE__ << " : continue --> bad 3D-DCA of mu+(-) with respect to mu-(+): " << mumuDCA << std::endl;
		      continue;
		    }


		  // ############################################
		  // # Cut on the dimuon inviariant mass and pT #
		  // ############################################
		  pT = sqrt((muTrackmTT.track().momentum().x() + muTrackpTT.track().momentum().x()) * (muTrackmTT.track().momentum().x() + muTrackpTT.track().momentum().x()) +
			    (muTrackmTT.track().momentum().y() + muTrackpTT.track().momentum().y()) * (muTrackmTT.track().momentum().y() + muTrackpTT.track().momentum().y()));
		  double MuMuInvMass = Utility->computeInvMass (muTrackmTT.track().momentum().x(),muTrackmTT.track().momentum().y(),muTrackmTT.track().momentum().z(),Utility->muonMass,
								muTrackpTT.track().momentum().x(),muTrackpTT.track().momentum().y(),muTrackpTT.track().momentum().z(),Utility->muonMass);
		  if ((pT < (MINMUMUPT*(1.0-MUVARTOLE))) || (MuMuInvMass < (MINMUMUINVMASS*(1.0-MUVARTOLE))) || (MuMuInvMass > (MAXMUMUINVMASS*(1.0+MUVARTOLE))))
		    {
		      if (printMsg == true) std::cout << __LINE__ << " : continue --> no good mumu pair pT: " << pT << "\tinv. mass: " << MuMuInvMass << std::endl;
		      continue;
		    }




		  // #######################################################
		  // # @@@ Make mu-mu and implement pre-selection cuts @@@ #
		  // #######################################################
		  if (printMsg == true) std::cout << "\n" << __LINE__ << " : @@@ I have 2 good oppositely-charged muons. I'm trying to vertex them @@@" << std::endl;


		  chi = 0.;
		  ndf = 0.;
		  // ####################################################
		  // # Try to vertex the two muons to get dimuon vertex #
		  // ####################################################
		  std::vector<RefCountedKinematicParticle> muonParticles;
		  muonParticles.push_back(partFactory.particle(muTrackmTT,muonMass,chi,ndf,muonMassErr));
		  muonParticles.push_back(partFactory.particle(muTrackpTT,muonMass,chi,ndf,muonMassErr));
	  
		  RefCountedKinematicTree mumuVertexFitTree = PartVtxFitter.fit(muonParticles);
		  if (mumuVertexFitTree->isValid() == false)
		    {
		      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid vertex from the mu+ mu- vertex fit" << std::endl;
		      continue; 
		    }
	  
		  mumuVertexFitTree->movePointerToTheTop();
		  RefCountedKinematicParticle mumu_KP = mumuVertexFitTree->currentParticle();
		  RefCountedKinematicVertex mumu_KV   = mumuVertexFitTree->currentDecayVertex();
		  if (mumu_KV->vertexIsValid() == false)
		    {
		      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid vertex from the mu+ mu- vertex fit" << std::endl;
		      continue;
		    }
	  
		  if (TMath::Prob(static_cast<double>(mumu_KV->chiSquared()), static_cast<int>(rint(mumu_KV->degreesOfFreedom()))) < CLMUMUVTX)
		    {
		      if (printMsg == true) std::cout << __LINE__ << " : continue --> bad vtx CL from mu+ mu- fit: " << TMath::Prob(static_cast<double>(mumu_KV->chiSquared()), static_cast<int>(rint(mumu_KV->degreesOfFreedom()))) << std::endl;
		      continue;
		    }


		  // #########################################################
		  // # Extract the re-fitted tracks after the dimuon vtx fit #
		  // #########################################################
		  mumuVertexFitTree->movePointerToTheTop();

		  mumuVertexFitTree->movePointerToTheFirstChild();
		  RefCountedKinematicParticle refitMum  = mumuVertexFitTree->currentParticle();
		  const reco::TransientTrack refitMumTT = refitMum->refittedTransientTrack();

		  mumuVertexFitTree->movePointerToTheNextChild();
		  RefCountedKinematicParticle refitMup  = mumuVertexFitTree->currentParticle();
		  const reco::TransientTrack refitMupTT = refitMup->refittedTransientTrack();


		  // ########################
		  // # Muon pT and eta cuts #
		  // ########################
		  pT  = sqrt(refitMupTT.track().momentum().x()*refitMupTT.track().momentum().x() + refitMupTT.track().momentum().y()*refitMupTT.track().momentum().y());
		  eta = Utility->computeEta (refitMupTT.track().momentum().x(),refitMupTT.track().momentum().y(),refitMupTT.track().momentum().z());
		  if ((pT < MUMINPT) || (fabs(eta) > MUMAXETA))
		    {
		      if (printMsg == true) std::cout << __LINE__ << " : continue --> too low pT of mu+ : " << pT << " or too high eta : " << eta << std::endl;
		      continue;
		    }

		  pT  = sqrt(refitMumTT.track().momentum().x()*refitMumTT.track().momentum().x() + refitMumTT.track().momentum().y()*refitMumTT.track().momentum().y());
		  eta = Utility->computeEta (refitMumTT.track().momentum().x(),refitMumTT.track().momentum().y(),refitMumTT.track().momentum().z());
		  if ((pT < MUMINPT) || (fabs(eta) > MUMAXETA))
		    {
		      if (printMsg == true) std::cout << __LINE__ << " : continue --> too low pT of mu- : " << pT << " or too high eta : " << eta << std::endl;
		      skip = true;
		      continue;
		    }


		  // ############################################
		  // # Cut on the dimuon inviariant mass and pT #
		  // ############################################
		  pT = sqrt((refitMumTT.track().momentum().x() + refitMupTT.track().momentum().x()) * (refitMumTT.track().momentum().x() + refitMupTT.track().momentum().x()) +
			    (refitMumTT.track().momentum().y() + refitMupTT.track().momentum().y()) * (refitMumTT.track().momentum().y() + refitMupTT.track().momentum().y()));
		  MuMuInvMass = mumu_KP->currentState().mass();
		  if ((pT < MINMUMUPT) || (MuMuInvMass < MINMUMUINVMASS) || (MuMuInvMass > MAXMUMUINVMASS))
		    {
		      if (printMsg == true) std::cout << __LINE__ << " : continue --> no good mumu pair pT: " << pT << "\tinv. mass: " << MuMuInvMass << std::endl;
		      continue;
		    }


		  // ######################################################
		  // # Compute the distance between mumu vtx and BeamSpot #
		  // ######################################################
		  double MuMuLSBS;
		  double MuMuLSBSErr;
		  Utility->computeLS (mumu_KV->position().x(),mumu_KV->position().y(),0.0,
				      beamSpot.position().x(),beamSpot.position().y(),0.0,
				      mumu_KV->error().cxx(),mumu_KV->error().cyy(),0.0,
				      mumu_KV->error().matrix()(0,1),0.0,0.0,
				      beamSpot.covariance()(0,0),beamSpot.covariance()(1,1),0.0,
				      beamSpot.covariance()(0,1),0.0,0.0,
				      &MuMuLSBS,&MuMuLSBSErr);
		  if (MuMuLSBS/MuMuLSBSErr < LSMUMUBS)
		    {
		      if (printMsg == true) std::cout << __LINE__ << " : continue --> bad mumu L/sigma with respect to BeamSpot: " << MuMuLSBS << "+/-" << MuMuLSBSErr << std::endl;
		      continue;
		    }
	  

		  // ###################################################################
		  // # Compute cos(alpha) between mumu momentum and mumuVtx - BeamSpot #
		  // ###################################################################
		  double MuMuCosAlphaBS;
		  double MuMuCosAlphaBSErr;
		  Utility->computeCosAlpha (mumu_KP->currentState().globalMomentum().x(),mumu_KP->currentState().globalMomentum().y(),0.0,
					    mumu_KV->position().x() - beamSpot.position().x(),mumu_KV->position().y() - beamSpot.position().y(),0.0,
					    mumu_KP->currentState().kinematicParametersError().matrix()(3,3),mumu_KP->currentState().kinematicParametersError().matrix()(4,4),0.0,
					    mumu_KP->currentState().kinematicParametersError().matrix()(3,4),0.0,0.0,
					    mumu_KV->error().cxx() + beamSpot.covariance()(0,0),mumu_KV->error().cyy() + beamSpot.covariance()(1,1),0.0,
					    mumu_KV->error().matrix()(0,1) + beamSpot.covariance()(0,1),0.0,0.0,
					    &MuMuCosAlphaBS,&MuMuCosAlphaBSErr);
		  if (MuMuCosAlphaBS < COSALPHAMUMUBS)
		    {
		      if (printMsg == true) std::cout << __LINE__ << " : continue --> bad mumu cos(alpha) with respect to BeamSpot: " << MuMuCosAlphaBS << "+/-" << MuMuCosAlphaBSErr << std::endl;
		      continue;
		    }


		  // ##############
		  // # Get Track- #
		  // ##############
		  for (std::vector<pat::GenericParticle>::const_iterator iTrackM = thePATTrackHandle->begin(); iTrackM != thePATTrackHandle->end(); iTrackM++)
		    {
		      bool skip = false;


		      // ###########################
		      // # Check Track- kinematics #
		      // ###########################
		      Trackm = iTrackM->track();
		      if ((Trackm.isNull() == true) || (Trackm->charge() != -1)) continue;

		      const reco::TransientTrack TrackmTT(Trackm, &(*bFieldHandle));


		      // ##########################
		      // # Track- pT and eta cuts #
		      // ##########################
		      pT = sqrt(TrackmTT.track().momentum().x()*TrackmTT.track().momentum().x() + TrackmTT.track().momentum().y()*TrackmTT.track().momentum().y());
		      if (pT < (MINHADPT*(1.0-HADVARTOLE)))
			{
			  if (printMsg == true) std::cout << __LINE__ << " : break --> too low pT of track- : " << pT << std::endl;
			  break;
			}
		      

		      // ######################################
		      // # Compute K*0 track- DCA to BeamSpot #
		      // ######################################
		      theDCAXBS = TrackmTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
		      if (theDCAXBS.isValid() == false)
			{
			  if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for track-" << std::endl;
			  continue;
			}
		      double DCAKstTrkmBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
		      double DCAKstTrkmBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
		      if (fabs(DCAKstTrkmBS/DCAKstTrkmBSErr) < HADDCASBS)
			{
			  if (printMsg == true) std::cout << __LINE__ << " : continue --> track- DCA/sigma with respect to BeamSpot is too small: " << DCAKstTrkmBS << "+/-" << DCAKstTrkmBSErr << std::endl;
			  continue;
			}
		      
		      
		      // ##############
		      // # Get Track+ #
		      // ##############
		      for (std::vector<pat::GenericParticle>::const_iterator iTrackP = thePATTrackHandle->begin(); iTrackP != thePATTrackHandle->end(); iTrackP++)
			{
			  // ###########################
			  // # Check Track+ kinematics #
			  // ###########################
			  Trackp = iTrackP->track();
			  if ((Trackp.isNull() == true) || (Trackp->charge() != 1)) continue;

			  const reco::TransientTrack TrackpTT(Trackp, &(*bFieldHandle));


			  // ##########################
			  // # Track+ pT and eta cuts #
			  // ##########################
			  pT = sqrt(TrackpTT.track().momentum().x()*TrackpTT.track().momentum().x() + TrackpTT.track().momentum().y()*TrackpTT.track().momentum().y());
			  if (pT < (MINHADPT*(1.0-HADVARTOLE)))
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : break --> too low pT of track+ : " << pT << std::endl;
			      break;
			    }


			  // ######################################
			  // # Compute K*0 track+ DCA to BeamSpot #
			  // ######################################
			  theDCAXBS = TrackpTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
			  if (theDCAXBS.isValid() == false)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for track+" << std::endl;
			      continue;
			    }
			  double DCAKstTrkpBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
			  double DCAKstTrkpBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
			  if (fabs(DCAKstTrkpBS/DCAKstTrkpBSErr) < HADDCASBS)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> track+ DCA/sigma with respect to BeamSpot is too small: " << DCAKstTrkpBS << "+/-" << DCAKstTrkpBSErr << std::endl;
			      continue;
			    }


			  // ##############################################
			  // # Check goodness of hadrons closest approach #
			  // ##############################################
			  ClosestApp.calculate(TrackpTT.initialFreeState(),TrackmTT.initialFreeState());
			  if (ClosestApp.status() == false)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> bad status of closest approach" << std::endl;
			      continue;
			    }
			  XingPoint = ClosestApp.crossingPoint();
			  if ((sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y()) > TRKMAXR) || (fabs(XingPoint.z()) > TRKMAXZ))
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> closest approach crossing point outside the tracker volume" << std::endl;
			      continue;
			    }


			  // ######################################################
			  // # Check if K*0 (OR K*0bar) mass is within acceptance #
			  // ######################################################
			  double kstInvMass    = Utility->computeInvMass (TrackmTT.track().momentum().x(),TrackmTT.track().momentum().y(),TrackmTT.track().momentum().z(),Utility->pionMass,
									  TrackpTT.track().momentum().x(),TrackpTT.track().momentum().y(),TrackpTT.track().momentum().z(),Utility->kaonMass);
			  double kstBarInvMass = Utility->computeInvMass (TrackmTT.track().momentum().x(),TrackmTT.track().momentum().y(),TrackmTT.track().momentum().z(),Utility->kaonMass,
									  TrackpTT.track().momentum().x(),TrackpTT.track().momentum().y(),TrackpTT.track().momentum().z(),Utility->pionMass);
			  if ((fabs(kstInvMass - Utility->kstMass)    > (KSTMASSWINDOW*Utility->kstSigma*(1.0+HADVARTOLE))) &&
			      (fabs(kstBarInvMass - Utility->kstMass) > (KSTMASSWINDOW*Utility->kstSigma*(1.0+HADVARTOLE))))
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> bad K*0 mass: " << kstInvMass << " AND K*0bar mass: " << kstBarInvMass << std::endl;
			      continue;
			    }




			  // ####################################################
			  // # @@@ Make K* and implement pre-selection cuts @@@ #
			  // ####################################################
			  if (printMsg == true) std::cout << "\n" << __LINE__ << " : @@@ I have 2 good oppositely-charged tracks. I'm trying to vertex them @@@" << std::endl;


			  chi = 0.;
			  ndf = 0.;
			  // ##############################################################################
			  // # Try to vertex the two Tracks to get K*0 vertex: pion = track- | k = track+ #
			  // ##############################################################################
			  std::vector<RefCountedKinematicParticle> kstParticles;
			  kstParticles.push_back(partFactory.particle(TrackmTT,pionMass,chi,ndf,pionMassErr));
			  kstParticles.push_back(partFactory.particle(TrackpTT,kaonMass,chi,ndf,kaonMassErr));
			  
			  RefCountedKinematicTree kstVertexFitTree = PartVtxFitter.fit(kstParticles);
			  if (kstVertexFitTree->isValid() == false)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid vertex from the K*0 vertex fit" << std::endl;
			      continue;
			    }
			  
			  kstVertexFitTree->movePointerToTheTop();
			  RefCountedKinematicParticle kst_KP = kstVertexFitTree->currentParticle();
			  RefCountedKinematicVertex kst_KV   = kstVertexFitTree->currentDecayVertex();
			  if (kst_KV->vertexIsValid() == false)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid vertex from the K*0 vertex fit" << std::endl;
			      continue;
			    }


			  chi = 0.;
			  ndf = 0.;
			  // #################################################################################
			  // # Try to vertex the two Tracks to get K*0bar vertex: pion = track+ | k = track- #
			  // #################################################################################
			  std::vector<RefCountedKinematicParticle> kstBarParticles;
			  kstBarParticles.push_back(partFactory.particle(TrackmTT,kaonMass,chi,ndf,kaonMassErr));
			  kstBarParticles.push_back(partFactory.particle(TrackpTT,pionMass,chi,ndf,pionMassErr));
			  
			  RefCountedKinematicTree kstBarVertexFitTree = PartVtxFitter.fit(kstBarParticles);
			  if (kstBarVertexFitTree->isValid() == false)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid vertex from the K*0bar vertex fit" << std::endl;
			      continue;
			    }
			  
			  kstBarVertexFitTree->movePointerToTheTop();
			  RefCountedKinematicParticle kstBar_KP = kstBarVertexFitTree->currentParticle();
			  RefCountedKinematicVertex kstBar_KV   = kstBarVertexFitTree->currentDecayVertex();
			  if (kstBar_KV->vertexIsValid() == false)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid vertex from the K*0bar vertex fit" << std::endl;
			      continue;
			    }


			  // ######################################################
			  // # Check if K*0 (OR K*0bar) mass is within acceptance #
			  // ######################################################
			  kstInvMass    = kst_KP->currentState().mass();
			  kstBarInvMass = kstBar_KP->currentState().mass();
			  if ((fabs(kstInvMass - Utility->kstMass) > KSTMASSWINDOW*Utility->kstSigma) && (fabs(kstBarInvMass - Utility->kstMass) > KSTMASSWINDOW*Utility->kstSigma))
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> bad K*0 mass: " << kstInvMass << " AND K*0bar mass: " << kstBarInvMass << std::endl;
			      continue;
			    }


			  // ###########################################################
			  // # Extract the re-fitted tracks after the dihadron vtx fit #
			  // ###########################################################
			  kstVertexFitTree->movePointerToTheTop();

			  kstVertexFitTree->movePointerToTheFirstChild();
			  RefCountedKinematicParticle refitTrkm  = kstVertexFitTree->currentParticle();
			  const reco::TransientTrack refitTrkmTT = refitTrkm->refittedTransientTrack();

			  kstVertexFitTree->movePointerToTheNextChild();
			  RefCountedKinematicParticle refitTrkp  = kstVertexFitTree->currentParticle();
			  const reco::TransientTrack refitTrkpTT = refitTrkp->refittedTransientTrack();


			  // ##########################
			  // # Hadron pT and eta cuts #
			  // ##########################
			  pT = sqrt(refitTrkpTT.track().momentum().x()*refitTrkpTT.track().momentum().x() + refitTrkpTT.track().momentum().y()*refitTrkpTT.track().momentum().y());
			  if (pT < MINHADPT)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : break --> too low pT of track+ : " << pT << std::endl;
			      continue;
			    }

			  pT = sqrt(refitTrkmTT.track().momentum().x()*refitTrkmTT.track().momentum().x() + refitTrkmTT.track().momentum().y()*refitTrkmTT.track().momentum().y());
			  if (pT < MINHADPT)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : break --> too low pT of track- : " << pT << std::endl;
			      skip = true;
			      continue;
			    }


			  // #################################################
			  // # Check if the hadron Track- is actually a muon #
			  // #################################################
			  MuMCat.clear();
			  MuMCat = "NotMatched";
			  for (std::vector<pat::Muon>::const_iterator iMuon = thePATMuonHandle->begin(); iMuon != thePATMuonHandle->end(); iMuon++)
			    {
			      muTrackTmp = iMuon->innerTrack();
			      if ((muTrackTmp.isNull() == true) || (muTrackTmp->charge() != -1)) continue;
			      if (Trackm == muTrackTmp)
				{
				  MuMCat.clear();
				  MuMCat.append(getMuCat(*iMuon));
				  if (printMsg == true) std::cout << __LINE__ << " : negative charged hadron is actually a muon (momentum: " << Trackm->p() << ") whose category is: " << MuMCat.c_str() << std::endl;
				  break;
				}
			    }


			  // #################################################
			  // # Check if the hadron Track+ is actually a muon #
			  // #################################################
			  MuPCat.clear();
			  MuPCat = "NotMatched";
			  for (std::vector<pat::Muon>::const_iterator iMuon = thePATMuonHandle->begin(); iMuon != thePATMuonHandle->end(); iMuon++)
			    {
			      muTrackTmp = iMuon->innerTrack();
			      if ((muTrackTmp.isNull() == true) || (muTrackTmp->charge() != 1)) continue;
			      if (Trackp == muTrackTmp)
				{
				  MuPCat.clear();
				  MuPCat.append(getMuCat(*iMuon));
				  if (printMsg == true) std::cout << __LINE__ << " : positive charged hadron is actually a muon (momentum: " << Trackp->p() << ") whose category is: " << MuPCat.c_str() << std::endl;
				  break;
				}
			    }




			  // ####################################################
			  // # @@@ Make B0 and implement pre-selection cuts @@@ #
			  // ####################################################
			  if (printMsg == true) std::cout << "\n" << __LINE__ << " : @@@ I have 4 good charged tracks. I'm trying to vertex them @@@" << std::endl;


			  chi = 0.;
			  ndf = 0.;
			  // #####################################
			  // # B0 vertex fit with vtx constraint #
			  // #####################################
			  std::vector<RefCountedKinematicParticle> bParticles;
			  bParticles.push_back(partFactory.particle(muTrackmTT,muonMass,chi,ndf,muonMassErr));
			  bParticles.push_back(partFactory.particle(muTrackpTT,muonMass,chi,ndf,muonMassErr));
			  bParticles.push_back(partFactory.particle(TrackmTT,pionMass,chi,ndf,pionMassErr));
			  bParticles.push_back(partFactory.particle(TrackpTT,kaonMass,chi,ndf,kaonMassErr));

			  RefCountedKinematicTree bVertexFitTree = PartVtxFitter.fit(bParticles);
			  if (bVertexFitTree->isValid() == false)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid vertex from the B0 vertex fit" << std::endl;
			      continue;
			    }

			  bVertexFitTree->movePointerToTheTop();
			  RefCountedKinematicParticle b_KP = bVertexFitTree->currentParticle();
			  RefCountedKinematicVertex b_KV   = bVertexFitTree->currentDecayVertex();
			  if (b_KV->vertexIsValid() == false)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid vertex from the B0 vertex fit" << std::endl;
			      continue;
			    }


			  chi = 0.;
			  ndf = 0.;
			  // ########################################
			  // # B0bar vertex fit with vtx constraint #
			  // ########################################
			  std::vector<RefCountedKinematicParticle> bBarParticles;
			  bBarParticles.push_back(partFactory.particle(muTrackmTT,muonMass,chi,ndf,muonMassErr));
			  bBarParticles.push_back(partFactory.particle(muTrackpTT,muonMass,chi,ndf,muonMassErr));
			  bBarParticles.push_back(partFactory.particle(TrackmTT,kaonMass,chi,ndf,kaonMassErr));
			  bBarParticles.push_back(partFactory.particle(TrackpTT,pionMass,chi,ndf,pionMassErr));
				      
			  RefCountedKinematicTree bBarVertexFitTree = PartVtxFitter.fit(bBarParticles);
			  if (bBarVertexFitTree->isValid() == false)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid vertex from the B0bar vertex fit" << std::endl;
			      continue;
			    }

			  bBarVertexFitTree->movePointerToTheTop();
			  RefCountedKinematicParticle bBar_KP = bBarVertexFitTree->currentParticle();
			  RefCountedKinematicVertex bBar_KV   = bBarVertexFitTree->currentDecayVertex();
			  if (bBar_KV->vertexIsValid() == false)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid vertex from the B0bar vertex fit" << std::endl;
			      continue;
			    }


			  // ########################
			  // # Cuts on B0 AND B0bar #
			  // ########################
			  if (((b_KP->currentState().mass() < B0MASSLOWLIMIT) || (b_KP->currentState().mass() > B0MASSUPLIMIT)) &&
			      ((bBar_KP->currentState().mass() < B0MASSLOWLIMIT) || (bBar_KP->currentState().mass() > B0MASSUPLIMIT)))
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> bad B0 mass: " << b_KP->currentState().mass() << " AND B0bar mass: " << bBar_KP->currentState().mass() << std::endl;
			      continue;
			    }
			  if ((TMath::Prob(static_cast<double>(b_KV->chiSquared()), static_cast<int>(rint(b_KV->degreesOfFreedom()))) < CLB0VTX) &&
			      (TMath::Prob(static_cast<double>(bBar_KV->chiSquared()), static_cast<int>(rint(bBar_KV->degreesOfFreedom()))) < CLB0VTX))
			    {
			      if (printMsg == true)
				{
				  std::cout << __LINE__ << " : continue --> bad vtx CL from B0 fit: " << TMath::Prob(static_cast<double>(b_KV->chiSquared()), static_cast<int>(rint(b_KV->degreesOfFreedom())));
				  std::cout << " AND bad vtx CL from B0bar fit: " << TMath::Prob(static_cast<double>(bBar_KV->chiSquared()), static_cast<int>(rint(bBar_KV->degreesOfFreedom()))) << std::endl;
				}
			      continue;
			    }


			  // ###############################
			  // # Cuts on B0 L/sigma BeamSpot #
			  // ###############################
			  Utility->computeLS (b_KV->position().x(),b_KV->position().y(),0.0,
					      beamSpot.position().x(),beamSpot.position().y(),0.0,
					      b_KV->error().cxx(),b_KV->error().cyy(),0.0,
					      b_KV->error().matrix()(0,1),0.0,0.0,
					      beamSpot.covariance()(0,0),beamSpot.covariance()(1,1),0.0,
					      beamSpot.covariance()(0,1),0.0,0.0,
					      &LSBS,&LSBSErr);


			  // ##############################
			  // # Compute B0 DCA to BeamSpot #
			  // ##############################
			  theDCAXBS = b_KP->refittedTransientTrack().trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
			  if (theDCAXBS.isValid() == false)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for B0" << std::endl;
			      continue;
			    }
			  double DCAB0BS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
			  double DCAB0BSErr = theDCAXBS.perigeeError().transverseImpactParameterError();

		
			  // #####################################
			  // # Compute B0 cos(alpha) to BeamSpot #
			  // #####################################
			  Utility->computeCosAlpha (b_KP->currentState().globalMomentum().x(),b_KP->currentState().globalMomentum().y(),0.0,
						    b_KV->position().x() - beamSpot.position().x(),b_KV->position().y() - beamSpot.position().y(),0.0,
						    b_KP->currentState().kinematicParametersError().matrix()(3,3),b_KP->currentState().kinematicParametersError().matrix()(4,4),0.0,
						    b_KP->currentState().kinematicParametersError().matrix()(3,4),0.0,0.0,
						    b_KV->error().cxx() + beamSpot.covariance()(0,0),b_KV->error().cyy() + beamSpot.covariance()(1,1),0.0,
						    b_KV->error().matrix()(0,1) + beamSpot.covariance()(0,1),0.0,0.0,
						    &cosAlphaBS,&cosAlphaBSErr);




			  // #############################################################################################################
			  // # If the tracks in the B reco candidate belongs to the primary vertex, then remove them and redo the Vertex #
			  // #############################################################################################################
			  std::vector<reco::TransientTrack> vertexTracks;
			  for (std::vector<reco::TrackBaseRef>::const_iterator iTrack = bestVtx.tracks_begin(); iTrack != bestVtx.tracks_end(); iTrack++)
			    {
			      // Compare primary tracks to check for matches with B candidate
			      reco::TrackRef trackRef = iTrack->castTo<reco::TrackRef>();
					  
			      if ((Trackm != trackRef) && (Trackp != trackRef) && (muTrackm != trackRef) && (muTrackp != trackRef))
				{
				  reco::TransientTrack TT(trackRef, &(*bFieldHandle));
				  vertexTracks.push_back(TT);
				}
			      else if (printMsg == true) std::cout << __LINE__ << " : some of the B0 reco candidate tracks belong to the Pri.Vtx" << std::endl;
			    }

			  if (vertexTracks.size() < 2)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> number of tracks of the new Pri.Vtx is too small: " << vertexTracks.size() << std::endl;
			      continue;
			    }
				      
			  TransientVertex bestTransVtx = theVtxFitter.vertex(vertexTracks);
			  if (bestTransVtx.isValid() == false)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid new Pri.Vtx" << std::endl;
			      continue;
			    }
			  bestVtxReFit = (reco::Vertex)bestTransVtx;




			  // ##########################
			  // # Compute mu- DCA to Vtx #
			  // ##########################
			  theDCAXVtx = IPTools::absoluteImpactParameter3D(muTrackmTT, bestVtxReFit);
			  if (theDCAXVtx.first == false)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 3D for mu-" << std::endl;
			      continue;
			    }
			  double DCAmumVtx    = theDCAXVtx.second.value();
			  double DCAmumVtxErr = theDCAXVtx.second.error();
			  
			  					  
			  // ##########################
			  // # Compute mu+ DCA to Vtx #
			  // ##########################
			  theDCAXVtx = IPTools::absoluteImpactParameter3D(muTrackpTT, bestVtxReFit);
			  if (theDCAXVtx.first == false)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 3D for mu+" << std::endl;
			      continue;
			    }
			  double DCAmupVtx    = theDCAXVtx.second.value();
			  double DCAmupVtxErr = theDCAXVtx.second.error();

			
			  // #################################
			  // # Compute K*0 track- DCA to Vtx #
			  // #################################
			  theDCAXVtx = IPTools::absoluteImpactParameter3D(TrackmTT, bestVtxReFit);
			  if (theDCAXVtx.first == false)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 3D for track-" << std::endl;
			      continue;
			    }
			  double DCAKstTrkmVtx    = theDCAXVtx.second.value();
			  double DCAKstTrkmVtxErr = theDCAXVtx.second.error();


			  // #################################
			  // # Compute K*0 track+ DCA to Vtx #
			  // #################################
			  theDCAXVtx = IPTools::absoluteImpactParameter3D(TrackpTT, bestVtxReFit);
			  if (theDCAXVtx.first == false)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 3D for track+" << std::endl;
			      continue;
			    }
			  double DCAKstTrkpVtx    = theDCAXVtx.second.value();
			  double DCAKstTrkpVtxErr = theDCAXVtx.second.error();


			  // #########################
			  // # Compute B0 DCA to Vtx #
			  // #########################
			  theDCAXVtx = IPTools::absoluteImpactParameter3D(b_KP->refittedTransientTrack(), bestVtxReFit);
			  if (theDCAXVtx.first == false)
			    {
			      if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 3D for B0" << std::endl;
			      continue;
			    }
			  double DCAB0Vtx    = theDCAXVtx.second.value();
			  double DCAB0VtxErr = theDCAXVtx.second.error();


			  // ################################
			  // # Compute B0 cos(alpha) to Vtx #
			  // ################################
			  Utility->computeCosAlpha (b_KP->currentState().globalMomentum().x(),b_KP->currentState().globalMomentum().y(),b_KP->currentState().globalMomentum().z(),
						    b_KV->position().x() - bestVtxReFit.x(),b_KV->position().y() - bestVtxReFit.y(), b_KV->position().z() - bestVtxReFit.z(),
						    b_KP->currentState().kinematicParametersError().matrix()(3,3),b_KP->currentState().kinematicParametersError().matrix()(4,4),b_KP->currentState().kinematicParametersError().matrix()(5,5),
						    b_KP->currentState().kinematicParametersError().matrix()(3,4),b_KP->currentState().kinematicParametersError().matrix()(3,5),b_KP->currentState().kinematicParametersError().matrix()(4,5),
						    b_KV->error().cxx() + bestVtxReFit.error()(0,0),b_KV->error().cyy() + bestVtxReFit.error()(1,1),b_KV->error().czz() + bestVtxReFit.error()(2,2),
						    b_KV->error().matrix()(0,1) + bestVtxReFit.error()(0,1),b_KV->error().matrix()(0,2) + bestVtxReFit.error()(0,2),b_KV->error().matrix()(1,2) + bestVtxReFit.error()(1,2),
						    &cosAlphaVtx,&cosAlphaVtxErr);


			  // #############################
			  // # Compute B0 L/sigma to Vtx #
			  // #############################
			  Utility->computeLS (b_KV->position().x(),b_KV->position().y(),b_KV->position().z(),
					      bestVtxReFit.x(),bestVtxReFit.y(),bestVtxReFit.z(),
					      b_KV->error().cxx(),b_KV->error().cyy(),b_KV->error().czz(),
					      b_KV->error().matrix()(0,1),b_KV->error().matrix()(0,2),b_KV->error().matrix()(1,2),
					      bestVtxReFit.error()(0,0),bestVtxReFit.error()(1,1),bestVtxReFit.error()(2,2),
					      bestVtxReFit.error()(0,1),bestVtxReFit.error()(0,2),bestVtxReFit.error()(1,2),
					      &LSVtx,&LSVtxErr);


			  // #####################################################
			  // # Compute ctau = mB * (BVtx - PVtx) dot pB / |pB|^2 #
			  // #####################################################
			  GlobalVector Bmomentum = b_KP->currentState().globalMomentum();
			  double Bmomentum2 = Bmomentum.dot(Bmomentum);
			  AlgebraicSymMatrix33 BmomentumErr2;
			  BmomentumErr2(0,0) = b_KP->currentState().kinematicParametersError().matrix()(3,3);
			  BmomentumErr2(0,1) = b_KP->currentState().kinematicParametersError().matrix()(3,4);
			  BmomentumErr2(0,2) = b_KP->currentState().kinematicParametersError().matrix()(3,5);
			  BmomentumErr2(1,0) = b_KP->currentState().kinematicParametersError().matrix()(4,3);
			  BmomentumErr2(1,1) = b_KP->currentState().kinematicParametersError().matrix()(4,4);
			  BmomentumErr2(1,2) = b_KP->currentState().kinematicParametersError().matrix()(4,5);
			  BmomentumErr2(2,0) = b_KP->currentState().kinematicParametersError().matrix()(5,3);
			  BmomentumErr2(2,1) = b_KP->currentState().kinematicParametersError().matrix()(5,4);
			  BmomentumErr2(2,2) = b_KP->currentState().kinematicParametersError().matrix()(5,5);

			  GlobalVector BVtxPVtxSep = GlobalPoint(b_KV->position()) - GlobalPoint(bestVtxReFit.position().x(), bestVtxReFit.position().y(), bestVtxReFit.position().z());
			  double BVtxPVtxSepDOTBmomentum = BVtxPVtxSep.dot(Bmomentum);

			  double bctauPVBS = Utility->B0Mass * BVtxPVtxSepDOTBmomentum / Bmomentum2;
			  double bctauPVBSErr = sqrt(Utility->B0Mass * Utility->B0Mass / (Bmomentum2 * Bmomentum2) *
					   
						     (Bmomentum.x() * Bmomentum.x() * bestVtxReFit.error()(0,0) +
						      Bmomentum.y() * Bmomentum.y() * bestVtxReFit.error()(1,1) +
						      Bmomentum.z() * Bmomentum.z() * bestVtxReFit.error()(2,2) +
					    
						      ((BVtxPVtxSep.x()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.x()) *
						       (BVtxPVtxSep.x()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.x()) *
						       BmomentumErr2(0,0) +
						       (BVtxPVtxSep.y()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.y()) *
						       (BVtxPVtxSep.y()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.y()) *
						       BmomentumErr2(1,1) +
						       (BVtxPVtxSep.z()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.z()) *
						       (BVtxPVtxSep.z()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.z()) *	
						       BmomentumErr2(2,2) +
					   
						       (BVtxPVtxSep.x()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.x()) * (BVtxPVtxSep.y()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.y()) * 2.*BmomentumErr2(0,1) +
						       (BVtxPVtxSep.x()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.x()) * (BVtxPVtxSep.z()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.z()) * 2.*BmomentumErr2(0,2) +
						       (BVtxPVtxSep.y()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.y()) * (BVtxPVtxSep.z()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.z()) * 2.*BmomentumErr2(1,2)) /

						      (Bmomentum2 * Bmomentum2)));




			  // #######################################
			  // # @@@ Fill B0-candidate variables @@@ #
			  // #######################################
			  if (printMsg == true) std::cout << __LINE__ << " : @@@ Filling B0 candidate variables @@@\n\n" << std::endl;


			  // ############
			  // # Save: B0 #
			  // ############
			  NTuple->nB++;

			  NTuple->bMass->push_back(b_KP->currentState().mass());
			  NTuple->bMassE->push_back(sqrt(b_KP->currentState().kinematicParametersError().matrix()(6,6)));
			  NTuple->bBarMass->push_back(bBar_KP->currentState().mass());
			  NTuple->bBarMassE->push_back(sqrt(bBar_KP->currentState().kinematicParametersError().matrix()(6,6)));

			  NTuple->bPx->push_back(b_KP->currentState().globalMomentum().x());
			  NTuple->bPy->push_back(b_KP->currentState().globalMomentum().y());
			  NTuple->bPz->push_back(b_KP->currentState().globalMomentum().z());	      

			  NTuple->bVtxCL->push_back(TMath::Prob(static_cast<double>(b_KV->chiSquared()), static_cast<int>(rint(b_KV->degreesOfFreedom()))));
			  NTuple->bVtxX->push_back(b_KV->position().x());
			  NTuple->bVtxY->push_back(b_KV->position().y());
			  NTuple->bVtxZ->push_back(b_KV->position().z());
	
			  NTuple->bCosAlphaVtx->push_back(cosAlphaVtx);
			  NTuple->bCosAlphaVtxE->push_back(cosAlphaVtxErr);
			  NTuple->bCosAlphaBS->push_back(cosAlphaBS);
			  NTuple->bCosAlphaBSE->push_back(cosAlphaBSErr);

			  NTuple->bLVtx->push_back(LSVtx);
			  NTuple->bLVtxE->push_back(LSVtxErr);
			  NTuple->bLBS->push_back(LSBS);
			  NTuple->bLBSE->push_back(LSBSErr);

			  NTuple->bDCAVtx->push_back(DCAB0Vtx);
			  NTuple->bDCAVtxE->push_back(DCAB0VtxErr);
			  NTuple->bDCABS->push_back(DCAB0BS);
			  NTuple->bDCABSE->push_back(DCAB0BSErr);
 
			  NTuple->bctauPVBS->push_back(bctauPVBS);
			  NTuple->bctauPVBSE->push_back(bctauPVBSErr);


			  // #############
			  // # Save: K*0 #
			  // #############
			  NTuple->kstMass->push_back(kstInvMass);
			  NTuple->kstMassE->push_back(sqrt(kst_KP->currentState().kinematicParametersError().matrix()(6,6)));
			  NTuple->kstBarMass->push_back(kstBarInvMass);
			  NTuple->kstBarMassE->push_back(sqrt(kstBar_KP->currentState().kinematicParametersError().matrix()(6,6)));

			  NTuple->kstPx->push_back(kst_KP->currentState().globalMomentum().x());
			  NTuple->kstPy->push_back(kst_KP->currentState().globalMomentum().y());
			  NTuple->kstPz->push_back(kst_KP->currentState().globalMomentum().z());

			  NTuple->kstVtxCL->push_back(TMath::Prob(static_cast<double>(kst_KV->chiSquared()), static_cast<int>(rint(kst_KV->degreesOfFreedom()))));
			  NTuple->kstVtxX->push_back(kst_KV->position().x());
			  NTuple->kstVtxY->push_back(kst_KV->position().y());
			  NTuple->kstVtxZ->push_back(kst_KV->position().z());


			  // #################
			  // # Save: mu+ mu- #
			  // #################
			  NTuple->mumuMass->push_back(MuMuInvMass);
			  NTuple->mumuMassE->push_back(sqrt(mumu_KP->currentState().kinematicParametersError().matrix()(6,6)));

			  NTuple->mumuPx->push_back(mumu_KP->currentState().globalMomentum().x());
			  NTuple->mumuPy->push_back(mumu_KP->currentState().globalMomentum().y());
			  NTuple->mumuPz->push_back(mumu_KP->currentState().globalMomentum().z());

			  NTuple->mumuVtxCL->push_back(TMath::Prob(static_cast<double>(mumu_KV->chiSquared()), static_cast<int>(rint(mumu_KV->degreesOfFreedom()))));
			  NTuple->mumuVtxX->push_back(mumu_KV->position().x());
			  NTuple->mumuVtxY->push_back(mumu_KV->position().y());
			  NTuple->mumuVtxZ->push_back(mumu_KV->position().z());

			  NTuple->mumuCosAlphaBS->push_back(MuMuCosAlphaBS);
			  NTuple->mumuCosAlphaBSE->push_back(MuMuCosAlphaBSErr);
			  NTuple->mumuLBS->push_back(MuMuLSBS);
			  NTuple->mumuLBSE->push_back(MuMuLSBSErr);
			  NTuple->mumuDCA->push_back(mumuDCA);


			  // #############
			  // # Save: mu- #
			  // #############
			  NTuple->mumHighPurity->push_back(muTrackm->quality(reco::Track::highPurity));
			  NTuple->mumCL->push_back(TMath::Prob(muTrackmTT.chi2(), static_cast<int>(rint(muTrackmTT.ndof()))));
			  NTuple->mumNormChi2->push_back(muTrackm->normalizedChi2());
			  NTuple->mumPx->push_back(refitMumTT.track().momentum().x());
			  NTuple->mumPy->push_back(refitMumTT.track().momentum().y());
			  NTuple->mumPz->push_back(refitMumTT.track().momentum().z());

			  NTuple->mumDCAVtx->push_back(DCAmumVtx);
			  NTuple->mumDCAVtxE->push_back(DCAmumVtxErr);
			  NTuple->mumDCABS->push_back(DCAmumBS);
			  NTuple->mumDCABSE->push_back(DCAmumBSErr);

			  NTuple->mumKinkChi2->push_back(iMuonM->combinedQuality().trkKink);
			  NTuple->mumFracHits->push_back(static_cast<double>(muTrackm->hitPattern().numberOfValidHits()) / static_cast<double>(muTrackm->hitPattern().numberOfValidHits() +
																	       muTrackm->hitPattern().numberOfLostHits() +
																	       muTrackm->trackerExpectedHitsInner().numberOfLostHits() +
																	       muTrackm->trackerExpectedHitsOuter().numberOfLostHits()));
			  theDCAXVtx = IPTools::absoluteTransverseImpactParameter(muTrackmTT, bestVtxReFit);
			  NTuple->mumdxyVtx->push_back(theDCAXVtx.second.value());
			  NTuple->mumdzVtx->push_back(muTrackmTT.track().dz(bestVtxReFit.position()));

			  NTuple->mumCat->push_back(getMuCat(*iMuonM));

			  NTuple->mumNPixHits->push_back(muTrackm->hitPattern().numberOfValidPixelHits());
			  NTuple->mumNPixLayers->push_back(muTrackm->hitPattern().pixelLayersWithMeasurement());  
			  NTuple->mumNTrkHits->push_back(muTrackm->hitPattern().numberOfValidTrackerHits());
			  NTuple->mumNTrkLayers->push_back(muTrackm->hitPattern().trackerLayersWithMeasurement());
 			  if (iMuonM->isGlobalMuon() == true) NTuple->mumNMuonHits->push_back(iMuonM->globalTrack()->hitPattern().numberOfValidMuonHits());
			  else NTuple->mumNMuonHits->push_back(0);
			  NTuple->mumNMatchStation->push_back(iMuonM->numberOfMatchedStations());

			  tmpString.clear();
			  const pat::Muon* patMuonM = &(*iMuonM);
			  for (unsigned int i = 0; i < TrigTable.size(); i++)
			    {
			      myString.clear(); myString.str(""); myString << TrigTable[i].c_str() << "*";
			      if (patMuonM->triggerObjectMatchesByPath(myString.str().c_str()).empty() == false) tmpString.append(TrigTable[i]+" ");
			    }
			  if (tmpString.size() == 0) tmpString.append("NotInTable");
			  NTuple->mumTrig->push_back(tmpString);


			  // #############
			  // # Save: mu+ #
			  // #############
			  NTuple->mupHighPurity->push_back(muTrackp->quality(reco::Track::highPurity));
			  NTuple->mupCL->push_back(TMath::Prob(muTrackpTT.chi2(), static_cast<int>(rint(muTrackpTT.ndof()))));
			  NTuple->mupNormChi2->push_back(muTrackp->normalizedChi2());
			  NTuple->mupPx->push_back(refitMupTT.track().momentum().x());
			  NTuple->mupPy->push_back(refitMupTT.track().momentum().y());
			  NTuple->mupPz->push_back(refitMupTT.track().momentum().z());

			  NTuple->mupDCAVtx->push_back(DCAmupVtx);
			  NTuple->mupDCAVtxE->push_back(DCAmupVtxErr);
			  NTuple->mupDCABS->push_back(DCAmupBS);
			  NTuple->mupDCABSE->push_back(DCAmupBSErr);
			  
			  NTuple->mupKinkChi2->push_back(iMuonP->combinedQuality().trkKink);
			  NTuple->mupFracHits->push_back(static_cast<double>(muTrackp->hitPattern().numberOfValidHits()) / static_cast<double>(muTrackp->hitPattern().numberOfValidHits() +
																	       muTrackp->hitPattern().numberOfLostHits() +
																	       muTrackp->trackerExpectedHitsInner().numberOfLostHits() +
																	       muTrackp->trackerExpectedHitsOuter().numberOfLostHits()));
			  theDCAXVtx = IPTools::absoluteTransverseImpactParameter(muTrackpTT, bestVtxReFit);
			  NTuple->mupdxyVtx->push_back(theDCAXVtx.second.value());
			  NTuple->mupdzVtx->push_back(muTrackpTT.track().dz(bestVtxReFit.position()));

			  NTuple->mupCat->push_back(getMuCat(*iMuonP));

			  NTuple->mupNPixHits->push_back(muTrackp->hitPattern().numberOfValidPixelHits());
			  NTuple->mupNPixLayers->push_back(muTrackp->hitPattern().pixelLayersWithMeasurement());  
			  NTuple->mupNTrkHits->push_back(muTrackp->hitPattern().numberOfValidTrackerHits());
			  NTuple->mupNTrkLayers->push_back(muTrackp->hitPattern().trackerLayersWithMeasurement());
			  if (iMuonP->isGlobalMuon() == true) NTuple->mupNMuonHits->push_back(iMuonP->globalTrack()->hitPattern().numberOfValidMuonHits());
			  else NTuple->mupNMuonHits->push_back(0);
			  NTuple->mupNMatchStation->push_back(iMuonP->numberOfMatchedStations());

 			  tmpString.clear();
			  const pat::Muon* patMuonP = &(*iMuonP);
			  for (unsigned int i = 0; i < TrigTable.size(); i++)
			    {
			      myString.clear(); myString.str(""); myString << TrigTable[i].c_str() << "*";
			      if (patMuonP->triggerObjectMatchesByPath(myString.str().c_str()).empty() == false) tmpString.append(TrigTable[i]+" ");
			    }
			  if (tmpString.size() == 0) tmpString.append("NotInTable");
			  NTuple->mupTrig->push_back(tmpString);


			  // ################
			  // # Save: Track- #
			  // ################
			  NTuple->kstTrkmHighPurity->push_back(Trackm->quality(reco::Track::highPurity));
			  NTuple->kstTrkmCL->push_back(TMath::Prob(TrackmTT.chi2(), static_cast<int>(rint(TrackmTT.ndof()))));
			  NTuple->kstTrkmNormChi2->push_back(Trackm->normalizedChi2());
			  NTuple->kstTrkmPx->push_back(refitTrkmTT.track().momentum().x());
			  NTuple->kstTrkmPy->push_back(refitTrkmTT.track().momentum().y());
			  NTuple->kstTrkmPz->push_back(refitTrkmTT.track().momentum().z());

			  NTuple->kstTrkmDCAVtx->push_back(DCAKstTrkmVtx);
			  NTuple->kstTrkmDCAVtxE->push_back(DCAKstTrkmVtxErr);
			  NTuple->kstTrkmDCABS->push_back(DCAKstTrkmBS);
			  NTuple->kstTrkmDCABSE->push_back(DCAKstTrkmBSErr);

			  NTuple->kstTrkmFracHits->push_back(static_cast<double>(Trackm->hitPattern().numberOfValidHits()) / static_cast<double>(Trackm->hitPattern().numberOfValidHits() +
																		 Trackm->hitPattern().numberOfLostHits() +
																		 Trackm->trackerExpectedHitsInner().numberOfLostHits()));
			  theDCAXVtx = IPTools::absoluteTransverseImpactParameter(TrackmTT, bestVtxReFit);
			  NTuple->kstTrkmdxyVtx->push_back(theDCAXVtx.second.value());
			  NTuple->kstTrkmdzVtx->push_back(TrackmTT.track().dz(bestVtxReFit.position()));

			  NTuple->kstTrkmMuMatch->push_back(MuMCat);

			  // I do NOT include the number of missing outer hits because the hadron might interact
			  NTuple->kstTrkmNPixHits->push_back(Trackm->hitPattern().numberOfValidPixelHits());
			  NTuple->kstTrkmNPixLayers->push_back(Trackm->hitPattern().pixelLayersWithMeasurement());  
			  NTuple->kstTrkmNTrkHits->push_back(Trackm->hitPattern().numberOfValidTrackerHits());
			  NTuple->kstTrkmNTrkLayers->push_back(Trackm->hitPattern().trackerLayersWithMeasurement());


			  // ################
			  // # Save: Track+ #
			  // ################
			  NTuple->kstTrkpHighPurity->push_back(Trackp->quality(reco::Track::highPurity));
			  NTuple->kstTrkpCL->push_back(TMath::Prob(TrackpTT.chi2(), static_cast<int>(rint(TrackpTT.ndof()))));
			  NTuple->kstTrkpNormChi2->push_back(Trackp->normalizedChi2());
			  NTuple->kstTrkpPx->push_back(refitTrkpTT.track().momentum().x());
			  NTuple->kstTrkpPy->push_back(refitTrkpTT.track().momentum().y());
			  NTuple->kstTrkpPz->push_back(refitTrkpTT.track().momentum().z());

			  NTuple->kstTrkpDCAVtx->push_back(DCAKstTrkpVtx);
			  NTuple->kstTrkpDCAVtxE->push_back(DCAKstTrkpVtxErr);
			  NTuple->kstTrkpDCABS->push_back(DCAKstTrkpBS);
			  NTuple->kstTrkpDCABSE->push_back(DCAKstTrkpBSErr);

			  NTuple->kstTrkpFracHits->push_back(static_cast<double>(Trackp->hitPattern().numberOfValidHits()) / static_cast<double>(Trackp->hitPattern().numberOfValidHits() +
																		 Trackp->hitPattern().numberOfLostHits() +
																		 Trackp->trackerExpectedHitsInner().numberOfLostHits()));
			  theDCAXVtx = IPTools::absoluteTransverseImpactParameter(TrackpTT, bestVtxReFit);
			  NTuple->kstTrkpdxyVtx->push_back(theDCAXVtx.second.value());
			  NTuple->kstTrkpdzVtx->push_back(TrackpTT.track().dz(bestVtxReFit.position()));

			  NTuple->kstTrkpMuMatch->push_back(MuPCat);

			  // I do NOT include the number of missing outer hits because the hadron might interact
			  NTuple->kstTrkpNPixHits->push_back(Trackp->hitPattern().numberOfValidPixelHits());
			  NTuple->kstTrkpNPixLayers->push_back(Trackp->hitPattern().pixelLayersWithMeasurement());  
			  NTuple->kstTrkpNTrkHits->push_back(Trackp->hitPattern().numberOfValidTrackerHits());
			  NTuple->kstTrkpNTrkLayers->push_back(Trackp->hitPattern().trackerLayersWithMeasurement());




			  // #####################
			  // # Clear all vectors #
			  // #####################
			  vertexTracks.clear();
			  bParticles.clear();
			  bBarParticles.clear();
			  kstParticles.clear();
			  kstBarParticles.clear();
			} // End for Track+
		      if (skip == true) continue;
		    } // End for Track-
		  muonParticles.clear(); 
		} // End for mu+
	      if (skip == true) continue;
	    } // End for mu-
	} // End if bestVtx is true
      else if (printMsg == true) std::cout << __LINE__ << " : continue --> invalid Pri.Vtx" << std::endl;
      if (B0MassConstraint != NULL) delete B0MassConstraint;


      if (NTuple->nB > 0)
	{
	  NTuple->runN   = iEvent.id().run();
	  NTuple->eventN = iEvent.id().event();
	  for (std::vector<reco::Vertex>::const_iterator iVertex = recVtx->begin(); iVertex != recVtx->end(); iVertex++)
	    {
	      if(iVertex->ndof() < PRIVTXNDOF)                 continue;
	      if(fabs(iVertex->z()) > PRIVTXMAXZ)              continue;
	      if(fabs(iVertex->position().rho()) > PRIVTXMAXR) continue;
	      NTuple->recoVtxN++;
	    }


	  // #####################################
	  // # Save: Primary Vertex and BeamSpot #
	  // #####################################
	  NTuple->priVtxCL = TMath::Prob(bestVtx.chi2(), static_cast<int>(rint(bestVtx.ndof())));
	  NTuple->priVtxX  = bestVtx.x();
	  NTuple->priVtxY  = bestVtx.y();
	  NTuple->priVtxZ  = bestVtx.z();
	  NTuple->bsX      = beamSpot.position().x();
	  NTuple->bsY      = beamSpot.position().y();
	}
      else if (printMsg == true) std::cout << __LINE__ << " : @@@ No B0 --> K*0 (K pi) mu+ mu- candidate were found in the event @@@" << std::endl;
    } // End if doGenReco_ == 1 || doGenReco_ == 2




  // ###########################################
  // # Get truth information from genParticles #
  // ###########################################


  // Particle IDs:
  //   11 = e-
  //   12 = nu_e
  //   13 = mu-
  //   14 = nu_mu
  //   16 = nu_tau
  //   22 = gamma
  //  130 = KL
  //  211 = pi+
  //  310 = Ks
  //  311 = K0
  //  313 = K*0
  //  321 = K+
  //  443 = J/psi; 100443 psi(2S)
  //  511 = B0
  //  521 = B+
  //  531 = Bs

  // 2212 = p
  // 5122 = Lambda_b


  if (doGenReco_ == 2)
    {
      // #################################
      // # Save pileup information in MC #
      // #################################
      edm::Handle< std::vector<PileupSummaryInfo> > PupInfo;
      iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
      for (std::vector<PileupSummaryInfo>::const_iterator PVI = PupInfo->begin(); PVI != PupInfo->end(); PVI++)
	{
	  NTuple->bunchXingMC->push_back(PVI->getBunchCrossing());
	  NTuple->numInteractionsMC->push_back(PVI->getPU_NumInteractions());
	  NTuple->trueNumInteractionsMC->push_back(PVI->getTrueNumInteractions());
	}
    }


  if ((doGenReco_ == 2) || (doGenReco_ == 3))
    {
      std::vector<std::vector<unsigned int>* > posDau;
      std::vector<std::vector<unsigned int>* > negDau;

      for (unsigned int itGen = 0; itGen < genParticles->size(); itGen++)
	{
	  const reco::Candidate* FirstPart = &(*genParticles)[itGen];


	  // ##########################
	  // # Check for:             #
	  // # B0 / B0bar             #
	  // # Bs / Bsbar             #
	  // # Lambda_b / Lambda_bbar #
	  // # B+ / B-                #
	  // ##########################
	  if ((abs(FirstPart->pdgId()) == 511) || (abs(FirstPart->pdgId()) == 531) || (abs(FirstPart->pdgId()) == 5122) || (abs(FirstPart->pdgId()) == 521))
	    {	     
	      if (printMsg == true) std::cout << "\n" << __LINE__ << " : @@@ Found B0/B0bar OR Bs/Bsbar OR Lambda_b/Lambda_bbar OR B+/B- in MC @@@" << std::endl;


	      // ########################################################
	      // # Search for oscillations                              #
	      // # If there are, then move to the non oscillating stage #
	      // ########################################################
	      FirstPart = skipOscillations(FirstPart);


	      int imum      = -1;
	      int imup      = -1;
	      int ikst      = -1;
	      int ikst_trkm = -1;
	      int ikst_trkp = -1;
	      int ipsi      = -1;
	      int ipsi_mum  = -1;
	      int ipsi_mup  = -1;

	      bool PhotonB0 = false;
	      bool PhotonKst = false;
	      bool PhotonPsi = false;

	      bool isSignal = true;
	      bool isPsi2SnotJPsi = false;


	      if (abs(FirstPart->pdgId()) == 511)
		{	
		  for (unsigned int i = 0; i < FirstPart->numberOfDaughters(); i++)
		    {
		      if (isSignal == false) break;
		      
		      const reco::Candidate* genDau = FirstPart->daughter(i);

		      if (genDau->pdgId() == 22)
			{
			  PhotonB0 = true;
			  if (printMsg == true) std::cout << __LINE__ << " : found B0/B0bar photon" << std::endl;
			  continue;
			}
		      if (genDau->pdgId() == 13)
			{
			  imum = i;
			  if (printMsg == true) std::cout << __LINE__ << " : found mu-" << std::endl;
			  continue;
			}
		      if (genDau->pdgId() == -13)
			{
			  imup = i;
			  if (printMsg == true) std::cout << __LINE__ << " : found mu+" << std::endl;
			  continue;
			}
		      if (abs(genDau->pdgId()) == 313)
			{
			  ikst = i;

			  if (printMsg == true) std::cout << __LINE__ << " : found K*0/K*0bar" << std::endl;

			  for (unsigned int j = 0; j < genDau->numberOfDaughters(); j++)
			    {
			      const reco::Candidate* genDau2 = genDau->daughter(j);
			  
			      if (genDau2->pdgId() == 22)
				{
				  PhotonKst = true;
				  if (printMsg == true) std::cout << __LINE__ << " : found K*0/K*0bar photon" << std::endl;
				  continue;
				}
			      if ((genDau2->pdgId() == -211) || (genDau2->pdgId() == -321))
				{
				  ikst_trkm = j;
				  if (printMsg == true) std::cout << __LINE__ << " : found K*0/K*0bar track-" << std::endl;
				  continue;
				}
			      if ((genDau2->pdgId() == 211) || (genDau2->pdgId() == 321))
				{
				  ikst_trkp = j;
				  if (printMsg == true) std::cout << __LINE__ << " : found K*0/K*0bar track+" << std::endl;
				  continue;
				}
			      isSignal = false;
			      break;
			    }
			  continue;
			}
		      if (abs(genDau->pdgId() == 443) || abs(genDau->pdgId() == 100443))
			{
			  if (abs(genDau->pdgId() == 443)) isPsi2SnotJPsi = false;
			  else isPsi2SnotJPsi = true;

			  ipsi = i;
		      
			  if (printMsg == true) std::cout << __LINE__ << " : found J/psi or psi(2S)" << std::endl;

			  for (unsigned int j = 0; j < genDau->numberOfDaughters(); j++)
			    {
			      const reco::Candidate* genDau2 = genDau->daughter(j);
			  
			      if (genDau2->pdgId() == 22)
				{
				  PhotonPsi = true;
				  if (printMsg == true) std::cout << __LINE__ << " : found J/psi or psi(2S) photon" << std::endl;
				  continue;
				}
			      if (genDau2->pdgId() == 13)
				{
				  ipsi_mum = j;
				  if (printMsg == true) std::cout << __LINE__ << " : found J/psi or psi(2S) mu-" << std::endl;
				  continue;
				}
			      if (genDau2->pdgId() == -13)
				{
				  ipsi_mup = j;
				  if (printMsg == true) std::cout << __LINE__ << " : found J/psi or psi(2S) mu+" << std::endl;
				  continue;
				}
			      isSignal = false;
			      break;
			    }
			  continue;
			}
		      isSignal = false;
		      break;
		    } // End for B0 / B0bar daughters
		} // End if B0 / B0bar


	      const reco::Candidate* genPsi = NULL;
	      const reco::Candidate* genMum = NULL;
	      const reco::Candidate* genMup = NULL;
	      const reco::Candidate* genKst = NULL;
	      const reco::Candidate* genKst_trkm = NULL;
	      const reco::Candidate* genKst_trkp = NULL;


	      bool foundSomething = false;
	      if ((abs(FirstPart->pdgId()) == 511) && (isSignal == true) &&
		  (((imum != -1) && (imup != -1) && (ikst != -1) && (ikst_trkm != -1) && (ikst_trkp != -1)) ||
		   ((NTuple->genSignal != 1) && (NTuple->genSignal != 2) &&
		    (ikst != -1) && (ikst_trkm != -1) && (ikst_trkp != -1) && (ipsi != -1) && (ipsi_mum != -1) && (ipsi_mup != -1))))
		{
		  // ################################################################
		  // # Found Signal B0/B0bar --> K*0/K*0bar (K pi) mu+ mu-          #
		  // # Found Signal B0/B0bar --> K*0/K*0bar (K pi) J/psi (mu+mu-)   #
		  // # Found Signal B0/B0bar --> K*0/K*0bar (K pi) psi(2S) (mu+mu-) #
		  // ################################################################

		  NTuple->ClearMonteCarlo();

		  if ((imum != -1) && (imup != -1) && (ikst != -1) && (ikst_trkm != -1) && (ikst_trkp != -1))
		    {
		      if (printMsg == true) std::cout << __LINE__ << " : @@@ Found B0/B0bar --> K*0/K*0bar (K pi) mu+ mu- @@@" << std::endl;
		      
		      if      (FirstPart->daughter(ikst)->pdgId() == 313)  NTuple->genSignal = 1;
		      else if (FirstPart->daughter(ikst)->pdgId() == -313) NTuple->genSignal = 2;
		      
		      genMum = FirstPart->daughter(imum);
		      genMup = FirstPart->daughter(imup);
		    }
		  else
		    {
		      if (isPsi2SnotJPsi == false)
			{
			  if (printMsg == true) std::cout << __LINE__ << " : @@@ Found B0/B0bar --> K*0/K*0bar (K pi) J/psi (mu+mu-) @@@" << std::endl;
		      
			  if      (FirstPart->daughter(ikst)->pdgId() == 313)  NTuple->genSignal = 3;
			  else if (FirstPart->daughter(ikst)->pdgId() == -313) NTuple->genSignal = 4;
			}
		      else
			{
			  if (printMsg == true) std::cout << __LINE__ << " : @@@ Found B0/B0bar --> K*0/K*0bar (K pi) psi(2S) (mu+mu-) @@@" << std::endl;
		      
			  if      (FirstPart->daughter(ikst)->pdgId() == 313)  NTuple->genSignal = 5;
			  else if (FirstPart->daughter(ikst)->pdgId() == -313) NTuple->genSignal = 6;
			}
		      
		      genPsi = FirstPart->daughter(ipsi);
		      genMum = genPsi->daughter(ipsi_mum);
		      genMup = genPsi->daughter(ipsi_mup);

		      NTuple->genSignPsiHasFSR = PhotonPsi;
		    }

		  genKst = FirstPart->daughter(ikst);
		  genKst_trkm = FirstPart->daughter(ikst)->daughter(ikst_trkm);
		  genKst_trkp = FirstPart->daughter(ikst)->daughter(ikst_trkp);

 
		  // ##################
		  // # Save info. FSR #
		  // ##################
		  NTuple->genSignHasFSR = PhotonB0;
		  NTuple->genSignKstHasFSR = PhotonKst;


		  foundSomething = true;
		} // End if B0/B0bar signal
	      else
		{
		  if (printMsg == true) std::cout << __LINE__ << " : @@@ Start particle decay-tree scan for opposite charged stable particles for background search @@@" << std::endl;
		  searchForStableChargedDaughters (FirstPart, itGen, &posDau, &negDau);


		  if ((NTuple->genSignal == 0) && (posDau.size() != 0) && (negDau.size() != 0))
		    {
		      // ###############################
		      // # Search for Background from: #
		      // # B0 / B0bar                  #
		      // # Bs / Bsbar                  #
		      // # Lambda_b / Lambda_bbar      #
		      // # B+ / B-                     #
		      // ###############################

		      for (unsigned int i = 0; i < negDau.size(); i++)
			{
			  genMum = findDaughter (genParticles, &negDau, i);
			  if (genMum->pdgId() == 13)
			    {
			      for (unsigned int j = 0; j < posDau.size(); j++)
				{
				  genMup = findDaughter (genParticles, &posDau, j);
				  if (genMup->pdgId() == -13)
				    {
				      if (printMsg == true) std::cout << __LINE__ << " : found dimuons for possible background" << std::endl;
				      foundSomething = true;
				      break;
				    }
				}
			      if (foundSomething == true) break;
			    }
			}

		      if ((foundSomething == true) && (posDau.size() == 2) && (negDau.size() == 1))
			{
			  foundSomething = false;
			  for (unsigned int i = 0; i < posDau.size(); i++)
			    {
			      genKst_trkp = findDaughter (genParticles, &posDau, i);
			      if (genKst_trkp != genMup)
				{
				  NTuple->ClearMonteCarlo();

				  NTuple->genMuMuBG = FirstPart->pdgId();
				  NTuple->genMuMuBGnTrksp = 1;
				  foundSomething = true;
				  if (printMsg == true) std::cout << __LINE__ << " : get background positive track: " << genKst_trkp->pdgId() << "\tfrom mother: " << genKst_trkp->mother()->pdgId() << std::endl;
				}
			    }
			}
		      else if ((foundSomething == true) && (posDau.size() == 1) && (negDau.size() == 2))
			{
			  foundSomething = false;
			  for (unsigned int i = 0; i < negDau.size(); i++)
			    {
			      genKst_trkm = findDaughter (genParticles, &negDau, i);
			      if (genKst_trkm != genMum)
				{
				  NTuple->ClearMonteCarlo();

				  NTuple->genMuMuBG = FirstPart->pdgId();			  
				  NTuple->genMuMuBGnTrksm = 1;
				  foundSomething = true;
				  if (printMsg == true) std::cout << __LINE__ << " : get background negative track: " << genKst_trkm->pdgId() << "\tfrom mother: " << genKst_trkm->mother()->pdgId() << std::endl;
				}
			    }
			}
		      else if ((foundSomething == true) && (posDau.size() == 2) && (negDau.size() == 2))
			{
			  foundSomething = false;
			  for (unsigned int i = 0; i < negDau.size(); i++)
			    {
			      genKst_trkm = findDaughter (genParticles, &negDau, i);
			      if (genKst_trkm != genMum)
				{
				  for (unsigned int j = 0; j < posDau.size(); j++)
				    {
				      genKst_trkp = findDaughter (genParticles, &posDau, j);
				      if (genKst_trkp != genMup)
					{
					  NTuple->ClearMonteCarlo();

					  NTuple->genMuMuBG = FirstPart->pdgId();
					  NTuple->genMuMuBGnTrksp = 1;
					  NTuple->genMuMuBGnTrksm = 1;
					  foundSomething = true;
					  if (printMsg == true)
					    {
					      std::cout << __LINE__ << " : get background negative track: " << genKst_trkm->pdgId() << "\tfrom mother: " << genKst_trkm->mother()->pdgId();
					      std::cout << "\tand positive track: " << genKst_trkp->pdgId() << "\tfrom mother: " << genKst_trkp->mother()->pdgId() << std::endl;
					    }
					}
				    }
				}
			    }
			}
		      else if ((foundSomething == true) && (posDau.size() >= 2) && (negDau.size() >= 2))
			{
			  foundSomething = false;

			  double bestMass = 0.;
			  unsigned int negDauIndx = 0;
			  unsigned int posDauIndx = 0;

			  for (unsigned int i = 0; i < negDau.size(); i++)
			    {
			      genKst_trkm = findDaughter (genParticles, &negDau, i);
			      if (genKst_trkm != genMum)
				{
				  for (unsigned int j = 0; j < posDau.size(); j++)
				    {
				      genKst_trkp = findDaughter (genParticles, &posDau, j);
				      if (genKst_trkp != genMup)
					{
					  double invMass = Utility->computeInvMass (genKst_trkm->px(),genKst_trkm->py(),genKst_trkm->pz(),genKst_trkm->mass(),
										    genKst_trkp->px(),genKst_trkp->py(),genKst_trkp->pz(),genKst_trkp->mass());
				      
					  if (fabs(invMass - Utility->kstMass) < fabs(bestMass - Utility->kstMass))
					    {
					      bestMass = invMass;
					      negDauIndx = i;
					      posDauIndx = j;
					    }
					}
				    }
				}
			    }

			  if (bestMass > 0.)
			    {
			      genKst_trkm = findDaughter (genParticles, &negDau, negDauIndx);
			      genKst_trkp = findDaughter (genParticles, &posDau, posDauIndx);

			      NTuple->ClearMonteCarlo();

			      NTuple->genMuMuBG = FirstPart->pdgId();
			      NTuple->genMuMuBGnTrksp = 1;
			      NTuple->genMuMuBGnTrksm = 1;
			      foundSomething = true;
			      if (printMsg == true)
				{
				  std::cout << __LINE__ << " : get background negative track: " << genKst_trkm->pdgId() << "\tfrom mother: " << genKst_trkm->mother()->pdgId();
				  std::cout << "\tand positive track: " << genKst_trkp->pdgId() << "\tfrom mother: " << genKst_trkp->mother()->pdgId() << std::endl;
				}
			    }
			}
		      else
			{
			  foundSomething = false;
			  if (printMsg == true) std::cout << __LINE__ << " : @@@ No background found @@@" << std::endl;
			}
		    }
		  else if (printMsg == true) std::cout << __LINE__ << " : @@@ No possible background found @@@" << std::endl;
		} // End else background


	      // ########################
	      // # Save gen information #
	      // ########################
	      if (foundSomething == true)
		{
		  if (printMsg == true) std::cout << __LINE__ << " : @@@ Saving signal OR background compatible with signal @@@" << std::endl;


		  // #############################
		  // # Search for primary vertex #
		  // #############################
		  const reco::Candidate* PVtx = FirstPart;
		  while (PVtx->numberOfMothers() > 0) PVtx = PVtx->mother(0);


		  // ############################
		  // # Generated Primary Vertex #
		  // ############################
		  NTuple->genPriVtxX = PVtx->vx();
		  NTuple->genPriVtxY = PVtx->vy();
		  NTuple->genPriVtxZ = PVtx->vz();

		  
		  // #####################
		  // # Generated B0 Mass #
		  // #####################
		  NTuple->genB0Mass = FirstPart->mass();
		  NTuple->genB0Px = FirstPart->px();
		  NTuple->genB0Py = FirstPart->py();
		  NTuple->genB0Pz = FirstPart->pz();


		  // ####################
		  // # Generated B0 Vtx #
		  // ####################
		  NTuple->genB0VtxX = FirstPart->vx();
		  NTuple->genB0VtxY = FirstPart->vy();
		  NTuple->genB0VtxZ = FirstPart->vz();
		  
	
		  if (NTuple->genSignal != 0)
		    {
		      // ######################
		      // # Generated K*0 Mass #
		      // ######################
		      NTuple->genKstMass = genKst->mass();
		      NTuple->genKstPx = genKst->px();
		      NTuple->genKstPy = genKst->py();
		      NTuple->genKstPz = genKst->pz();
		      

		      // #####################
		      // # Generated K*0 Vtx #
		      // #####################
		      NTuple->genKstVtxX = genKst->vx();
		      NTuple->genKstVtxY = genKst->vy();
		      NTuple->genKstVtxZ = genKst->vz();
		    }
		  else if ((NTuple->genMuMuBGnTrksm != 0) && (NTuple->genMuMuBGnTrksp != 0) &&
			   (genKst_trkm->vx() == genKst_trkp->vx()) &&
			   (genKst_trkm->vy() == genKst_trkp->vy()) &&
			   (genKst_trkm->vz() == genKst_trkp->vz()))
		    {
		      double invMass = Utility->computeInvMass (genKst_trkm->momentum().x(),genKst_trkm->momentum().y(),genKst_trkm->momentum().z(),genKst_trkm->mass(),
								genKst_trkp->momentum().x(),genKst_trkp->momentum().y(),genKst_trkp->momentum().z(),genKst_trkp->mass());
		      

		      // ######################
		      // # Generated K*0 Mass #
		      // ######################
		      NTuple->genKstMass = invMass;
		      NTuple->genKstPx = genKst_trkm->momentum().x() + genKst_trkp->momentum().x();
		      NTuple->genKstPy = genKst_trkm->momentum().y() + genKst_trkp->momentum().y();
		      NTuple->genKstPz = genKst_trkm->momentum().z() + genKst_trkp->momentum().z();


		      // #####################
		      // # Generated K*0 Vtx #
		      // #####################
		      NTuple->genKstVtxX = genKst_trkm->vx();
		      NTuple->genKstVtxY = genKst_trkm->vy();
		      NTuple->genKstVtxZ = genKst_trkm->vz();
		    }
		  
		  
		  // ###########################################
		  // # Generated J/psi or psi(2S) Mass and Vtx #
		  // ###########################################
		  if ((NTuple->genSignal == 3) || (NTuple->genSignal == 4) || (NTuple->genSignal == 5) || (NTuple->genSignal == 6))
		    {
		      NTuple->genPsiMass = genPsi->mass();
		      NTuple->genPsiVtxX = genPsi->vx();
		      NTuple->genPsiVtxY = genPsi->vy();
		      NTuple->genPsiVtxZ = genPsi->vz();
		    }
		  else if ((NTuple->genSignal == 0) &&
			   (genMum->vx() == genMup->vx()) &&
			   (genMum->vy() == genMup->vy()) &&
			   (genMum->vz() == genMup->vz()))
		    {
		      double invMass = Utility->computeInvMass (genMum->momentum().x(),genMum->momentum().y(),genMum->momentum().z(),genMum->mass(),
								genMup->momentum().x(),genMup->momentum().y(),genMup->momentum().z(),genMup->mass());
		      

		      NTuple->genPsiMass = invMass;
		      NTuple->genPsiVtxX = genMum->vx();
		      NTuple->genPsiVtxY = genMum->vy();
		      NTuple->genPsiVtxZ = genMum->vz();
		    }


		  // #################
		  // # Generated mu- #
		  // #################
		  NTuple->genMumMother = genMum->mother()->pdgId();
		  NTuple->genMumPx = genMum->px();
		  NTuple->genMumPy = genMum->py();
		  NTuple->genMumPz = genMum->pz();
		  

		  // #################
		  // # Generated mu+ #
		  // #################
		  NTuple->genMupMother = genMup->mother()->pdgId();
		  NTuple->genMupPx = genMup->px();
		  NTuple->genMupPy = genMup->py();
		  NTuple->genMupPz = genMup->pz();


		  // ########################
		  // # Generated K*0 track- #
		  // ########################
		  if ((NTuple->genSignal != 0) || (NTuple->genMuMuBGnTrksm != 0))
		    {
		      NTuple->genKstTrkmMother = genKst_trkm->mother()->pdgId();
		      NTuple->genKstTrkmID = genKst_trkm->pdgId();
		      NTuple->genKstTrkmPx = genKst_trkm->px();
		      NTuple->genKstTrkmPy = genKst_trkm->py();
		      NTuple->genKstTrkmPz = genKst_trkm->pz();
		    }


		  // ########################
		  // # Generated K*0 track+ #
		  // ########################
		  if ((NTuple->genSignal != 0) || (NTuple->genMuMuBGnTrksp != 0))
		    {
		      NTuple->genKstTrkpMother = genKst_trkp->mother()->pdgId();
		      NTuple->genKstTrkpID = genKst_trkp->pdgId();
		      NTuple->genKstTrkpPx = genKst_trkp->px();
		      NTuple->genKstTrkpPy = genKst_trkp->py();
		      NTuple->genKstTrkpPz = genKst_trkp->pz();
		    }
		} // End if foundSomething
	    } // End if B0/B0bar OR Bs/Bsbar OR Lambda_b/Lambda_bbar OR B+/B-


	  // ##############################################
	  // # Check to see if J/psi or psi(2S) is prompt #
	  // ##############################################
	  bool isPrompt = false;
	  const reco::Candidate& PsiCand = (*genParticles)[itGen];
	  
	  if ((abs(PsiCand.pdgId()) == 443) || (abs(PsiCand.pdgId()) == 100443))
	    {
	      isPrompt = true;

	      for (unsigned int i = 0; i < PsiCand.numberOfMothers(); i++)
		{
		  const reco::Candidate* psiMom = PsiCand.mother(i);
		  
		  if (((abs(psiMom->pdgId()) < 600) && (abs(psiMom->pdgId()) > 500)) || ((abs(psiMom->pdgId()) < 6000) && (abs(psiMom->pdgId()) > 5000)))
		    {
		      isPrompt = false;
		      continue;
		    }
		  else
		    {
		      for (unsigned int i = 0; i < psiMom->numberOfMothers(); i++)
			{
			  const reco::Candidate* psiMom2 = psiMom->mother(i);

			  if (((abs(psiMom2->pdgId()) < 600) && (abs(psiMom2->pdgId()) > 500)) || ((abs(psiMom2->pdgId()) < 6000) && (abs(psiMom2->pdgId()) > 5000)))
			    {
			      isPrompt = false;
			      continue;
			    }
			  else
			    {
			      for (unsigned int i = 0; i < psiMom2->numberOfMothers(); i++)
				{
				  const reco::Candidate* psiMom3 = psiMom2->mother(i);

				  if (((abs(psiMom3->pdgId()) < 600) && (abs(psiMom3->pdgId()) > 500)) || ((abs(psiMom3->pdgId()) < 6000) && (abs(psiMom3->pdgId()) > 5000)))
				    {
				      isPrompt = false;
				      continue;
				    }
				}
			    }
			}
		    }
		}
	      
	      if (isPrompt == true)
		{
		  NTuple->genPsiPrompt = 1;
		  if (printMsg == true) std::cout << __LINE__ << " : found prompt J/psi or psi(2S)" << std::endl;
		}
	      continue;
	    }

	  for (unsigned int i = 0; i < posDau.size(); i++)
	    {
	      posDau[i]->clear();
	      delete posDau[i];
	    }
	  posDau.clear();
	  for (unsigned int i = 0; i < negDau.size(); i++)
	    {
	      negDau[i]->clear();
	      delete negDau[i];
	    }
	  negDau.clear();
	} // End for genParticles
    } // End if doGenReco_ == 2 || doGenReco_ == 3




  // ####################################
  // # Perform matching with candidates #
  // ####################################
  if ((NTuple->genSignal != 0) || (NTuple->genMuMuBG != 0))
    for (unsigned int i = 0; i < NTuple->nB; i++)
      {
	// ###########################
	// # Check matching with mu- #
	// ###########################
	deltaEtaPhi = Utility->computeEtaPhiDistance (NTuple->genMumPx,NTuple->genMumPy,NTuple->genMumPz, NTuple->mumPx->at(i),NTuple->mumPy->at(i),NTuple->mumPz->at(i));
	NTuple->mumDeltaRwithMC->push_back(deltaEtaPhi);
	if (deltaEtaPhi < RCUTMU)
	  {
	    NTuple->truthMatchMum->push_back(true);
	    if (printMsg == true) std::cout << __LINE__ << " : found matched mu-" << std::endl;
	  }
	else NTuple->truthMatchMum->push_back(false);
		      

	// ###########################
	// # Check matching with mu+ #
	// ###########################
	deltaEtaPhi = Utility->computeEtaPhiDistance (NTuple->genMupPx,NTuple->genMupPy,NTuple->genMupPz, NTuple->mupPx->at(i),NTuple->mupPy->at(i),NTuple->mupPz->at(i));
	NTuple->mupDeltaRwithMC->push_back(deltaEtaPhi);
	if (deltaEtaPhi < RCUTMU)
	  {
	    NTuple->truthMatchMup->push_back(true);
	    if (printMsg == true) std::cout << __LINE__ << " : found matched mu+" << std::endl;
	  }
	else NTuple->truthMatchMup->push_back(false);
		      

	// ##################################
	// # Check matching with K*0 track- #
	// ##################################
	if ((NTuple->genSignal != 0) || (NTuple->genMuMuBGnTrksm != 0))
	  {
	    deltaEtaPhi = Utility->computeEtaPhiDistance (NTuple->genKstTrkmPx,NTuple->genKstTrkmPy,NTuple->genKstTrkmPz, NTuple->kstTrkmPx->at(i),NTuple->kstTrkmPy->at(i),NTuple->kstTrkmPz->at(i));
	    NTuple->kstTrkmDeltaRwithMC->push_back(deltaEtaPhi);
	    if (deltaEtaPhi < RCUTTRK)
	      {
		NTuple->truthMatchTrkm->push_back(true);
		if (printMsg == true) std::cout << __LINE__ << " : found matched track-" << std::endl;
	      }
	    else NTuple->truthMatchTrkm->push_back(false);
	  }
	else
	  {
	    NTuple->truthMatchTrkm->push_back(false);
	    NTuple->kstTrkmDeltaRwithMC->push_back(-1.0);
	  }


	// ##################################
	// # Check matching with K*0 track+ #
	// ##################################
	if ((NTuple->genSignal != 0) || (NTuple->genMuMuBGnTrksp != 0))
	  {
	    deltaEtaPhi = Utility->computeEtaPhiDistance (NTuple->genKstTrkpPx,NTuple->genKstTrkpPy,NTuple->genKstTrkpPz, NTuple->kstTrkpPx->at(i),NTuple->kstTrkpPy->at(i),NTuple->kstTrkpPz->at(i));
	    NTuple->kstTrkpDeltaRwithMC->push_back(deltaEtaPhi);
	    if (deltaEtaPhi < RCUTTRK) 
	      {
		NTuple->truthMatchTrkp->push_back(true);
		if (printMsg == true) std::cout << __LINE__ << " : found matched track+" << std::endl;
	      }
	    else NTuple->truthMatchTrkp->push_back(false);
	  }
	else
	  {
	    NTuple->truthMatchTrkp->push_back(false);
	    NTuple->kstTrkpDeltaRwithMC->push_back(-1.0);
	  }


	// ####################################################
	// # Check matching with B0 --> track+ track- mu+ mu- #
	// ####################################################
	if ((NTuple->truthMatchTrkm->back() == true) && (NTuple->truthMatchTrkp->back() == true) &&
	    (NTuple->truthMatchMum->back() == true) && (NTuple->truthMatchMup->back() == true))
	  {
	    NTuple->truthMatchSignal->push_back(true);
	    if (printMsg == true) std::cout << __LINE__ << " : @@@ Found matched B0 --> track+ track- mu+ mu- @@@" << std::endl;
	  }
	else NTuple->truthMatchSignal->push_back(false);
      }
  else
    for (unsigned int i = 0; i < NTuple->nB; i++)
      {
	NTuple->mumDeltaRwithMC->push_back(-1.0);
	NTuple->mupDeltaRwithMC->push_back(-1.0);
	NTuple->kstTrkmDeltaRwithMC->push_back(-1.0);
	NTuple->kstTrkpDeltaRwithMC->push_back(-1.0);
	
	NTuple->truthMatchMum->push_back(false);
	NTuple->truthMatchMup->push_back(false);
	NTuple->truthMatchTrkm->push_back(false);
	NTuple->truthMatchTrkp->push_back(false);
	NTuple->truthMatchSignal->push_back(false);
      }




  // ###################
  // # Fill the ntuple #
  // ###################
  if (printMsg == true) std::cout << __LINE__ << " : @@@ Filling the tree @@@" << std::endl;
  theTree->Fill();
  NTuple->ClearNTuple();
}


std::string B0KstMuMu::getMuCat (reco::Muon const& muon)
{
  std::stringstream muCat;
  muCat.str("");

  if (muon.isGlobalMuon() == true)
    {
      muCat << " GlobalMuon";
      if (muon::isGoodMuon(muon, muon::GlobalMuonPromptTight) == true) muCat << " GlobalMuonPromptTight";
    }
  if (muon.isTrackerMuon() == true)
    {
      muCat << " TrackerMuon";
      if (muon::isGoodMuon(muon, muon::TrackerMuonArbitrated) == true) muCat << " TrackerMuonArbitrated";

      if (muon::isGoodMuon(muon, muon::TMLastStationTight)     == true) muCat << " TMLastStationTight";
      if (muon::isGoodMuon(muon, muon::TMLastStationLoose)     == true) muCat << " TMLastStationLoose";
      if (muon::isGoodMuon(muon, muon::TM2DCompatibilityTight) == true) muCat << " TM2DCompatibilityTight";
      if (muon::isGoodMuon(muon, muon::TM2DCompatibilityLoose) == true) muCat << " TM2DCompatibilityLoose";
      if (muon::isGoodMuon(muon, muon::TMOneStationTight)      == true) muCat << " TMOneStationTight";
      if (muon::isGoodMuon(muon, muon::TMOneStationLoose)      == true) muCat << " TMOneStationLoose";
      if (muon::isGoodMuon(muon, muon::TMLastStationAngTight)  == true) muCat << " TMLastStationAngTight";
      if (muon::isGoodMuon(muon, muon::TMLastStationAngLoose)  == true) muCat << " TMLastStationAngLoose";
      if (muon::isGoodMuon(muon, muon::TMOneStationAngTight)   == true) muCat << " TMOneStationAngTight";
      if (muon::isGoodMuon(muon, muon::TMOneStationAngLoose)   == true) muCat << " TMOneStationAngLoose";
    }
  if (muon.isStandAloneMuon() == true) muCat << " StandAloneMuon";
  if (muon.isCaloMuon()       == true) muCat << " CaloMuon";
  if ((muon.isGlobalMuon() == false) && (muon.isTrackerMuon() == false) && (muon.isStandAloneMuon() == false) && (muon.isCaloMuon() == false)) muCat << " NotInTable";

  return muCat.str();
}


const reco::Candidate* B0KstMuMu::skipOscillations (const reco::Candidate* Mother)
{
  if (abs(Mother->pdgId()) != 521)
    for (unsigned int i = 0; i < Mother->numberOfDaughters(); i++)
      if ((abs(Mother->daughter(i)->pdgId()) == 511) || (abs(Mother->daughter(i)->pdgId()) == 531) || (abs(Mother->daughter(i)->pdgId()) == 5122))
	{
	  if (printMsg == true) std::cout << __LINE__ << " : @@@ Found oscillating B0/B0bar OR Bs/Bsbar OR Lambda_b/Lambda_bbar @@@" << std::endl;
	  Mother = Mother->daughter(i);
	}

  return Mother;
}


const reco::Candidate* B0KstMuMu::findDaughter (edm::Handle<reco::GenParticleCollection> genParticles,
						std::vector<std::vector<unsigned int>* >* Dau,
						unsigned int it)
{
  const reco::Candidate* gen;
  if (Dau != NULL)
    {
      gen = skipOscillations(&(*genParticles)[(*Dau)[it]->operator[](0)]);
      for (unsigned int i = 1; i < (*Dau)[it]->size(); i++) gen = gen->daughter((*Dau)[it]->operator[](i));
      return gen;
    }
  else return NULL;
}


void B0KstMuMu::searchForStableChargedDaughters (const reco::Candidate* Mother,
						 unsigned int motherIndex,
						 std::vector<std::vector<unsigned int>* >* posDau,
						 std::vector<std::vector<unsigned int>* >* negDau)
{
  const reco::Candidate* genDau;
  unsigned int sizeNegVec = negDau->size();
  unsigned int sizePosVec = posDau->size();

  for (unsigned int i = 0; i < Mother->numberOfDaughters(); i++)
    {
      genDau = Mother->daughter(i);

      if (printMsg == true) std::cout << __LINE__ << " : start exploring daughter: " << genDau->pdgId() << "\tindex: " << i << "\tof mother: " << Mother->pdgId() << std::endl;

      if ((genDau->pdgId() == 11) ||
	  (genDau->pdgId() == 13) ||
	  (genDau->pdgId() == -211) ||
	  (genDau->pdgId() == -321) ||
	  (genDau->pdgId() == -2212))
	{
	  std::vector<unsigned int>* vecDau = new std::vector<unsigned int>;
	  vecDau->push_back(i);
	  negDau->push_back(vecDau);
	  if (printMsg == true) std::cout << __LINE__ << " : found possible background negative track: " << genDau->pdgId() << "\tfrom mother: " << genDau->mother()->pdgId() << std::endl;
	}
      else if ((genDau->pdgId() == -11) ||
	       (genDau->pdgId() == -13) ||
	       (genDau->pdgId() == 211) ||
	       (genDau->pdgId() == 321) ||
	       (genDau->pdgId() == 2212))
	{
	  std::vector<unsigned int>* vecDau = new std::vector<unsigned int>;
	  vecDau->push_back(i);
	  posDau->push_back(vecDau);
	  if (printMsg == true) std::cout << __LINE__ << " : found possible background positive track: " << genDau->pdgId() << "\tfrom mother: " << genDau->mother()->pdgId() << std::endl;
	}
      else if ((abs(genDau->pdgId()) != 311) &&
	       (abs(genDau->pdgId()) != 310) &&
	       (abs(genDau->pdgId()) != 130) &&
	       (abs(genDau->pdgId()) != 22) &&
	       (abs(genDau->pdgId()) != 12) &&
	       (abs(genDau->pdgId()) != 14) &&
	       (abs(genDau->pdgId()) != 16))
	searchForStableChargedDaughters (genDau, i, posDau, negDau);
      else if (printMsg == true) std::cout << __LINE__ << " : found long living neutral particle: " << genDau->pdgId() << std::endl;
    }

  for (unsigned int it = negDau->size(); it > sizeNegVec; it--) negDau->operator[](it-1)->insert(negDau->operator[](it-1)->begin(), motherIndex);
  for (unsigned int it = posDau->size(); it > sizePosVec; it--) posDau->operator[](it-1)->insert(posDau->operator[](it-1)->begin(), motherIndex);
}


void B0KstMuMu::beginJob ()
{
  edm::Service<TFileService> fs;
  theTree = fs->make<TTree>("B0KstMuMuNTuple","B0KstMuMuNTuple");
  NTuple->MakeTreeBranches(theTree);


  // #####################
  // # Loading HLT-paths #
  // #####################
  Utility->ReadHLTpaths(parameterFile_.c_str(),&TrigTable);


  // ############################
  // # Loading HLT-trigger cuts #
  // ############################
  Utility->ReadPreselectionCut(parameterFile_.c_str());
  std::cout << "\n@@@ Pre-selection cuts @@@" << std::endl;
  CLMUMUVTX      = Utility->GetPreCut("MuMuVtxCL");      std::cout << __LINE__ << " : CLMUMUVTX      = " << CLMUMUVTX << std::endl;
  LSMUMUBS       = Utility->GetPreCut("MuMuLsBS");       std::cout << __LINE__ << " : LSMUMUBS       = " << LSMUMUBS << std::endl;
  DCAMUMU        = Utility->GetPreCut("DCAMuMu");        std::cout << __LINE__ << " : DCAMUMU        = " << DCAMUMU << std::endl;
  DCAMUBS        = Utility->GetPreCut("DCAMuBS");        std::cout << __LINE__ << " : DCAMUBS        = " << DCAMUBS << std::endl;
  COSALPHAMUMUBS = Utility->GetPreCut("cosAlphaMuMuBS"); std::cout << __LINE__ << " : COSALPHAMUMUBS = " << COSALPHAMUMUBS << std::endl;
  MUMINPT        = Utility->GetPreCut("MinMupT");        std::cout << __LINE__ << " : MUMINPT        = " << MUMINPT << std::endl;
  MUMAXETA       = Utility->GetPreCut("MuEta");          std::cout << __LINE__ << " : MUMAXETA       = " << MUMAXETA << std::endl;
  MINMUMUPT      = Utility->GetPreCut("MuMupT");         std::cout << __LINE__ << " : MINMUMUPT      = " << MINMUMUPT << std::endl;
  MINMUMUINVMASS = Utility->GetPreCut("MinMuMuMass");    std::cout << __LINE__ << " : MINMUMUINVMASS = " << MINMUMUINVMASS << std::endl;
  MAXMUMUINVMASS = Utility->GetPreCut("MaxMuMuMass");    std::cout << __LINE__ << " : MAXMUMUINVMASS = " << MAXMUMUINVMASS << std::endl;
  
  // ##############################
  // # Loading pre-selection cuts #
  // ##############################
  B0MASSLOWLIMIT = Utility->GetPreCut("MinB0Mass"); std::cout << __LINE__ << " : B0MASSLOWLIMIT = " << B0MASSLOWLIMIT << std::endl;
  B0MASSUPLIMIT  = Utility->GetPreCut("MaxB0Mass"); std::cout << __LINE__ << " : B0MASSUPLIMIT  = " << B0MASSUPLIMIT << std::endl;
  CLB0VTX        = Utility->GetPreCut("B0VtxCL");   std::cout << __LINE__ << " : CLB0VTX        = " << CLB0VTX << std::endl;
  KSTMASSWINDOW  = Utility->GetPreCut("KstMass");   std::cout << __LINE__ << " : KSTMASSWINDOW  = " << KSTMASSWINDOW << std::endl;
  HADDCASBS      = Utility->GetPreCut("HadDCASBS"); std::cout << __LINE__ << " : HADDCASBS      = " << HADDCASBS << std::endl;
  MINHADPT       = Utility->GetPreCut("HadpT");     std::cout << __LINE__ << " : MINHADPT       = " << MINHADPT << std::endl;

  std::cout << "\n@@@ Global constants @@@" << std::endl;
  std::cout << __LINE__ << " : TRKMAXR    = " << TRKMAXR << std::endl;
  
  std::cout << __LINE__ << " : PRIVTXNDOF = " << PRIVTXNDOF << std::endl;
  std::cout << __LINE__ << " : PRIVTXMAXZ = " << PRIVTXMAXZ << std::endl;
  std::cout << __LINE__ << " : PRIVTXMAXR = " << PRIVTXMAXR << std::endl;
  
  std::cout << __LINE__ << " : MUVARTOLE  = " << MUVARTOLE << std::endl;
  std::cout << __LINE__ << " : HADVARTOLE = " << HADVARTOLE << std::endl;

  std::cout << __LINE__ << " : RCUTMU     = " << RCUTMU << std::endl;
  std::cout << __LINE__ << " : RCUTTRK    = " << RCUTTRK << std::endl;
}


void B0KstMuMu::endJob ()
{
  theTree->GetDirectory()->cd();
  theTree->Write();
}


void B0KstMuMu::endLuminosityBlock(const edm::LuminosityBlock& lumiBlock, const edm::EventSetup& iSetup)
{
  if ((doGenReco_ == 2) || (doGenReco_ == 3))
    {
      edm::Handle<GenFilterInfo> genFilter;
      lumiBlock.getByLabel("genFilterEfficiencyProducer", genFilter);

      NTuple->numEventsTried  = genFilter->numEventsTried();
      NTuple->numEventsPassed = genFilter->numEventsPassed();

      if (printMsg == true)
	{
	  std::cout << "\n@@@ End of a luminosity block @@@" << std::endl;
	  std::cout << __LINE__ << " : number of events tried  = " << NTuple->numEventsTried << std::endl;
	  std::cout << __LINE__ << " : number of events passed = " << NTuple->numEventsPassed << std::endl;
	}
    }
}


// Define this as a plug-in
DEFINE_FWK_MODULE(B0KstMuMu);
