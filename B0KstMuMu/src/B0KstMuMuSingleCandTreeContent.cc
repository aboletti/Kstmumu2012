#include "../interface/B0KstMuMuSingleCandTreeContent.h"

B0KstMuMuSingleCandTreeContent::B0KstMuMuSingleCandTreeContent ()
{
  B0KstMuMuTreeContent();

  // ### B0 variables ###
  B0MassArb      = 0;
  B0pT           = 0;
  B0Eta          = 0;
  B0Phi          = 0;
  B0notB0bar     = true;
  kstMassArb     = 0;
  rightFlavorTag = true;

  // ### Trigger category ###
  TrigCat = 0;

  // ### Angular anaylsis ###
  CosThetaKArb       = 0;
  CosThetaMuArb      = 0;
  PhiKstMuMuPlaneArb = 0;
}

void B0KstMuMuSingleCandTreeContent::Init ()
{
  B0KstMuMuTreeContent::Init();
}

void B0KstMuMuSingleCandTreeContent::ClearNTuple ()
{
  B0KstMuMuTreeContent::ClearNTuple();

  // ### B0 variables ###
  B0MassArb      = 0;
  B0pT           = 0;
  B0Eta          = 0;
  B0Phi          = 0;
  B0notB0bar     = true;
  kstMassArb     = 0;
  rightFlavorTag = true;

  // ### Trigger category ###
  TrigCat = 0;

  // ### Angular anaylsis ###
  CosThetaKArb       = 0;
  CosThetaMuArb      = 0;
  PhiKstMuMuPlaneArb = 0;
}

void B0KstMuMuSingleCandTreeContent::MakeTreeBranches (TTree* theTree)
{
  B0KstMuMuTreeContent::MakeTreeBranches(theTree);

  // ### B0 variables ###
  theTree->Branch("B0MassArb",      &B0MassArb,      "B0MassArb/D");
  theTree->Branch("B0pT",           &B0pT,           "B0pT/D");
  theTree->Branch("B0Eta",          &B0Eta,          "B0Eta/D");
  theTree->Branch("B0Phi",          &B0Phi,          "B0Phi/D");
  theTree->Branch("B0notB0bar",     &B0notB0bar,     "B0notB0bar/O");
  theTree->Branch("kstMassArb",     &kstMassArb,     "kstMassArb/D");
  theTree->Branch("rightFlavorTag", &rightFlavorTag, "rightFlavorTag/O");

  // ### Trigger category ###
  theTree->Branch("TrigCat", &TrigCat, "TrigCat/I");

  // ### Angular anaylsis ###
  theTree->Branch("CosThetaKArb",       &CosThetaKArb,       "CosThetaKArb/D");
  theTree->Branch("CosThetaMuArb",      &CosThetaMuArb,      "CosThetaMuArb/D");
  theTree->Branch("PhiKstMuMuPlaneArb", &PhiKstMuMuPlaneArb, "PhiKstMuMuPlaneArb/D");
}

void B0KstMuMuSingleCandTreeContent::SetBranchAddresses (TTree* theTree)
{
  B0KstMuMuTreeContent::SetBranchAddresses(theTree);

  // ### B0 variables ###
  theTree->SetBranchAddress("B0MassArb",      &B0MassArb);
  theTree->SetBranchAddress("B0pT",           &B0pT);
  theTree->SetBranchAddress("B0Eta",          &B0Eta);
  theTree->SetBranchAddress("B0Phi",          &B0Phi);
  theTree->SetBranchAddress("B0notB0bar",     &B0notB0bar);
  theTree->SetBranchAddress("kstMassArb",     &kstMassArb);
  theTree->SetBranchAddress("rightFlavorTag", &rightFlavorTag);

  // ### Trigger category ###
  theTree->SetBranchAddress("TrigCat", &TrigCat);

  // ### Angular anaylsis ###
  theTree->SetBranchAddress("CosThetaKArb",       &CosThetaKArb);
  theTree->SetBranchAddress("CosThetaMuArb",      &CosThetaMuArb);
  theTree->SetBranchAddress("PhiKstMuMuPlaneArb", &PhiKstMuMuPlaneArb);
}

void B0KstMuMuSingleCandTreeContent::CopyCandidate (B0KstMuMuSingleCandTreeContent* NTupleIn, int index)
{
  CopyScalars(NTupleIn);
  B0KstMuMuTreeContent::CopyCandidate(NTupleIn,index);
}

void B0KstMuMuSingleCandTreeContent::CopyAllCandidates (B0KstMuMuSingleCandTreeContent* NTupleIn)
{
  CopyScalars(NTupleIn);
  B0KstMuMuTreeContent::CopyAllCandidates(NTupleIn);
}

void B0KstMuMuSingleCandTreeContent::CopyScalars (B0KstMuMuSingleCandTreeContent* NTupleIn)
{
  // ### B0 variables ###
  B0MassArb      = NTupleIn->B0MassArb;
  B0pT           = NTupleIn->B0pT;
  B0Eta          = NTupleIn->B0Eta;
  B0Phi          = NTupleIn->B0Phi;
  B0notB0bar     = NTupleIn->B0notB0bar;
  kstMassArb     = NTupleIn->kstMassArb;
  rightFlavorTag = NTupleIn->rightFlavorTag;

  // ### Trigger category ###
  TrigCat = NTupleIn->TrigCat;

  // ### Angular anaylsis ###
  CosThetaKArb       = NTupleIn->CosThetaKArb;
  CosThetaMuArb      = NTupleIn->CosThetaMuArb;
  PhiKstMuMuPlaneArb = NTupleIn->PhiKstMuMuPlaneArb;
}
