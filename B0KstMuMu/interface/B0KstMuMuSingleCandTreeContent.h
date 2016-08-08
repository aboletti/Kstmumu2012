#ifndef B0KSTMUMUSINGLECANDTREECONTENT_H
#define B0KSTMUMUSINGLECANDTREECONTENT_H

#include "TTree.h"
#include "../interface/B0KstMuMuTreeContent.h"


class B0KstMuMuSingleCandTreeContent : public B0KstMuMuTreeContent
{
 public:
  
  B0KstMuMuSingleCandTreeContent ();

  void Init ();
  void ClearNTuple ();
  void MakeTreeBranches (TTree* theTree);
  void SetBranchAddresses (TTree* theTree);
  void CopyCandidate (B0KstMuMuSingleCandTreeContent* NTupleIn, int index);
  void CopyAllCandidates (B0KstMuMuSingleCandTreeContent* NTupleIn);


  // ################
  // # B0 variables #
  // ################
  double B0MassArb;
  double B0pT, B0Eta, B0Phi;
  bool   B0notB0bar;
  double kstMassArb;
  bool   rightFlavorTag;

  // ################################
  // # Trigger category:            #
  // # 0   = not in trigger list    #
  // # > 0 = trigger category index #
  // ################################
  int TrigCat;

  // ####################
  // # Angular anaylsis #
  // ####################
  double CosThetaKArb, CosThetaMuArb, PhiKstMuMuPlaneArb;


 private:

  void CopyScalars (B0KstMuMuSingleCandTreeContent* NTupleIn);
};

#endif
