#!/bin/bash

cd /lustre/cmswork/boletti/Kstmumu/CMSSW_5_3_28/src/Stefano/B0KstMuMu/plugins/

source  /cvmfs/cms.cern.ch/cmsset_default.sh
pwd
eval `scram runtime -sh`

mkdir Data_scan_hprof2
cd Data_scan_hprof2
mkdir ${2}
cd ${2}
ln -s ../../../efficiency/effKEpdf_out_RT.root . 
ln -s ../../../efficiency/effKEpdf_out_WT.root .
export ANALYPATH=/lustre/cmswork/boletti/Kstmumu/CMSSW_5_3_28/src/Stefano/B0KstMuMu/
export HOME=/home/boletti/
# mv effKEpdf_out.root effKEpdf_out_RT.root
# ln -s ../../Closure_sys2_wt_1.0/${2}/effKEpdf_out.root effKEpdf_out_WT.root

# ../../ExtractYield 106 dcap://t2-srm-02.lnl.infn.it////pnfs/lnl.infn.it/data/cms/store/user/slacapra/B0KsMuMu/Data2012B0KstMuMuResults/MonteCarlo2012/SingleCand/singleCand_B0ToKstMuMu_MC_NTuple.root yesEffCorr ${1}

# ../../ExtractYield 116 /lustre/cmswork/boletti/Kstmumu/CMSSW_5_3_28/src/B0Analysis/B0KstMuMu/plugins/toyDatasets.root yesEffCorr ${1}

# ../../ExtractYield 6 dcap://t2-srm-02.lnl.infn.it////pnfs/lnl.infn.it/data/cms/store/user/slacapra/B0KsMuMu/Data2012B0KstMuMuResults/singleCand_B0ToKstMuMu_Data2012ABCD_NTuples.root yesEffCorr ${1} 0 "" ${2}

../../ExtractYield 6 /lustre/cmswork/boletti/Kstmumu/CMSSW_5_3_28/src/Stefano/B0KstMuMu/plugins/singleCand_B0ToKstMuMu_Data2012ABCD_NTuples.root yesEffCorr ${1} 0 "" ${2}
