#!/bin/bash

cd /lustre/cmswork/boletti/Kstmumu/CMSSW_5_3_28/src/Stefano/B0KstMuMu/plugins/

source  /cvmfs/cms.cern.ch/cmsset_default.sh
pwd
eval `scram runtime -sh`

mkdir TestP1_freefit
cd TestP1_freefit
ln -s ../../efficiency/effKEpdf_out_RT.root . 
ln -s ../../efficiency/effKEpdf_out_WT.root .
export ANALYPATH=/lustre/cmswork/boletti/Kstmumu/CMSSW_5_3_28/src/Stefano/B0KstMuMu/
export HOME=/home/boletti/

../ExtractYieldTest 16 /lustre/cmswork/boletti/Kstmumu/CMSSW_5_3_28/src/Stefano/B0KstMuMu/plugins/Data_test yesEffCorr ${1} 0 "" 0

cd ..

