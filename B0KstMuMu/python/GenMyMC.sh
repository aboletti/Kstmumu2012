################
# Instructions #
################

# 1. Copy GEN-Fragment.py in Configuration/Generator/python/
# 2. If not in repository copy name.dec in GeneratorInterface/ExternalDecays/data/
# 3. scram b
# 4. source GenMyMC.sh


### Create just GEN-MC ###
# Remove only the mumugenfilter from GEN-Fragment.py
cmsDriver.py Configuration/Generator/python/$1 -s GEN --conditions START53_V19F::All --datatier GEN-SIM --eventcontent GENRAW -n 100 --no_exec


### Create just GEN-SIM-MC for official production tests ###
#cmsDriver.py Configuration/Generator/python/$1 -s GEN,SIM --conditions START53_V19F::All --datatier GEN-SIM --eventcontent RAWSIM -n 100 --no_exec


### Create MC up to HLT ###
#cmsDriver.py Configuration/Generator/python/$1 -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT --conditions START53_V19F::All --datatier GEN-SIM-RAW --eventcontent RAWSIM -n 100 --no_exec


### Create RECO-MC from MC up to HLT (output RECO + SIM) ###
#cmsDriver.py Configuration/Generator/python/$1 -s RAW2DIGI,RECO --conditions START53_V19F::All --datatier GEN-SIM-RECO --eventcontent RECOSIM --filein file:name.root -n 100 --no_exec


### Create RECO-MC from MC up to HLT (output AOD + SIM) ###
#cmsDriver.py Configuration/Generator/python/$1 -s RAW2DIGI,RECO --conditions START53_V19F::All --datatier GEN-SIM-RECO --eventcontent AODSIM --filein file:name.root -n 100 --no_exec
