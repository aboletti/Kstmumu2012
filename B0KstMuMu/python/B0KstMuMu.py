#################
### Variables ###
#################
runDataMC          = 1 # 1 = Data; 2 = MC (Reco + Gen); 3 = MC (Gen)
useJSON            = True
printMsg           = False
run2012not2011     = False
triggerProcessName = 'HLT' # 'GEN' or 'HLT' or 'RECO' or 'TEST' or ...

print "\n@@@ CMSSW run configuration flags @@@"
print "runDataMC          : ", runDataMC
print "useJSON            : ", useJSON
print "printMsg           : ", printMsg
print "run2012not2011     : ", run2012not2011
print "triggerProcessName : ", triggerProcessName


#####################
### CMSSW configs ###
#####################
import FWCore.ParameterSet.Config as cms
process = cms.Process('B0KSTMUMUNTUPLIZER')

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryIdeal_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.MessageLogger.suppressWarning = cms.untracked.vstring('B0KstMuMu')


#################
### GlobalTag ###
#################
### Auto detect GlobalTag ###
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:com10_8E33v2', '')
if (runDataMC != 1):
    if (run2012not2011 == True):
        process.GlobalTag.globaltag = cms.string('START53_V19F::All') # Run dep. MC: START53_V19F; J/psi X MC: START53_V7G
    else:
        process.GlobalTag.globaltag = cms.string('START53_LV6::All')
else:
    if (run2012not2011 == True):
        process.GlobalTag.globaltag = cms.string('FT53_V21A_AN6::All')
    else:
        process.GlobalTag.globaltag = cms.string('FT_53_LV5_AN1::All')


#######################
### Read input file ###
#######################
import sys
if (len(sys.argv) > 2):
    readFiles = sys.argv[2]
### GEN MC ###
    path = 'file:/mnt/hadoop/store/user/dinardo/PYTHIA6_Bd2KstarMuMu_EtaPtFilter_TuneZ2star_8TeV_GEN_NoFilter/'
#    path = 'file:/mnt/hadoop/store/user/dinardo/PYTHIA6_Bd2JpsiKstar_EtaPtFilter_TuneZ2star_8TeV_GEN_NoFilter/'
#    path = 'file:/mnt/hadoop/store/user/dinardo/PYTHIA6_Bd2Psi2SKstar_EtaPtFilter_TuneZ2star_8TeV_GEN_NoFilter/'
    file = readFiles.replace(path, '')
else:
### GEN MC ###
    from B0ToKstMuMu_GEN_NoFilter_MC_cff import readFiles
#    from B0ToJPsiKst_GEN_NoFilter_MC_cff import readFiles
#    from B0ToPsi2SKst_GEN_NoFilter_MC_cff import readFiles
process.source = cms.Source('PoolSource',
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            skipEvents = cms.untracked.uint32(0),
                            fileNames = cms.untracked.vstring(readFiles))


##################################
### Use JSON file interatively ###
##################################
if (runDataMC == 1 and useJSON == True):
    import FWCore.PythonUtilities.LumiList as LumiList
    import FWCore.ParameterSet.Types as CfgTypes
    if (run2012not2011 == True):
        myLumis = LumiList.LumiList(filename = 'Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON_MuonPhys.txt').getCMSSWString().split(',')
    else:
        myLumis = LumiList.LumiList(filename = 'Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON_MuonPhys_v2.txt').getCMSSWString().split(',')
    process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
    process.source.lumisToProcess.extend(myLumis)


###################
### Output file ###
###################
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.maxLuminosityBlocks = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.TFileService = cms.Service('TFileService', fileName = cms.string(
    'B0ToKstMuMu.root'
### GEN MC ###
#    '/nfs/data36/cms/dinardo/B0ToKstMuMu_GEN_NoFilter_MC_NTuples/' + file
#    '/nfs/data36/cms/dinardo/B0ToJPsiKst_GEN_NoFilter_MC_NTuples/' + file
#    '/nfs/data36/cms/dinardo/B0ToPsi2SKst_GEN_NoFilter_MC_NTuples/' + file
    ), closeFileFast = cms.untracked.bool(True))


###################################
### Import PAT skeleton process ###
###################################
process.load("PhysicsTools.PatAlgos.patSequences_cff")


############################
### Add track candidates ###
############################
process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
# Require NOT to check overlap with muons and electrons
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps     = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)


#########################
### Set up PAT tracks ###
#########################
from PhysicsTools.PatAlgos.tools.trackTools import *
from PhysicsTools.PatAlgos.tools.coreTools  import *
makeTrackCandidates(process,
                    label        = 'TrackCands',                  # output collection
                    tracks       = cms.InputTag('generalTracks'), # input track collection
                    particleType = 'pi+',                         # particle type (for assigning a mass)
                    preselection = 'pt > 0.1',                    # preselection cut on candidates
                    selection    = 'pt > 0.1',                    # selection cut on candidates
                    isolation    = {},                            # isolations to use (set to {} for None)
                    isoDeposits  = [],
                    mcAs         = None)                          # replicate MC match as the one used for Muons
removeMCMatching(process, ['All'], outputModules = [])


#####################################
### Do trigger matching for muons ###
#####################################

if (run2012not2011 == True):
    process.cleanMuonTriggerMatchHLT = cms.EDProducer(
        'PATTriggerMatcherDRLessByR',
        src                   = cms.InputTag('cleanPatMuons'),
        matched               = cms.InputTag('patTrigger'),
        matchedCuts           = cms.string('path("HLT_DoubleMu3p5_LowMass_Displaced_v*")'),
        maxDeltaR             = cms.double(0.1),
        resolveAmbiguities    = cms.bool(True),
        resolveByMatchQuality = cms.bool(True))
else:
    process.cleanMuonTriggerMatchHLT1 = cms.EDProducer(
        'PATTriggerMatcherDRLessByR',
        src                   = cms.InputTag('cleanPatMuons'),
        matched               = cms.InputTag('patTrigger'),
        matchedCuts           = cms.string('path("HLT_Dimuon7_LowMass_Displaced_v*")'),
        maxDeltaR             = cms.double(0.1),
        resolveAmbiguities    = cms.bool(True),
        resolveByMatchQuality = cms.bool(True))
    process.cleanMuonTriggerMatchHLT2 = cms.EDProducer(
        'PATTriggerMatcherDRLessByR',
        src                   = cms.InputTag('cleanPatMuons'),
        matched               = cms.InputTag('patTrigger'),
        matchedCuts           = cms.string('path("HLT_DoubleMu4_LowMass_Displaced_v*")'),
        maxDeltaR             = cms.double(0.1),
        resolveAmbiguities    = cms.bool(True),
        resolveByMatchQuality = cms.bool(True))
    process.cleanMuonTriggerMatchHLT3 = cms.EDProducer(
        'PATTriggerMatcherDRLessByR',
        src                   = cms.InputTag('cleanPatMuons'),
        matched               = cms.InputTag('patTrigger'),
        matchedCuts           = cms.string('path("HLT_DoubleMu4p5_LowMass_Displaced_v*")'),
        maxDeltaR             = cms.double(0.1),
        resolveAmbiguities    = cms.bool(True),
        resolveByMatchQuality = cms.bool(True))
    process.cleanMuonTriggerMatchHLT4 = cms.EDProducer(
        'PATTriggerMatcherDRLessByR',
        src                   = cms.InputTag('cleanPatMuons'),
        matched               = cms.InputTag('patTrigger'),
        matchedCuts           = cms.string('path("HLT_DoubleMu5_LowMass_Displaced_v*")'),
        maxDeltaR             = cms.double(0.1),
        resolveAmbiguities    = cms.bool(True),
        resolveByMatchQuality = cms.bool(True))


##################################
### Switch on PAT trigger info ###
##################################
from PhysicsTools.PatAlgos.tools.trigTools import *
if (run2012not2011 == True):
    switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT'], hltProcess = triggerProcessName, outputModule = '')
else:
    switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT1','cleanMuonTriggerMatchHLT2','cleanMuonTriggerMatchHLT3','cleanMuonTriggerMatchHLT4'], hltProcess = triggerProcessName, outputModule = '')


#########################
### Remove unused PAT ###
#########################
process.patDefaultSequence.remove(process.muonMatch)
process.patDefaultSequence.remove(process.electronMatch)
process.patDefaultSequence.remove(process.photonMatch)
process.patDefaultSequence.remove(process.tauMatch)
process.patDefaultSequence.remove(process.tauGenJetMatch)

process.patDefaultSequence.remove(process.tauGenJets)
process.patDefaultSequence.remove(process.tauGenJetsSelectorAllHadrons)

process.patDefaultSequence.remove(process.patJetPartonMatch)
process.patDefaultSequence.remove(process.patJetGenJetMatch)
process.patDefaultSequence.remove(process.patJetCorrFactors)
process.patDefaultSequence.remove(process.patJetCharge)
process.patDefaultSequence.remove(process.patJetPartons)
process.patDefaultSequence.remove(process.patJetPartonAssociation)
process.patDefaultSequence.remove(process.patJetFlavourAssociation)
process.patDefaultSequence.remove(process.patJets)

process.patDefaultSequence.remove(process.selectedPatJets)
process.patDefaultSequence.remove(process.cleanPatJets)
process.patDefaultSequence.remove(process.countPatJets)

process.patDefaultSequence.remove(process.patCandidateSummary)
process.patDefaultSequence.remove(process.selectedPatCandidateSummary)
process.patDefaultSequence.remove(process.cleanPatCandidateSummary)

process.patDefaultSequence.remove(process.patMETs)


##########################
### B0 --> K*0 mu+ mu- ###
##########################
process.B0KstMuMu = cms.EDAnalyzer('B0KstMuMu',
                                   HLTriggerResults = cms.untracked.string(triggerProcessName),
                                   VtxSample        = cms.untracked.string('offlinePrimaryVertices'),
                                   BeamSpot         = cms.untracked.string('offlineBeamSpot'),
                                   GenParticles     = cms.untracked.string('genParticles'),
                                   MuonType         = cms.untracked.string('cleanPatMuonsTriggerMatch'),
                                   TrackType        = cms.untracked.string('cleanPatTrackCands'),
                                   ParameterFile    = cms.untracked.string('ParameterFile.txt'),
                                   doGenReco        = cms.untracked.uint32(runDataMC),
                                   printMsg         = cms.untracked.bool(printMsg))


###########
### RUN ###
###########
### Run unscheduled = create and read just what I need ###
#process.options = cms.untracked.PSet(allowUnscheduled = cms.untracked.bool(True))

process.patPath  = cms.Path(process.patDefaultSequence)
process.ntupPath = cms.Path(process.B0KstMuMu)

if (runDataMC == 3):
    process.schedule = cms.Schedule(process.ntupPath)
else:
    process.schedule = cms.Schedule(process.patPath, process.ntupPath)
