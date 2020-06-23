
#DATA_TAG = "ReReco" # Change to PromptReco for Run2016 period H
LEPTON_SETUP = 2018  # current default = 2018 
#ELECORRTYPE = "None" # "None" to switch off
#ELEREGRESSION = "None" # "None" to switch off
#APPLYMUCORR = False  # Switch off muon scale corrections
#APPLYJEC = False     #
#APPLYJER = False     #
#RECORRECTMET = False #
KINREFIT = False    # control KinZFitter (very slow)
PROCESS_CR = False   # Uncomment to run CR paths and trees
#ADDLOOSEELE = True  # Run paths for loose electrons
#APPLYTRIG = False    # Skip events failing required triggers. They are stored with sel<0 if set to False 
#KEEPLOOSECOMB = True # Do not skip loose lepton ZZ combinations (for debugging)
ADDZTREE = False      # Add tree for Z analysis
BESTCANDCOMPARATOR = "byBestZ1bestZ2"

# tau parameters
TAUCUT = "pt>15"
APPLYTESCORRECTION = False # shift the central value of the tau energy scale before computing up/down variations
TAUDISCRIMINATOR = "byIsolationMVA3oldDMwoLTraw"


PD = ""
MCFILTER = ""

## ****************************************
## *** choose the xsec ( xsec = xsec*BR )
## *** XSEC HH = 0.03105 pb  # https://twiki.cern.ch/twiki/bin/view/LHCPhysics/LHCHXSWGHH

#XSEC = 0.03105 * 0.00014     # HH->4lbb
#XSEC = 0.03105 * 0.000015    # HH->4ltautau
#XSEC = 0.03105 * 0.0000023   # HH->4lww         #FIXME
#XSEC = 0.03105 * 0.00000054  # HH->4lgammagamma

#GENXSEC = XSEC
#GENBR   = 1
## ****************************************

# ***************************
# xsec for sync
#XSEC=48.58 * 0.0002745 #ggH
# ***************************

#For DATA: 
IsMC = False
PD = "DoubleMu"

# var parsing
import FWCore.ParameterSet.VarParsing as VarParsing
import sys

options = VarParsing.VarParsing()

options.register('inputFile',
'/store/data/Run2018B/DoubleMuon/MINIAOD/17Sep2018-v1/60000/FF098453-BB5E-5049-BE8E-6F144AC12FBA.root', #data2018B
#'/store/mc/RunIIAutumn18MiniAOD/GluGluHToZZTo4L_M125_13TeV_powheg2_JHUGenV7011_pythia8/MINIAODSIM/102X_upgrade2018_realistic_v15-v2/80000/D65A4D51-2E80-AD41-B50D-E4083BA2A668.root', #ggH 2018
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 "inputFile")

options.parseArguments()




# Get absolute path
import os
PyFilePath = os.environ['CMSSW_BASE'] + "/src/ZZXAnalysis/AnalysisStep/test/"

### ----------------------------------------------------------------------
### Standard sequence
### ----------------------------------------------------------------------

execfile(PyFilePath + "analyzer.py")
execfile(PyFilePath + "prod/pyFragments/RecoProbabilities.py")

if not IsMC:
	process.source.inputCommands = cms.untracked.vstring("keep *", "drop LHERunInfoProduct_*_*_*", "drop LHEEventProduct_*_*_*") ###FIXME In 9X this removes all collections for MC


### ----------------------------------------------------------------------
### Replace parameters
### ----------------------------------------------------------------------

# keep events even if they have same runNumber:lumiNumber:eventNumber
# NOT TO BE USED
# process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')


# --- input files
process.source.fileNames = cms.untracked.vstring( options.inputFile

### HH files
    
    #'root://eoscms//eos/cms/store/user/covarell/HH/SM/4lbb/testMINIAOD_HHSM_4lbb_1.root'
    #'root://eoscms//eos/cms/store/user/covarell/HH/SM/4lgammagamma/testMINIAOD_HHSM_4lgammagamma_3.root'
    #'root://eoscms//eos/cms/store/user/covarell/HH/SM/4ltautau/testMINIAOD_HHSM_4ltautau_1.root' 
    #'root://eoscms//eos/cms/store/user/covarell/HH/SM/6l2nu/testMINIAOD_HHSM_4lWW_1.root'

    )


#process.calibratedPatElectrons.isSynchronization = cms.bool(True)
process.calibratedMuons.isSynchronization = cms.bool(True)

#process.maxEvents.input = -1
process.maxEvents.input = 5000
#process.source.skipEvents = cms.untracked.uint32(5750)

# Silence output
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000


### ----------------------------------------------------------------------
### Analyzer for Plots
### ----------------------------------------------------------------------


#process.source.eventsToProcess = cms.untracked.VEventRange("1:8670")

# Debug
process.dumpUserData =  cms.EDAnalyzer("dumpUserData",
     dumpTrigger = cms.untracked.bool(True),
     muonSrcs = cms.PSet(
#       slimmedMuons = cms.InputTag("slimmedMuons"),
        muons = cms.InputTag("appendPhotons:muons"),
     ),
     electronSrcs = cms.PSet(
#       slimmedElectron = cms.InputTag("slimmedElectrons"),
        electrons = cms.InputTag("appendPhotons:electrons"),
#        RSE = cms.InputTag("appendPhotons:looseElectrons"),
#        TLE = cms.InputTag("appendPhotons:electronstle"), 
     ),
     candidateSrcs = cms.PSet(
        Z     = cms.InputTag("ZCand"),
#        ZRSE     = cms.InputTag("ZCandlooseEle"),
#        ZTLE     = cms.InputTag("ZCandtle"),
        ZZ  = cms.InputTag("ZZCand"),
#        ZZRSE     = cms.InputTag("ZZCandlooseEle"),
#        ZZTLE     = cms.InputTag("ZZCandtle"),
        ZLL  = cms.InputTag("ZLLCand"),
        ZL  = cms.InputTag("ZlCand"),
     ),
     jetSrc = cms.InputTag("cleanJets"),
#     photonSrc = cms.InputTag("pikaPhotons"),
)

# Create lepton sync file
#process.PlotsZZ.dumpForSync = True;
#process.p = cms.EndPath( process.PlotsZZ)

# Keep all events in the tree, even if no candidate is selected
#process.ZZTree.skipEmptyEvents = False

# replace the paths in analyzer.py
#process.trees = cms.EndPath(process.ZZTree)

#Dump reconstructed variables
#process.appendPhotons.debug = cms.untracked.bool(True)
#process.fsrPhotons.debug = cms.untracked.bool(True)
#process.dump = cms.Path(process.dumpUserData)

#Print MC history
#process.mch = cms.EndPath(process.printTree)


#Monitor memory usage
#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
#    ignoreTotal = cms.untracked.int32(1)
#)
