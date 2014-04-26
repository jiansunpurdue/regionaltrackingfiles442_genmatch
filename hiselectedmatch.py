import FWCore.ParameterSet.Config as cms
isMC = True
#isMC = False
hiReco = True
#hiReco = False
#reReco = False
reReco = True
#hasSimInfo = False
hasSimInfo = True
#genTag = "'hiSignal','generator'"
genTag = "hiSignal"
#genTag = "generator"
hltFilter = "*"
#hltFilter = "HLT_HIL2Mu*_v*"
#hltFilter = "HLT_HIJet80_v*"
trigResults = 'TriggerResults::RECO'
#trigResults = 'TriggerResults::HISIGNAL'
gTag = 'STARTHI44_V12::All'
#gTag = 'GR_R_44_V15::All'
hiMode = True
redoPFJets = False

# some important triggers:  HLT_Jet40_v1, HLT_HIL2Mu7_v1'

if hiReco:
    svTracks = "hiSecondaryVertexSelectedTracks"
#    svTracks = "hiGeneralTracks"
    pvProducer = cms.untracked.vstring("offlinePrimaryVertices")
    #print "hacked to look at hiGeneralTracks"
    #svTracks = "hiGeneralTracks"
    #pvProducer = "hiSelectedVertex"
else:
    svTracks = "generalTracks"
    pvProducer = cms.vstring("offlinePrimaryVertices")
    
print "Reco'ing SV's w/ ", svTracks, ", PV w/ ", pvProducer 

#process.SimpleMemoryCheck = cms.Service("SimpleMemoryCheck")
#process.Timing = cms.Service("Timing")
process = cms.Process("Analysis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')

if hiReco:
    process.load('Configuration.StandardSequences.ReconstructionHeavyIons_cff')
    process.load('Configuration.EventContent.EventContentHeavyIons_cff')
else:
    process.load('Configuration.StandardSequences.Reconstruction_cff')
    process.load('Configuration.EventContent.EventContent_cff')
    
if isMC:
    process.load('SimGeneral.MixingModule.mixNoPU_cfi')
    process.load('Configuration.StandardSequences.RawToDigi_cff')
else:
    process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
    
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')


process.GlobalTag.globaltag = "STARTHI44_V12::All"
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(15))
process.source = cms.Source( "PoolSource",
                             fileNames = cms.untracked.vstring(
 'file:/afs/cern.ch/work/j/jisun/public/020C2B72-5A10-E311-9997-848F69FD2484_HiFall11_DR44X_NoPileUp_STARTHI44_V12-v1.root'
#                'file:./BsMM_7TeV_cfi_py_GEN_SIM_DIGI.root'
#                             'file:./608BDABE-93E4-E211-864F-7845C4FC37AF.root'
                 )
                           )

# FileService is mandatory, as the following analyzer module
# will want it, to create output histogram file
#
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("regionaltracking.root")
)

######################
# Hi specific reco
######################
# load centrality
from CmsHi.Analysis2010.CommonFunctions_cff import *
overrideCentrality(process)

if isMC:
    process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowers"),
    nonDefaultGlauberModel = cms.string("Hydjet_Drum"),
    centralitySrc = cms.InputTag("hiCentrality")
    )
else:
    process.HeavyIonGlobalParameters = cms.PSet(
    centralityVariable = cms.string("HFtowers"),
    nonDefaultGlauberModel = cms.string(""),
    centralitySrc = cms.InputTag("hiCentrality")
    )


# need to produce centrality for other analyzers
process.load("RecoHI.HiCentralityAlgos.HiCentrality_cfi")
process.hiCentrality.producePixelTracks = cms.bool(False)
process.hiCentrality.produceETmidRapidity = cms.bool(False)
process.hiCentrality.producePixelhits = cms.bool(False)
process.hiCentrality.produceEcalhits = cms.bool(False)
process.hiCentrality.produceZDChits = cms.bool(False)
process.hiCentrality.produceBasicClusters = cms.bool(False)
process.hiCentrality.produceHFhits = cms.bool(True)
process.hiCentrality.produceTracks = cms.bool(False)



# add gen step for HI

if isMC:
    process.load('CmsHi.JetAnalysis.ExtraGenReco_cff')
    process.HiGenParticleAna = cms.EDAnalyzer("HiGenAnalyzer")
    process.HiGenParticleAna.src= cms.untracked.InputTag(genTag)   
    process.hiGenParticles.useCrossingFrame = cms.untracked.bool(False)
    process.hiGenParticles.srcVector = cms.vstring(genTag)
    process.higen_step          = cms.Path(     
        process.hiGenParticles * process.hiGenParticlesForJets * process.genPartons * process.hiPartons * process.hiRecoGenJets #* process.HiGenParticleAna
        )
    
########  ADD EXTRA RECO WITH REGIT

#######################
#   Analyzers
########################


#  Track Analyzers

process.load("MitHig.PixelTrackletAnalyzer.trackAnalyzer_cff")

#process.anaTrack.trackSrc = svTracks
process.anaTrack.trackSrc = "hiSelectedTracks"
process.anaTrack.vertexSrc = cms.vstring("hiSelectedVertex")
process.anaTrack.qualityString = "highPurity"

process.anaTrack.trackPtMin = 0
process.anaTrack.doTrack = cms.untracked.bool(False)
process.anaTrack.useQuality = False
process.anaTrack.doPFMatching = True
process.anaTrack.doSimTrack = hasSimInfo and isMC
process.anaTrack.fillSimTrack = cms.untracked.bool(hasSimInfo)
if hiMode:process.anaTrack.useCentrality = cms.untracked.bool(True)
else: process.anaTrack.useCentrality = cms.untracked.bool(False)


process.load("edwenger.HiTrkEffAnalyzer.hitrkEffAnalyzer_cff")

process.trackAnalyzers = cms.Sequence(
        process.anaTrack
        )            

#
process.ana_step          = cms.Path(         
    process.trackAnalyzers
    )
process.dcharged = cms.EDAnalyzer( "Dcharged_ccbar_Analyzer" )
process.dstar = cms.EDAnalyzer( "Dstar_ccbar_Analyzer" )
process.d0 = cms.EDAnalyzer( "D0_ccbar_Analyzer" )
process.dsphi = cms.EDAnalyzer( "Ds_ccbar_Analyzer" )
process.dsks = cms.EDAnalyzer( "DsKs_ccbar_Analyzer" )
process.dmeson = cms.Path(process.d0*process.dcharged*process.dstar*process.dsphi*process.dsks)


process.schedule = cms.Schedule(
#    process.higen_step,
	process.dmeson,
	process.ana_step
)

