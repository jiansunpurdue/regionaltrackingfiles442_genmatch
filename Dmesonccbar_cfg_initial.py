import FWCore.ParameterSet.Config as cms

process = cms.Process("Analysis")

# run the input file through the end;
# for a limited number of events, replace -1 with the desired number 
#
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load( "SimGeneral.HepPDTESSource.pythiapdt_cfi" )
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
process.load("Configuration.StandardSequences.Geometry_cff")
process.load("RecoVertex.BeamSpotProducer.BeamSpot_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = "STARTHI44_V12::All"

process.source = cms.Source( "PoolSource",
                             fileNames = cms.untracked.vstring(
 'file:/afs/cern.ch/work/j/jisun/public/020C2B72-5A10-E311-9997-848F69FD2484_HiFall11_DR44X_NoPileUp_STARTHI44_V12-v1.root'
#'root://xrootd.unl.edu//store/himc/HiFall11/Pyquen_DiJet_Pt50_TuneZ2_Unquenched_2760GeV/GEN-SIM-RECODEBUG/HiFall11_DR44X_NoPileUp_STARTHI44_V12-v1/10000/001C9742-6F10-E311-94BC-008CFA000F5C.root'
#			     'file:./BsMM_7TeV_cfi_py_GEN_SIM_DIGI.root'
#                             'file:./608BDABE-93E4-E211-864F-7845C4FC37AF.root'
			     )
                           )
	      
# FileService is mandatory, as the following analyzer module 
# will want it, to create output histogram file
# 
process.TFileService = cms.Service("TFileService",
        fileName = cms.string("dmeson.root")
)


#process.load("PhysicsTools.HepMCCandAlgos.genParticles_cfi")
#process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
#genParticles = cms.EDProducer("GenParticleProducer",
#                  src = cms.InputTag("source")
#)
#process.load( "" )
# the analyzer itself - empty parameter set 
#
#process.TestHepMCEvt = cms.EDAnalyzer( "HZZ4muExampleAnalyzer" )

process.dcharged = cms.EDAnalyzer( "Dcharged_ccbar_Analyzer" )
process.dstar = cms.EDAnalyzer( "Dstar_ccbar_Analyzer" )
process.d0 = cms.EDAnalyzer( "D0_ccbar_Analyzer" )
process.dsphi = cms.EDAnalyzer( "Ds_ccbar_Analyzer" )
process.dsks = cms.EDAnalyzer( "DsKs_ccbar_Analyzer" )
process.p1 = cms.Path(process.dcharged * process.dstar * process.d0 * process.dsphi * process.dsks)
#process.p1 = cms.Path(process.d0)
#process.p1 = cms.Path(process.genParticles * process.d0)
