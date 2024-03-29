import FWCore.ParameterSet.Config as cms

# include  particle flow local reconstruction
from RecoParticleFlow.PFClusterProducer.particleFlowCluster_cff import *

from RecoParticleFlow.PFTracking.pfTrack_cfi import *
pfTrack.UseQuality = cms.bool(True)
pfTrack.TrackQuality = cms.string('highPurity')
pfTrack.TkColList = cms.VInputTag("hiSelectedTracks")
pfTrack.PrimaryVertexLabel = cms.InputTag("hiSelectedVertex")
pfTrack.MuColl = cms.InputTag("muons")

# run a trimmed down PF sequence with heavy-ion vertex, no conversions, nucl int, etc.
from RecoParticleFlow.Configuration.RecoParticleFlow_cff import *
particleFlowBlock.useConvBremPFRecTracks = cms.bool(False)
particleFlowBlock.useIterTracking = cms.bool(False)
particleFlowBlock.useNuclear = cms.bool(False)
particleFlowBlock.useConversions = cms.bool(False)

particleFlowTmp.vertexCollection = cms.InputTag("hiSelectedVertex")
particleFlowTmp.usePFElectrons = cms.bool(True)
particleFlowTmp.muons = cms.InputTag("muons")
particleFlowTmp.usePFConversions = cms.bool(False)

from RecoParticleFlow.PFTracking.pfTrackElec_cfi import *
pfTrackElec.applyGsfTrackCleaning = cms.bool(True)
pfTrackElec.PrimaryVertexLabel = cms.InputTag("hiSelectedVertex")

electronsCiCLoose.verticesCollection = cms.InputTag("hiSelectedVertex")

# local reco must run before electrons (RecoHI/HiEgammaAlgos), due to PF integration
HiParticleFlowLocalReco = cms.Sequence(particleFlowCluster
                                       * pfTrack
                                       * pfTrackElec
                                       )

#PF Reco runs after electrons
HiParticleFlowReco = cms.Sequence(pfGsfElectronCiCSelectionSequence
                                  * particleFlowBlock
                                  * particleFlowTmp
                                  )
