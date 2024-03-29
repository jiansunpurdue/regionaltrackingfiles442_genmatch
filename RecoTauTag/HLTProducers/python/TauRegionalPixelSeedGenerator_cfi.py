import FWCore.ParameterSet.Config as cms

tauRegionalPixelSeedGenerator = cms.EDProducer("SeedGeneratorFromRegionHitsEDProducer",
    OrderedHitsFactoryPSet = cms.PSet(
        ComponentName = cms.string('StandardHitPairGenerator'),
        SeedingLayers = cms.string('PixelLayerPairs')
    ),
    SeedComparitorPSet = cms.PSet(
        ComponentName = cms.string('none')
    ),
    RegionFactoryPSet = cms.PSet(
        ComponentName = cms.string('TauRegionalPixelSeedGenerator'),
        RegionPSet = cms.PSet(
            precise = cms.bool(True),
            deltaPhiRegion = cms.double(0.4),
            originHalfLength = cms.double(0.2),
            originRadius = cms.double(0.2),
            deltaEtaRegion = cms.double(0.3),
            ptMin = cms.double(5.0),
            JetSrc = cms.InputTag("icone5Tau1"),
            originZPos = cms.double(0.0),
            vertexSrc = cms.InputTag("pixelVertices")
        )
    ),
    TTRHBuilder = cms.string('WithTrackAngle')
)


