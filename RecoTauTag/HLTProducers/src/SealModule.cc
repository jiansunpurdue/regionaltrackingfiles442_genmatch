#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoTauTag/HLTProducers/interface/IsolatedTauJetsSelector.h"
#include "RecoTauTag/HLTProducers/interface/PFTauToJetProducer.h"
#include "RecoTauTag/HLTProducers/interface/PFJetToCaloProducer.h"
#include "RecoTauTag/HLTProducers/interface/EMIsolatedTauJetsSelector.h"
#include "RecoTauTag/HLTProducers/interface/L2TauJetsProvider.h"
#include "RecoTauTag/HLTProducers/interface/L1HLTJetsMatching.h"
#include "RecoTauTag/HLTProducers/interface/L2TauJetsMerger.h"
#include "RecoTauTag/HLTProducers/interface/CaloTowerCreatorForTauHLT.h"
#include "RecoTauTag/HLTProducers/interface/HLTTauProducer.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionProducerFactory.h" 	 
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionProducer.h" 	 
#include "TauRegionalPixelSeedGenerator.h" 	 
#include "RecoTauTag/HLTProducers/interface/L2TauIsolationSelector.h"
#include "RecoTauTag/HLTProducers/interface/L2TauRelaxingIsolationSelector.h"
#include "RecoTauTag/HLTProducers/interface/L2TauIsolationProducer.h"
#include "RecoTauTag/HLTProducers/interface/L2TauNarrowConeIsolationProducer.h"
#include "RecoTauTag/HLTProducers/interface/L2TauModularIsolationProducer.h"
#include "RecoTauTag/HLTProducers/interface/L2TauModularIsolationSelector.h"
#include "RecoTauTag/HLTProducers/interface/TauJetSelectorForHLTTrackSeeding.h"
//#include "RecoTauTag/HLTProducers/interface/PFJetIsolator.h"



 	 
DEFINE_EDM_PLUGIN(TrackingRegionProducerFactory, TauRegionalPixelSeedGenerator, "TauRegionalPixelSeedGenerator"); 	 

//
DEFINE_FWK_MODULE(IsolatedTauJetsSelector);
DEFINE_FWK_MODULE(EMIsolatedTauJetsSelector);
DEFINE_FWK_MODULE(L2TauJetsProvider);
DEFINE_FWK_MODULE(L2TauJetsMerger);
DEFINE_FWK_MODULE(L1HLTJetsMatching);
DEFINE_FWK_MODULE(CaloTowerCreatorForTauHLT);
DEFINE_FWK_MODULE(HLTTauProducer);
DEFINE_FWK_MODULE(PFTauToJetProducer);
DEFINE_FWK_MODULE(PFJetToCaloProducer);
DEFINE_FWK_MODULE(L2TauIsolationProducer);
DEFINE_FWK_MODULE(L2TauNarrowConeIsolationProducer);
DEFINE_FWK_MODULE(L2TauModularIsolationProducer);
DEFINE_FWK_MODULE(L2TauModularIsolationSelector);
DEFINE_FWK_MODULE(L2TauIsolationSelector);
DEFINE_FWK_MODULE(L2TauRelaxingIsolationSelector);
DEFINE_FWK_MODULE(TauJetSelectorForHLTTrackSeeding);
//DEFINE_FWK_MODULE(PFJetIsolator);
