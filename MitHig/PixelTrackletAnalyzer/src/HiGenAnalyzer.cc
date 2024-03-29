// -*- C++ -*-
//
// Package:    HiGenAnalyzer
// Class:      HiGenAnalyzer
// 
/**\class HiGenAnalyzer HiGenAnalyzer.cc

Description: Analyzer that studies (HI) gen event info

Implementation:
<Notes on implementation>
 */
//
// Original Author:  Yetkin Yilmaz, Frank Ma
//         Created:  Tue Dec 18 09:44:41 EST 2007
// $Id: HiGenAnalyzer.cc,v 1.7 2012/06/04 22:19:34 yilmaz Exp $
//
//


// system include files
#include <memory>
#include <iostream>
#include <string>
#include <fstream>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/HiGenData/interface/GenHIEvent.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "HepMC/GenEvent.h"
#include "HepMC/HeavyIon.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

// root include file
#include "TFile.h"
#include "TNtuple.h"

using namespace std;

static const Int_t MAXPARTICLES = 50000;
static const Int_t MAXVTX = 1000;
static const Int_t ETABINS = 3; // Fix also in branch string

//
// class decleration
//

struct HydjetEvent{

  Int_t event;
  Float_t b;
  Float_t npart;
  Float_t ncoll;
  Float_t nhard;
  Float_t phi0;
  Float_t scale;

  Int_t n[ETABINS];
  Float_t ptav[ETABINS];

  Int_t mult;
  Float_t pt[MAXPARTICLES];
  Float_t eta[MAXPARTICLES];
  Float_t phi[MAXPARTICLES];
  Int_t pdg[MAXPARTICLES];
  Int_t chg[MAXPARTICLES];
  Int_t sube[MAXPARTICLES];

  Float_t vx;
  Float_t vy;
  Float_t vz;
  Float_t vr;

};

class HiGenAnalyzer : public edm::EDAnalyzer {
  public:
    explicit HiGenAnalyzer(const edm::ParameterSet&);
    ~HiGenAnalyzer();


  private:
    virtual void beginRun(const edm::Run&, const edm::EventSetup&) ;
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    // ----------member data ---------------------------

    std::ofstream out_b;
    std::string fBFileName;

    std::ofstream out_n;
    std::string fNFileName;

    std::ofstream out_m;
    std::string fMFileName;


    TTree* hydjetTree_;
    HydjetEvent hev_;

    TNtuple *nt;

    std::string output;           // Output filename

    Bool_t doAnalysis_;
    Bool_t printLists_;
    Bool_t doCF_;
    Bool_t doVertex_;
    Bool_t useHepMCProduct_;
    Bool_t doHI_;
    Bool_t doParticles_;

    Double_t etaMax_;
    Double_t ptMin_;
    Bool_t chargedOnly_;
    edm::InputTag src_;
    edm::InputTag genParticleSrc_;
    edm::InputTag genHIsrc_;

    edm::ESHandle < ParticleDataTable > pdt;
    edm::Service<TFileService> f;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HiGenAnalyzer::HiGenAnalyzer(const edm::ParameterSet& iConfig)
{
  //now do what ever initialization is needed
  fBFileName = iConfig.getUntrackedParameter<std::string>("output_b", "b_values.txt");
  fNFileName = iConfig.getUntrackedParameter<std::string>("output_n", "n_values.txt");
  fMFileName = iConfig.getUntrackedParameter<std::string>("output_m", "m_values.txt");
  doAnalysis_ = iConfig.getUntrackedParameter<Bool_t>("doAnalysis", true);
  useHepMCProduct_ = iConfig.getUntrackedParameter<Bool_t>("useHepMCProduct", false);
  printLists_ = iConfig.getUntrackedParameter<Bool_t>("printLists", false);
  doCF_ = iConfig.getUntrackedParameter<Bool_t>("doMixed", false);
  doVertex_ = iConfig.getUntrackedParameter<Bool_t>("doVertex", false);
  etaMax_ = iConfig.getUntrackedParameter<Double_t>("etaMax", 2);
  ptMin_ = iConfig.getUntrackedParameter<Double_t>("ptMin", 0);
  chargedOnly_ = iConfig.getUntrackedParameter<Bool_t>("chargedOnly", true);
  src_ = iConfig.getUntrackedParameter<edm::InputTag>("src",edm::InputTag("generator"));
  genParticleSrc_ = iConfig.getUntrackedParameter<edm::InputTag>("genpSrc",edm::InputTag("hiGenParticles"));
  genHIsrc_ = iConfig.getUntrackedParameter<edm::InputTag>("genHiSrc",edm::InputTag("heavyIon"));
  doParticles_ = iConfig.getUntrackedParameter<Bool_t>("doParticles", true);
}


HiGenAnalyzer::~HiGenAnalyzer()
{
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called to for each event  ------------
  void
HiGenAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace HepMC;

  hev_.event = iEvent.id().event();
  for(Int_t ieta = 0; ieta < ETABINS; ++ieta){
    hev_.n[ieta] = 0;
    hev_.ptav[ieta] = 0;
  }
  hev_.mult = 0;

  Double_t phi0 = 0;
  Double_t b = -1;
  Double_t scale = -1;
  Int_t npart = -1;
  Int_t ncoll = -1;
  Int_t nhard = -1;
  Double_t vx = -99;
  Double_t vy = -99;
  Double_t vz = -99;
  Double_t vr = -99;
  const GenEvent* evt;

  Int_t nmix = -1;
  Int_t np = 0;
  Int_t sig = -1;
  Int_t src = -1;

  if(useHepMCProduct_){
    if(doCF_){
      Handle<CrossingFrame<HepMCProduct> > cf;
      iEvent.getByLabel(InputTag("mix","source"),cf);
      MixCollection<HepMCProduct> mix(cf.product());
      nmix = mix.size();
      cout<<"Mix Collection Size: "<<mix<<endl;

      MixCollection<HepMCProduct>::iterator mbegin = mix.begin();
      MixCollection<HepMCProduct>::iterator mend = mix.end();

      for(MixCollection<HepMCProduct>::iterator mixit = mbegin; mixit != mend; ++mixit){
	const GenEvent* subevt = (*mixit).GetEvent();
	Int_t all = subevt->particles_size();
	np += all;
	HepMC::GenEvent::particle_const_iterator begin = subevt->particles_begin();
	HepMC::GenEvent::particle_const_iterator end = subevt->particles_end();
	for(HepMC::GenEvent::particle_const_iterator it = begin; it != end; ++it){
	  if ((*it)->momentum().perp()<ptMin_) continue;
	  if((*it)->status() == 1){
	    Int_t pdg_id = (*it)->pdg_id();
	    Float_t eta = (*it)->momentum().eta();
	    Float_t phi = (*it)->momentum().phi();
	    Float_t pt = (*it)->momentum().perp();
	    const ParticleData * part = pdt->particle(pdg_id );
	    Int_t charge = static_cast<Int_t>(part->charge());
	    if (chargedOnly_&&charge==0) continue;

	    hev_.pt[hev_.mult] = pt;
	    hev_.eta[hev_.mult] = eta;
	    hev_.phi[hev_.mult] = phi;
	    hev_.pdg[hev_.mult] = pdg_id;
	    hev_.chg[hev_.mult] = charge;

	    eta = fabs(eta);
	    Int_t etabin = 0;
	    if(eta > 0.5) etabin = 1;
	    if(eta > 1.) etabin = 2;
	    if(eta < 2.){
	      hev_.ptav[etabin] += pt;
	      ++(hev_.n[etabin]);
	    }
	    ++(hev_.mult);
	  }
	}
      }
    }else{

      Handle<HepMCProduct> mc;
      iEvent.getByLabel(src_,mc);
      evt = mc->GetEvent();
      scale = evt->event_scale();

      const HeavyIon* hi = evt->heavy_ion();
      if(hi){
	b = hi->impact_parameter();
	npart = hi->Npart_proj()+hi->Npart_targ();
	ncoll = hi->Ncoll();
	nhard = hi->Ncoll_hard();
	phi0 = hi->event_plane_angle();

	if(printLists_){
	  out_b<<b<<endl;
	  out_n<<npart<<endl;
	}
      }

      src = evt->particles_size();

      HepMC::GenEvent::particle_const_iterator begin = evt->particles_begin();
      HepMC::GenEvent::particle_const_iterator end = evt->particles_end();
      for(HepMC::GenEvent::particle_const_iterator it = begin; it != end; ++it){
	if ((*it)->momentum().perp()<ptMin_) continue;
	if((*it)->status() == 1){
	  Int_t pdg_id = (*it)->pdg_id();
	  Float_t eta = (*it)->momentum().eta();
	  Float_t phi = (*it)->momentum().phi();
	  Float_t pt = (*it)->momentum().perp();
	  const ParticleData * part = pdt->particle(pdg_id );
	  Int_t charge = static_cast<Int_t>(part->charge());
	  if (chargedOnly_&&charge==0) continue;

	  hev_.pt[hev_.mult] = pt;
	  hev_.eta[hev_.mult] = eta;
	  hev_.phi[hev_.mult] = phi;
	  hev_.pdg[hev_.mult] = pdg_id;
	  hev_.chg[hev_.mult] = charge;

	  eta = fabs(eta);
	  Int_t etabin = 0;
	  if(eta > 0.5) etabin = 1; 
	  if(eta > 1.) etabin = 2;
	  if(eta < 2.){
	    hev_.ptav[etabin] += pt;
	    ++(hev_.n[etabin]);
	  }
	  ++(hev_.mult);
	}
      }
    }
  }else{
    edm::Handle<reco::GenParticleCollection> parts;
    iEvent.getByLabel(genParticleSrc_,parts);
    for(UInt_t i = 0; i < parts->size(); ++i){
      const reco::GenParticle& p = (*parts)[i];
      if (p.status()!=1) continue;
      if (p.pt()<ptMin_) continue;
      if (chargedOnly_&&p.charge()==0) continue;
      hev_.pt[hev_.mult] = p.pt();
      hev_.eta[hev_.mult] = p.eta();
      hev_.phi[hev_.mult] = p.phi();
      hev_.pdg[hev_.mult] = p.pdgId();
      hev_.chg[hev_.mult] = p.charge();
      hev_.sube[hev_.mult] = p.collisionId();
      Double_t eta = fabs(p.eta());

      Int_t etabin = 0;
      if(eta > 0.5) etabin = 1;
      if(eta > 1.) etabin = 2;
      if(eta < 2.){
	hev_.ptav[etabin] += p.pt();
	++(hev_.n[etabin]);
      }
      ++(hev_.mult);
    }
    if(doHI_){
      edm::Handle<GenHIEvent> higen;
      iEvent.getByLabel(genHIsrc_,higen);
    }
  }

  if(doVertex_){
    edm::Handle<edm::SimVertexContainer> simVertices;
    iEvent.getByType<edm::SimVertexContainer>(simVertices);

    if (! simVertices.isValid() ) throw cms::Exception("FatalError") << "No vertices found\n";
    Int_t inum = 0;

    edm::SimVertexContainer::const_iterator it=simVertices->begin();
    if(it != simVertices->end()){
       SimVertex vertex = (*it);
       cout<<" Vertex position "<< inum <<" " << vertex.position().rho()<<" "<<vertex.position().z()<<endl;
       vx = vertex.position().x();
       vy = vertex.position().y();
       vz = vertex.position().z();
       vr = vertex.position().rho();
    }
  }

  for(Int_t i = 0; i<3; ++i){
    hev_.ptav[i] = hev_.ptav[i]/hev_.n[i];
  }

  hev_.b = b;
  hev_.scale = scale;
  hev_.npart = npart;
  hev_.ncoll = ncoll;
  hev_.nhard = nhard;
  hev_.phi0 = phi0;
  hev_.vx = vx;
  hev_.vy = vy;
  hev_.vz = vz;
  hev_.vr = vr;

  nt->Fill(nmix,np,src,sig);

  hydjetTree_->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
  void
HiGenAnalyzer::beginRun(const edm::Run&, const edm::EventSetup& iSetup) 
{
  iSetup.getData(pdt);
}

  void 
HiGenAnalyzer::beginJob()
{

  if(printLists_){
    out_b.open(fBFileName.c_str());
    if(out_b.good() == false)
      throw cms::Exception("BadFile") << "Can\'t open file " << fBFileName;
    out_n.open(fNFileName.c_str());
    if(out_n.good() == false)
      throw cms::Exception("BadFile") << "Can\'t open file " << fNFileName;
    out_m.open(fMFileName.c_str());
    if(out_m.good() == false)
      throw cms::Exception("BadFile") << "Can\'t open file " << fMFileName;
  }   

  if(doAnalysis_){
    nt = f->make<TNtuple>("nt","Mixing Analysis","mix:np:src:sig");

    hydjetTree_ = f->make<TTree>("hi","Tree of Hi gen Event");
    hydjetTree_->Branch("event",&hev_.event,"event/I");
    hydjetTree_->Branch("b",&hev_.b,"b/F");
    hydjetTree_->Branch("npart",&hev_.npart,"npart/F");
    hydjetTree_->Branch("ncoll",&hev_.ncoll,"ncoll/F");
    hydjetTree_->Branch("nhard",&hev_.nhard,"nhard/F");
    hydjetTree_->Branch("phi0",&hev_.phi0,"phi0/F");
    hydjetTree_->Branch("scale",&hev_.scale,"scale/F");

    hydjetTree_->Branch("n",hev_.n,"n[3]/I");
    hydjetTree_->Branch("ptav",hev_.ptav,"ptav[3]/F");

    if(doParticles_){

      hydjetTree_->Branch("mult",&hev_.mult,"mult/I");
      hydjetTree_->Branch("pt",hev_.pt,"pt[mult]/F");
      hydjetTree_->Branch("eta",hev_.eta,"eta[mult]/F");
      hydjetTree_->Branch("phi",hev_.phi,"phi[mult]/F");
      hydjetTree_->Branch("pdg",hev_.pdg,"pdg[mult]/I");
      hydjetTree_->Branch("chg",hev_.chg,"chg[mult]/I");
      hydjetTree_->Branch("sube",hev_.sube,"sube[mult]/I");

      hydjetTree_->Branch("vx",&hev_.vx,"vx/F");
      hydjetTree_->Branch("vy",&hev_.vy,"vy/F");
      hydjetTree_->Branch("vz",&hev_.vz,"vz/F");
      hydjetTree_->Branch("vr",&hev_.vr,"vr/F");
    }
  }

}

// ------------ method called once each job just after ending the event loop  ------------
void 
HiGenAnalyzer::endJob() {
}

//define this as a plug-in
DEFINE_FWK_MODULE(HiGenAnalyzer);
