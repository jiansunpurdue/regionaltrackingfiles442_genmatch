#include <iostream>
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"
#include "DataFormats/TrackReco/interface/Track.h"
// essentials !!!
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ServiceRegistry/interface/Service.h" 
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1.h"
#include "TTree.h"
#include "TObject.h"
#include "UserCode/OpenHF/container/daughterparticle.h"
#include "UserCode/OpenHF/container/motherquark.h"
#include "UserCode/OpenHF/container/middleparticle.h"
#include "UserCode/OpenHF/container/dmeson.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
using namespace edm;
using namespace std;
using namespace reco;



class D0_ccbar : public edm::EDAnalyzer
{

   public:
   
      //
      explicit D0_ccbar( const edm::ParameterSet& ) ;
      virtual ~D0_ccbar() {} // no need to delete ROOT stuff
                                   // as it'll be deleted upon closing TFile
      
      virtual void analyze( const edm::Event&, const edm::EventSetup& ) ;
      virtual void beginJob() ;
      virtual void endRun( const edm::Run&, const edm::EventSetup& ) ;
      virtual void endJob() ;
      void daughterfill(daughterparticle * daughter,const Candidate * candidate);
      void motherfill(motherquark * mother, const Candidate * candidate);
      void middlefill(middleparticle * middle, const Candidate * candidate);
      void dmesonfill(dmeson * d, const Candidate * candidate);    

   private:
     
     TH1D*       d0pt ;
     TTree *     d0tree;
     double      pt_dmeson;
     double      eta_paifromd0;
//     daughterparticle *paionfromd0;
     dmeson *    dmeson_d0;
     daughterparticle * particle1fromdmeson;
     daughterparticle * particle2fromdmeson;
     motherquark *      dmesonmotherquark;
     
}; 


D0_ccbar::D0_ccbar( const ParameterSet& pset )
{
    particle1fromdmeson = new daughterparticle();
    particle2fromdmeson = new daughterparticle();
    dmesonmotherquark =  new motherquark();
    dmeson_d0 = new dmeson();
// actually, pset is NOT in use - we keep it here just for illustratory putposes

//  d0tree = new TTree("d0tree", "d0tree");
//  pai = new daughterparticle();
//  d0tree->Branch("pai","daughterparticle",&pai,32000,1);

}

void D0_ccbar::beginJob()
{
  Service<TFileService> fs; 
  d0tree = fs->make<TTree>("d0tree","d0tree");
//  d0tree->Branch("pai","daughterparticle",&paionfromd0,32000,0);
  d0tree->Branch("d0","dmeson",&dmeson_d0,4000,99);  
  
  return ;
  
}

void D0_ccbar::daughterfill(daughterparticle *daughter, const Candidate * candidate)
{  
       daughter->pt = candidate->pt();
       daughter->charge = candidate->charge(); 
       daughter->motherid = (candidate->mother())->pdgId(); 
       daughter->pid = candidate->pdgId(); 
       daughter->phi = candidate->phi(); 
       daughter->eta = candidate->eta(); 
       daughter->mass = candidate->mass();
  return;
}

void D0_ccbar::motherfill(motherquark *mother, const Candidate * candidate)
{
       mother->pt = candidate->pt();
       mother->charge = candidate->charge();
       mother->pid = candidate->pdgId();
       mother->phi = candidate->phi();
       mother->eta = candidate->eta();
       mother->mass = candidate->mass();
  return;
}

void D0_ccbar::middlefill(middleparticle *middle, const Candidate * candidate)
{
       middle->pt = candidate->pt();
       middle->charge = candidate->charge();
       middle->motherid = (candidate->mother())->pdgId();
       middle->pid = candidate->pdgId();
       middle->phi = candidate->phi();
       middle->eta = candidate->eta();
       middle->mass = candidate->mass();
       middle->numberofdaughter = candidate->numberOfDaughters();
  return;
}

void D0_ccbar::dmesonfill(dmeson *d, const Candidate * candidate)
{
       d->pt = candidate->pt();
       d->charge = candidate->charge();
       d->motherid = (candidate->mother())->pdgId();
       d->pid = candidate->pdgId();
       d->phi = candidate->phi();
       d->eta = candidate->eta();
       d->mass = candidate->mass();
       d->numberofdaughter = candidate->numberOfDaughters();
  return;
}


void D0_ccbar::analyze( const Event& e, const EventSetup& )
{
  
  // here's an example of accessing GenEventInfoProduct
  Handle< GenEventInfoProduct > GenInfoHandle;
  e.getByLabel( "generator", GenInfoHandle );
  double qScale = GenInfoHandle->qScale();
  double pthat = ( GenInfoHandle->hasBinningValues() ? 
                  (GenInfoHandle->binningValues())[0] : 0.0);
//  cout << " qScale = " << qScale << " pthat = " << pthat << endl;
  //
  // this (commented out) code below just exemplifies how to access certain info 
  //
  //double evt_weight1 = GenInfoHandle->weights()[0]; // this is "stanrd Py6 evt weight;
                                                    // corresponds to PYINT1/VINT(97)
  //double evt_weight2 = GenInfoHandle->weights()[1]; // in case you run in CSA mode or otherwise
                                                    // use PYEVWT routine, this will be weight
						    // as returned by PYEVWT, i.e. PYINT1/VINT(99)
  //std::cout << " evt_weight1 = " << evt_weight1 << std::endl;
  //std::cout << " evt_weight2 = " << evt_weight2 << std::endl;
  //double weight = GenInfoHandle->weight();
  //std::cout << " as returned by the weight() method, integrated event weight = " << weight << std::endl;
  
  // here's an example of accessing particles in the event record (HepMCProduct)
  //
  Handle< HepMCProduct > EvtHandle ;
  
  // find initial (unsmeared, unfiltered,...) HepMCProduct
  //
  e.getByLabel( "generator", EvtHandle ) ;
  
  const HepMC::GenEvent* Evt = EvtHandle->GetEvent() ;
  
  // this a pointer - and not an array/vector/... 
  // because this example explicitely assumes
  // that there one and only Higgs in the record
  //
  HepMC::GenVertex* D0DecVtx = 0 ;

  int numberofd0 = 0;
  int numberofd0tokaipai = 0;
 
  
  Handle<GenParticleCollection> genParticles;
//  e.getByLabel("genParticles", genParticles);
  e.getByLabel("hiGenParticles", genParticles);




  for(size_t i = 0; i < genParticles->size(); ++ i) {
     const GenParticle & p = (*genParticles)[i];
     int id = p.pdgId();
     int st = p.status();  
     if ( id != 421 || abs(st) != 2)   continue;
     numberofd0++;
     const Candidate * mom1 = p.mother();
     const Candidate * mom = mom1->mother();
//     cout << "pid of mother's mother = " << (mom->mother())->pdgId() << endl;
//     cout << "pid of mother's mother number = " << mom1->numberOfMothers() << endl;
//     const Candidate * ancester = p.();
//     double pt = p.pt(), eta = p.eta(), phi = p.phi(), mass = p.mass();
//     pt_dmeson = pt;
//     double vx = p.vx(), vy = p.vy(), vz = p.vz();
//     int charge = p.charge();
     int n = p.numberOfDaughters();
     bool target = false;
     if ( n != 2) continue;
     const Candidate * d0 =  &p;
     const Candidate * paionfromd0 = NULL;
     const Candidate * kaionfromd0 = NULL;
  //    cout << "D0 to pai d0" << endl;
  //    cout << "D0 daughter number" << d0->numberOfDaughters()<<endl;
     if (!((abs(d0->daughter(0)->pdgId())==321&&abs(d0->daughter(1)->pdgId())==211)||(abs(d0->daughter(0)->pdgId())==211&&abs(d0->daughter(1)->pdgId())==321)))
           continue;
     if(d0->daughter(0)->charge() * d0->daughter(1)->charge() > 0)   continue;
     if(abs(d0->daughter(0)->pdgId())==321&&abs(d0->daughter(1)->pdgId())==211)   {  paionfromd0 = d0->daughter(1); kaionfromd0 = d0->daughter(0) ;}
     if(abs(d0->daughter(0)->pdgId())==211&&abs(d0->daughter(1)->pdgId())==321)   {  paionfromd0 = d0->daughter(0); kaionfromd0 = d0->daughter(1) ;}
     cout << "d0 information pt =" << d0->pt() << " eta="<< d0->eta()  <<endl;//<< " x=" << d0->Point().x()<< " y=" << d0->Point().y() <<endl;
 //    cout << "D0 to kai pai" << endl;
     target =  true;
  
   
 //    if(target){
     daughterfill(particle1fromdmeson,paionfromd0);
     motherfill(dmesonmotherquark,mom);
     daughterfill(particle2fromdmeson,kaionfromd0);
     dmesonfill (dmeson_d0,d0);
     (dmeson_d0->momquark).push_back(*dmesonmotherquark);
     (dmeson_d0->daughters).push_back(*particle1fromdmeson);
     (dmeson_d0->daughters).push_back(*particle2fromdmeson);
  //   paionfromd0->CandidatePtr();
     d0tree->Fill();
//    }
//   particle1fromdmeson->(simtks);   
    
   dmeson_d0->clear();
   particle1fromdmeson->clear();
   particle2fromdmeson->clear();
   dmesonmotherquark->clear();

   }
 

   cout << "number of D0 =====" << numberofd0 << endl;
   
//   dmeson_d0->clear();
//   particle1fromdmeson->clear();
//   particle2fromdmeson->clear();
//   dmesonmotherquark->clear();
// 
   return ;
   
}

void D0_ccbar::endRun( const edm::Run& r, const edm::EventSetup& )
{

   return;

}


void D0_ccbar::endJob()
{
   
   return ;
}
 
typedef D0_ccbar D0_ccbar_Analyzer;
DEFINE_FWK_MODULE(D0_ccbar_Analyzer);
