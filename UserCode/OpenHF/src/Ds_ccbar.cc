#include <iostream>
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

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



class Ds_ccbar : public edm::EDAnalyzer
{

   public:
   
      //
      explicit Ds_ccbar( const edm::ParameterSet& ) ;
      virtual ~Ds_ccbar() {} // no need to delete ROOT stuff
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
     
     TH1D*       dspt ;
     TTree *     dstree;
     double      pt_dmeson;
     double      eta_paifromds;
//     daughterparticle *paionfromds;
     dmeson *    dmeson_ds;
     daughterparticle * particlefromdmeson;
     daughterparticle * particle1frommiddle;
     daughterparticle * particle2frommiddle;
     motherquark *      dmesonmotherquark;
     middleparticle *   particledecayed;
     
}; 


Ds_ccbar::Ds_ccbar( const ParameterSet& pset )
{
    particlefromdmeson = new daughterparticle();
    particle1frommiddle = new daughterparticle();
    particle2frommiddle = new daughterparticle();
    dmesonmotherquark =  new motherquark();
    particledecayed = new middleparticle();
    dmeson_ds = new dmeson();
// actually, pset is NOT in use - we keep it here just for illustratory putposes

//  dstree = new TTree("dstree", "dstree");
//  pai = new daughterparticle();
//  dstree->Branch("pai","daughterparticle",&pai,32000,1);

}

void Ds_ccbar::beginJob()
{
  Service<TFileService> fs; 
  dstree = fs->make<TTree>("dstree","dstree");
//  dstree->Branch("pai","daughterparticle",&paionfromds,32000,0);
  dstree->Branch("ds","dmeson",&dmeson_ds,32000,1);  
  
  return ;
  
}

void Ds_ccbar::daughterfill(daughterparticle *daughter, const Candidate * candidate)
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

void Ds_ccbar::motherfill(motherquark *mother, const Candidate * candidate)
{
       mother->pt = candidate->pt();
       mother->charge = candidate->charge();
       mother->pid = candidate->pdgId();
       mother->phi = candidate->phi();
       mother->eta = candidate->eta();
       mother->mass = candidate->mass();
  return;
}

void Ds_ccbar::middlefill(middleparticle *middle, const Candidate * candidate)
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

void Ds_ccbar::dmesonfill(dmeson *d, const Candidate * candidate)
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


void Ds_ccbar::analyze( const Event& e, const EventSetup& )
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

  int numberofds = 0;
  int numberofphitokaipai = 0;
 
  
  Handle<GenParticleCollection> genParticles;
  e.getByLabel("hiGenParticles", genParticles);
//  cout << "particle number = " << genParticles->size() << endl;


  for(size_t i = 0; i < genParticles->size(); ++ i) {
     const GenParticle & p = (*genParticles)[i];
     int id = p.pdgId();
     int st = p.status();  
     if (abs(id) != 431 || abs(st) != 2)   continue;
     numberofds++;
     const Candidate * mom1 = p.mother();
     const Candidate * mom = mom1->mother();
//     cout << "pid of mother's mother = " << (mom->mother())->pdgId() << endl;
//     const Candidate * ancester = p.();
//     double pt = p.pt(), eta = p.eta(), phi = p.phi(), mass = p.mass();
//     pt_dmeson = pt;
//     double vx = p.vx(), vy = p.vy(), vz = p.vz();
//     int charge = p.charge();
     int n = p.numberOfDaughters();
     if ( n != 2) continue;
     const Candidate * phi =  NULL;
     const Candidate * pai =  NULL;
     const Candidate * kaion1fromphi = NULL;
     const Candidate * kaion2fromphi = NULL;
     if (!((p.daughter(0)->pdgId() == 333 && abs( p.daughter(1)->pdgId()) == 211) ||( abs(p.daughter(0)->pdgId()) == 211 &&  p.daughter(1)->pdgId() == 333)))   
           continue;
     if ( p.daughter(0)->pdgId() == 333 && abs( p.daughter(1)->pdgId()) == 211 )   { phi = p.daughter(0); pai = p.daughter(1); }
     if ( abs(p.daughter(0)->pdgId()) == 211 &&  p.daughter(1)->pdgId() == 333 )   { phi = p.daughter(1); pai = p.daughter(0); }
  //    cout << "Ds to pai phi" << endl;
  //    cout << "D0 daughter number" << phi->numberOfDaughters()<<endl;
     if (!(phi->status() == 2 && phi->numberOfDaughters() == 2)) continue;
     if (!(abs(phi->daughter(0)->pdgId())==321&&abs(phi->daughter(1)->pdgId())==321))         continue;

     if(phi->daughter(0)->charge() * phi->daughter(1)->charge() > 0)   continue;
     if(abs(phi->daughter(0)->pdgId())==321&&abs(phi->daughter(1)->pdgId())==321)   {  kaion1fromphi = phi->daughter(1); kaion2fromphi = phi->daughter(0) ;}

     const Candidate * dds = &p;
 //    cout << "Ds to D0 pai to pai kai pai" << endl;
     
     daughterfill(particlefromdmeson,pai);
     motherfill(dmesonmotherquark,mom);
     middlefill(particledecayed,phi);
     daughterfill(particle1frommiddle,kaion1fromphi);
     daughterfill(particle2frommiddle,kaion2fromphi);
     dmesonfill (dmeson_ds,dds);
     (dmeson_ds->momquark).push_back(*dmesonmotherquark);
     (dmeson_ds->middle).push_back(*particledecayed);
     (dmeson_ds->daughters).push_back(*particlefromdmeson);
     (dmeson_ds->daughters).push_back(*particle1frommiddle);
     (dmeson_ds->daughters).push_back(*particle2frommiddle); 

     dstree->Fill();
     dmeson_ds->clear();
     dmesonmotherquark->clear();
     particledecayed->clear();
     particle1frommiddle->clear();
     particle2frommiddle->clear();
   }
 

//   cout << "number of Ds =====" << numberofds << endl;
  
   
   return ;
   
}

void Ds_ccbar::endRun( const edm::Run& r, const edm::EventSetup& )
{

   return;

}


void Ds_ccbar::endJob()
{
   
   return ;
}
 
typedef Ds_ccbar Ds_ccbar_Analyzer;
DEFINE_FWK_MODULE(Ds_ccbar_Analyzer);
