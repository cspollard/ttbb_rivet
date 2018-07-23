#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Tools/Cutflow.hh"

namespace Rivet {


  class ATLAS_XX : public Analysis {
  public:

    DEFAULT_RIVET_ANA_CONSTRUCTOR(ATLAS_XX);

    void init() {
      // Eta ranges
      Cut eta_full = (Cuts::abseta < 5.0);
      //      Cut eta_full = Cuts::abseta < 5.0 && Cuts::pT >= 1.0*MeV;
      // Lepton cuts
      Cut lep_cuts25 = (Cuts::abseta < 2.5) && (Cuts::pT >= 25*GeV);
      // All final state particles
      FinalState fs(eta_full);

      // Get photons to dress leptons
      IdentifiedFinalState pho_id(fs);
      pho_id.acceptIdPair(PID::PHOTON);

      PromptFinalState photons(pho_id);
      photons.acceptTauDecays(true);
      // photons.acceptMuonDecays(true);

      // Projection to find the electrons
      IdentifiedFinalState el_id(fs);
      el_id.acceptIdPair(PID::ELECTRON);

      PromptFinalState electrons(el_id);
      electrons.acceptTauDecays(true);

      //      DressedLeptons dressedelectrons25(photons, electrons, 0.1, lep_cuts25, true, true);
      DressedLeptons dressedelectrons25(photons, electrons, 0.1, lep_cuts25, true, true);
      addProjection(dressedelectrons25, "DressedElectrons25");

      // Projection to find the muons
      IdentifiedFinalState mu_id(fs);
      mu_id.acceptIdPair(PID::MUON);
      PromptFinalState muons(mu_id);
      muons.acceptTauDecays(true);

      DressedLeptons dressedmuons25(photons, muons, 0.1, lep_cuts25, true, true);
      //      DressedLeptons dressedmuons25(photons, muons, 0.1, lep_cuts25, true, true);
      addProjection(dressedmuons25, "DressedMuons25");

      // From here on we are just setting up the jet clustering

      IdentifiedFinalState nu_id;
      nu_id.acceptNeutrinos();
      PromptFinalState neutrinos(nu_id);
      neutrinos.acceptTauDecays(true);

      VetoedFinalState jet_photons(photons);
      jet_photons.addDecayProductsVeto(15);
      jet_photons.addDecayProductsVeto(-15);

      VetoedFinalState jet_electrons(electrons);
      jet_electrons.addDecayProductsVeto(22);
      jet_electrons.addDecayProductsVeto(111);
      jet_electrons.addDecayProductsVeto(-111);

      DressedLeptons all_dressed_electrons(jet_photons, jet_electrons, 0.1, eta_full, true);
      DressedLeptons all_dressed_muons(jet_photons, muons, 0.1, eta_full, true);

      VetoedFinalState vfs(fs);
      vfs.addVetoOnThisFinalState(all_dressed_electrons);
      vfs.addVetoOnThisFinalState(all_dressed_muons);
      vfs.addVetoOnThisFinalState(neutrinos);
      addProjection(vfs, "vfs");

      FastJets jets(vfs, FastJets::ANTIKT, 0.4);
      jets.useInvisibles(JetAlg::DECAY_INVISIBLES);
      jets.useMuons(JetAlg::DECAY_MUONS);
      declare(jets, "jets");

      bin_ptlead_SL= { 30.0, 90.0, 120.0, 150.0, 185.0, 2025.0 };
      bin_ptsublead_SL= { 25.0, 65.0, 85.0, 120.0, 1510.0 };
      bin_pt3rdjet_SL = { 25.0, 45.0, 65.0, 90.0, 670.0 };
      bin_pt4thjet_SL = { 25.0, 30.0, 45.0, 480.0 };
      bin_hthad_SL = { 185.0, 390.0, 470.0, 575.0, 710.0, 5925.0 };
      bin_ht_SL = { 230.0, 450.0, 535.0, 640.0, 765.0, 6000.0 };
      bin_htb_SL = { 110.0, 240.0, 310.0, 400.0, 3665.0 };
      bin_mbb_SL = { 20.0, 125.0, 175.0, 230.0, 295.0, 390.0, 4410.0 };
      bin_ptbb_SL = { 0.0, 50.0, 85.0, 120.0, 160.0, 220.0, 1720.0 };
      bin_drbb_SL = { .35, 1.25, 1.85, 2.25, 2.55, 2.80, 3.05, 5.8 };
      bin_mbbclosest_SL ={ 10.0, 40.0, 60.0, 85.0, 120.0, 1460.0 };
      bin_ptbbclosest_SL = { 15.0, 90.0, 120.0, 160.0, 210.0, 2210.0 };
      bin_drbbclosest_SL = { .35, .60, .75, .90, 1.10, 1.35, 3.15 };
      bin_ptlightjet_SL = { 25, 60, 80, 105, 140, 195, 2190 };

      bin_ptlead_DIL = { 25, 65, 100, 150, 200, 1000 };
      bin_ptsublead_DIL = { 25, 50, 90, 150, 1000 };
      bin_pt3rdjet_DIL = { 25, 50, 80, 1000 };
      bin_pt4thjet_DIL = { 25, 30, 35, 1000 };
      bin_hthad_DIL = { 100, 225, 350, 600, 1400 };
      bin_ht_DIL = { 100, 270, 350, 450, 700, 1400 };
      bin_mbb_DIL = { 20, 50, 100, 200, 400, 2000 };
      bin_ptbb_DIL = { 10, 65, 100, 140, 1000 };
      bin_drbb_DIL = { .4, 1.3, 1.8, 2.3, 2.8, 8.0 };
      bin_mbbclosest_DIL = { 10, 55, 100, 145, 250, 1000 };
      bin_ptbbclosest_DIL = { 10, 65, 100, 150, 250, 1000 };
      bin_drbbclosest_DIL = { .4, .7, 1.0, 1.3, 3.0 };


      _hmbb_SL                       = bookHisto1D( "mbb_SL", bin_mbb_SL);
      _hpt_leadingBJet_SL            = bookHisto1D("leadBjetpt_SL", bin_ptlead_SL);                                                                                                       
      _heta_leadingBJet_SL           = bookHisto1D( "leadBjeteta_SL", 5 , 0, 2.5);                                                                                                                      
      _hpt_2ndleadingBJet_SL         = bookHisto1D( "secondleadingBJetpt_SL", bin_ptsublead_SL);                                                                                                        
      _heta_2ndleadingBJet_SL        = bookHisto1D( "secondleadingBJeteta_SL",5 ,0 ,2.5 );                                                                                                       
      _hnjets_SL                     = bookHisto1D( "njets_SL",6 ,2 ,8);                                                                                                                     
      _hnbjets_SL                    = bookHisto1D( "nbjets_SL",6 ,2 ,8);                                                                                                               
      _hptbb_SL                      = bookHisto1D( "ptbb_SL", bin_ptbb_SL);                                                                                                                   
      _hmbb_close_SL                 = bookHisto1D( "mbbclose_SL", bin_mbbclosest_SL);                                                                                                        
      _hptbb_close_SL                = bookHisto1D( "ptbbclose_SL", bin_ptbbclosest_SL);                                                                                                      
      _hdrbb_close_SL                = bookHisto1D( "drbbclose_SL", bin_drbbclosest_SL);                                                                                                            
      _hdrbb_SL                      = bookHisto1D( "drbb_SL", bin_drbb_SL);                                                                                                                  
      _hpt_3rdleadingBJet_SL         = bookHisto1D( "thirdleadingBJetpt_SL", bin_pt3rdjet_SL);                                                                                                      
      _heta_3rdleadingBJet_SL        = bookHisto1D( "thirdleadingBJeteta_SL",5 ,0, 2.5);                                                                                                               
      _hpt_4thleadingBJet_SL         = bookHisto1D( "fourthleadingBJetpt_SL", bin_pt4thjet_SL);
      _hHT_SL                        = bookHisto1D( "HT_SL", bin_ht_SL);       
      _hHTHAD_SL                     = bookHisto1D( "HTHAD_SL", bin_hthad_SL);      
      _hHTB_SL                       = bookHisto1D( "HTB_SL", bin_htb_SL);
      _hpt_leadingLightJet_SL        = bookHisto1D( "leadingLightJetpt_SL", bin_ptlightjet_SL);      


      _hmbb_DIL                      = bookHisto1D( "mbb_DIL", bin_mbb_DIL);
      _hpt_leadingBJet_DIL           = bookHisto1D("leadBjetpt_DIL", bin_ptlead_DIL);
      _heta_leadingBJet_DIL          = bookHisto1D( "leadBjeteta_DIL", 5 ,0 ,2.5);
      _hpt_2ndleadingBJet_DIL        = bookHisto1D( "secondleadingBJetpt_DIL", bin_ptsublead_DIL);
      _heta_2ndleadingBJet_DIL       = bookHisto1D( "secondleadingBJeteta_DIL", 5 ,0 ,2.5 );
      _hnjets_DIL                    = bookHisto1D( "njets_DIL",6 ,2 ,8);
      _hnbjets_DIL                   = bookHisto1D( "nbjets_DIL",6 ,2 ,8);
      _hptbb_DIL                     = bookHisto1D( "ptbb_DIL", bin_ptbb_DIL);
      _hmbb_close_DIL                = bookHisto1D( "mbbclose_DIL", bin_mbbclosest_DIL);
      _hptbb_close_DIL               = bookHisto1D( "ptbbclose_DIL", bin_ptbbclosest_DIL);
      _hdrbb_close_DIL               = bookHisto1D( "drbbclose_DIL", bin_drbbclosest_DIL);
      _hdrbb_DIL                     = bookHisto1D( "drbb_DIL", bin_drbb_DIL);
      _hpt_3rdleadingBJet_DIL        = bookHisto1D( "thirdleadingBJetpt_DIL", bin_pt3rdjet_DIL);
      _heta_3rdleadingBJet_DIL       = bookHisto1D( "thirdleadingBJeteta_DIL",5 ,0, 2.5);
      _hpt_4thleadingBJet_DIL        = bookHisto1D( "fourthleadingBJetpt_DIL", bin_pt4thjet_DIL);
      _hHT_DIL                       = bookHisto1D( "HT_DIL", bin_ht_DIL);
      _hHTHAD_DIL                    = bookHisto1D( "HTHAD_DIL", bin_hthad_DIL);
      iEvent = 0;
    }


    void analyze(const Event& event) {
      iEvent++;
      std::cout << iEvent << std::endl;
      vector<DressedLepton> electrons = applyProjection<DressedLeptons>(event, "DressedElectrons25").dressedLeptons();
      vector<DressedLepton> muons     = applyProjection<DressedLeptons>(event, "DressedMuons25").dressedLeptons();
      weight = event.weight();
      vector<DressedLepton> kept_electrons25; vector<DressedLepton> kept_electrons27;
      vector<DressedLepton> kept_muons25; vector<DressedLepton> kept_muons27;
      const Jets jets = applyProjection<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
      variableClear();

      //overlap removal muons, then electrons
      doOR(jets,electrons,kept_electrons25,kept_electrons27);
      doOR(jets,muons,kept_muons25,kept_muons27);
      fillJets(jets);
      // Evaluate basic event selection
	

      bool no25 = kept_muons25.size() > 0 || kept_electrons25.size() > 0;
      bool passSL = ((kept_electrons27.size() == 1 && kept_muons27.size() == 0) || (kept_electrons27.size() == 0 && kept_muons27.size() == 1)) && !no25;
      bool pass_emu = kept_electrons27.size() == 1 && kept_muons27.size() == 1 && !no25;	
      //      if(njets >= 2 && nbjets >= 2 && pass_emu)      std::cout << njets << " "<< nbjets << " " << kept_muons27.size() << " " <<kept_electrons27.size() <<" " << kept_electrons27[0].pT()*1000 <<" "<< kept_electrons27[0].eta()  << " " <<kept_muons27[0].pT()*1000 <<" "<< kept_muons27[0].eta() << endl;
      //   else cout << " " << endl;

      DressedLepton *lepton0 = kept_muons27.empty() ? &kept_electrons27[0] : &kept_muons27[0];
      DressedLepton *lepton1 = NULL;
      if(!pass_emu && !passSL) vetoEvent;
      


      if(pass_emu) {
	//	lepton0 = &kept_muons27[0];
	lepton1 = &kept_electrons27[0];
	if(lepton1->charge() == lepton0->charge()) vetoEvent;
      }

      if(nbjets >= 2 && passSL){
	_hnjets_SL->fill(njets, weight);
        _hnbjets_SL->fill(nbjets, weight);
      }
      if(nbjets >= 2 && pass_emu) {
	_hnjets_DIL->fill(njets, weight);
        _hnbjets_DIL->fill(nbjets, weight);
      }

      if(passSL && (nbjets < 4 || njets < 6)) vetoEvent;
      if(pass_emu && (nbjets < 3 || njets < 3)) vetoEvent;

      getMindr();

      ht = 0; hthad = 0; htb = 0;
      ht += lepton0->pT();
      if(pass_emu) ht += lepton1->pT();
      hthad = sum(jets,pT,0.0);
      htb = sum(bjets,pT,0.0);
      ht += hthad;
     
      bjet0 = bjets[0].momentum();
      bjet1 = bjets[1].momentum();
      bjet2 = bjets[2].momentum();
      jsum = bjet0 + bjet1;

      if(nbjets >= 4) bjet3 = bjets[3].momentum();
      if(passSL) fillSLHists();
      if(pass_emu) fillEMUHists();
    }


      void finalize() {
	normalize({_hpt_leadingBJet_SL, _hpt_2ndleadingBJet_SL, _hmbb_SL, _hptbb_SL,
	      _hpt_leadingBJet_DIL, _hpt_2ndleadingBJet_DIL, _hmbb_DIL, _hptbb_DIL,
	      _hptbb_close_SL, _hmbb_close_SL, _hdrbb_close_SL, _hptbb_close_DIL, _hmbb_close_DIL, _hdrbb_close_DIL,
	      _hdrbb_SL, _hdrbb_DIL, _heta_3rdleadingBJet_SL , _heta_3rdleadingBJet_DIL , _hpt_3rdleadingBJet_SL , _hpt_3rdleadingBJet_DIL,
	      _heta_leadingBJet_SL , _heta_2ndleadingBJet_SL, _heta_leadingBJet_DIL, _heta_2ndleadingBJet_DIL,
	      _hHT_DIL, _hHT_SL, _hHTHAD_DIL, _hHTHAD_SL, _hpt_4thleadingBJet_SL, _hpt_4thleadingBJet_DIL,
	      _hpt_leadingLightJet_SL, _hHTB_SL},1);
	      
      }


   private:
    Histo1DPtr _hpt_leadingBJet_SL, _hpt_2ndleadingBJet_SL, _hmbb_SL, _hptbb_SL, _hnjets_SL, _hnbjets_SL;
    Histo1DPtr _hpt_leadingBJet_DIL, _hpt_2ndleadingBJet_DIL, _hmbb_DIL, _hptbb_DIL, _hnjets_DIL, _hnbjets_DIL;
    Histo1DPtr _hptbb_close_SL, _hmbb_close_SL, _hdrbb_close_SL, _hptbb_close_DIL, _hmbb_close_DIL, _hdrbb_close_DIL;
    Histo1DPtr _hdrbb_SL, _hdrbb_DIL, _heta_3rdleadingBJet_SL , _heta_3rdleadingBJet_DIL , _hpt_3rdleadingBJet_SL , _hpt_3rdleadingBJet_DIL ;
    Histo1DPtr _heta_leadingBJet_SL , _heta_2ndleadingBJet_SL, _heta_leadingBJet_DIL, _heta_2ndleadingBJet_DIL;
    Histo1DPtr _hHT_DIL, _hHT_SL, _hHTHAD_DIL, _hHTHAD_SL, _hpt_4thleadingBJet_SL, _hpt_4thleadingBJet_DIL;
    Histo1DPtr _hpt_leadingLightJet_SL, _hHTB_SL;

    int iEvent;
    int nbjets;
    int njets;
    double mindrOld;
    double mbbclose;
    double ptbbclose;
    double ht; double hthad; double htb;
    Jets lightjets;
    Jets bjets;

    FourMomentum bjet0;
    FourMomentum bjet1;
    FourMomentum bjet2;
    FourMomentum jsum;
    FourMomentum bjet3;
    double weight;

    std::vector<double> bin_ptlead_SL;
    std::vector<double> bin_ptsublead_SL;
    std::vector<double> bin_pt3rdjet_SL;
    std::vector<double> bin_pt4thjet_SL;
    std::vector<double> bin_hthad_SL;
    std::vector<double> bin_ht_SL;
    std::vector<double> bin_mbb_SL;
    std::vector<double> bin_ptbb_SL;
    std::vector<double> bin_drbb_SL;
    std::vector<double> bin_mbbclosest_SL;
    std::vector<double> bin_ptbbclosest_SL;
    std::vector<double> bin_drbbclosest_SL;
    std::vector<double> bin_ptlightjet_SL;
    std::vector<double> bin_htb_SL;

    std::vector<double> bin_ptlead_DIL;
    std::vector<double> bin_ptsublead_DIL;
    std::vector<double> bin_pt3rdjet_DIL;
    std::vector<double> bin_pt4thjet_DIL;
    std::vector<double> bin_hthad_DIL;
    std::vector<double> bin_ht_DIL;
    std::vector<double> bin_mbb_DIL;
    std::vector<double> bin_ptbb_DIL;
    std::vector<double> bin_drbb_DIL;
    std::vector<double> bin_mbbclosest_DIL;
    std::vector<double> bin_ptbbclosest_DIL;
    std::vector<double> bin_drbbclosest_DIL;
    
    void variableClear() {
      lightjets.clear();
      bjets.clear();
      mindrOld = 9999;
      mbbclose = -1111;
      ptbbclose = -1111;
      nbjets = 0;
      njets = 0;
    }

    void getMindr()
    {
      double mindr = 999;
      for (size_t i = 0; i < bjets.size(); i++) {
        for (size_t j = 0; j < bjets.size(); j++) {
          if(i == j) continue;
	  mindr = deltaR(bjets[i], bjets[j]);
          if(mindr < mindrOld){
            FourMomentum closefour = bjets[i].momentum() + bjets[j].momentum();
            mbbclose = closefour.mass()*GeV;
            ptbbclose = closefour.pT()*GeV;
            mindrOld = mindr;
          }
        }
      }
    }
    
    void doOR(const Jets jets,vector<DressedLepton> & initial_Lepton, vector<DressedLepton> & saved_Lepton, vector<DressedLepton> & saved_Lepton2){
      foreach (DressedLepton& leptons,initial_Lepton) {
	bool dRPass = true;
	foreach (const Jet& jet, jets) {
	  double deltaRjets = deltaR(leptons.momentum(), jet.momentum());
	  //	  std::cout << deltaRjets << std::endl;
	  if (deltaRjets < 0.4) dRPass = false;
	}
	if(dRPass && leptons.pT()*GeV >= 25 && leptons.pT()*GeV < 27) saved_Lepton.push_back(leptons);
	if(dRPass && leptons.pT()*GeV >= 27) saved_Lepton2.push_back(leptons);
      }
    }
    
    void fillSLHists(){
      double drbb = deltaR(bjet0,bjet1);
      if(lightjets.size() > 0) _hpt_leadingLightJet_SL->fill(lightjets[0].pT()*GeV,weight);
      _hHT_SL->fill(ht*GeV,weight);
      _hHTHAD_SL->fill(hthad*GeV,weight);
      _hHTB_SL->fill(htb*GeV,weight);
      _hmbb_close_SL->fill(mbbclose*GeV,weight);
      _hptbb_close_SL->fill(ptbbclose*GeV,weight);
      _hdrbb_close_SL->fill(mindrOld,weight);
      _hptbb_SL->fill(jsum.pT()*GeV,weight);
      _hmbb_SL->fill(jsum.mass()*GeV,weight);
      _hdrbb_SL->fill(drbb,weight);
      _hpt_leadingBJet_SL->fill(bjet0.pT()*GeV,weight);
      _hpt_2ndleadingBJet_SL->fill(bjet1.pT()*GeV,weight);
      _hpt_3rdleadingBJet_SL->fill(bjet2.pT()*GeV,weight);
      _heta_leadingBJet_SL->fill(bjet0.abseta(), weight);
      _heta_2ndleadingBJet_SL->fill(bjet1.abseta(), weight);
      _heta_3rdleadingBJet_SL->fill(bjet2.abseta(), weight);
      if(nbjets >= 4) _hpt_4thleadingBJet_SL->fill(bjet3.pT()*GeV,weight);
    }
    void fillEMUHists(){
      double drbb = deltaR(bjet0,bjet1);
      _hHT_DIL->fill(ht*GeV,weight);
      _hHTHAD_DIL->fill(hthad*GeV,weight);
      _hptbb_DIL->fill(jsum.pT()*GeV,weight);
      _hmbb_DIL->fill(jsum.mass()*GeV,weight);
      _hdrbb_DIL->fill(drbb,weight);
      _hpt_leadingBJet_DIL->fill(bjet0.pT()*GeV,weight);
      _hpt_2ndleadingBJet_DIL->fill(bjet1.pT()*GeV,weight);
      _hpt_3rdleadingBJet_DIL->fill(bjet2.pT()*GeV,weight);
      _heta_leadingBJet_DIL->fill(bjet0.abseta(), weight);
      _heta_2ndleadingBJet_DIL->fill(bjet1.abseta(), weight);
      _heta_3rdleadingBJet_DIL->fill(bjet2.abseta(), weight);
      _hptbb_close_DIL->fill(ptbbclose*GeV,weight);
      _hmbb_close_DIL->fill(mbbclose*GeV,weight);
      _hdrbb_close_DIL->fill(mindrOld,weight);
      if(nbjets >= 4) _hpt_4thleadingBJet_DIL->fill(bjet3.pT()*GeV,weight);
    }    
    
    void fillJets(const Jets& jets){
      foreach (const Jet& jet, jets){
        njets++;
        bool isBtagged = jet.bTagged(Cuts::pT >= 5.0*GeV);
        if(!isBtagged) lightjets += jet;
        if(isBtagged){
          bjets += jet;
          nbjets++;
	}
      }
    }
    
    
  };
  

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_XX);

}
