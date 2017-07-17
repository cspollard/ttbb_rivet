// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Analyses/MC_JetAnalysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/MissingMomentum.hh"



namespace Rivet {


  /// @brief Add a short analysis description here
  class MCTTB : public Analysis {
  public:

    /// Constructor
    MCTTB(const string& name="MCTTB", size_t n_ljet=4, size_t n_bjet=2, size_t n_Bjet=0, double ptcut=15.0, size_t bdef=2, size_t cparton=0)
      : Analysis(name), _h_pT_ljet(n_ljet), _h_pT_bjet(n_bjet),_h_pT_Bjet(n_Bjet),_h_Bjet_btags(n_Bjet)
    {  num_Bjets = n_Bjet;
       num_bjets = n_bjet;
       num_ljets = n_ljet;
       m_ptcut=ptcut;
       m_Bdef=bdef;
       m_parton=cparton;
       MSG_INFO("constructor, bool cparton:" << cparton);
    }


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {
      const double sqrts = sqrtS() ? sqrtS() : 8000.*GeV;

      ChargedLeptons lfs(FinalState(-5., 5., 30*GeV));
      //declare(lfs, "LFS");

      VetoedFinalState fs(FinalState(-5., 5., 0*GeV));
      fs.addVetoOnThisFinalState(lfs);
      declare(FastJets(fs, FastJets::ANTIKT, 0.4), "Jets");
      //declare(MissingMomentum(fs), "MissingET");

      const int nbins_pT = 100;


      for (size_t i = 0; i < num_ljets; ++i) {
          const string pTname = "light_jet_pT_" + to_str(i+1);
          const double pTmax = 1.0/(double(i)+3.0) * sqrts/GeV/2.0;
          _h_pT_ljet[i] = bookHisto1D(pTname, logspace(nbins_pT, m_ptcut, pTmax));

        }

      for (size_t i = 0; i < num_bjets; ++i) {
          const string pTname = "b_jet_pT_" + to_str(i+1);
          const double pTmax = 1.0/(double(i)+3.0) * sqrts/GeV/2.0;
          _h_pT_bjet[i] = bookHisto1D(pTname, logspace(nbins_pT, m_ptcut, pTmax));
        }

      for (size_t i = 0; i < num_Bjets; ++i) {
          const string pTname = "B_jet_pT_" + to_str(i+1);
          const double pTmax = 1.0/(double(i)+3.0) * sqrts/GeV/2.0;
          _h_Bjet_btags[i] = bookHisto1D("btags_Bjet"+to_str(i+1),10,0 -0.5 ,10 -0.5);
          _h_pT_Bjet[i] = bookHisto1D(pTname, logspace(nbins_pT, m_ptcut, pTmax));
        }

      _h_ljet_multi = bookHisto1D("ljet_multi", 9, num_ljets-0.5, num_ljets + 9 -0.5);
      _h_bjet_multi = bookHisto1D("bjet_multi", 5, num_bjets-0.5, num_bjets + 5 -0.5);
      _h_Bjet_multi = bookHisto1D("Bjet_multi", 3, num_Bjets-0.5, num_Bjets + 3 -0.5);


    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      double weight = event.weight();

      const FastJets& jetpro = apply<FastJets>(event, "Jets");
      const Jets alljets = jetpro.jetsByPt(Cuts::pT>m_ptcut*GeV);
      if (alljets.size() < 4) {
          MSG_DEBUG("Event failed jet multiplicity cut");
          vetoEvent;
        }

      Jets bjets,Bjets,ljets;
      MSG_INFO("all jets:" << alljets.size());
      MSG_INFO("m_parton:" << m_parton);

      foreach (Jet jet, alljets) {
          if (m_parton==0){  //hadron level
              MSG_INFO("in hadron level");
              if(jet.bTagged()){
                  Particles btags = jet.bTags();
                  if(btags.size() >= m_Bdef) {   //is this sufficient to have a fat bjet?
                      Bjets.push_back(jet);
                    }else{
                      bjets.push_back(jet);
                    }
                }else ljets.push_back(jet);

            }else{ //parton level
              MSG_INFO("in parton level");
              size_t nb(0);
              Particles jet_content= jet.constituents();
              foreach(Particle part, jet_content){
                  if (part.hasBottom()) nb++;
                }
              if (nb>0){
                  if (nb >= m_Bdef){
                      MSG_INFO("B jet");
                      Bjets.push_back(jet);
                    }else {
                      MSG_INFO("b jet");
                      bjets.push_back(jet);
                    }
                }else{
                  MSG_INFO("l jet");
                  ljets.push_back(jet);
                }
            }
        }
      if (ljets.size() < num_ljets) vetoEvent;
      if (bjets.size() < num_bjets) vetoEvent;
      if (Bjets.size() < num_Bjets) vetoEvent;


      //fill histograms:
      for (size_t i=0; i<num_ljets;i++){
          _h_pT_ljet[i]->fill(ljets[i].pT(), weight);
        }
      for (size_t i=0; i<num_bjets;i++){
          _h_pT_bjet[i]->fill(bjets[i].pT(), weight);
        }
      for (size_t i=0; i<num_Bjets;i++){
          _h_pT_Bjet[i]->fill(Bjets[i].pT(), weight);
          _h_Bjet_btags[i] ->fill(Bjets[i].bTags().size(), weight);
        }
      _h_ljet_multi->fill(ljets.size(), weight);
      _h_bjet_multi->fill(bjets.size(), weight);
      _h_Bjet_multi->fill(Bjets.size(), weight);


    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const double xs = crossSection()/picobarn;

      for (size_t i=0; i<num_ljets;i++){
          scale(_h_pT_ljet[i], xs/sumOfWeights());
        }
      for (size_t i=0; i<num_bjets;i++){
          scale(_h_pT_bjet[i], xs/sumOfWeights());
        }
      for (size_t i=0; i<num_Bjets;i++){
          scale(_h_pT_Bjet[i], xs/sumOfWeights());
          scale(_h_Bjet_btags[i], xs/sumOfWeights());
        }
      scale(_h_ljet_multi, xs/sumOfWeights());
      scale(_h_bjet_multi, xs/sumOfWeights());
      scale(_h_Bjet_multi, xs/sumOfWeights());

    }

    //@}

  private:
    size_t num_ljets;
    size_t num_bjets;
    size_t num_Bjets;
    double m_ptcut;
    size_t m_Bdef;
    size_t m_parton;

    /// @name Histograms
    // vector histograms: create with required size of jets
    std::vector<Histo1DPtr> _h_pT_ljet;
    std::vector<Histo1DPtr> _h_pT_bjet;
    std::vector<Histo1DPtr> _h_pT_Bjet;
    std::vector<Histo1DPtr> _h_Bjet_btags;


    Histo1DPtr _h_ljet_multi;
    Histo1DPtr _h_bjet_multi;
    Histo1DPtr _h_Bjet_multi;



    //@{

    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(MCTTB);


  //hadron level analyses

  class MCTTB_L2_b2_B0_PT30 : public MCTTB {
  public:
    MCTTB_L2_b2_B0_PT30()
      :MCTTB("MCTTB_L2_b2_B0_PT30",2,2,0,30)
    {   }
  };

  class MCTTB_L0_b2_B0_PT30 : public MCTTB {
  public:
    MCTTB_L0_b2_B0_PT30()
      :MCTTB("MCTTB_L0_b2_B0_PT30",0,2,0,30)
    {   }
  };

  class MCTTB_L0_b4_B0_PT30 : public MCTTB {
  public:
    MCTTB_L0_b4_B0_PT30()
      :MCTTB("MCTTB_L0_b4_B0_PT30",0,4,0,30)
    {   }
  };

  class MCTTB_L0_b3_B1_PT30 : public MCTTB {
  public:
    MCTTB_L0_b3_B1_PT30()
      :MCTTB("MCTTB_L0_b3_B1_PT30",0,3,1,30)
    {   }
  };

  class MCTTB_L0_b2_B1_PT30 : public MCTTB {
  public:
    MCTTB_L0_b2_B1_PT30()
      :MCTTB("MCTTB_L0_b2_B1_PT30",0,2,1,30)
    {   }
  };

  DECLARE_RIVET_PLUGIN(MCTTB_L2_b2_B0_PT30);
  DECLARE_RIVET_PLUGIN(MCTTB_L0_b2_B0_PT30);
  DECLARE_RIVET_PLUGIN(MCTTB_L0_b4_B0_PT30);
  DECLARE_RIVET_PLUGIN(MCTTB_L0_b3_B1_PT30);
  DECLARE_RIVET_PLUGIN(MCTTB_L0_b2_B1_PT30);


  //parton level analyses

  class MCTTB_L2_b2_B0_PT30_parton : public MCTTB {
  public:
    MCTTB_L2_b2_B0_PT30_parton()
      :MCTTB("MCTTB_L2_b2_B0_PT30_parton",2,2,0,30,2,1)
    {   }
  };

  class MCTTB_L0_b2_B0_PT30_parton : public MCTTB {
  public:
    MCTTB_L0_b2_B0_PT30_parton()
      :MCTTB("MCTTB_L0_b2_B0_PT30_parton",0,2,0,30,2,1)
    {   }
  };

  class MCTTB_L0_b4_B0_PT30_parton : public MCTTB {
  public:
    MCTTB_L0_b4_B0_PT30_parton()
      :MCTTB("MCTTB_L0_b4_B0_PT30_parton",0,4,0,30,2,1)
    {   }
  };

  class MCTTB_L0_b3_B1_PT30_parton : public MCTTB {
  public:
    MCTTB_L0_b3_B1_PT30_parton()
      :MCTTB("MCTTB_L0_b3_B1_PT30_parton",0,3,1,30,2,1)
    {   }
  };

  class MCTTB_L0_b2_B1_PT30_parton : public MCTTB {
  public:
    MCTTB_L0_b2_B1_PT30_parton()
      :MCTTB("MCTTB_L0_b2_B1_PT30_parton",0,2,1,30,2,1)
    {   }
  };

  DECLARE_RIVET_PLUGIN(MCTTB_L2_b2_B0_PT30_parton);
  DECLARE_RIVET_PLUGIN(MCTTB_L0_b2_B0_PT30_parton);
  DECLARE_RIVET_PLUGIN(MCTTB_L0_b4_B0_PT30_parton);
  DECLARE_RIVET_PLUGIN(MCTTB_L0_b3_B1_PT30_parton);
  DECLARE_RIVET_PLUGIN(MCTTB_L0_b2_B1_PT30_parton);





}
