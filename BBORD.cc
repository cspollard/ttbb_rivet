// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class BBORD : public Analysis {
  public:

    /// Constructor
    BBORD(const string& name="BBORD", double dr=0.4, double ptcut=15.0, bool cparton=true)
      : Analysis(name)
    {  m_ptcut=ptcut;
       m_dr=dr;
       m_parton=cparton;
    }



    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      ChargedLeptons lfs(FinalState(-5., 5., 30*GeV));
      //declare(lfs, "LFS");

      VetoedFinalState fs(FinalState(-5., 5., 0*GeV));
      fs.addVetoOnThisFinalState(lfs);
      declare(FastJets(fs, FastJets::ANTIKT, m_dr), "Jets");

      _h_bjet_multi = bookHisto1D("bjet_multi", 7, -0.5,  7 -0.5);


    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      double weight = event.weight();

      const FastJets& jetpro = apply<FastJets>(event, "Jets");
      const Jets alljets = jetpro.jetsByPt(Cuts::pT>m_ptcut*GeV);
      if (alljets.size() < 2) {
          MSG_DEBUG("Event failed jet multiplicity cut");
          vetoEvent;
        }

      Jets bjets;
      size_t count(0);
      foreach (Jet jet, alljets) {
          if (bjets.size()!=count) break;
          if (m_parton==0){  //hadron level
              if(jet.bTagged()) {
                bjets.push_back(jet);
                }

            }else{ //parton level
              size_t nb(0);
              Particles jet_content= jet.constituents();
              //MSG_INFO("jet content" << jet_content);
              foreach(Particle part, jet_content){
                  if (part.hasBottom()) nb++;
                }
              if (nb>0) {
                  bjets.push_back(jet);
              }

            }
          count++;
        }
      if (bjets.size()>0)  MSG_INFO("bjet size: " << bjets.size());
      _h_bjet_multi->fill(bjets.size(), weight);


    }


    /// Normalise histograms etc., after the run
    void finalize() {

      scale(_h_bjet_multi, crossSection()/picobarn/sumOfWeights());


    }

    //@}

   private:

    double m_dr;
    double m_ptcut;
    bool m_parton;

    /// @name Histograms
    //@{
    Histo1DPtr _h_bjet_multi;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(BBORD);



  class BBORD_DR04_PT30_parton : public BBORD {
  public:
    BBORD_DR04_PT30_parton()
      :BBORD("BBORD_DR04_PT30_parton",0.4,30,true)
    {   }
  };
  DECLARE_RIVET_PLUGIN(BBORD_DR04_PT30_parton);



  class BBORD_DR04_PT15_parton : public BBORD {
  public:
    BBORD_DR04_PT15_parton()
      :BBORD("BBORD_DR04_PT15_parton",0.4,15,true)
    {   }
  };
  DECLARE_RIVET_PLUGIN(BBORD_DR04_PT15_parton);



  class BBORD_DR04_PT80_parton : public BBORD {
  public:
    BBORD_DR04_PT80_parton()
      :BBORD("BBORD_DR04_PT80_parton",0.4,80,true)
    {   }
  };
  DECLARE_RIVET_PLUGIN(BBORD_DR04_PT80_parton);



  class BBORD_DR02_PT30_parton : public BBORD {
  public:
    BBORD_DR02_PT30_parton()
      :BBORD("BBORD_DR02_PT30_parton",0.2,30,true)
    {   }
  };
  DECLARE_RIVET_PLUGIN(BBORD_DR02_PT30_parton);

  class BBORD_DR08_PT30_parton : public BBORD {
  public:
    BBORD_DR08_PT30_parton()
      :BBORD("BBORD_DR08_PT30_parton",0.8,30,true)
    {   }
  };
  DECLARE_RIVET_PLUGIN(BBORD_DR08_PT30_parton);

  class BBORD_DR04_PT10_parton : public BBORD {
  public:
    BBORD_DR04_PT10_parton()
      :BBORD("BBORD_DR04_PT15_parton",0.4,10,true)
    {   }
  };
  DECLARE_RIVET_PLUGIN(BBORD_DR04_PT10_parton);



}



