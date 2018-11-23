#include <unordered_map>

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
        int mode = 1;
        //mode implmentation for unfolded and background subtraction
        Cut eta_full = (Cuts::abseta < 5.0);
        // Lepton cuts
        Cut lep_cuts25 = (Cuts::abseta < 2.5) && (Cuts::pT >= 25*GeV);
        // All final state particles
        FinalState fs(eta_full);

        // Get photons to dress leptons
        IdentifiedFinalState pho_id(fs);
        pho_id.acceptIdPair(PID::PHOTON);

        PromptFinalState photons(pho_id);
        photons.acceptTauDecays(true);

        // Projection to find the electrons
        IdentifiedFinalState el_id(fs);
        el_id.acceptIdPair(PID::ELECTRON);
        PromptFinalState electrons(el_id);
        electrons.acceptTauDecays(true);

        // Projection to find the muons
        IdentifiedFinalState mu_id(fs);
        mu_id.acceptIdPair(PID::MUON);
        PromptFinalState muons(mu_id);
        muons.acceptTauDecays(true);

        DressedLeptons dressedelectrons25(photons, electrons, 0.1, lep_cuts25, true, true);
        DressedLeptons dressedmuons25(photons, muons, 0.1, lep_cuts25, true, true);

        addProjection(dressedelectrons25, "electrons");
        addProjection(dressedmuons25, "muons");

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

        // Eventually this will all be taken from the data file anyway...
        std::vector<double> bin_ptlead_SL = {30.0, 90.0, 120.0, 150.0, 190.0, 2025.0};
        std::vector<double> bin_ptsublead_SL = {25.0, 65.0, 85.0, 120.0, 1510.0};
        std::vector<double> bin_pt3rdjet_SL = {25.0, 45.0, 65.0, 90.0, 670.0};
        std::vector<double> bin_pt4thjet_SL = {25.0, 30.0, 45.0, 480.0};
        std::vector<double> bin_hthad_SL = {185.0, 390.0, 470.0, 575.0, 710.0, 5925.0};
        std::vector<double> bin_ht_SL = {230.0, 450.0, 535.0, 640.0, 765.0, 6000.0};
        std::vector<double> bin_htb_SL = {110.0, 240.0, 310.0, 400.0, 3665.0};
        std::vector<double> bin_mbb_SL = {20.0, 125.0, 175.0, 230.0, 295.0, 390.0, 4410.0};
        std::vector<double> bin_ptbb_SL = {0.0, 50.0, 85.0, 120.0, 160.0, 220.0, 1720.0};
        std::vector<double> bin_drbb_SL = {.35, 1.25, 1.85, 2.25, 2.55, 2.80, 3.05, 5.8};
        std::vector<double> bin_mbbclosest_SL ={10.0, 40.0, 60.0, 85.0, 120.0, 1460.0};
        std::vector<double> bin_ptbbclosest_SL = {15.0, 90.0, 120.0, 160.0, 210.0, 2210.0};
        std::vector<double> bin_drbbclosest_SL = {.35, .60, .75, .90, 1.10, 1.35, 3.15};
        std::vector<double> bin_ptlightjet_SL = {25, 60, 80, 105, 140, 195, 2190};

        std::vector<double> bin_ptlead_DIL = {25, 65, 100, 150, 200, 1000};  // Checked
        std::vector<double> bin_ptsublead_DIL = {25, 50, 90, 150, 1000};  // Checked
        std::vector<double> bin_pt3rdjet_DIL = {25, 50, 80, 1000};  // Checked
        std::vector<double> bin_hthad_DIL = {100, 225, 350, 600, 2000};  // Checked
        std::vector<double> bin_ht_DIL = {100, 270, 350, 450, 700, 2000};  // Checked
        std::vector<double> bin_mbb_DIL = {0, 50, 100, 200, 400, 2000};  // Checked
        std::vector<double> bin_ptbb_DIL = {0, 65, 100, 140, 1000};  // Checked
        std::vector<double> bin_drbb_DIL = {0.4, 1.3, 1.8, 2.3, 2.8, 8.0};   // Checked
        std::vector<double> bin_mbbclosest_DIL = {0, 55, 100, 145, 250, 1000};  // Checked
        std::vector<double> bin_ptbbclosest_DIL = {0, 65, 100, 150, 250, 1000};  // Checked
        std::vector<double> bin_drbbclosest_DIL = {0.4, 0.7, 1.0, 1.3, 3.0};  // Checked


        // We are going to book the histograms in the order that they 
        // appear in the paper...

        // fiducial cross-section histogram
        _histograms["fid_xsec"] = bookHisto1D(1, 1, 1);

        _histograms["nbjets_emu"] = bookHisto1D(2, 1, 1);

        // HT
        book_hist_emu("ht", 3, 1, 1);
        book_hist_emu("ht_had", 4, 1, 1);

        book_hist_ljets("ht", 5, 1, 1);
        book_hist_ljets("ht_had", 6, 1, 1);

        // b-jet pTs
        book_hist_emu("lead_bjet_pt", 7, 1, 1);
        book_hist_emu("sublead_bjet_pt", 8, 1, 1);
        book_hist_emu("third_bjet_pt", 9, 1, 1);

        book_hist_ljets("lead_bjet_pt", 10, 1, 1);
        book_hist_ljets("sublead_bjet_pt", 11, 1, 1);
        book_hist_ljets("third_bjet_pt", 12, 1, 1);
        book_hist_ljets("fourth_bjet_pt", 13, 1, 1);

        // leading bb pair
        book_hist_emu("m_bb_leading", 14, 1, 1);
        book_hist_emu("pt_bb_leading", 15, 1, 1);
        book_hist_emu("dR_bb_leading", 16, 1, 1);

        book_hist_ljets("m_bb_leading", 17, 1, 1);
        book_hist_ljets("pt_bb_leading", 18, 1, 1);
        book_hist_ljets("dR_bb_leading", 19, 1, 1);

        // closest bb pair
        book_hist_emu("m_bb_closest", 20, 1, 1);
        book_hist_emu("pt_bb_closest", 21, 1, 1);
        book_hist_emu("dR_bb_closest", 22, 1, 1);

        book_hist_ljets("m_bb_closest", 23, 1, 1);
        book_hist_ljets("pt_bb_closest", 24, 1, 1);
        book_hist_ljets("dR_bb_closest", 25, 1, 1);

        // // fiducial cross-section histogram
        // _histograms["fid_xsec"] = bookHisto1D("fid_xsec", 4, 0, 4);

        // _histograms["nbjets_emu"] = bookHisto1D("nbjets_emu", 3, 2, 5);

        // // HT
        // book_hist_emu("ht", bin_ht_DIL);
        // book_hist_emu("ht_had", bin_hthad_DIL);

        // book_hist_ljets("ht", bin_ht_SL);
        // book_hist_ljets("ht_had", bin_hthad_SL);

        // // b-jet pTs
        // book_hist_emu("lead_bjet_pt", bin_ptlead_DIL);
        // book_hist_emu("sublead_bjet_pt", bin_ptsublead_DIL);
        // book_hist_emu("third_bjet_pt", bin_pt3rdjet_DIL);

        // book_hist_ljets("lead_bjet_pt", bin_ptlead_SL);
        // book_hist_ljets("sublead_bjet_pt", bin_ptsublead_SL);
        // book_hist_ljets("third_bjet_pt", bin_pt3rdjet_SL);
        // book_hist_ljets("fourth_bjet_pt", bin_pt4thjet_SL);

        // // leading bb pair
        // book_hist_emu("m_bb_leading", bin_mbb_DIL);
        // book_hist_emu("pt_bb_leading", bin_ptbb_DIL);
        // book_hist_emu("dR_bb_leading", bin_drbb_DIL);

        // book_hist_ljets("m_bb_leading", bin_mbb_SL);
        // book_hist_ljets("pt_bb_leading", bin_ptbb_SL);
        // book_hist_ljets("dR_bb_leading", bin_drbb_SL);

        // // closest bb pair
        // book_hist_emu("m_bb_closest", bin_mbbclosest_DIL);
        // book_hist_emu("pt_bb_closest", bin_ptbbclosest_DIL);
        // book_hist_emu("dR_bb_closest", bin_drbbclosest_DIL);

        // book_hist_ljets("m_bb_closest", bin_mbbclosest_SL);
        // book_hist_ljets("pt_bb_closest", bin_ptbbclosest_SL);
        // book_hist_ljets("dR_bb_closest", bin_drbbclosest_SL);
      }


      void analyze(const Event& event) {
        double weight = event.weight();

        vector<DressedLepton> electrons = applyProjection<DressedLeptons>(event, "electrons").dressedLeptons();
        vector<DressedLepton> muons = applyProjection<DressedLeptons>(event, "muons").dressedLeptons();
        const Jets jets = applyProjection<FastJets>(event, "jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);

        vector<DressedLepton> leptons = overlap_removal_leptons(jets, electrons, muons);

        Jets bjets;
        foreach(const Jet& jet, jets) {
          if (jet.bTagged(Cuts::pT >= 5.0*GeV)) {
            bjets += jet;
          }
        }

        int njets = jets.size();
        int nbjets = bjets.size();

        // Evaluate basic event selection
        bool pass_ljets = (leptons.size() == 1 && leptons[0].pT() * GeV > 27);

        bool pass_emu =
          // 2 leptons > 27 GeV
          (leptons.size() == 2) &&
          (leptons[0].pT() * GeV > 27 && leptons[1].pT() * GeV > 27) &&
          // emu events
          ((leptons[0].abspid() == 11 && leptons[1].abspid() == 13) ||
           (leptons[0].abspid() == 13 && leptons[1].abspid() == 11)) &&
          // opposite charge
          (leptons[0].charge() != leptons[1].charge());

        // If we don't have exactly 1 or 2 leptons then veto the event
        if (!pass_emu && !pass_ljets) vetoEvent;

        if (pass_emu) {
          if (nbjets >= 2) {
            _histograms["nbjets_emu"]->fill(nbjets, weight);
          }
          if (nbjets >= 3) {
            _histograms["fid_xsec"]->fill(1.5, weight);
          }
          if (nbjets >= 4) {
            _histograms["fid_xsec"]->fill(0.5, weight);
          }
        }

        if (pass_ljets) {
          if (nbjets >= 3 && njets >= 5) {
            _histograms["fid_xsec"]->fill(3.5, weight);
          }
          if (nbjets >= 4 && njets >= 6) {
            _histograms["fid_xsec"]->fill(2.5, weight);
          }
        }

        if (pass_emu && (nbjets < 3 || njets < 3)) vetoEvent;
        if (pass_ljets && (nbjets < 4 || njets < 6)) vetoEvent;

        double hthad = sum(jets, pT, 0.0);
        double ht = sum(leptons, pT, hthad);

        FourMomentum jsum = bjets[0].momentum() + bjets[1].momentum();
        double dr_leading = deltaR(bjets[0].momentum(), bjets[1].momentum());

        std::pair<size_t, size_t> indices = get_min_dr(bjets);
        FourMomentum bb_closest =
          bjets[indices.first].momentum() + bjets[indices.second].momentum();
        double dr_closest =
          deltaR(bjets[indices.first].momentum(), bjets[indices.second].momentum());

        if (pass_ljets) {
          // b-jet pTs
          _histograms["lead_bjet_pt_ljets"]->fill(bjets[0].pT() * GeV, weight);
          _histograms["sublead_bjet_pt_ljets"]->fill(bjets[1].pT() * GeV, weight);
          _histograms["third_bjet_pt_ljets"]->fill(bjets[2].pT() * GeV, weight);
          if (nbjets >= 4) {
            _histograms["fourth_bjet_pt_ljets"]->fill(bjets[3].pT() * GeV, weight);
          }

          // HT
          _histograms["ht_ljets"]->fill(ht * GeV, weight);
          _histograms["ht_had_ljets"]->fill(hthad * GeV, weight);

          // leading bb pair
          _histograms["m_bb_leading_ljets"]->fill(jsum.mass() * GeV, weight);
          _histograms["pt_bb_leading_ljets"]->fill(jsum.pT() * GeV, weight);
          _histograms["dR_bb_leading_ljets"]->fill(dr_leading, weight);

          // closest bb pair
          _histograms["m_bb_closest_ljets"]->fill(bb_closest.mass() * GeV, weight);
          _histograms["pt_bb_closest_ljets"]->fill(bb_closest.pT() * GeV, weight);
          _histograms["dR_bb_closest_ljets"]->fill(dr_closest, weight);
        }
        if (pass_emu) {
          // b-jet pTs
          _histograms["lead_bjet_pt_emu"]->fill(bjets[0].pT() * GeV, weight);
          _histograms["sublead_bjet_pt_emu"]->fill(bjets[1].pT() * GeV, weight);
          _histograms["third_bjet_pt_emu"]->fill(bjets[2].pT() * GeV, weight);

          // HT
          _histograms["ht_emu"]->fill(ht * GeV, weight);
          _histograms["ht_had_emu"]->fill(hthad * GeV, weight);

          // leading bb pair
          _histograms["m_bb_leading_emu"]->fill(jsum.mass() * GeV, weight);
          _histograms["pt_bb_leading_emu"]->fill(jsum.pT() * GeV, weight);
          _histograms["dR_bb_leading_emu"]->fill(dr_leading, weight);

          // closest bb pair
          _histograms["m_bb_closest_emu"]->fill(bb_closest.mass() * GeV, weight);
          _histograms["pt_bb_closest_emu"]->fill(bb_closest.pT() * GeV, weight);
          _histograms["dR_bb_closest_emu"]->fill(dr_closest, weight);
        }
      }


      void finalize() {
        // Normalise all histograms
        for (auto const& h : _histograms) {
          if (h.first == "fid_xsec") continue;
          normalize(h.second, 1.0);
        }
        const double sf = crossSection() / femtobarn / sumOfWeights();
        scale(_histograms["fid_xsec"], sf);
        _histograms["fid_xsec"]->setAnnotation("bin1", "emu_4b");
        _histograms["fid_xsec"]->setAnnotation("bin2", "emu_3b");
        _histograms["fid_xsec"]->setAnnotation("bin3", "ljets_4b");
        _histograms["fid_xsec"]->setAnnotation("bin4", "ljets_3b");
      }

    private:
      std::unordered_map<std::string, Histo1DPtr> _histograms;

      std::pair<size_t, size_t> get_min_dr(const Jets& bjets) {
        double mindr = 999;
        std::pair<size_t, size_t> indices;
        for (size_t i = 0; i < bjets.size(); i++) {
          for (size_t j = 0; j < bjets.size(); j++) {
            if (i == j) continue;
            double dr = deltaR(bjets[i], bjets[j]);
            if (dr < mindr) {
              indices.first = i;
              indices.second = j;
              mindr = dr;
            }
          }
        }
        return indices;
      }

      vector<DressedLepton>
        overlap_removal_leptons(const Jets& jets,
            const vector<DressedLepton>& electrons,
            const vector<DressedLepton>& muons) {
          std::vector<DressedLepton> leptons;
          foreach(const DressedLepton& electron, electrons) {
            if (keep_lepton(jets, electron)) leptons.push_back(electron);
          }
          foreach(const DressedLepton& muon, muons) {
            if (keep_lepton(jets, muon)) leptons.push_back(muon);
          }
          return leptons;
        }

      bool keep_lepton(const Jets& jets, const DressedLepton& lepton) {
        foreach(const Jet& jet, jets) {
          if (deltaR(jet.momentum(), lepton.momentum()) < 0.4) return false;
        }
        return true;
      }

      inline void book_hist_ljets(const std::string& name, std::vector<double> bins) {
        _histograms[name + "_ljets"] = bookHisto1D(name + "_ljets", bins);
      }

      inline void book_hist_ljets(const std::string& name, int d, int x=1, int y=1) {
        _histograms[name + "_ljets"] = bookHisto1D(d, x, y);
      }

      inline void book_hist_emu(const std::string& name, std::vector<double> bins) {
        _histograms[name + "_emu"] = bookHisto1D(name + "_emu", bins);
      }
      inline void book_hist_emu(const std::string& name, int d, int x=1, int y=1) {
        _histograms[name + "_emu"] = bookHisto1D(d, x, y);
      }
  };

  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(ATLAS_XX);

}  // namespace Rivet

