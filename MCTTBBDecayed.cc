// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Tools/JetSmearingFunctions.hh"
#include "TTBBHist.hh"

namespace Rivet {

  using V4s = vector<FourMomentum>;

  class TTBBDecayedHists {
    public:

      TTBBDecayedHists() { };

      TTBBDecayedHists(const string& prefix) {
        h_njl = make_shared<TTBBHist>("h_" + prefix + "_njl", 10, -0.5, 9.5, "", "light-jet multiplicity", dsdx("n", "1"));
        h_njb = make_shared<TTBBHist>("h_" + prefix + "_njb", 8, -0.5, 7.5, "", "$b$-jet mulitplicity", dsdx("n", "1"));
        h_nj2b = make_shared<TTBBHist>("h_" + prefix + "_nj2b", 5,-0.5, 4.5, "", "$bb$-jet multiplicity", dsdx("n", "1"));
        h_nj1b = make_shared<TTBBHist>("h_" + prefix + "_nj1b", 5, -0.5, 4.5, "", "$b1$-jet multiplicity", dsdx("n", "1"));

        p_fracj1b = make_shared<Profile1D>(7, 0.5, 7.5, "p_" + prefix + "_fracj1b", "");
        p_fracj1b->setAnnotation("XTitle", "$b1$-jet fraction");
        p_fracj1b->setAnnotation("YTitle", dsdx("frac", "1"));

        p_fracj2b = make_shared<Profile1D>(7, 0.5, 7.5, "p_" + prefix + "_fracj2b", "");
        p_fracj2b->setAnnotation("XTitle", "$bb$-jet fraction");
        p_fracj2b->setAnnotation("YTitle", dsdx("frac", "1"));

        h_jl1pt = make_shared<TTBBHist>("h_" + prefix + "_jl1pt", 50, 0, 500*GeV, "", "leading light-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_jl2pt = make_shared<TTBBHist>("h_" + prefix + "_jl2pt", 50, 0, 500*GeV, "", "subleading light-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));

        h_jb1pt = make_shared<TTBBHist>("h_" + prefix + "_jb1pt", 50, 0, 500*GeV, "", "leading $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_jb2pt = make_shared<TTBBHist>("h_" + prefix + "_jb2pt", 50, 0, 500*GeV, "", "subleading $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_jb3pt = make_shared<TTBBHist>("h_" + prefix + "_jb3pt", 50, 0, 500*GeV, "", "3rd $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_jb4pt = make_shared<TTBBHist>("h_" + prefix + "_jb4pt", 50, 0, 500*GeV, "", "4th $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));

        h_j2b1pt = make_shared<TTBBHist>("h_" + prefix + "_j2b1pt", 50, 0, 500*GeV, "", "leading $bb$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_j2b2pt = make_shared<TTBBHist>("h_" + prefix + "_j2b2pt", 50, 0, 500*GeV, "", "subleading $bb$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_j2b3pt = make_shared<TTBBHist>("h_" + prefix + "_j2b3pt", 50, 0, 500*GeV, "", "3rd $bb$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_j2b4pt = make_shared<TTBBHist>("h_" + prefix + "_j2b4pt", 50, 0, 500*GeV, "", "4th $bb$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));

        h_j1b1pt = make_shared<TTBBHist>("h_" + prefix + "_j1b1pt", 50, 0, 500*GeV, "", "leading $b1$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_j1b2pt = make_shared<TTBBHist>("h_" + prefix + "_j1b2pt", 50, 0, 500*GeV, "", "subleading $b1$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_j1b3pt = make_shared<TTBBHist>("h_" + prefix + "_j1b3pt", 50, 0, 500*GeV, "", "3rd $b1$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_j1b4pt = make_shared<TTBBHist>("h_" + prefix + "_j1b4pt", 50, 0, 500*GeV, "", "4th $b1$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));


        h_jl1pt_log = make_shared<TTBBHist>("h_" + prefix + "_jl1pt_log", logspace(50, 10, 1000*GeV), "", "leading light-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_jl2pt_log = make_shared<TTBBHist>("h_" + prefix + "_jl2pt_log", logspace(50, 10, 1000*GeV), "", "subleading light-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_jb1pt_log = make_shared<TTBBHist>("h_" + prefix + "_jb1pt_log", logspace(50, 10, 1000*GeV), "", "leading $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_jb2pt_log = make_shared<TTBBHist>("h_" + prefix + "_jb2pt_log", logspace(50, 10, 1000*GeV), "", "subleading $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_jb3pt_log = make_shared<TTBBHist>("h_" + prefix + "_jb3pt_log", logspace(50, 10, 1000*GeV), "", "3rd $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_jb4pt_log = make_shared<TTBBHist>("h_" + prefix + "_jb4pt_log", logspace(50, 10, 1000*GeV), "", "4th $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));

        h_j2b1pt_log = make_shared<TTBBHist>("h_" + prefix + "_j2b1pt_log", logspace(50, 10, 1000*GeV), "", "leading $bb$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_j2b2pt_log = make_shared<TTBBHist>("h_" + prefix + "_j2b2pt_log", logspace(50, 10, 1000*GeV), "", "subleading $bb$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_j2b3pt_log = make_shared<TTBBHist>("h_" + prefix + "_j2b3pt_log", logspace(50, 10, 1000*GeV), "", "3rd $bb$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_j2b4pt_log = make_shared<TTBBHist>("h_" + prefix + "_j2b4pt_log", logspace(50, 10, 1000*GeV), "", "4th $bb$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));

        h_j1b1pt_log = make_shared<TTBBHist>("h_" + prefix + "_j1b1pt_log", logspace(50, 10, 1000*GeV), "", "leading $b1$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_j1b2pt_log = make_shared<TTBBHist>("h_" + prefix + "_j1b2pt_log", logspace(50, 10, 1000*GeV), "", "subleading $b1$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_j1b3pt_log = make_shared<TTBBHist>("h_" + prefix + "_j1b3pt_log", logspace(50, 10, 1000*GeV), "", "3rd $b1$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_j1b4pt_log = make_shared<TTBBHist>("h_" + prefix + "_j1b4pt_log", logspace(50, 10, 1000*GeV), "", "4th $b1$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));

        h_jl1eta = make_shared<TTBBHist>("h_" + prefix + "_jl1eta", 30, -3, 3, "", "leading light-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        h_jl2eta = make_shared<TTBBHist>("h_" + prefix + "_jl2eta", 30, -3, 3, "", "subleading light-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        h_jb1eta = make_shared<TTBBHist>("h_" + prefix + "_jb1eta", 30, -3, 3, "", "leading $b$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        h_jb2eta = make_shared<TTBBHist>("h_" + prefix + "_jb2eta", 30, -3, 3, "", "subleading $b$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        h_j2b1eta = make_shared<TTBBHist>("h_" + prefix + "_j2b1eta", 30, -3, 3, "", "leading $bb$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        h_j2b2eta = make_shared<TTBBHist>("h_" + prefix + "_j2b2eta", 30, -3, 3, "", "subleading $bb$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        h_j1b1eta = make_shared<TTBBHist>("h_" + prefix + "_j1b1eta", 30, -3, 3, "", "leading $b1$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        h_j1b2eta = make_shared<TTBBHist>("h_" + prefix + "_j1b2eta", 30, -3, 3, "", "subleading $b1$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));

        h_lep1pt = make_shared<TTBBHist>("h_" + prefix + "_lep1pt", 50, 0, 500*GeV, "", "leading lepton " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_lep2pt = make_shared<TTBBHist>("h_" + prefix + "_lep2pt", 50, 0, 500*GeV, "", "subleading lepton " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_lep1eta = make_shared<TTBBHist>("h_" + prefix + "_lep1eta", 30, -3, 3, "", "leading lepton " + seta, dsdx(seta, "\\mathrm{GeV}"));
        h_lep2eta = make_shared<TTBBHist>("h_" + prefix + "_lep2eta", 30, -3, 3, "", "subleading lepton " + seta, dsdx(seta, "\\mathrm{GeV}"));

        h_dphibb = make_shared<TTBBHist>("h_" + prefix + "_dphibb", 50, 0, 4, "", sdphibb, dsdx(sdphibb, "1"));
        h_drbb = make_shared<TTBBHist>("h_" + prefix + "_drbb", 50, 0, 5, "", sdrbb, dsdx(sdrbb, "1"));

        h_mbb_avg = make_shared<TTBBHist>("h_" + prefix + "_mbb_avg", 50, 0, 500*GeV, "", "average " + smbb + " [GeV]", dsdx(smbb, "\\mathrm{GeV}"));
        h_mbb_wind = make_shared<TTBBHist>("h_" + prefix + "_mbb_wind", 5, -0.5, 4.5, "", "number of " + smbb + " combinations in (100, 150) GeV", dsdx("n", "1"));

        h_mbb = make_shared<TTBBHist>("h_" + prefix + "_mbb", 50, 0, 500*GeV, "", smbb + " [GeV]", dsdx(smbb, "\\mathrm{GeV}"));
        h_ptbb = make_shared<TTBBHist>("h_" + prefix + "_ptbb", 50, 0, 500*GeV, "", sptbb + " [GeV]", dsdx(sptbb, "\\mathrm{GeV}"));
        h_ht = make_shared<TTBBHist>("h_" + prefix + "_ht", 50, 0, 2000*GeV, "", sht + " [GeV]", dsdx(sht, "\\mathrm{GeV}"));

        h_mbb_log = make_shared<TTBBHist>("h_" + prefix + "_mbb_log", logspace(50, 10, 1000*GeV), "", smbb + " [GeV]", dsdx(smbb, "\\mathrm{GeV}"));
        h_ptbb_log = make_shared<TTBBHist>("h_" + prefix + "_ptbb_log", logspace(50, 10, 1000*GeV), "", sptbb + " [GeV]", dsdx(sptbb, "\\mathrm{GeV}"));
        h_ht_log = make_shared<TTBBHist>("h_" + prefix + "_ht_log", logspace(50, 10, 2000*GeV), "", sht + " [GeV]", dsdx(sht, "\\mathrm{GeV}"));

      }

      void fill(double weight, const V4s& jls, const V4s& jbs, const vector<DressedLepton>& leps) {
        size_t njbs = jbs.size();

        h_njl->fill(jls.size(), weight);
        h_njb->fill(njbs, weight);

        if (jls.size() >= 1) {
          h_jl1pt->fill(jls[0].pt(), weight);
          h_jl1pt_log->fill(jls[0].pt(), weight);
          h_jl1eta->fill(jls[0].eta(), weight);
        }

        if (jls.size() >= 2) {
          h_jl2pt->fill(jls[1].pt(), weight);
          h_jl2pt_log->fill(jls[1].pt(), weight);
          h_jl2eta->fill(jls[1].eta(), weight);
        }

        if (jbs.size() >= 1) {
          h_jb1pt->fill(jbs[0].pt(), weight);
          h_jb1pt_log->fill(jbs[0].pt(), weight);
          h_jb1eta->fill(jbs[0].eta(), weight);
        }

        if (jbs.size() >= 2) {
          h_jb2pt->fill(jbs[1].pt(), weight);
          h_jb2pt_log->fill(jbs[1].pt(), weight);
          h_jb2eta->fill(jbs[1].eta(), weight);
        }

        if (jbs.size() >= 3) {
          h_jb3pt->fill(jbs[2].pt(), weight);
          h_jb3pt_log->fill(jbs[2].pt(), weight);
        }

        if (jbs.size() >= 4) {
          h_jb4pt->fill(jbs[3].pt(), weight);
          h_jb4pt_log->fill(jbs[3].pt(), weight);
        }

        if (leps.size() >= 1) {
          h_lep1pt->fill(leps[0].pt(), weight);
          h_lep1eta->fill(leps[0].eta(), weight);
        }

        if (leps.size() >= 2) {
          h_lep2pt->fill(leps[1].pt(), weight);
          h_lep2eta->fill(leps[1].eta(), weight);
        }

        // we define the ht as the scalar sum of all the light, b, and B jet pts
        double ht = 0.0;
        for (const Jet& lj: jls)
          ht += lj.pt();

        for (const Jet& bj: jbs)
          ht += bj.pt();

        h_ht->fill(ht, weight);
        h_ht_log->fill(ht, weight);

        // find the two leading b-jets
        FourMomentum b1, b2;
        if (jbs.size() >= 2) {
          b1 = jbs[0];
          b2 = jbs[1];
        } else
          return;

        h_dphibb->fill(abs(deltaPhi(b1, b2)), weight);
        h_drbb->fill(deltaR(b1, b2), weight);

        h_mbb->fill((b1 + b2).mass(), weight);
        h_ptbb->fill((b1 + b2).pt(), weight);
        h_mbb_log->fill((b1 + b2).mass(), weight);
        h_ptbb_log->fill((b1 + b2).pt(), weight);

        size_t ncomb = jbs.size()*(jbs.size() - 1);
        size_t nwind = 0;
        double mbbsum = 0;
        for (size_t ib = 0; ib < jbs.size(); ib++) {
          for (size_t jb = ib+1; jb < jbs.size(); jb++) {
            double mbb = (jbs[ib] + jbs[jb]).mass();
            mbbsum += mbb;

            if (100*GeV < mbb && mbb < 150*GeV)
              nwind++;
          }
        }

        h_mbb_avg->fill(mbbsum/ncomb, weight);
        h_mbb_wind->fill(nwind, weight);

      };

      // this has to be a vector<Histo1DPtr> rather than vector<Histo1D> because
      // the Histo1D copy constructor loses all annotations?!?!?!?
      vector<shared_ptr<TTBBHist>> histograms() {
        return
          { h_njl, h_njb, h_nj2b, h_nj1b
          , h_jl1pt, h_jl2pt
          , h_jb1pt, h_jb2pt, h_jb3pt, h_jb4pt
          , h_j2b1pt, h_j2b2pt, h_j2b3pt, h_j2b4pt
          , h_j1b1pt, h_j1b2pt, h_j1b3pt, h_j1b4pt
          , h_jl1pt_log, h_jl2pt_log
          , h_jb1pt_log, h_jb2pt_log, h_jb3pt_log, h_jb4pt_log
          , h_j2b1pt_log, h_j2b2pt_log, h_j2b3pt_log, h_j2b4pt_log
          , h_j1b1pt_log, h_j1b2pt_log, h_j1b3pt_log, h_j1b4pt_log
          , h_jl1eta, h_jl2eta, h_jb1eta, h_jb2eta
          , h_j2b1eta, h_j2b2eta, h_j1b1eta, h_j1b2eta
          , h_lep1eta, h_lep2eta, h_lep1pt, h_lep2pt
          , h_mbb, h_mbb_log, h_mbb_avg, h_mbb_wind
          , h_dphibb, h_drbb
          , h_ptbb, h_ht
          , h_ptbb_log, h_ht_log
          };
      }

      vector<Profile1DPtr> profiles() {
        return { p_fracj1b, p_fracj2b };
      }

    private:
      string dsdx(const string& x, const string& xunit) {
        return "\\ensuremath{\\frac{d\\sigma}{d" + x + "} \\frac{{pb}}{" + xunit + "}}";
      }


      shared_ptr<TTBBHist>
          h_njl, h_njb, h_nj2b, h_nj1b
          , h_jl1pt, h_jl2pt
          , h_jb1pt, h_jb2pt, h_jb3pt, h_jb4pt
          , h_j2b1pt, h_j2b2pt, h_j2b3pt, h_j2b4pt
          , h_j1b1pt, h_j1b2pt, h_j1b3pt, h_j1b4pt
          , h_jl1pt_log, h_jl2pt_log
          , h_jb1pt_log, h_jb2pt_log, h_jb3pt_log, h_jb4pt_log
          , h_j2b1pt_log, h_j2b2pt_log, h_j2b3pt_log, h_j2b4pt_log
          , h_j1b1pt_log, h_j1b2pt_log, h_j1b3pt_log, h_j1b4pt_log
          , h_jl1eta, h_jl2eta, h_jb1eta, h_jb2eta
          , h_j2b1eta, h_j2b2eta, h_j1b1eta, h_j1b2eta
          , h_lep1eta, h_lep2eta, h_lep1pt, h_lep2pt
          , h_mbb, h_mbb_log, h_mbb_avg, h_mbb_wind
          , h_dphibb, h_drbb
          , h_ptbb, h_ht
          , h_ptbb_log, h_ht_log;

      Profile1DPtr p_fracj1b, p_fracj2b;

  };

  class MCTTBBDecayed : public Analysis {
    public:

      DEFAULT_RIVET_ANALYSIS_CTOR(MCTTBBDecayed);

      void init() {
        FinalState fps(Cuts::abseta < 5);
        PromptFinalState pfps(fps);

        IdentifiedFinalState ys(fps);
        ys.acceptId(PID::PHOTON);

        IdentifiedFinalState bareleps(pfps);
        bareleps.acceptChLeptons();

        DressedLeptons dressedleps(ys, bareleps, 0.1
          , Cuts::abseta < 2.5 && Cuts::pT > 25*GeV, true, true);
        declare(dressedleps, "Leptons");

        VetoedFinalState jfs(fps);
        jfs.addVetoOnThisFinalState(dressedleps);

        declare(FastJets(jfs, FastJets::ANTIKT, 0.4), "Jets");

        histmap = map< string , TTBBDecayedHists >();

        histmap["onelep_inclusive"] = TTBBDecayedHists("onelep_inclusive");
        histmap["onelep_ge6j_eq2jb"] = TTBBDecayedHists("onelep_ge6j_eq2jb");
        histmap["onelep_ge6j_eq3jb"] = TTBBDecayedHists("onelep_ge6j_eq3jb");
        histmap["onelep_ge6j_ge4jb"] = TTBBDecayedHists("onelep_ge6j_ge4jb");

        histmap["dilep_inclusive"] = TTBBDecayedHists("dilep_inclusive");
        histmap["dilep_ge4j_eq2jb"] = TTBBDecayedHists("dilep_ge4j_eq2jb");
        histmap["dilep_ge4j_eq3jb"] = TTBBDecayedHists("dilep_ge4j_eq3jb");
        histmap["dilep_ge4j_ge4jb"] = TTBBDecayedHists("dilep_ge4j_ge4jb");


      }


      void analyze(const Event& event) {

        // 70% b-tagging efficiency, 5% charm mistag rate
        JET_BTAG_EFFS tagger(0.85, 0.05, 0.0);


        double weight = event.weight();

        const Jets& jets =
          apply<FastJets>(event, "Jets").jetsByPt(
              Cuts::pT > 25*GeV && Cuts::abseta < 2.5);

        const vector<DressedLepton>& leps =
          apply<DressedLeptons>(event, "Leptons").dressedLeptons();

        // find the light, b, and B jets by checking the number of b-hadrons
        // associated to the jet.
        // jls = light jets
        // jbs = jets with at least one associated b-hadron
        // j1bs = jets with exactly one associated b-hadron
        // j2bs = jets with at least two associated b-hadrons

        vector<V4s> j25cats = jet_categories(jets, 1.0, tagger);
        V4s jls = j25cats[0];
        V4s jbs = j25cats[1];
        int njets = jls.size() + jbs.size();

        if (leps.size() == 1) {
          histmap["onelep_inclusive"].fill(weight, jls, jbs, leps);

          if (njets >= 6) {
            if (jbs.size() == 2) histmap["onelep_ge6j_eq2jb"].fill(weight, jls, jbs, leps);
            else if (jbs.size() == 3) histmap["onelep_ge6j_eq3jb"].fill(weight, jls, jbs, leps);
            else if (jbs.size() >= 4) histmap["onelep_ge6j_ge4jb"].fill(weight, jls, jbs, leps);
          }

        } else if (leps.size() == 2) {
          histmap["dilep_inclusive"].fill(weight, jls, jbs, leps);

          if (njets >= 4) {
            if (jbs.size() == 2) histmap["dilep_ge4j_eq2jb"].fill(weight, jls, jbs, leps);
            else if (jbs.size() == 3) histmap["dilep_ge4j_eq3jb"].fill(weight, jls, jbs, leps);
            else if (jbs.size() >= 4) histmap["dilep_ge4j_ge4jb"].fill(weight, jls, jbs, leps);
          }

        }

        return;
      }

      void finalize() {
        for (pair<const string, TTBBDecayedHists>& keyval: histmap) {
          TTBBDecayedHists& ttbbhists = keyval.second;

          for (shared_ptr<TTBBHist>& ttbbhist: ttbbhists.histograms()) {
            Histo1DPtr ph = ttbbhist->nominal();
            ph->setPath(histoDir() + ph->path());
            scale(ph, crossSection()/picobarn/sumOfWeights());
            addAnalysisObject(ph);

            Histo1DPtr neg = ttbbhist->negative();
            neg->setPath(histoDir() + neg->path());
            scale(neg, crossSection()/picobarn/sumOfWeights());
            addAnalysisObject(neg);

            Histo1DPtr pos = ttbbhist->positive();
            pos->setPath(histoDir() + pos->path());
            scale(pos, crossSection()/picobarn/sumOfWeights());
            addAnalysisObject(pos);

          }

          for (Profile1DPtr& ph: ttbbhists.profiles()) {
            ph->setPath(histoDir() + ph->path());
            ph->scaleW(crossSection()/picobarn/sumOfWeights());
            addAnalysisObject(ph);
          }
        }
      }

    private:
      map< string , TTBBDecayedHists > histmap;

      vector<V4s> jet_categories(
          const Jets& jets
        , float etscale
        , JET_BTAG_EFFS tagrates
        ) {

        V4s jls, jbs;
        for (const Jet& j: jets) {
          FourMomentum mom = j.mom()*etscale;

          if (efffilt(j, tagrates)) 
            jbs.push_back(mom);
          else
            jls.push_back(mom);
        }

        return {jls, jbs};
      }

  };

  DECLARE_RIVET_PLUGIN(MCTTBBDecayed);
}

