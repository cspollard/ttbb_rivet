// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FinalPartons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"

namespace Rivet {


  class MCTTBB2 : public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(MCTTBB2);

    void init() {
      FinalPartons fps(Cuts::abseta < 5 && Cuts::abspid != 6);
      declare(PartonicTops(), "Tops");

      declare(FastJets(fps, FastJets::ANTIKT, 0.4), "Jets");

      h_nlj = bookHisto1D("h_nlj", 10, 0, 10);
      h_nbj = bookHisto1D("h_nbj", 5, 0, 5);
      h_nBj = bookHisto1D("h_nBj", 5, 0, 5);

      h_ptl1 = bookHisto1D("h_ptl1", 50, 0, 500*GeV);
      h_ptl2 = bookHisto1D("h_ptl2", 50, 0, 500*GeV);
      h_ptb1 = bookHisto1D("h_ptb1", 50, 0, 500*GeV);
      h_ptb2 = bookHisto1D("h_ptb2", 50, 0, 500*GeV);
      h_ptB1 = bookHisto1D("h_ptB1", 50, 0, 500*GeV);
      h_ptB2 = bookHisto1D("h_ptB2", 50, 0, 500*GeV);

      h_etal1 = bookHisto1D("h_etal1", 50, -2.5, 2.5);
      h_etal2 = bookHisto1D("h_etal2", 50, -2.5, 2.5);
      h_etab1 = bookHisto1D("h_etab1", 50, -2.5, 2.5);
      h_etab2 = bookHisto1D("h_etab2", 50, -2.5, 2.5);
      h_etaB1 = bookHisto1D("h_etaB1", 50, -2.5, 2.5);
      h_etaB2 = bookHisto1D("h_etaB2", 50, -2.5, 2.5);

      h_mbb = bookHisto1D("h_mbb", 50, 0, 500*GeV);
      h_dphibb = bookHisto1D("h_dphibb", 50, -4, 4);
      h_drbb = bookHisto1D("h_drbb", 50, -5, 5);
      h_ptbb = bookHisto1D("h_ptbb", 50, 0, 500*GeV);
      h_ht = bookHisto1D("h_ht", 50, 0, 1000*GeV);

      h_mbb_mbbgt100 = bookHisto1D("h_mbb_mbbgt100", 50, 0, 500*GeV);
      h_dphibb_mbbgt100 = bookHisto1D("h_dphibb_mbbgt100", 50, -4, 4);
      h_drbb_mbbgt100 = bookHisto1D("h_drbb_mbbgt100", 50, -5, 5);
      h_ptbb_mbbgt100 = bookHisto1D("h_ptbb_mbbgt100", 50, 0, 500*GeV);
      h_ht_mbbgt100 = bookHisto1D("h_ht_mbbgt100", 50, 0, 1000*GeV);

      hists =
        { h_nlj, h_nbj, h_nBj, h_ptl1, h_ptl2, h_ptb1, h_ptb2, h_ptB1, h_ptB2
        , h_etal1, h_etal2, h_etab1, h_etab2, h_etaB1, h_etaB2, h_mbb, h_dphibb
        , h_drbb, h_ptbb, h_ht, h_mbb_mbbgt100, h_dphibb_mbbgt100
        , h_drbb_mbbgt100, h_ptbb_mbbgt100, h_ht_mbbgt100
        };
    }


    void analyze(const Event& event) {

      double weight = event.weight();

      const Jets& jets =
        apply<FastJets>(event, "Jets").jetsByPt(
          (Cuts::mass > 150*GeV || Cuts::pT > 25*GeV) && Cuts::abseta < 2.5);
      const Particles& tquarks = apply<PartonicTops>(event, "Tops").tops();

      // find the light, b, and B jets by checking the number of b-quark
      // constituents in the jet.
      Jets ljs, bjs, Bjs;
      for (const Jet& j: jets) {

        size_t nbs = 0;
        for (const Particle& c: j.constituents()) {
          if (abspid(c) == 5)
            nbs++;
        }

        if (nbs == 0)
          ljs.push_back(j);
        else if (nbs == 1)
          bjs.push_back(j);
        else
          Bjs.push_back(j);
      }


      // fill the jet multiplicity histograms
      h_nlj->fill(ljs.size());
      h_nbj->fill(bjs.size());
      h_nBj->fill(Bjs.size());

      if (ljs.size() >= 1) {
        h_ptl1->fill(ljs[0].pt(), weight);
        h_etal1->fill(ljs[0].eta(), weight);
      }

      if (ljs.size() >= 2) {
        h_ptl2->fill(ljs[1].pt(), weight);
        h_etal2->fill(ljs[1].eta(), weight);
      }

      if (bjs.size() >= 1) {
        h_ptb1->fill(bjs[0].pt(), weight);
        h_etab1->fill(bjs[0].eta(), weight);
      }

      if (bjs.size() >= 2) {
        h_ptb2->fill(bjs[1].pt(), weight);
        h_etab2->fill(bjs[1].eta(), weight);
      }

      if (Bjs.size() >= 1) {
        h_ptB1->fill(Bjs[0].pt(), weight);
        h_etaB1->fill(Bjs[0].eta(), weight);
      }

      if (Bjs.size() >= 2) {
        h_ptB2->fill(Bjs[1].pt(), weight);
        h_etaB2->fill(Bjs[1].eta(), weight);
      }

      // we define the ht as the scalar sum of all the light, b, and B jet pts
      // plus the masses of all top quarks in the event.
      double ht = 0;
      for (const Particle& tq: tquarks)
        ht += tq.mass();

      for (const Jet& lj: ljs)
        ht += lj.pt();

      for (const Jet& bj: bjs)
        ht += bj.pt();

      for (const Jet& Bj: Bjs)
        ht += Bj.pt();

      h_ht->fill(ht, weight);

      // find the two leading "b-jets" in the event, where "b-jet" here means
      // the two leading jets with at least one b-quark constituent, with jets
      // with exactly one b-quark constituent taking priority.
      FourMomentum b1, b2;
      if (bjs.size() >= 2) {
        b1 = bjs[0].mom();
        b2 = bjs[1].mom();
      } else if (bjs.size() >= 1 && Bjs.size() >= 1) {
        b1 = bjs[0].mom();
        b2 = Bjs[0].mom();
      } else if (Bjs.size() >= 2) {
        b1 = Bjs[0].mom();
        b2 = Bjs[1].mom();
      } else
        return;

      h_mbb->fill((b1 + b2).mass(), weight);
      h_dphibb->fill(deltaPhi(b1, b2), weight);
      h_drbb->fill(deltaR(b1, b2), weight);

      if ((b1 + b2).mass() < 100*GeV)
        return;

      h_ht_mbbgt100->fill(deltaR(b1, b2), weight);
      h_mbb_mbbgt100->fill((b1 + b2).mass(), weight);
      h_dphibb_mbbgt100->fill(deltaPhi(b1, b2), weight);
      h_drbb_mbbgt100->fill(deltaR(b1, b2), weight);

      return;
    }


      void finalize() {
        for (Histo1DPtr& h: hists)
          scale(h, crossSection()/picobarn/sumOfWeights());
      }

    private:
      vector<Histo1DPtr> hists;
      Histo1DPtr h_nlj;
      Histo1DPtr h_nbj;
      Histo1DPtr h_nBj;

      Histo1DPtr h_ptl1;
      Histo1DPtr h_ptl2;
      Histo1DPtr h_ptb1;
      Histo1DPtr h_ptb2;
      Histo1DPtr h_ptB1;
      Histo1DPtr h_ptB2;

      Histo1DPtr h_etal1;
      Histo1DPtr h_etal2;
      Histo1DPtr h_etab1;
      Histo1DPtr h_etab2;
      Histo1DPtr h_etaB1;
      Histo1DPtr h_etaB2;

      Histo1DPtr h_mbb;
      Histo1DPtr h_dphibb;
      Histo1DPtr h_drbb;
      Histo1DPtr h_ptbb;
      Histo1DPtr h_ht;

      Histo1DPtr h_mbb_mbbgt100;
      Histo1DPtr h_dphibb_mbbgt100;
      Histo1DPtr h_drbb_mbbgt100;
      Histo1DPtr h_ptbb_mbbgt100;
      Histo1DPtr h_ht_mbbgt100;
    };

    DECLARE_RIVET_PLUGIN(MCTTBB2);

  }
