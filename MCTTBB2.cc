// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FinalPartons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"

namespace Rivet {


  string dsdx(const string& x, const string& xunit) {
    return "\\ensuremath{\\frac{d\\sigma}{d" + x + "} \\frac{{pb}}{" + xunit + "}}";
  }


  class MCTTBB2 : public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(MCTTBB2);

    void init() {
      FinalPartons fps(Cuts::abseta < 5 && Cuts::abspid != 6);
      declare(PartonicTops(), "Tops");

      string pt = "\\ensuremath{p_\\mathrm{T}}";
      string eta = "\\ensuremath{\\eta}";
      string dphibb = "\\ensuremath{\\Delta\\phi(b,b)}";
      string drbb = "\\ensuremath{\\Delta R(b,b)}";
      string ptbb = "\\ensuremath{p_{\\mathrm{T}, bb}}";
      string mbb = "\\ensuremath{m_{bb}}";
      string ht = "\\ensuremath{h_\\mathrm{T}}";

      declare(FastJets(fps, FastJets::ANTIKT, 0.4), "Jets");

      h_nlj = bookHisto1D("h_nlj", 10, 0, 10, "", "light-jet multiplicity", dsdx("n", "1"));
      h_nbj = bookHisto1D("h_nbj", 5, 0, 5, "", "$b$-jet mulitplicity", dsdx("n", "1"));
      h_nBj = bookHisto1D("h_nBj", 5, 0, 5, "", "$bb$-jet multiplicity", dsdx("n", "1"));

      h_ptl1 = bookHisto1D("h_ptl1", 50, 0, 500*GeV, "", "leading light-jet " + pt + " [GeV]", dsdx(pt, "\\mathrm{GeV}"));

      h_ptl2 = bookHisto1D("h_ptl2", 50, 0, 500*GeV, "", "subleading light-jet " + pt + " [GeV]", dsdx(pt, "\\mathrm{GeV}"));

      h_ptb1 = bookHisto1D("h_ptb1", 50, 0, 500*GeV, "", "leading $b$-jet " + pt + " [GeV]", dsdx(pt, "\\mathrm{GeV}"));

      h_ptb2 = bookHisto1D("h_ptb2", 50, 0, 500*GeV, "", "subleading $b$-jet " + pt + " [GeV]", dsdx(pt, "\\mathrm{GeV}"));

      h_ptB1 = bookHisto1D("h_ptB1", 50, 0, 500*GeV, "", "leading $bb$-jet " + pt + " [GeV]", dsdx(pt, "\\mathrm{GeV}"));

      h_ptB2 = bookHisto1D("h_ptB2", 50, 0, 500*GeV, "", "subleading $bb$-jet " + pt + " [GeV]", dsdx(pt, "\\mathrm{GeV}"));


      h_etal1 = bookHisto1D("h_etal1", 30, -3, 3, "", "leading light-jet " + eta, dsdx(eta, "\\mathrm{GeV}"));

      h_etal2 = bookHisto1D("h_etal2", 30, -3, 3, "", "subleading light-jet " + eta, dsdx(eta, "\\mathrm{GeV}"));

      h_etab1 = bookHisto1D("h_etab1", 30, -3, 3, "", "leading $b$-jet " + eta, dsdx(eta, "\\mathrm{GeV}"));

      h_etab2 = bookHisto1D("h_etab2", 30, -3, 3, "", "subleading $b$-jet " + eta, dsdx(eta, "\\mathrm{GeV}"));

      h_etaB1 = bookHisto1D("h_etaB1", 30, -3, 3, "", "leading $bb$-jet " + eta, dsdx(eta, "\\mathrm{GeV}"));

      h_etaB2 = bookHisto1D("h_etaB2", 30, -3, 3, "", "subleading $bb$-jet " + eta, dsdx(eta, "\\mathrm{GeV}"));

      h_mbb = bookHisto1D("h_mbb", 50, 0, 500*GeV, "", mbb + " [GeV]", dsdx(mbb, "\\mathrm{GeV}"));
      h_dphibb = bookHisto1D("h_dphibb", 50, 0, 4, "", dphibb, dsdx(dphibb, "1"));
      h_drbb = bookHisto1D("h_drbb", 50, 0, 5, "", drbb, dsdx(drbb, "1"));
      h_ptbb = bookHisto1D("h_ptbb", 50, 0, 500*GeV, "", ptbb + " [GeV]", dsdx(ptbb, "\\mathrm{GeV}"));
      h_ht = bookHisto1D("h_ht", 50, 0, 2000*GeV, "", ht + " [GeV]", dsdx(ht, "\\mathrm{GeV}"));

      h_mbbgt100_mbb = bookHisto1D("h_mbbgt100_mbb", 50, 0, 500*GeV, "", mbb + " [GeV]", dsdx(mbb, "\\mathrm{GeV}"));
      h_mbbgt100_dphibb = bookHisto1D("h_mbbgt100_dphibb", 50, 0, 4, "", dphibb, dsdx(dphibb, "1"));
      h_mbbgt100_drbb = bookHisto1D("h_mbbgt100_drbb", 50, 0, 5, "", drbb, dsdx(drbb, "1"));
      h_mbbgt100_ptbb = bookHisto1D("h_mbbgt100_ptbb", 50, 0, 500*GeV, "", ptbb + " [GeV]", dsdx(ptbb, "\\mathrm{GeV}"));
      h_mbbgt100_ht = bookHisto1D("h_mbbgt100_ht", 50, 0, 2000*GeV, "", ht + " [GeV]", dsdx(ht, "\\mathrm{GeV}"));

      hists =
        { h_nlj, h_nbj, h_nBj, h_ptl1, h_ptl2, h_ptb1, h_ptb2, h_ptB1, h_ptB2
        , h_etal1, h_etal2, h_etab1, h_etab2, h_etaB1, h_etaB2, h_mbb, h_dphibb
        , h_drbb, h_ptbb, h_ht, h_mbbgt100_mbb, h_mbbgt100_dphibb
        , h_mbbgt100_drbb, h_mbbgt100_ptbb, h_mbbgt100_ht
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
      h_nlj->fill(ljs.size(), weight);
      h_nbj->fill(bjs.size(), weight);
      h_nBj->fill(Bjs.size(), weight);

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
      h_ptbb->fill((b1 + b2).pt(), weight);

      if ((b1 + b2).mass() < 100*GeV)
        return;

      h_mbbgt100_ht->fill(ht, weight);
      h_mbbgt100_mbb->fill((b1 + b2).mass(), weight);
      h_mbbgt100_dphibb->fill(deltaPhi(b1, b2), weight);
      h_mbbgt100_drbb->fill(deltaR(b1, b2), weight);
      h_mbbgt100_ptbb->fill((b1 + b2).pt(), weight);

      return;
    }


      void finalize() {
        for (Histo1DPtr& h: hists) {
          scale(h, crossSection()/picobarn/sumOfWeights());
          Histo1DPtr hnorm = make_shared<YODA::Histo1D>(*h);
          hnorm->normalize();
          hnorm->setPath(h->path() + "_norm");
          hnorm->setAnnotation("YLabel", "\\ensuremath{\\frac{1}{\\sigma}}" + hnorm->annotation("YLabel"));
          addAnalysisObject(hnorm);
        }
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

      Histo1DPtr h_mbbgt100_mbb;
      Histo1DPtr h_mbbgt100_dphibb;
      Histo1DPtr h_mbbgt100_drbb;
      Histo1DPtr h_mbbgt100_ptbb;
      Histo1DPtr h_mbbgt100_ht;
    };

    DECLARE_RIVET_PLUGIN(MCTTBB2);

  }
