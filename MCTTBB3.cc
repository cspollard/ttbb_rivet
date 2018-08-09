// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FinalPartons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "TTBBHist.hh"

namespace Rivet {

  class TTBBHists {
    public:

      TTBBHists() { };

      TTBBHists(const string& prefix) {
        h_njl = make_shared<TTBBHist>("h_" + prefix + "_njl", 10, -0.5, 9.5, "", "light-jet multiplicity", dsdx("n", "1"));
        h_njb = make_shared<TTBBHist>("h_" + prefix + "_njb", 5, -0.5, 4.5, "", "$b$-jet mulitplicity", dsdx("n", "1"));
        h_nj0b = make_shared<TTBBHist>("h_" + prefix + "_nj0b", 5, -0.5, 4.5, "", "$bb$-jet multiplicity", dsdx("n", "1"));
        h_nj1b = make_shared<TTBBHist>("h_" + prefix + "_nj1b", 5, -0.5, 4.5, "", "$b1$-jet multiplicity", dsdx("n", "1"));

        h_jl1pt = make_shared<TTBBHist>("h_" + prefix + "_jl1pt", 50, 0, 500*GeV, "", "leading light-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_jl2pt = make_shared<TTBBHist>("h_" + prefix + "_jl2pt", 50, 0, 500*GeV, "", "subleading light-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_jb1pt = make_shared<TTBBHist>("h_" + prefix + "_jb1pt", 50, 0, 500*GeV, "", "leading $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_jb2pt = make_shared<TTBBHist>("h_" + prefix + "_jb2pt", 50, 0, 500*GeV, "", "subleading $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_j0b1pt = make_shared<TTBBHist>("h_" + prefix + "_j0b1pt", 50, 0, 500*GeV, "", "leading $bb$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_j0b2pt = make_shared<TTBBHist>("h_" + prefix + "_j0b2pt", 50, 0, 500*GeV, "", "subleading $bb$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_j1b1pt = make_shared<TTBBHist>("h_" + prefix + "_j1b1pt", 50, 0, 500*GeV, "", "leading $b1$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        h_j1b2pt = make_shared<TTBBHist>("h_" + prefix + "_j1b2pt", 50, 0, 500*GeV, "", "subleading $b1$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));

        h_jl1eta = make_shared<TTBBHist>("h_" + prefix + "_jl1eta", 30, -3, 3, "", "leading light-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        h_jl2eta = make_shared<TTBBHist>("h_" + prefix + "_jl2eta", 30, -3, 3, "", "subleading light-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        h_jb1eta = make_shared<TTBBHist>("h_" + prefix + "_jb1eta", 30, -3, 3, "", "leading $b$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        h_jb2eta = make_shared<TTBBHist>("h_" + prefix + "_jb2eta", 30, -3, 3, "", "subleading $b$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        h_j0b1eta = make_shared<TTBBHist>("h_" + prefix + "_j0b1eta", 30, -3, 3, "", "leading $bb$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        h_j0b2eta = make_shared<TTBBHist>("h_" + prefix + "_j0b2eta", 30, -3, 3, "", "subleading $bb$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        h_j1b1eta = make_shared<TTBBHist>("h_" + prefix + "_j1b1eta", 30, -3, 3, "", "leading $b1$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        h_j1b2eta = make_shared<TTBBHist>("h_" + prefix + "_j1b2eta", 30, -3, 3, "", "subleading $b1$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));

        h_dphitb = make_shared<TTBBHist>("h_" + prefix + "_dphitb", 50, 0, 4, "", sdphitb, dsdx(sdphitb, "\\mathrm{1}"));
        h_detatb = make_shared<TTBBHist>("h_" + prefix + "_detatb", 50, 0, 4, "", sdetatb, dsdx(sdetatb, "\\mathrm{1}"));
        h_drtb = make_shared<TTBBHist>("h_" + prefix + "_drtb", 50, 0, 4, "", sdrtb, dsdx(sdrtb, "\\mathrm{1}"));

        h_mbb = make_shared<TTBBHist>("h_" + prefix + "_mbb", 50, 0, 500*GeV, "", smbb + " [GeV]", dsdx(smbb, "\\mathrm{GeV}"));
        h_dphibb = make_shared<TTBBHist>("h_" + prefix + "_dphibb", 50, 0, 4, "", sdphibb, dsdx(sdphibb, "1"));
        h_drbb = make_shared<TTBBHist>("h_" + prefix + "_drbb", 50, 0, 5, "", sdrbb, dsdx(sdrbb, "1"));
        h_ptbb = make_shared<TTBBHist>("h_" + prefix + "_ptbb", 50, 0, 500*GeV, "", sptbb + " [GeV]", dsdx(sptbb, "\\mathrm{GeV}"));
        h_ht = make_shared<TTBBHist>("h_" + prefix + "_ht", 50, 0, 2000*GeV, "", sht + " [GeV]", dsdx(sht, "\\mathrm{GeV}"));
      }

      void fill(double weight, const Jets& jls, const Jets& jbs, const Jets& j0bs, const Jets& j1bs, const Particles& tops) {
        h_njl->fill(jls.size(), weight);
        h_njb->fill(j0bs.size()+j1bs.size(), weight);
        h_nj0b->fill(j0bs.size(), weight);
        h_nj1b->fill(j1bs.size(), weight);

        if (jls.size() >= 1) {
          h_jl1pt->fill(jls[0].pt(), weight);
          h_jl1eta->fill(jls[0].eta(), weight);
        }

        if (jls.size() >= 2) {
          h_jl2pt->fill(jls[1].pt(), weight);
          h_jl2eta->fill(jls[1].eta(), weight);
        }

        if (jbs.size() >= 1) {
          h_jb1pt->fill(jbs[0].pt(), weight);
          h_jb1eta->fill(jbs[0].eta(), weight);
        }

        if (jbs.size() >= 2) {
          h_jb2pt->fill(jbs[1].pt(), weight);
          h_jb2eta->fill(jbs[1].eta(), weight);
        }

        if (j0bs.size() >= 1) {
          h_j0b1pt->fill(j0bs[0].pt(), weight);
          h_j0b1eta->fill(j0bs[0].eta(), weight);
        }

        if (j0bs.size() >= 2) {
          h_j0b2pt->fill(j0bs[1].pt(), weight);
          h_j0b2eta->fill(j0bs[1].eta(), weight);
        }

        if (j1bs.size() >= 1) {
          h_j1b1pt->fill(j1bs[0].pt(), weight);
          h_j1b1eta->fill(j1bs[0].eta(), weight);
        }

        if (j1bs.size() >= 2) {
          h_j1b2pt->fill(j1bs[1].pt(), weight);
          h_j1b2eta->fill(j1bs[1].eta(), weight);
        }

        // we define the ht as the scalar sum of all the light, b, and B jet pts
        // plus the masses of all top quarks in the event.
        double ht = 0.0;
        for (const Particle& t: tops)
          ht += t.mass();

        for (const Jet& lj: jls)
          ht += lj.pt();

        for (const Jet& bj: jbs)
          ht += bj.pt();

        h_ht->fill(ht, weight);

        // find the two leading "b-jets" in the event, where "b-jet" here means
        // the two leading jets with at least one b-quark constituent, with jets
        // with exactly one b-quark constituent taking priority.
        FourMomentum b1, b2;
        if (j1bs.size() >= 2) {
          b1 = j1bs[0].mom();
          b2 = j1bs[1].mom();
        } else if (j1bs.size() >= 1 && j0bs.size() >= 1) {
          b1 = j1bs[0].mom();
          b2 = j0bs[0].mom();
        } else if (j0bs.size() >= 2) {
          b1 = j0bs[0].mom();
          b2 = j0bs[1].mom();
        } else
          return;

        // find the closest (in dR) top to a jet with at least one b
        // and store relative angles
        double mindr = -1;
        double mindrdeta = -1;
        double mindrdphi = -1;
        for (const Jet& jb: jbs) {
          for (const Particle& t: tops) {
            double dr = deltaR(jb, t);
            if (mindr < 1 || dr < mindr) {
              mindr = dr;
              mindrdphi = abs(deltaPhi(jb, t));
              mindrdeta = abs(deltaEta(jb, t));
            }
          }
        }

        if (mindr >= 0) {
          h_drtb->fill(mindr, weight);
          h_dphitb->fill(mindrdphi, weight);
          h_detatb->fill(mindrdeta, weight);
        }

        h_mbb->fill((b1 + b2).mass(), weight);
        h_dphibb->fill(abs(deltaPhi(b1, b2)), weight);
        h_drbb->fill(deltaR(b1, b2), weight);
        h_ptbb->fill((b1 + b2).pt(), weight);

      };

      // this has to be a vector<Histo1DPtr> rather than vector<Histo1D> because
      // the Histo1D copy constructor loses all annotations?!?!?!?
      vector<shared_ptr<TTBBHist>> histograms() {
        return
        { h_njl, h_njb, h_nj0b, h_nj1b
          , h_jl1pt, h_jl2pt, h_jb1pt, h_jb2pt
            , h_j0b1pt, h_j0b2pt, h_j1b1pt, h_j1b2pt
            , h_jl1eta, h_jl2eta, h_jb1eta, h_jb2eta
            , h_j0b1eta, h_j0b2eta, h_j1b1eta, h_j1b2eta
            , h_dphitb, h_detatb, h_drtb
            , h_mbb, h_dphibb, h_drbb, h_ptbb, h_ht
        };
      }

    private:
      string dsdx(const string& x, const string& xunit) {
        return "\\ensuremath{\\frac{d\\sigma}{d" + x + "} \\frac{{pb}}{" + xunit + "}}";
      }


      shared_ptr<TTBBHist>
        h_njl, h_njb, h_nj0b, h_nj1b
        , h_jl1pt, h_jl2pt, h_jb1pt, h_jb2pt
        , h_j0b1pt, h_j0b2pt, h_j1b1pt, h_j1b2pt
        , h_jl1eta, h_jl2eta, h_jb1eta, h_jb2eta
        , h_j0b1eta, h_j0b2eta, h_j1b1eta, h_j1b2eta
        , h_dphitb, h_detatb, h_drtb
        , h_mbb, h_dphibb, h_drbb, h_ptbb, h_ht;

  };



  class MCTTBB3 : public Analysis {
    public:

      DEFAULT_RIVET_ANALYSIS_CTOR(MCTTBB3);

      void init() {
        FinalState fps(Cuts::abseta < 5 && Cuts::abspid != 6);
        declare(FinalState(Cuts::abspid==6), "Tops");

        declare(FastJets(fps, FastJets::ANTIKT, 0.4), "Jets");

        h_inclusive = TTBBHists("inclusive");
        h_inclusive_j15 = TTBBHists("inclusive_j15");
        h_inclusive_j30 = TTBBHists("inclusive_j30");
        h_inclusive_j40 = TTBBHists("inclusive_j40");

        h_zerob = TTBBHists("zerob");
        h_atleastoneb = TTBBHists("atleastoneb");
        h_onej1b = TTBBHists("onej1b");
        h_twoj1b = TTBBHists("twoj1b");
        h_atleasttwoj1b = TTBBHists("atleasttwoj1b");
        h_threej1b = TTBBHists("threej1b");
        h_fourj1b = TTBBHists("fourj1b");
        h_onej0b = TTBBHists("onej0b");
        h_mbbgt100 = TTBBHists("mbbgt100");

      }


      void analyze(const Event& event) {

        double weight = event.weight();

        const Particles& tquarks = apply<FinalState>(event, "Tops").particles();

        const Jets& jets =
          apply<FastJets>(event, "Jets").jetsByPt(
              Cuts::pT > 25*GeV && Cuts::abseta < 2.5);

        const Jets& jets15 =
          apply<FastJets>(event, "Jets").jetsByPt(
              Cuts::pT > 15*GeV && Cuts::abseta < 2.5);

        const Jets& jets30 =
          apply<FastJets>(event, "Jets").jetsByPt(
              Cuts::pT > 30*GeV && Cuts::abseta < 2.5);

        const Jets& jets40 =
          apply<FastJets>(event, "Jets").jetsByPt(
              Cuts::pT > 40*GeV && Cuts::abseta < 2.5);



        vector<Jets> j25cats = jet_categories(jets);
        vector<Jets> j15cats = jet_categories(jets15);
        vector<Jets> j30cats = jet_categories(jets30);
        vector<Jets> j40cats = jet_categories(jets40);

        // find the light, b, and B jets by checking the number of b-quark
        // constituents in the jet.
        // jls = light jets
        // jbs = jets with at least one b-quark inside
        // j0bs = jets with at least one b-quark, but a corresponding
        //        anti-b-quark for each b-quark
        // j1bs = jets with an imbalance between b-quarks and anti-b-quarks

        Jets jls = j25cats[0];
        Jets jbs = j25cats[1];
        Jets j0bs = j25cats[2];
        Jets j1bs = j25cats[3];

        h_inclusive_j15.fill(weight, j15cats[0], j15cats[1], j15cats[2], j15cats[3], tquarks);
        h_inclusive_j30.fill(weight, j30cats[0], j30cats[1], j30cats[2], j30cats[3], tquarks);
        h_inclusive_j40.fill(weight, j40cats[0], j40cats[1], j40cats[2], j40cats[3], tquarks);

        // inclusive selection
        h_inclusive.fill(weight, jls, jbs, j0bs, j1bs, tquarks);

        // no accepted jets with b-quarks inside
        if (jbs.size() == 0) h_zerob.fill(weight, jls, jbs, j0bs, j1bs, tquarks);
        // at least one jet with b-quarks inside
        else h_atleastoneb.fill(weight, jls, jbs, j0bs, j1bs, tquarks);

        // exactly one jet with an imbalance of b and anti-b
        if (j1bs.size() == 1) h_onej1b.fill(weight, jls, jbs, j0bs, j1bs, tquarks);
        // exactly two jets with an imbalance of b and anti-b
        else if (j1bs.size() == 2) h_twoj1b.fill(weight, jls, jbs, j0bs, j1bs, tquarks);
        // exactly three jets with an imbalance of b and anti-b
        else if (j1bs.size() == 3) h_threej1b.fill(weight, jls, jbs, j0bs, j1bs, tquarks);
        // exactly four jets with an imbalance of b and anti-b
        else if (j1bs.size() == 4) h_fourj1b.fill(weight, jls, jbs, j0bs, j1bs, tquarks);

        // at least two jets with an imbalance of b and anti-b
        if (j1bs.size() >= 2) h_atleasttwoj1b.fill(weight, jls, jbs, j0bs, j1bs, tquarks);

        // exactly one jet with balance of b and anti-b
        if (j0bs.size() == 1) h_onej0b.fill(weight, jls, jbs, j0bs, j1bs, tquarks);

        // find the two leading "b-jets" in the event, where "b-jet" here means
        // the two leading jets with at least one b-quark constituent, with jets
        // with exactly one b-quark constituent taking priority.
        FourMomentum b1, b2;
        if (j1bs.size() >= 2) {
          b1 = j1bs[0].mom();
          b2 = j1bs[1].mom();
        } else if (j1bs.size() >= 1 && j0bs.size() >= 1) {
          b1 = j1bs[0].mom();
          b2 = j0bs[0].mom();
        } else if (j0bs.size() >= 2) {
          b1 = j0bs[0].mom();
          b2 = j0bs[1].mom();
        } else
          return;

        if ((b1+b2).mass() > 100*GeV) h_mbbgt100.fill(weight, jls, jbs, j0bs, j1bs, tquarks);

        return;
      }

      void finalize() {
        vector<TTBBHists> hists =
        { h_inclusive, h_inclusive_j15, h_inclusive_j30, h_inclusive_j40
          , h_zerob, h_atleastoneb, h_onej0b
          , h_onej1b, h_twoj1b, h_atleasttwoj1b, h_threej1b
          , h_fourj1b, h_mbbgt100
        };

        for (TTBBHists& ttbbhists: hists) {
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

            /*
               Scatter2DPtr ps = ttbbhist->posfrac();
               ps->setPath(histoDir() + ps->path());
               addAnalysisObject(ps);
               */
          }
        }
      }

    private:
      TTBBHists h_inclusive, h_inclusive_j15, h_inclusive_j30, h_inclusive_j40,
                h_zerob, h_atleastoneb, h_onej1b, h_twoj1b, h_atleasttwoj1b, h_threej1b, h_fourj1b, h_onej0b, h_mbbgt100;

      vector<Jets> jet_categories(const Jets& jets) {
        Jets jls, jbs, j0bs, j1bs;
        for (const Jet& j: jets) {

          size_t absnbs = 0;
          int nbs = 0;
          for (const Particle& c: j.constituents()) {
            int thispid = pid(c);
            if (abs(thispid) == 5) {
              absnbs++;

              if (thispid > 0) nbs++;

              else nbs--;
            }
          }

          if (absnbs == 0) jls.push_back(j);

          else {
            jbs.push_back(j);
            if (nbs == 0) j0bs.push_back(j);
            else j1bs.push_back(j);
          }
        }

        return {jls, jbs, j0bs, j1bs};
      }


  };

  DECLARE_RIVET_PLUGIN(MCTTBB3);

}
