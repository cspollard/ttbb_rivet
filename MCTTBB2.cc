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

  string spt = "\\ensuremath{p_\\mathrm{T}}";
  string seta = "\\ensuremath{\\eta}";
  string sdphibb = "\\ensuremath{|\\Delta\\phi(b,b)|}";
  string sdrbb = "\\ensuremath{\\Delta R(b,b)}";
  string sptbb = "\\ensuremath{p_{\\mathrm{T}, bb}}";
  string smbb = "\\ensuremath{m_{bb}}";
  string sht = "\\ensuremath{h_\\mathrm{T}}";


  Histo1D histo1D(const string& name, size_t nb, double low, double high, const string& title, const string& xlab, const string& ylab) {
    Histo1D h = Histo1D(nb, low, high, name, title);
    h.setAnnotation("XLabel", xlab);
    h.setAnnotation("YLabel", ylab);
    return h;
  }


  class TTBBHists {
  public:

    TTBBHists() { };

    TTBBHists(const string& prefix) {
      h_njl = histo1D("h_" + prefix + "_njl", 10, 0, 10, "", "light-jet multiplicity", dsdx("n", "1"));
      h_njb = histo1D("h_" + prefix + "_njb", 5, 0, 5, "", "$b$-jet mulitplicity", dsdx("n", "1"));
      h_nj0b = histo1D("h_" + prefix + "_nj0b", 5, 0, 5, "", "$bb$-jet multiplicity", dsdx("n", "1"));
      h_nj1b = histo1D("h_" + prefix + "_nj1b", 5, 0, 5, "", "$b1$-jet multiplicity", dsdx("n", "1"));

      h_jl1pt = histo1D("h_" + prefix + "_jl1pt", 50, 0, 500*GeV, "", "leading light-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
      h_jl2pt = histo1D("h_" + prefix + "_jl2pt", 50, 0, 500*GeV, "", "subleading light-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
      h_jb1pt = histo1D("h_" + prefix + "_jb1pt", 50, 0, 500*GeV, "", "leading $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
      h_jb2pt = histo1D("h_" + prefix + "_jb2pt", 50, 0, 500*GeV, "", "subleading $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
      h_j0b1pt = histo1D("h_" + prefix + "_j0b1pt", 50, 0, 500*GeV, "", "leading $bb$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
      h_j0b2pt = histo1D("h_" + prefix + "_j0b2pt", 50, 0, 500*GeV, "", "subleading $bb$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
      h_j1b1pt = histo1D("h_" + prefix + "_j1b1pt", 50, 0, 500*GeV, "", "leading $b1$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
      h_j1b2pt = histo1D("h_" + prefix + "_j1b2pt", 50, 0, 500*GeV, "", "subleading $b1$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));

      h_jl1eta = histo1D("h_" + prefix + "_jl1eta", 30, -3, 3, "", "leading light-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
      h_jl2eta = histo1D("h_" + prefix + "_jl2eta", 30, -3, 3, "", "subleading light-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
      h_jb1eta = histo1D("h_" + prefix + "_jb1eta", 30, -3, 3, "", "leading $b$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
      h_jb2eta = histo1D("h_" + prefix + "_jb2eta", 30, -3, 3, "", "subleading $b$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
      h_j0b1eta = histo1D("h_" + prefix + "_j0b1eta", 30, -3, 3, "", "leading $bb$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
      h_j0b2eta = histo1D("h_" + prefix + "_j0b2eta", 30, -3, 3, "", "subleading $bb$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
      h_j1b1eta = histo1D("h_" + prefix + "_j1b1eta", 30, -3, 3, "", "leading $b1$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
      h_j1b2eta = histo1D("h_" + prefix + "_j1b2eta", 30, -3, 3, "", "subleading $b1$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));

      h_mbb = histo1D("h_" + prefix + "_mbb", 50, 0, 500*GeV, "", smbb + " [GeV]", dsdx(smbb, "\\mathrm{GeV}"));
      h_dphibb = histo1D("h_" + prefix + "_dphibb", 50, 0, 4, "", sdphibb, dsdx(sdphibb, "1"));
      h_drbb = histo1D("h_" + prefix + "_drbb", 50, 0, 5, "", sdrbb, dsdx(sdrbb, "1"));
      h_ptbb = histo1D("h_" + prefix + "_ptbb", 50, 0, 500*GeV, "", sptbb + " [GeV]", dsdx(sptbb, "\\mathrm{GeV}"));
      h_ht = histo1D("h_" + prefix + "_ht", 50, 0, 2000*GeV, "", sht + " [GeV]", dsdx(sht, "\\mathrm{GeV}"));

    }

    void fill(double weight, const Jets& jls, const Jets& jbs, const Jets& j0bs, const Jets& j1bs, double ht_other) {
        h_njl.fill(jls.size(), weight);
        h_njb.fill(j0bs.size()+j1bs.size(), weight);
        h_nj0b.fill(j0bs.size(), weight);
        h_nj1b.fill(j1bs.size(), weight);

        if (jls.size() >= 1) {
          h_jl1pt.fill(jls[0].pt(), weight);
          h_jl1eta.fill(jls[0].eta(), weight);
        }

        if (jls.size() >= 2) {
          h_jl2pt.fill(jls[1].pt(), weight);
          h_jl2eta.fill(jls[1].eta(), weight);
        }

        if (jbs.size() >= 1) {
          h_jb1pt.fill(jbs[0].pt(), weight);
          h_jb1eta.fill(jbs[0].eta(), weight);
        }

        if (jbs.size() >= 2) {
          h_jb2pt.fill(jbs[1].pt(), weight);
          h_jb2eta.fill(jbs[1].eta(), weight);
        }

        if (j0bs.size() >= 1) {
          h_j0b1pt.fill(j0bs[0].pt(), weight);
          h_j0b1eta.fill(j0bs[0].eta(), weight);
        }

        if (j0bs.size() >= 2) {
          h_j0b2pt.fill(j0bs[1].pt(), weight);
          h_j0b2eta.fill(j0bs[1].eta(), weight);
        }

        if (j1bs.size() >= 1) {
          h_j1b1pt.fill(j1bs[0].pt(), weight);
          h_j1b1eta.fill(j1bs[0].eta(), weight);
        }

        if (j1bs.size() >= 2) {
          h_j1b2pt.fill(j1bs[1].pt(), weight);
          h_j1b2eta.fill(j1bs[1].eta(), weight);
        }

        // we define the ht as the scalar sum of all the light, b, and B jet pts
        // plus the masses of all top quarks in the event.
        double ht = ht_other;
        for (const Jet& lj: jls)
          ht += lj.pt();

        for (const Jet& bj: jbs)
          ht += bj.pt();

        h_ht.fill(ht, weight);

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

        h_mbb.fill((b1 + b2).mass(), weight);
        h_dphibb.fill(abs(deltaPhi(b1, b2)), weight);
        h_drbb.fill(deltaR(b1, b2), weight);
        h_ptbb.fill((b1 + b2).pt(), weight);

      };

      const vector<Histo1D> histograms() const {
        return
          { h_njl, h_njb, h_nj0b, h_nj1b
          , h_jl1pt, h_jl2pt, h_jb1pt, h_jb2pt
          , h_j0b1pt, h_j0b2pt, h_j1b1pt, h_j1b2pt
          , h_jl1eta, h_jl2eta, h_jb1eta, h_jb2eta
          , h_j0b1eta, h_j0b2eta, h_j1b1eta, h_j1b2eta
          , h_mbb, h_dphibb, h_drbb, h_ptbb, h_ht
          };
      }

    private:

      Histo1D
          h_njl, h_njb, h_nj0b, h_nj1b
        , h_jl1pt, h_jl2pt, h_jb1pt, h_jb2pt
        , h_j0b1pt, h_j0b2pt, h_j1b1pt, h_j1b2pt
        , h_jl1eta, h_jl2eta, h_jb1eta, h_jb2eta
        , h_j0b1eta, h_j0b2eta, h_j1b1eta, h_j1b2eta
        , h_mbb, h_dphibb, h_drbb, h_ptbb, h_ht;

    };



  class MCTTBB2 : public Analysis {
  public:

    DEFAULT_RIVET_ANALYSIS_CTOR(MCTTBB2);

    void init() {
      FinalState fps(Cuts::abseta < 5 && Cuts::abspid != 6);
      declare(FinalState(Cuts::abspid==6), "Tops");

      declare(FastJets(fps, FastJets::ANTIKT, 0.4), "Jets");

      h_inclusive = TTBBHists("inclusive");
      h_zerob = TTBBHists("zerob");
      h_atleastoneb = TTBBHists("atleastoneb");
      h_onej1b = TTBBHists("onej1b");
      h_twoj1b = TTBBHists("twoj1b");
      h_onej0b = TTBBHists("onej0b");
      h_mbbgt100 = TTBBHists("mbbgt100");

    }


    void analyze(const Event& event) {

      double weight = event.weight();

      const Jets& jets =
      apply<FastJets>(event, "Jets").jetsByPt(
        Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
        const Particles& tquarks = apply<FinalState>(event, "Tops").particles();

        // find the light, b, and B jets by checking the number of b-quark
        // constituents in the jet.
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

        // we define the ht as the scalar sum of all the light, b, and B jet pts
        // plus the masses of all top quarks in the event.
        double ht = 0;
        for (const Particle& tq: tquarks)
          ht += tq.mass();

        h_inclusive.fill(weight, jls, jbs, j0bs, j1bs, ht);

        if (jbs.size() == 0) h_zerob.fill(weight, jls, jbs, j0bs, j1bs, ht);

        else h_atleastoneb.fill(weight, jls, jbs, j0bs, j1bs, ht);

        if (j1bs.size() == 1) h_onej1b.fill(weight, jls, jbs, j0bs, j1bs, ht);
        else if (j1bs.size() == 2) h_twoj1b.fill(weight, jls, jbs, j0bs, j1bs, ht);

        if (j0bs.size() == 1) h_onej0b.fill(weight, jls, jbs, j0bs, j1bs, ht);

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

        if ((b1+b2).mass() > 100*GeV) h_mbbgt100.fill(weight, jls, jbs, j0bs, j1bs, ht);


        return;
      }


        void finalize() {
          vector<TTBBHists> hists =
            {h_inclusive, h_zerob, h_atleastoneb, h_onej0b, h_onej1b, h_twoj1b, h_mbbgt100};

          for (const TTBBHists& hist: hists) {
            for (const Histo1D& h: hist.histograms()) {
              Histo1DPtr ph = make_shared<Histo1D>(h);
              ph->setPath(histoDir() + ph->path());
              scale(ph, crossSection()/picobarn/sumOfWeights());
              addAnalysisObject(ph);
          }
        }
      }

    private:
      TTBBHists h_inclusive, h_zerob, h_atleastoneb, h_onej1b, h_twoj1b, h_onej0b, h_mbbgt100;

    };

    DECLARE_RIVET_PLUGIN(MCTTBB2);

  }
