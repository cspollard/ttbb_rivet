// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FinalPartons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"

namespace Rivet {

  using V4s = vector<FourMomentum>;
  using Region = map< string , Histo1DPtr >;
  using RegMap = map< string , Region >;

  using StringMap = map< string , string >;

  string spt = "\\ensuremath{p_\\mathrm{T}}";
  string seta = "\\ensuremath{\\eta}";
  string sdphitb = "\\ensuremath{|\\Delta\\phi(t,b)|}";
  string sdetatb = "\\ensuremath{|\\Delta\\eta(t,b)|}";
  string sdrtb = "\\ensuremath{|\\Delta R(t,b)|}";
  string sdphibb = "\\ensuremath{|\\Delta\\phi(b,b)|}";
  string sdrbb = "\\ensuremath{\\Delta R(b,b)}";
  string sptbb = "\\ensuremath{p_{\\mathrm{T}, bb}}";
  string smbb = "\\ensuremath{m_{bb}}";
  string sht = "\\ensuremath{h_\\mathrm{T}}";
  string shtlep = "\\ensuremath{h_\\mathrm{T} ^\\mathrm{lep}}";

  string dsdx(const string& x, const string& xunit) {
    return "\\ensuremath{\\frac{d\\sigma}{d" + x + "} \\frac{{pb}}{" + xunit + "}}";
  }


  class MCTTBB3 : public Analysis {
    private:

      RegMap regions;

      StringMap xlabs, ylabs, titles;

      Histo1DPtr histo1D(
            const string& name, size_t nb, float low, float high
          , const string& title, const string& xlab, const string& ylab) {
      
        Histo1DPtr h;
        h = book(h, name, nb, low, high);

        // unfortunately I can no longer set annotations here
        // and have to store them for finalize :@

        xlabs[name] = xlab;
        ylabs[name] = ylab;
        titles[name] = title;
        return h;
      }


      Histo1DPtr histo1D(
            const string& name, const vector<double>& binning
          , const string& title, const string& xlab, const string& ylab) {

        Histo1DPtr h;
        h = book(h, name, binning);

        xlabs[name] = xlab;
        ylabs[name] = ylab;
        titles[name] = title;
        return h;
      }

      Region initRegion(const string& prefix) {
        Region hm;

        hm["njl"] = histo1D(prefix + "_njl", 10, -0.5, 9.5, "", "light-jet multiplicity", dsdx("n", "1"));
        hm["njb"] = histo1D(prefix + "_njb", 5, -0.5, 4.5, "", "$b$-jet mulitplicity", dsdx("n", "1"));
        hm["nj0b"] = histo1D(prefix + "_nj0b", 5, -0.5, 4.5, "", "$bb$-jet multiplicity", dsdx("n", "1"));
        hm["nj1b"] = histo1D(prefix + "_nj1b", 5, -0.5, 4.5, "", "$b1$-jet multiplicity", dsdx("n", "1"));

        hm["jl1pt"] = histo1D(prefix + "_jl1pt", 50, 0, 500*GeV, "", "leading light-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["jl2pt"] = histo1D(prefix + "_jl2pt", 50, 0, 500*GeV, "", "subleading light-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["jb1pt"] = histo1D(prefix + "_jb1pt", 50, 0, 500*GeV, "", "leading $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["jb2pt"] = histo1D(prefix + "_jb2pt", 50, 0, 500*GeV, "", "subleading $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["j0b1pt"] = histo1D(prefix + "_j0b1pt", 50, 0, 500*GeV, "", "leading $bb$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["j0b2pt"] = histo1D(prefix + "_j0b2pt", 50, 0, 500*GeV, "", "subleading $bb$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["j1b1pt"] = histo1D(prefix + "_j1b1pt", 50, 0, 500*GeV, "", "leading $b1$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["j1b2pt"] = histo1D(prefix + "_j1b2pt", 50, 0, 500*GeV, "", "subleading $b1$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));

        hm["jl1pt_log"] = histo1D(prefix + "_jl1pt_log", logspace(50, 10, 1000*GeV), "", "leading light-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["jl2pt_log"] = histo1D(prefix + "_jl2pt_log", logspace(50, 10, 1000*GeV), "", "subleading light-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["jb1pt_log"] = histo1D(prefix + "_jb1pt_log", logspace(50, 10, 1000*GeV), "", "leading $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["jb2pt_log"] = histo1D(prefix + "_jb2pt_log", logspace(50, 10, 1000*GeV), "", "subleading $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["j0b1pt_log"] = histo1D(prefix + "_j0b1pt_log", logspace(50, 10, 1000*GeV), "", "leading $bb$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["j0b2pt_log"] = histo1D(prefix + "_j0b2pt_log", logspace(50, 10, 1000*GeV), "", "subleading $bb$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["j1b1pt_log"] = histo1D(prefix + "_j1b1pt_log", logspace(50, 10, 1000*GeV), "", "leading $b1$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["j1b2pt_log"] = histo1D(prefix + "_j1b2pt_log", logspace(50, 10, 1000*GeV), "", "subleading $b1$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));

        hm["jl1eta"] = histo1D(prefix + "_jl1eta", 30, -3, 3, "", "leading light-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        hm["jl2eta"] = histo1D(prefix + "_jl2eta", 30, -3, 3, "", "subleading light-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        hm["jb1eta"] = histo1D(prefix + "_jb1eta", 30, -3, 3, "", "leading $b$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        hm["jb2eta"] = histo1D(prefix + "_jb2eta", 30, -3, 3, "", "subleading $b$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        hm["j0b1eta"] = histo1D(prefix + "_j0b1eta", 30, -3, 3, "", "leading $bb$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        hm["j0b2eta"] = histo1D(prefix + "_j0b2eta", 30, -3, 3, "", "subleading $bb$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        hm["j1b1eta"] = histo1D(prefix + "_j1b1eta", 30, -3, 3, "", "leading $b1$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        hm["j1b2eta"] = histo1D(prefix + "_j1b2eta", 30, -3, 3, "", "subleading $b1$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));

        hm["dphitb"] = histo1D(prefix + "_dphitb", 50, 0, 4, "", sdphitb, dsdx(sdphitb, "\\mathrm{1}"));
        hm["detatb"] = histo1D(prefix + "_detatb", 50, 0, 4, "", sdetatb, dsdx(sdetatb, "\\mathrm{1}"));
        hm["drtb"] = histo1D(prefix + "_drtb", 50, 0, 4, "", sdrtb, dsdx(sdrtb, "\\mathrm{1}"));

        hm["dphibb"] = histo1D(prefix + "_dphibb", 50, 0, 4, "", sdphibb, dsdx(sdphibb, "1"));
        hm["drbb"] = histo1D(prefix + "_drbb", 50, 0, 5, "", sdrbb, dsdx(sdrbb, "1"));

        hm["mbb"] = histo1D(prefix + "_mbb", 50, 0, 500*GeV, "", smbb + " [GeV]", dsdx(smbb, "\\mathrm{GeV}"));
        hm["ptbb"] = histo1D(prefix + "_ptbb", 50, 0, 500*GeV, "", sptbb + " [GeV]", dsdx(sptbb, "\\mathrm{GeV}"));
        hm["ht"] = histo1D(prefix + "_ht", 50, 0, 2000*GeV, "", sht + " [GeV]", dsdx(sht, "\\mathrm{GeV}"));

        hm["mbb_log"] = histo1D(prefix + "_mbb_log", logspace(50, 10, 1000*GeV), "", smbb + " [GeV]", dsdx(smbb, "\\mathrm{GeV}"));
        hm["ptbb_log"] = histo1D(prefix + "_ptbb_log", logspace(50, 10, 1000*GeV), "", sptbb + " [GeV]", dsdx(sptbb, "\\mathrm{GeV}"));
        hm["ht_log"] = histo1D(prefix + "_ht_log", logspace(50, 10, 2000*GeV), "", sht + " [GeV]", dsdx(sht, "\\mathrm{GeV}"));

        return hm;
      }


      RegMap initRegions() {
       
        RegMap rm;

        rm["inclusive"] = initRegion("inclusive");
        rm["inclusive_j15"] = initRegion("inclusive_j15");
        rm["inclusive_j30"] = initRegion("inclusive_j30");
        rm["inclusive_j40"] = initRegion("inclusive_j40");

        rm["zerob"] = initRegion("zerob");
        rm["atleastoneb"] = initRegion("atleastoneb");
        rm["onej1b"] = initRegion("onej1b");
        rm["twoj1b"] = initRegion("twoj1b");
        rm["atleasttwoj1b"] = initRegion("atleasttwoj1b");
        rm["threej1b"] = initRegion("threej1b");
        rm["fourj1b"] = initRegion("fourj1b");
        rm["onej0b"] = initRegion("onej0b");
        rm["mbbgt100"] = initRegion("mbbgt100");
        
        return rm;
      }

      void fillRegion(Region& region, const V4s& jls, const V4s& jbs, const V4s& j0bs, const V4s& j1bs, const V4s& tops) {
        region["njl"]->fill(jls.size());
        region["njb"]->fill(j0bs.size()+j1bs.size());
        region["nj0b"]->fill(j0bs.size());
        region["nj1b"]->fill(j1bs.size());

        if (jls.size() >= 1) {
          region["jl1pt"]->fill(jls[0].pt());
          region["jl1pt_log"]->fill(jls[0].pt());
          region["jl1eta"]->fill(jls[0].eta());
        }

        if (jls.size() >= 2) {
          region["jl2pt"]->fill(jls[1].pt());
          region["jl2pt_log"]->fill(jls[1].pt());
          region["jl2eta"]->fill(jls[1].eta());
        }

        if (jbs.size() >= 1) {
          region["jb1pt"]->fill(jbs[0].pt());
          region["jb1pt_log"]->fill(jbs[0].pt());
          region["jb1eta"]->fill(jbs[0].eta());
        }

        if (jbs.size() >= 2) {
          region["jb2pt"]->fill(jbs[1].pt());
          region["jb2pt_log"]->fill(jbs[1].pt());
          region["jb2eta"]->fill(jbs[1].eta());
        }

        if (j0bs.size() >= 1) {
          region["j0b1pt"]->fill(j0bs[0].pt());
          region["j0b1pt_log"]->fill(j0bs[0].pt());
          region["j0b1eta"]->fill(j0bs[0].eta());
        }

        if (j0bs.size() >= 2) {
          region["j0b2pt"]->fill(j0bs[1].pt());
          region["j0b2pt_log"]->fill(j0bs[1].pt());
          region["j0b2eta"]->fill(j0bs[1].eta());
        }

        if (j1bs.size() >= 1) {
          region["j1b1pt"]->fill(j1bs[0].pt());
          region["j1b1pt_log"]->fill(j1bs[0].pt());
          region["j1b1eta"]->fill(j1bs[0].eta());
        }

        if (j1bs.size() >= 2) {
          region["j1b2pt"]->fill(j1bs[1].pt());
          region["j1b2pt_log"]->fill(j1bs[1].pt());
          region["j1b2eta"]->fill(j1bs[1].eta());
        }

        // we define the ht as the scalar sum of all the light, b, and B jet pts
        // plus the masses of all top quarks in the event.
        double ht = 0.0;
        for (const FourMomentum& t: tops)
          ht += t.mass();

        for (const FourMomentum& lj: jls)
          ht += lj.pt();

        for (const FourMomentum& bj: jbs)
          ht += bj.pt();

        region["ht"]->fill(ht);
        region["ht_log"]->fill(ht);

        // find the two leading "b-jets" in the event, where "b-jet" here means
        // the two leading jets with at least one b-quark constituent, with jets
        // with exactly one b-quark constituent taking priority.
        FourMomentum b1, b2;
        if (j1bs.size() >= 2) {
          b1 = j1bs[0];
          b2 = j1bs[1];
        } else if (j1bs.size() >= 1 && j0bs.size() >= 1) {
          b1 = j1bs[0];
          b2 = j0bs[0];
        } else if (j0bs.size() >= 2) {
          b1 = j0bs[0];
          b2 = j0bs[1];
        } else
          return;

        // find the closest (in dR) top to a jet with at least one b
        // and store relative angles
        double mindr = -1;
        double mindrdeta = -1;
        double mindrdphi = -1;
        for (const FourMomentum& jb: jbs) {
          for (const FourMomentum& t: tops) {
            double dr = deltaR(jb, t);
            if (mindr < 1 || dr < mindr) {
              mindr = dr;
              mindrdphi = abs(deltaPhi(jb, t));
              mindrdeta = abs(deltaEta(jb, t));
            }
          }
        }

        if (mindr >= 0) {
          region["drtb"]->fill(mindr);
          region["dphitb"]->fill(mindrdphi);
          region["detatb"]->fill(mindrdeta);
        }

        region["dphibb"]->fill(abs(deltaPhi(b1, b2)));
        region["drbb"]->fill(deltaR(b1, b2));

        region["mbb"]->fill((b1 + b2).mass());
        region["ptbb"]->fill((b1 + b2).pt());

        region["mbb_log"]->fill((b1 + b2).mass());
        region["ptbb_log"]->fill((b1 + b2).pt());

      };


    public:

      DEFAULT_RIVET_ANALYSIS_CTOR(MCTTBB3);

      void init() {
        FinalState fps(Cuts::abseta < 5 && Cuts::abspid != 6);
        declare(FinalState(Cuts::abspid==6), "Tops");

        declare(FastJets(fps, FastJets::ANTIKT, 0.4), "Jets");

        regions = initRegions();
      }


      void analyze(const Event& event) {

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



        vector<V4s> j25cats = jet_categories(jets);
        vector<V4s> j15cats = jet_categories(jets15);
        vector<V4s> j30cats = jet_categories(jets30);
        vector<V4s> j40cats = jet_categories(jets40);

        // find the light, b, and B jets by checking the number of b-quark
        // constituents in the jet.
        // jls = light jets
        // jbs = jets with at least one b-quark inside
        // j0bs = jets with at least one b-quark, but a corresponding
        //        anti-b-quark for each b-quark
        // j1bs = jets with an imbalance between b-quarks and anti-b-quarks

        V4s jls = j25cats[0];
        V4s jbs = j25cats[1];
        V4s j0bs = j25cats[2];
        V4s j1bs = j25cats[3];

        fillRegion(regions["inclusive_j15"], j15cats[0], j15cats[1], j15cats[2], j15cats[3], tquarks);
        fillRegion(regions["inclusive_j30"], j30cats[0], j30cats[1], j30cats[2], j30cats[3], tquarks);
        fillRegion(regions["inclusive_j40"], j40cats[0], j40cats[1], j40cats[2], j40cats[3], tquarks);

        // inclusive selection
        fillRegion(regions["inclusive"], jls, jbs, j0bs, j1bs, tquarks);

        // no accepted jets with b-quarks inside
        if (jbs.size() == 0) fillRegion(regions["zerob"], jls, jbs, j0bs, j1bs, tquarks);
        // at least one jet with b-quarks inside
        else fillRegion(regions["atleastoneb"], jls, jbs, j0bs, j1bs, tquarks);

        // exactly one jet with an imbalance of b and anti-b
        if (j1bs.size() == 1) fillRegion(regions["onej1b"], jls, jbs, j0bs, j1bs, tquarks);
        // exactly two jets with an imbalance of b and anti-b
        else if (j1bs.size() == 2) fillRegion(regions["twoj1b"], jls, jbs, j0bs, j1bs, tquarks);
        // exactly three jets with an imbalance of b and anti-b
        else if (j1bs.size() == 3) fillRegion(regions["threej1b"], jls, jbs, j0bs, j1bs, tquarks);
        // exactly four jets with an imbalance of b and anti-b
        else if (j1bs.size() == 4) fillRegion(regions["fourj1b"], jls, jbs, j0bs, j1bs, tquarks);

        // at least two jets with an imbalance of b and anti-b
        if (j1bs.size() >= 2) fillRegion(regions["atleasttwoj1b"], jls, jbs, j0bs, j1bs, tquarks);

        // exactly one jet with balance of b and anti-b
        if (j0bs.size() == 1) fillRegion(regions["onej0b"], jls, jbs, j0bs, j1bs, tquarks);

        // find the two leading "b-jets" in the event, where "b-jet" here means
        // the two leading jets with at least one b-quark constituent, with jets
        // with exactly one b-quark constituent taking priority.
        FourMomentum b1, b2;
        if (j1bs.size() >= 2) {
          b1 = j1bs[0];
          b2 = j1bs[1];
        } else if (j1bs.size() >= 1 && j0bs.size() >= 1) {
          b1 = j1bs[0];
          b2 = j0bs[0];
        } else if (j0bs.size() >= 2) {
          b1 = j0bs[0];
          b2 = j0bs[1];
        } else
          return;

        if ((b1+b2).mass() > 100*GeV) fillRegion(regions["mbbgt100"], jls, jbs, j0bs, j1bs, tquarks);

        return;
      }

      void finalize() {
          for (pair<const string, Region>& keyval1: regions) {
            Region& region = keyval1.second;

            for (pair<const string , Histo1DPtr>& keval2: region) {
              const string& name = keyval1.first + "/" + keval2.first;
              Histo1DPtr& ph = keval2.second;

              scale(ph, crossSection()/picobarn/sumOfWeights());

              ph->setAnnotation("XLabel", xlabs[name]);
              ph->setAnnotation("YLabel", ylabs[name]);
              ph->setAnnotation("Title", titles[name]);
            }
          }
        }

    private:
      vector<V4s> jet_categories(const Jets& jets) {
        V4s jls, jbs, j0bs, j1bs;
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
