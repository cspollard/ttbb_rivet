// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/IdentifiedFinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/DressedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/PartonicTops.hh"
#include "Rivet/Tools/JetSmearingFunctions.hh"

namespace Rivet {

  using V4s = vector<FourMomentum>;
  using Region = map< string , Histo1DPtr >;
  using Variation = map< string , Region >;
  using VarMap = map< string , Variation >;

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


  class MCTTBBDecayed : public Analysis {
    public:

      DEFAULT_RIVET_ANALYSIS_CTOR(MCTTBBDecayed);


    private:
      string dsdx(const string& x, const string& xunit) {
        return "\\ensuremath{\\frac{d\\sigma}{d" + x + "} \\frac{{pb}}{" + xunit + "}}";
      }


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

      VarMap initVars() {
       
        VarMap vm;
        vm["nominal"] = initVariation("nominal");
        vm["jes"] = initVariation("jes");
        vm["btagrate"] = initVariation("btagrate");
        vm["ctagrate"] = initVariation("ctagrate");
        vm["nodet"] = initVariation("nodet");

        return vm;
      }


      Variation initVariation(const string& prefix) {
       
        Variation rm;

        rm["onelep_inclusive"] = initRegion(prefix + "/onelep_inclusive");
        rm["onelep_eq5j_eq2jb"] = initRegion(prefix + "/onelep_eq5j_eq2jb");
        rm["onelep_eq5j_eq3jb"] = initRegion(prefix + "/onelep_eq5j_eq3jb");
        rm["onelep_eq5j_ge4jb"] = initRegion(prefix + "/onelep_eq5j_ge4jb");
        rm["onelep_ge6j_eq2jb"] = initRegion(prefix + "/onelep_ge6j_eq2jb");
        rm["onelep_ge6j_eq3jb"] = initRegion(prefix + "/onelep_ge6j_eq3jb");
        rm["onelep_ge6j_ge4jb"] = initRegion(prefix + "/onelep_ge6j_ge4jb");

        rm["dilep_inclusive"] = initRegion(prefix + "/dilep_inclusive");
        rm["dilep_eq3j_eq2jb"] = initRegion(prefix + "/dilep_eq3j_eq2jb");
        rm["dilep_eq3j_eq3jb"] = initRegion(prefix + "/dilep_eq3j_eq3jb");
        rm["dilep_eq3j_eq4jb"] = initRegion(prefix + "/dilep_eq3j_eq4jb");
        rm["dilep_ge4j_eq2jb"] = initRegion(prefix + "/dilep_ge4j_eq2jb");
        rm["dilep_ge4j_eq3jb"] = initRegion(prefix + "/dilep_ge4j_eq3jb");
        rm["dilep_ge4j_ge4jb"] = initRegion(prefix + "/dilep_ge4j_ge4jb");
        
        return rm;
      }

      Region initRegion(const string& prefix) {

        Region hm;

        hm["nj"] = histo1D(prefix + "/nj", 10, -0.5, 9.5, "", "jet multiplicity", dsdx("n", "1"));
        hm["njl"] = histo1D(prefix + "/njl", 10, -0.5, 9.5, "", "light-jet multiplicity", dsdx("n", "1"));
        hm["njb"] = histo1D(prefix + "/njb", 8, -0.5, 7.5, "", "$b$-jet mulitplicity", dsdx("n", "1"));

        hm["jl1pt"] = histo1D(prefix + "/jl1pt", 50, 0, 500*GeV, "", "leading light-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["jl2pt"] = histo1D(prefix + "/jl2pt", 50, 0, 500*GeV, "", "subleading light-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));

        hm["jb1pt"] = histo1D(prefix + "/jb1pt", 50, 0, 500*GeV, "", "leading $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["jb2pt"] = histo1D(prefix + "/jb2pt", 50, 0, 500*GeV, "", "subleading $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["jb3pt"] = histo1D(prefix + "/jb3pt", 50, 0, 500*GeV, "", "3rd $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["jb4pt"] = histo1D(prefix + "/jb4pt", 50, 0, 500*GeV, "", "4th $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));

        hm["jl1pt_log"] = histo1D(prefix + "/jl1pt_log", logspace(50, 10, 1000*GeV), "", "leading light-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["jl2pt_log"] = histo1D(prefix + "/jl2pt_log", logspace(50, 10, 1000*GeV), "", "subleading light-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["jb1pt_log"] = histo1D(prefix + "/jb1pt_log", logspace(50, 10, 1000*GeV), "", "leading $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["jb2pt_log"] = histo1D(prefix + "/jb2pt_log", logspace(50, 10, 1000*GeV), "", "subleading $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["jb3pt_log"] = histo1D(prefix + "/jb3pt_log", logspace(50, 10, 1000*GeV), "", "3rd $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["jb4pt_log"] = histo1D(prefix + "/jb4pt_log", logspace(50, 10, 1000*GeV), "", "4th $b$-jet " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));

        hm["jl1eta"] = histo1D(prefix + "/jl1eta", 30, -3, 3, "", "leading light-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        hm["jl2eta"] = histo1D(prefix + "/jl2eta", 30, -3, 3, "", "subleading light-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        hm["jb1eta"] = histo1D(prefix + "/jb1eta", 30, -3, 3, "", "leading $b$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));
        hm["jb2eta"] = histo1D(prefix + "/jb2eta", 30, -3, 3, "", "subleading $b$-jet " + seta, dsdx(seta, "\\mathrm{GeV}"));

        hm["lep1pt"] = histo1D(prefix + "/lep1pt", 50, 0, 500*GeV, "", "leading lepton " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["lep2pt"] = histo1D(prefix + "/lep2pt", 50, 0, 500*GeV, "", "subleading lepton " + spt + " [GeV]", dsdx(spt, "\\mathrm{GeV}"));
        hm["lep1eta"] = histo1D(prefix + "/lep1eta", 30, -3, 3, "", "leading lepton " + seta, dsdx(seta, "\\mathrm{GeV}"));
        hm["lep2eta"] = histo1D(prefix + "/lep2eta", 30, -3, 3, "", "subleading lepton " + seta, dsdx(seta, "\\mathrm{GeV}"));

        hm["dphibb"] = histo1D(prefix + "/dphibb", 50, 0, 4, "", sdphibb, dsdx(sdphibb, "1"));
        hm["drbb"] = histo1D(prefix + "/drbb", 50, 0, 5, "", sdrbb, dsdx(sdrbb, "1"));

        hm["mbb_avg"] = histo1D(prefix + "/mbb_avg", 50, 0, 500*GeV, "", "average " + smbb + " [GeV]", dsdx(smbb, "\\mathrm{GeV}"));
        hm["mbb_wind"] = histo1D(prefix + "/mbb_wind", 5, -0.5, 4.5, "", "number of " + smbb + " combinations in (100, 150) GeV", dsdx("n", "1"));

        hm["mbb"] = histo1D(prefix + "/mbb", 50, 0, 500*GeV, "", smbb + " [GeV]", dsdx(smbb, "\\mathrm{GeV}"));
        hm["ptbb"] = histo1D(prefix + "/ptbb", 50, 0, 500*GeV, "", sptbb + " [GeV]", dsdx(sptbb, "\\mathrm{GeV}"));
        hm["ht"] = histo1D(prefix + "/ht", 50, 0, 2000*GeV, "", sht + " [GeV]", dsdx(sht, "\\mathrm{GeV}"));
        hm["htlep"] = histo1D(prefix + "/htlep", 50, 0, 2000*GeV, "", shtlep + " [GeV]", dsdx(shtlep, "\\mathrm{GeV}"));

        hm["mbb_log"] = histo1D(prefix + "/mbb_log", logspace(50, 10, 1000*GeV), "", smbb + " [GeV]", dsdx(smbb, "\\mathrm{GeV}"));
        hm["ptbb_log"] = histo1D(prefix + "/ptbb_log", logspace(50, 10, 1000*GeV), "", sptbb + " [GeV]", dsdx(sptbb, "\\mathrm{GeV}"));
        hm["ht_log"] = histo1D(prefix + "/ht_log", logspace(50, 10, 2000*GeV), "", sht + " [GeV]", dsdx(sht, "\\mathrm{GeV}"));
        hm["htlep_log"] = histo1D(prefix + "/htlep_log", logspace(50, 10, 2000*GeV), "", shtlep + " [GeV]", dsdx(shtlep, "\\mathrm{GeV}"));

        return hm;
      }


      void fillRegion(Region& region, const V4s& jls, const V4s& jbs, const V4s& leps) {
        size_t njbs = jbs.size();
        size_t njls = jls.size();

        region["nj"]->fill(njbs + njls);
        region["njl"]->fill(njls);
        region["njb"]->fill(njbs);

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

        if (jbs.size() >= 3) {
          region["jb3pt"]->fill(jbs[2].pt());
          region["jb3pt_log"]->fill(jbs[2].pt());
        }

        if (jbs.size() >= 4) {
          region["jb4pt"]->fill(jbs[3].pt());
          region["jb4pt_log"]->fill(jbs[3].pt());
        }

        if (leps.size() >= 1) {
          region["lep1pt"]->fill(leps[0].pt());
          region["lep1eta"]->fill(leps[0].eta());
        }

        if (leps.size() >= 2) {
          region["lep2pt"]->fill(leps[1].pt());
          region["lep2eta"]->fill(leps[1].eta());
        }

        float ht = 0.0;
        for (const FourMomentum& lj: jls)
          ht += lj.pt();

        for (const FourMomentum& bj: jbs)
          ht += bj.pt();

        region["ht"]->fill(ht);
        region["ht_log"]->fill(ht);

        float htlep = ht;
        for (const FourMomentum& l: leps)
          htlep += l.pt();

        region["htlep"]->fill(htlep);
        region["htlep_log"]->fill(htlep);

        if (jbs.size() < 2)
          return;

        FourMomentum b1, b2;
        float bestmbb = -1;
        for (unsigned int i = 0; i < jbs.size(); i++) {
          FourMomentum pi = jbs.at(i);
          for (unsigned int k = i+1; k < jbs.size(); k++) {
            FourMomentum pk = jbs.at(k);

            float mbb = (pi + pk).mass();

            if (bestmbb < 0 || (fabs(mbb - 125*GeV) < fabs(bestmbb - 125*GeV))) {
              b1 = jbs.at(i);
              b2 = jbs.at(k);
              bestmbb = mbb;
            }
          }
        }

        region["dphibb"]->fill(abs(deltaPhi(b1, b2)));
        region["drbb"]->fill(deltaR(b1, b2));

        region["mbb"]->fill((b1 + b2).mass());
        region["ptbb"]->fill((b1 + b2).pt());
        region["mbb_log"]->fill((b1 + b2).mass());
        region["ptbb_log"]->fill((b1 + b2).pt());

        size_t ncomb = jbs.size()*(jbs.size() - 1);
        size_t nwind = 0;
        float mbbsum = 0;
        for (size_t ib = 0; ib < jbs.size(); ib++) {
          for (size_t jb = ib+1; jb < jbs.size(); jb++) {
            float mbb = (jbs[ib] + jbs[jb]).mass();
            mbbsum += mbb;

            if (100*GeV < mbb && mbb < 150*GeV)
              nwind++;
          }
        }

        region["mbb_avg"]->fill(mbbsum/ncomb);
        region["mbb_wind"]->fill(nwind);

      };

      void fillVar(Variation& variation, float jescale, JET_BTAG_EFFS tagger, const Jets& js, const vector<DressedLepton>& dls) {

        V4s leps;
        for (const DressedLepton& dl : dls)
          leps.push_back(dl.mom());
        
        vector<V4s> jcats = jet_categories(js, 25*GeV, jescale, tagger);
        V4s jls = jcats[0];
        V4s jbs = jcats[1];

        int njets = jls.size() + jbs.size();

        if (leps.size() == 1) {
          fillRegion(variation["onelep_inclusive"], jls, jbs, leps);

          if (njets == 5) {
            if (jbs.size() == 2) fillRegion(variation["onelep_eq5j_eq2jb"], jls, jbs, leps);
            else if (jbs.size() == 3) fillRegion(variation["onelep_eq5j_eq3jb"], jls, jbs, leps);
            else if (jbs.size() >= 4) fillRegion(variation["onelep_eq5j_ge4jb"], jls, jbs, leps);
          } else if (njets >= 6) {
            if (jbs.size() == 2) fillRegion(variation["onelep_ge6j_eq2jb"], jls, jbs, leps);
            else if (jbs.size() == 3) fillRegion(variation["onelep_ge6j_eq3jb"], jls, jbs, leps);
            else if (jbs.size() >= 4) fillRegion(variation["onelep_ge6j_ge4jb"], jls, jbs, leps);
          }

        } else if (leps.size() == 2) {
          fillRegion(variation["dilep_inclusive"], jls, jbs, leps);

          if (njets == 3) {
            if (jbs.size() == 2) fillRegion(variation["dilep_eq3j_eq2jb"], jls, jbs, leps);
            else if (jbs.size() == 3) fillRegion(variation["dilep_eq3j_eq3jb"], jls, jbs, leps);
            else if (jbs.size() >= 4) fillRegion(variation["dilep_eq3j_eq3jb"], jls, jbs, leps);
          } else if (njets >= 4) {
            if (jbs.size() == 2) fillRegion(variation["dilep_ge4j_eq2jb"], jls, jbs, leps);
            else if (jbs.size() == 3) fillRegion(variation["dilep_ge4j_eq3jb"], jls, jbs, leps);
            else if (jbs.size() >= 4) fillRegion(variation["dilep_ge4j_ge4jb"], jls, jbs, leps);
          }

        }

        return;
      }


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


        variations = initVars();
      }


      void analyze(const Event& event) {

        const Jets& jets =
          apply<FastJets>(event, "Jets").jetsByPt(
              Cuts::pT > 20*GeV && Cuts::abseta < 2.5);

        const vector<DressedLepton>& leps =
          apply<DressedLeptons>(event, "Leptons").dressedLeptons();


        // nominal: 70% b-tagging efficiency, 20% charm mistag rate

        float btagrate = 0.70;
        float ctagrate = 0.20;
        float jescale = 1.0;

        JET_BTAG_EFFS tagger(btagrate, ctagrate, 0.0);

        fillVar(variations["nominal"], jescale, tagger, jets, leps);


        // btag eff variation: 72% b-tagging efficiency, 20% charm mistag rate

        btagrate = 0.72;
        ctagrate = 0.20;
        jescale = 1.0;

        tagger = JET_BTAG_EFFS(btagrate, ctagrate, 0.0);

        fillVar(variations["btagrate"], jescale, tagger, jets, leps);


        // ctag eff variation: 70% b-tagging efficiency, 25% charm mistag rate

        btagrate = 0.70;
        ctagrate = 0.25;
        jescale = 1.0;

        tagger = JET_BTAG_EFFS(btagrate, ctagrate, 0.0);

        fillVar(variations["ctagrate"], jescale, tagger, jets, leps);


        // jes variation: jes -> 1.05

        btagrate = 0.70;
        ctagrate = 0.20;
        jescale = 1.05;

        tagger = JET_BTAG_EFFS(btagrate, ctagrate, 0.0);

        fillVar(variations["jes"], jescale, tagger, jets, leps);

        // no detector
        btagrate = 1.0;
        ctagrate = 0.0;
        jescale = 1.0;

        tagger = JET_BTAG_EFFS(btagrate, ctagrate, 0.0);

        fillVar(variations["nodet"], jescale, tagger, jets, leps);

        return;
      }

      void finalize() {
        for (pair<const string, Variation>& keyval: variations) {
          Variation& variation = keyval.second;

          for (pair<const string, Region>& keyval1: variation) {
            Region& region = keyval1.second;

            for (pair<const string , Histo1DPtr>& keval2: region) {
              const string& name = keyval.first + "/" + keyval1.first + "/" + keval2.first;
              Histo1DPtr& ph = keval2.second;

              scale(ph, crossSection()/picobarn/sumOfWeights());

              ph->setAnnotation("XLabel", xlabs[name]);
              ph->setAnnotation("YLabel", ylabs[name]);
              ph->setAnnotation("Title", titles[name]);
            }
          }
        }
      }

    private:
      VarMap variations;
      StringMap xlabs, ylabs, titles;

      vector<V4s> jet_categories(
          const Jets& jets
        , float etcut
        , float jescale
        , JET_BTAG_EFFS tagrates
        ) {

        V4s jls, jbs;
        for (const Jet& j: jets) {
          FourMomentum mom = j.mom()*jescale;
          if (mom.Et() < etcut)
              continue;

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