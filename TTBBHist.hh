namespace Rivet {

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

Histo1DPtr histo1D(
  const string& name, size_t nb, double low, double high
  , const string& title, const string& xlab, const string& ylab) {

    Histo1DPtr h = make_shared<Histo1D>(nb, low, high, name, title);
    h->setAnnotation("XLabel", xlab);
    h->setAnnotation("YLabel", ylab);
    return h;
  }


  // a simple wrapper around a histogram that stores histograms corresponding to
  // positive- and negative-weighted events.

  class TTBBHist {
  public:
    TTBBHist() { };

    TTBBHist(
      const string& name, size_t nb, double low, double high
      , const string& title, const string& xlab, const string& ylab) {
        hnom = histo1D(name, nb, low, high, title, xlab, ylab);
        hpos = histo1D(name + "_onlypos", nb, low, high, title, xlab, ylab);
        hneg = histo1D(name + "_onlyneg", nb, low, high, title, xlab, ylab);
      }

      void fill(double xval, double weight) {
        hnom->fill(xval, weight);
        if (weight >= 0) hpos->fill(xval, weight);
        else hneg->fill(xval, -weight);
        return;
      }

      Histo1DPtr nominal() {
        return hnom;
      }

      Histo1DPtr positive() {
        return hpos;
      }

      Histo1DPtr negative() {
        return hneg;
      }

      vector<Histo1DPtr> histograms() {
        return { hnom, hpos, hneg };
      }


    private:
      Histo1DPtr hnom;
      Histo1DPtr hpos;
      Histo1DPtr hneg;
    };

}
