/*
 15/05/2015
 Class to store one PDF
*/

// avoid multiple def of the same class
#ifndef PDF_H  // check that this keyword has been set already
#define PDF_H  // define the keyword to be checked

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

/// Class object implementing a probability density function
class PDF {
 private:
  bool scaleLinear = 1;
  double scaleStep = -99, invScaleStep = -99;
  double scaleMin = 0, scaleMax = 0;

  size_t vsize;

 public:
  vector<double> vPDF;
  vector<double> chi2, xaxis, secondX, secondP;
  vector<int> ind, secondInd;

  // Prototype
  PDF() { ; }  // necessary for the pdfmap in the onesource constructor
  PDF(const double min, const double max,
      const size_t size);  ///< initialize the pdf with minimal value, maximal
                           ///< value, and total number of points.
  ~PDF() {
    vPDF.clear();
    chi2.clear();
    xaxis.clear();
  }

  void normalization();  ///< normalize the pdf to 1.
  void chi2toPDF();      ///< compute \f$\exp(-\chi^2/2)\f$ from the vector of
                         ///< \f$\chi^2\f$ values
  int chi2mini();        ///< index of the minimum \f$\chi^2\f$ value
  double
  int_parab();  ///< compute the local \f$\chi^2\f$ by parabolic interpolation
  pair<double, double> uncMin(
      double dchi);  ///< compute the bounds of the confidence interval
                     ///< defined by \f$\chi^2_{min} + dchi\f$
  double probaBay(
      double xvalue);  ///< return the inverse cumulative function at xvalue
  size_t index(
      const double inVal) const;     ///< return the vector index where the
                                     ///< pdf vector is closest to inVal
  void secondMax(const double win);  ///< search for high peaks in The ML
                                     ///< function vs xaxis, and sort them from
                                     ///< highest to smallest peaks in ML
  double levelCumu2x(float xval);    // find the xaxis value corresponding to a
                                     // level in the cumulative function
  pair<double, double> credible_interval(
      float xval,
      double redLoc);  // Look for the confidence interval around redLoc
  inline size_t get_maxid() const {
    return max_element(vPDF.begin(), vPDF.end()) - vPDF.begin();
  }
  inline double get_max() const { return xaxis[get_maxid()]; }
  inline double linear_interp(double x, double x1, double y1, double x2,
                              double y2) {
    if (x2 > x1)
      return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
    else
      return y1;
  }
  // for marginalized pdf use
  double int_parabL() {
    for (size_t k = 0; k < vPDF.size(); k++) {
      chi2[k] = -vPDF[k];
    }
    return int_parab();
  }

  inline size_t size() const { return vsize; }
};

#endif
