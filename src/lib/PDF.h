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

//! @class PDF
/*!
Implement a probability density function.
In order to avoid duplicating code, this class is also able to
store the chi2 curve. By construction the xaxis of the pdf is linearly
sampled.
*/
class PDF {
 private:
  double scaleStep = -99, invScaleStep = -99;
  double scaleMin = 0, scaleMax = 0;

  size_t vsize;

 public:
  vector<double> vPDF;
  vector<double> chi2, xaxis, secondX, secondP;
  vector<int> ind, secondInd;

  //! Empty constructor, needed by \ref onesource constructor
  PDF() { ; }

  //! Standard constructor
  /*!
    \param min start value for the x axis
    \param max end value for the x axis
    \param size number of points on the x axis grid
    The resulting points on the grid are thus : \f$x[k] = min +
    k\frac{max-min}{size-1}\f$
  */
  PDF(const double min, const double max, const size_t size);

  //! Destructor. Clear internal vectors xaxis, vPDF, chi2.
  ~PDF() {
    vPDF.clear();
    chi2.clear();
    xaxis.clear();
  }

  //! Normalize the vPDF vector to 1.
  /*!
    A simple trapezoidal rule is used.
   */
  double normalization();

  vector<double> cumulant();

  /*!
    Compute \f$e^{-\chi^2/2}\f$ from the vector
    of \f$\chi^2\f$ values chi2.
  */
  void chi2toPDF();

  int chi2mini();  ///< index of the minimum \f$\chi^2\f$ value
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

  pair<double, double> improve_extremum(bool is_chi2) const;

  pair<double, double> confidence_interval(float level);

  pair<double, double> credible_interval(float level, double val);

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

pair<double, double> quadratic_extremum(double x1, double x2, double x3,
                                        double y1, double y2, double y3);
#endif
