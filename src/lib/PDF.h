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

  inline double linear_interp(double x, double x1, double y1, double x2,
                              double y2) {
    if (x2 > x1)
      return y1 + (x - x1) * (y2 - y1) / (x2 - x1);
    else
      return y1;
  }

 public:
  vector<double> vPDF;  ///< probability density value at each #xaxis point
  vector<double> chi2,  ///< chi2 value at each #xaxis point (alternative
                        ///< storage to #vPDF for the chi2 curve)
      xaxis,            ///< linearly-sampled grid the PDF/chi2 is stored on
      secondX,          ///< x position of each secondary peak found by
                        ///< secondMax(), sorted from highest to smallest
      secondP;          ///< probability of each secondary peak in #secondX
  /// @warning intended (per the constructor's comment) to hold, for each
  /// #xaxis bin, the index of the best-matching SED template at that
  /// redshift/parameter value, but it is initialised to 0 and never
  /// populated anywhere in the current codebase. As a result secondMax()'s
  /// #secondInd output is always 0, which onesource::secondpeak() then uses
  /// to index the SED library: every secondary-solution attribute it
  /// reports other than the redshift itself and its probability (i.e.
  /// zsecEbv, zsecExtlaw, zsecScale, zsecMod, zsecAge, indminSec) is
  /// silently wrong.
  vector<int> ind,
      secondInd;  ///< rank-sorted #ind value at each secondary peak found
                  ///< by secondMax(); see the @warning on #ind above

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

  /// Compute the cumulative distribution of #vPDF over #xaxis (trapezoidal
  /// rule), un-normalized (starts at 0, ends at the total integral)
  /// @return the cumulative distribution, one value per #xaxis point
  vector<double> cumulant();

  /*!
    Compute \f$e^{-\chi^2/2}\f$ from the vector
    of \f$\chi^2\f$ values chi2.
  */
  void chi2toPDF();

  int chi2mini();  ///< index of the minimum \f$\chi^2\f$ value

  pair<double, double> uncMin(
      double dchi);  ///< compute the bounds of the confidence interval
                     ///< defined by \f$\chi^2_{min} + dchi\f$
  double probaBay(
      double xvalue);  ///< return the inverse cumulative function at xvalue

  /*! return the vector index where the pdf vector is closest to inVal
   *
   * @param inVal: input value
   * @return : return i so that pdf[i]<=inVal<pdf[i+1] if possible
   * else return i=0 if inVal<pdf[0] or i=size-1 if inVal>pdf[size-1].
   * In other words, overflows go to boundary bins.
   */
  size_t index(
      const double inVal) const;     ///< return the vector index where the
                                     ///< pdf vector is closest to inVal
  void secondMax(const double win);  ///< search for high peaks in The ML
                                     ///< function vs xaxis, and sort them from
                                     ///< highest to smallest peaks in ML
  /*! Find the #xaxis value corresponding to a given level of the normalized
   * cumulative distribution
   * @param xval: target cumulative-probability level, in [0,1] (or in
   * percent, i.e. up to 100, which is converted internally)
   * @return the interpolated #xaxis value at that cumulative level, or
   * -99.9 if it could not be bracketed
   */
  double levelCumu2x(float xval);

  /*!
   * Improve the the grid extremum by quadratic approximation around it
   *
   * @param is_chi2: is the extremum to be found for the chi2 yaxis (a minimum
   * or the PDF yaxis (a maximum)
   *
   * @return (xmin, ymin): the pair x,y at the found minimum unless :
   *  - for chi2 fit one of the 3 y values set to HIGH_CHI2
   *  - for the PDF fit one of the 3 y values set to 0
   *  - if the distance between the midpoint in x (the grid minimum)
   * and any of the other 2 points is more than 0.5 in x.
   *
   * In all these cases, the returned pair is current min and its y value.
   */
  pair<double, double> improve_extremum(bool is_chi2) const;

  /*! Compute the confidence interval around the chi2min
   *  @param level: confidence level (value to add to the min chi2
   * to define the interval)
   * @return pair: the bounds of the interval
   */
  pair<double, double> confidence_interval(float level);

  /*! Compute a credible interval on the PDF
   * @param level: confidence level
   * @param val: the value around which the interval is computed
   * The interval is symmetric as built : level/2 is the mass on each
   * side of `val`, unless the boundaries are met. In this case
   * the mass is maximised to the boundary, and the rest is allocated on
   * the other side
   */
  pair<double, double> credible_interval(float level, double val);

  //! Get the index of the PDF max value
  inline size_t get_maxid() const {
    return max_element(vPDF.begin(), vPDF.end()) - vPDF.begin();
  }

  //! Get the x value at PDF maximum (as returned by `get_maxid`)
  inline double get_max() const { return xaxis[get_maxid()]; }

  //! Return the size of the x array
  inline size_t size() const { return vsize; }
};

/*! extremum finding by quadratic approximation
 *
 * @param x1, y1, x2, y2,x3,y3: are the three points used to
 * obtain a parabolic approximation around them
 *
 * @return (xmin, ymin): the pair of the approximate minimum
 */
pair<double, double> quadratic_extremum(double x1, double x2, double x3,
                                        double y1, double y2, double y3);
#endif
