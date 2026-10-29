/*
 15/05/2015
 Class to store one PDF
*/

// avoid multiple def of the same class
#ifndef PDF_H  // check that this keyword has been set already
#define PDF_H  // define the keyword to be checked

#include <algorithm>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

/// Version of the definition of the redshift quality flag. 1 is the historical
/// Python definition (Gaussian least-squares fit for sigma, #PDF is not
/// involved), 2 is the definition implemented in zflag() below.
constexpr int ZFLAG_DEFINITION = 2;
/// Value of the flag for a source whose flag was not computed
constexpr int ZFLAG_NOT_COMPUTED = -999;
/// Value of the flag when it cannot be computed (null or invalid PDF)
constexpr int ZFLAG_INVALID = -99;

//! Parameters controlling the components of the redshift quality flag
/*!
  The flag is a sum of bits, a higher value meaning a lower quality PDF:
  - 1: error (pseudo standard deviation of the PDF around the estimate)
    above #error, but below half of the redshift range
  - 2: ratio mean(PDF)/max(PDF) above #peak_ratio
  - 4: probability mass outside +- #n_window sigma around the estimate above
    #tail, where sigma is derived from the FWHM of the main peak
  - 8: at least #nb_peak significant peaks
  - 16: error above half of the redshift range
  The parameters can be set from a string "name=value,name=value,...", see
  from_string(). Unspecified parameters keep their default value.
*/
struct FlagParams {
  int nb_peak = 2;          ///< number of peaks from which bit 8 is set
  double height = 0.43;     ///< peak detection threshold (fraction of the max)
  double peak_dz = 0.05;    ///< minimal separation of two peaks (in x units)
  double tail = 0.23;       ///< tail mass threshold for bit 4
  double n_window = 2.0;    ///< half-width of the core, in sigma, for the tail
  double peak_ratio = 0.1;  ///< mean/max threshold for bit 2
  double error = 0.1;       ///< error threshold for bit 1
  double fwhm_level = 0.5;  ///< level (fraction of the max) used to measure
                            ///< the width of the main peak: 0.5 is the FWHM

  /*! Build the parameters from "name=value,name=value,...".
   * Names: nb_peak, height, peak_dz, tail, n_window, peak_ratio, error,
   * fwhm_level (case insensitive). An empty string, or NONE, gives the
   * defaults. Throws std::invalid_argument on an unknown name, a malformed
   * item or an out of range value.
   */
  static FlagParams from_string(const string& str);
  /// Canonical "name=value,..." string, readable by from_string()
  string to_string() const;
};

//! Result of the computation of the quality flag of a PDF
struct ZFlagResult {
  int flag = ZFLAG_INVALID;  ///< the flag (sum of bits), -99 if invalid
  double estimate = -99.;    ///< estimate around which statistics are made
  double error = -99.;       ///< pseudo standard deviation around estimate
  double peak_ratio = -99.;  ///< mean/max of the PDF
  double tail_mass = -99.;   ///< probability mass outside the core
  int n_peaks = -99;         ///< number of significant peaks
  double sigma = -99.;       ///< Gaussian equivalent width of the main peak
};

/*! Compute the quality flag of a PDF sampled on a regular grid
 * @param x: grid, increasing and regularly sampled
 * @param v: PDF values on the grid (not necessarily normalized)
 * @param n: number of points
 * @param params: parameters of the flag
 * @param estimate: central value for the statistics. If NaN, the position of
 * the PDF maximum is used.
 * @return the flag and the quantities it is made of. For a null PDF, or one
 * with negative or non finite values, flag is #ZFLAG_INVALID.
 */
ZFlagResult zflag(const double* x, const double* v, size_t n,
                  const FlagParams& params = FlagParams(),
                  double estimate = numeric_limits<double>::quiet_NaN());

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
  /// Index, in the SED library, of the best-matching (lowest chi2) template
  /// at each #xaxis bin. Populated for the GAL/QSO redshift PDFs
  /// (pdfmap[9]/pdfmap[10]) by onesource::generatePDF() as the fit loop
  /// runs; not populated for the other PDF types (physical parameters,
  /// Bayesian redshift) since secondMax()/secondInd is only used for the
  /// GAL redshift PDF by onesource::secondpeak().
  vector<int> ind,
      secondInd;  ///< rank-sorted #ind value at each secondary peak found
                  ///< by secondMax()

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

  //! Quality flag of this PDF, see zflag() and FlagParams
  inline ZFlagResult quality_flag(
      const FlagParams& params = FlagParams(),
      double estimate = numeric_limits<double>::quiet_NaN()) const {
    return zflag(xaxis.data(), vPDF.data(), min(xaxis.size(), vPDF.size()),
                 params, estimate);
  }

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
 * @return (xmin, ymin): the pair of the approximate minimum. Return (x2, y2) if
 * the points are aligned.
 */
pair<double, double> quadratic_extremum(double x1, double x2, double x3,
                                        double y1, double y2, double y3);
#endif
