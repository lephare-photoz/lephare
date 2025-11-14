/*
10/11/2014
Class to store one filter
*/

// avoid multiple def of the same class
#ifndef FLT_H  // check that this keyword has been set already
#define FLT_H  // define the keyword to be checked

#include <fstream>   // print output file
#include <iostream>  // print standard file
#include <string>
#include <vector>

#include "globals.h"
#include "oneElLambda.h"

using namespace std;

/// Class to store filter information. A filter class instance is in essence a
/// vector of oneElLambda elements.
class flt {
 private:
  /*! \brief transform a transmission curve in photon units into a transmission
   * curve in energy units.
   *
   * If TRANS_TYPE=1, the transmission curve is provided in photon units and
   * needs to be converted as \f$T(\lambda)\rightarrow
   * T(\lambda)\cdot\lambda/\lambda_m\f$ where \f$\lambda_m\f$ is the mean
   * wavelength of the filter transmission curve, see flt::lambdaMean.
   * This transformation leaves the transmission curve integral invariant.
   */
  void trans();
  /// Sort the filter in lambda, remove the values with a low transmission on
  /// the edge, be sure that the extreme points are ending with 0
  void clean();

 public:
  vector<oneElLambda> lamb_trans;
  int id;
  string name;
  int transtyp, calibtyp;
  double leff, lmean, dwidth, ab, tg, veg, msun, fcorr, tpeak;

  flt() {
    leff = -999999.;
    lmean = -999999.;
    dwidth = -999999.;
    ab = -99.;
    tg = -99.;
    veg = -99.;
    msun = -99.;
    fcorr = -99.;
  }

  /// \brief generic constructor, with all internals set to unphysical defaults
  ///
  /// @param k index of the filter in the list of filters
  /// @param cname filter filename
  /// @param transt configuration keyword TRANS_TYPE
  /// @param calibt configuration keyword FILTER_CALIB
  flt(const int k, string cname, const int transt, const int calibt) : flt() {
    id = k;
    transtyp = transt;
    calibtyp = calibt;
    // if cname is a path
    if (cname.find_first_of("/") != string::npos) {
      // path is relative
      if (cname.at(0) != '/' && cname.at(0) != '.') {
        cname = lepharedir + "/filt/" + cname;
      }
    }
    name = cname;
    read(cname);
    compute_all();
  }

  /// \brief generic constructor, with all internals set to unphysical defaults
  ///
  /// @param k index of the filter in the list of filters
  /// @param filestream filter file stream
  /// @param transt configuration keyword TRANS_TYPE
  /// @param calibt configuration keyword FILTER_CALIB
  flt(const int k, ifstream &cname, const int transt, const int calibt)
      : flt() {
    id = k;
    transtyp = transt;
    calibtyp = calibt;
    read(cname);
    compute_all();
  }

  /// \brief constructor for a heavyside filter.
  ///
  /// The filter transmission is 1 between lmin and lmax, and 0 outside (a value
  /// of O is added at lmin-1 and lmax+1 to enforce it).
  /// @param lmin lower bound of the Heaviside function
  /// @param lmax upper bound of the Heaviside function
  /// @param nsteps number of steps in \f$\lambda\f$ from lmin to lmax
  flt(const double lmin, const double lmax, const int nsteps) : flt() {
    id = 0;
    transtyp = 0;
    calibtyp = 0;
    // Generate a heavyside filter. Transmission at 1. Not renormalized.
    // Delta lambda
    double dlamb = (lmax - lmin) / double(nsteps);
    // First element at T=0
    oneElLambda litBeg(lmin - 1, 0, 0);
    lamb_trans.push_back(litBeg);
    for (int k = 0; k <= nsteps; k++) {
      double lamb = lmin + double(k) * dlamb;
      oneElLambda litOne(lamb, 1, 0);
      lamb_trans.push_back(litOne);
    }
    // Last element at T=0
    oneElLambda litFin(lmax + 1, 0, 0);
    lamb_trans.push_back(litFin);
  }

  /* MP: erase all entries in lamb_trans */
  ~flt() { lamb_trans.clear(); }

  // Prototype of the functions declared in flt.cpp
  /// read a filter ascii file and store its content as a vector of oneElLambda
  /// elements
  void read(const string &fltFile);
  /// build the vector of oneElLambda elements out of a filter stream directly
  void read(ifstream &sfiltIn);
  /// Mean wavelength of the filter: \f$\lambda_{mean} = \frac{\int
  /// T(\lambda)\lambda d\lambda}{\int T(\lambda) d\lambda}\f$.
  double lambdaMean();
  /// full width at half maximum
  double width();
  /// maximum transmission value of the filter
  double peak();
  /// \brief effective wavelength based on the Vega spectrum as calibration SED.
  ///
  /// The Vega spectrum \f$V(\lambda)\f$ is included in LePhare
  /// ($LEPHAREDIR/vega/VegaLCB.sed). Then \f$\lambda_{eff} = \frac{\int
  /// V(\lambda)\, T(\lambda)\, \lambda\, d\lambda}{\int V(\lambda)\,
  /// T(\lambda)\, d\lambda}\f$.
  double lambdaEff();
  /// \brief effective wavelength based on a specific calibration SED.
  ///
  /// If \f$C(\lambda)\f$
  /// is a specific calibration SED defined by keyword FILTER_CALIB (see
  /// SED#generateCalib), then \f$\lambda_{eff2} = \frac{\int C(\lambda)\,
  /// T(\lambda)\, \lambda\, d\lambda}{\int C(\lambda)\, T(\lambda)\,
  /// d\lambda}\f$.
  double lambdaEff2();
  /// absolute magnitude of the Sun
  double magsun();
  /// Vega magnitude in this filter: \f$mag(Vega)
  /// = 2.5\cdot\log10\left(\frac{\int Vega(\lambda)\, T(\lambda)\,
  /// d\lambda}{\int T(\lambda)\, d\lambda}\right)\f$ where \f$Vega(\lambda)\f$
  /// is the Vega SED.
  double vega();

  /*! \brief compute the flux correction based on keyword FILTER_CALIB
   *
   * \f$fcorrec = B_0\cdot\lambda_{eff2}^2\cdot\frac{\int
   * T(\lambda)\,\lambda^{-2}\, d\lambda}{\int C(\lambda)\, T(\lambda)\,
   * d\lambda}\f$, where \f$C(\lambda)\f$ is a calibration function defined
   * through the value of FILTER_CALIB (see SED::generateCalib),
   * \f$\lambda_{eff2}\f$ is the alternative effective wavelength using
   * \f$C(\lambda)\f$ (see flt::lambdaEff2), and
   * - FILTER_CALIB=0 \f$B_0 = \lambda_{eff}^{-2}\f$
   * - FILTER_CALIB=1 \f$B_0 = \lambda_{eff}^{-1}\f$
   * - FILTER_CALIB=2 \f$B_0 = \lambda_{eff}^{-3}\f$
   * - FILTER_CALIB=3 \f$B_0 = Blackbody(10000K,\lambda_{eff})\f$
   * - FILTER_CALIB=4 \f$B_0 = Blackbody(10000K,\lambda_{eff})\f$
   * - FILTER_CALIB=5 \f$B_0 = \lambda_{eff}^{-3}\f$
   *
   * Pour les cas 4 et 5, \f$\lambda_{eff2}\f$ est évaluée en utilisant le
   * calibrator \f$C(\lambda)\f$ sélectionné par FILTER_CALIB=1.
   */
  double fcorrec();

  /// AB-Vega correction \f$ =  -2.5\,\log10\left(\frac{\int
  /// F(Vega)\,T(\lambda)\,d\lambda}{\int T(\lambda)\,
  /// c\lambda^{-2}\,d\lambda}\right)\f$
  double abcorr();

  /// Thuan Gunn correction \f$ =  2.5\,\log10\left(\frac{\int
  /// F(BD+17o4708)\,T(\lambda)\,d\lambda}{\int Vega(\lambda)\, T(\lambda)\,
  /// d\lambda}\right) + 9.5 - 0.03\f$
  double tgcorr();

  /// lowest stored lambda value
  double lmin() const { return lamb_trans.front().lamb; }
  /// highest stored lambda value
  double lmax() const { return lamb_trans.back().lamb; }

  void compute_all();
};

void write_output_filter(string &filtfile, string &filtdoc, vector<flt> vecFlt);
vector<flt> read_doc_filters(const string filtFile);

#endif
