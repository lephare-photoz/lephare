/*
  05/05/2015
  Class to store one source and store all its attributes
*/

// avoid multiple def of the same class
#ifndef SOURCE_H  // check that this keyword has been set already
#define SOURCE_H  // define the keyword to be checked

#include <array>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

#include "PDF.h"
#include "cosmology.h"
#include "flt.h"  // filter class
#include "globals.h"

using namespace std;

class SED;
struct SEDlight;

static vector<string> phys_par_names = {"AGE",  "LDUST", "LIR",  "MASS", "SFR",
                                        "SSFR", "COL1",  "COL2", "MREF"};
static vector<string> photoz_par_names = {"MIN_ZG", "MIN_ZQ", "BAY_ZG",
                                          "BAY_ZQ"};
// How to associate a name to the corresponding interger
static unordered_map<string, int> maptype = {
    {"MASS", 0},    {"SFR", 1},     {"SSFR", 2},   {"LDUST", 3}, {"LIR", 4},
    {"AGE", 5},     {"COL1", 6},    {"COL2", 7},   {"MREF", 8},  {"MIN_ZG", 9},
    {"MIN_ZQ", 10}, {"BAY_ZG", 11}, {"BAY_ZQ", 12}};

//! @class onesource
/*!
represents an object from a catalogue, and manages its fitting to an SED.
*/
class onesource {
 private:
  bool verbose;

 public:
  /// Best-fit physical parameters for this source (see phys_par_names for
  /// the keys), initialised to #INVALID_PHYS and filled by
  /// compute_best_fit_physical_quantities()
  unordered_map<string, double> results = {
      {"MASS_BEST", INVALID_PHYS},    {"SFR_BEST", INVALID_PHYS},
      {"SSFR_BEST", INVALID_PHYS},    {"LDUST_BEST", INVALID_PHYS},
      {"LUM_TIR_BEST", INVALID_PHYS}, {"AGE_BEST", INVALID_PHYS},
      {"EBV_BEST", INVALID_PHYS},     {"EXTLAW_BEST", INVALID_PHYS},
      {"LUM_NUV_BEST", INVALID_PHYS}, {"LUM_R_BEST", INVALID_PHYS},
      {"LUM_K_BEST", INVALID_PHYS},
  };
  /// Emission line flux (`EM_FLUX_<line>`) and equivalent width
  /// (`EM_EW_<line>`) for the best-fit galaxy template, filled by
  /// computeEmFlux()
  unordered_map<string, double> results_emission_lines;

      dm,           ///< unused legacy member (kept for ABI/output-format
                    ///< stability; see SED::dm for the per-template scaling
                    ///< actually used in the fit)
  long cont,     ///< input context: bitmask of bands used, as read from the
                 ///< catalogue (see globals.h CHECK_CONTEXT_BIT)
      new_cont;  ///< context rebuilt after masking bands with invalid flux
                 ///< or error, see fltUsed()
  vector<double> ab,     ///< observed flux, one value per band
      sab,               ///< uncertainty on #ab
      mab,               ///< observed AB magnitude, one value per band
      msab,              ///< uncertainty on #mab
      magm,              ///< predicted (model) apparent magnitude for the
                         ///< best-fit template, one value per band
      magm0,             ///< predicted magnitude at z=0 for the best-fit
                         ///< template, one value per band
      absmagPred,        ///< predicted absolute magnitude, one value per
                         ///< band
      magPred,           ///< predicted apparent magnitude, one value per
                         ///< band
      kap,               ///< k-correction, one value per band
      mabs,              ///< absolute magnitude, one value per band (see
                         ///< absmag())
      emabs,             ///< uncertainty on #mabs, one value per band
      ab_ori,            ///< copy of #ab before the MW-dust / zero-point
                         ///< corrections, see keepOri()
      sab_ori,           ///< copy of #sab before the corrections applied
                         ///< to #ab_ori, see keepOri()
      mab_ori,           ///< copy of #mab before the corrections applied to
                         ///< #ab_ori, see keepOri()
      abIR,              ///< observed flux used for the FIR fit, i.e. #ab
                         ///< minus the stellar predicted flux where
                         ///< applicable (see subtract_stellar_component())
      sabIR;             ///< uncertainty on #abIR
  vector<int> busnorma,  ///< per-band flag: 1 if the band is used in the
                         ///< main fit, 0 otherwise (see fltUsed())
      busul,             ///< per-band flag: 1 if the band is treated as an
                         ///< upper limit, 0 otherwise
      busfir,            ///< per-band flag: 1 if the band is used in the
                         ///< FIR fit, 0 otherwise (see fltUsedIR())
      bscfir,            ///< per-band flag: 1 if the band is used to scale
                         ///< the FIR fit, 0 otherwise
      absfilt;           ///< index of the filter used to compute the
                         ///< absolute magnitude, one value per element of
                         ///< #mabs (see absmag())
  string spec,           ///< identifier of the source (catalogue "spec" column)
      str_inp;    ///< additional free-form input string carried through to
                  ///< the output (see readsource())
  int pos,        ///< position (row index) of the source in the catalogue
      nbused,     ///< number of bands used in the main fit
      nbul,       ///< number of bands treated as upper limits
      nbusIR,     ///< number of bands used to scale the FIR fit
      indminSec,  ///< index in the library of the secondary chi2 minimum
                  ///< (see secondpeak())
      indminIR,   ///< index in the FIR library of the chi2 minimum
      imasminIR;  ///< model number of the FIR fit minimum
  double zs,      ///< spectroscopic redshift, as read from the catalogue
      consiz;     ///< redshift adopted for this source when computing
                  ///< derived/rest-frame quantities (e.g. absolute
                  ///< magnitudes, k-corrections); typically the best-fit
                  ///< or median photo-z depending on the calling context
  array<double, 3> zmin,  ///< redshift of the chi2 minimum, for GAL/QSO/STAR
                          ///< (indices 0/1/2 respectively)
      chimin,             ///< chi2 of the minimum, for GAL/QSO/STAR
      dmmin;              ///< template scaling of the minimum, for
                          ///< GAL/QSO/STAR
  array<int, 3> indmin,   ///< index in the library of the chi2 minimum, for
                          ///< GAL/QSO/STAR
      imasmin;            ///< model (template) number of the chi2 minimum,
                          ///< for GAL/QSO/STAR
  double zminIR,          ///< redshift of the FIR fit chi2 minimum
      chiminIR,           ///< chi2 of the FIR fit minimum
      dmminIR;            ///< template scaling of the FIR fit minimum
  array<double, 4>
      priorLib;  ///< absolute magnitude prior range, as [bright,faint] for
                 ///< the galaxy library followed by [bright,faint] for the
                 ///< AGN library (see setPriors())

  vector<double> chibay;
  vector<double> gridzg, gridLdustIR, gridEbv, gridLIR;
  PDF PDFebv;

  /// Marginalized redshift PDF summary for galaxy templates: median (resp.
  /// chi2-minimum, PDF mode) at index 0, followed by the 68/90/99% credible
  /// interval [low,high] bounds at indices 1-6 (see generatePDF())
  vector<double> zgmed,  ///< median of the marginalized PDF, and CI bounds
      zgmin,             ///< chi2-minimum solution, and CI bounds
      zgmode;            ///< mode of the marginalized PDF, and CI bounds
  /// Same as zgmed/zgmin/zgmode, for the QSO/AGN library
  vector<double> zqmed,  ///< median of the marginalized PDF, and CI bounds
      zqmin,             ///< chi2-minimum solution, and CI bounds
      zqmode;            ///< mode of the marginalized PDF, and CI bounds
  /// Marginalized PDF summary (median, then 68/90/99% credible interval
  /// bounds, same layout as #zgmed) of, respectively, log stellar mass,
  /// log SFR, log sSFR, log age, log dust luminosity, and the first and
  /// second rest-frame colors (see generatePDF())
  vector<double> massmed,  ///< log stellar mass, and CI bounds
      SFRmed,              ///< log SFR, and CI bounds
      sSFRmed,             ///< log sSFR, and CI bounds
      agemed,              ///< log age, and CI bounds
      Ldustmed,            ///< log dust luminosity, and CI bounds
      col1med,             ///< first rest-frame color, and CI bounds
      col2med,             ///< second rest-frame color, and CI bounds
      ebvmed,              ///< E(B-V), same layout as #massmed
      Mrefmed;  ///< reference absolute magnitude, same layout as #massmed
  /// log IR luminosity marginalized PDF summary (median, then 68/90/99%
  /// credible interval bounds), filled by uncertaintiesBayIR()
  array<double, 7> LIRmed;

  /// Emission-line flux of the best-fit template rescaled to this source,
  /// indexed as in SED::fac_line (see computeEmFlux())
  array<double, 65> fluxEL_SED = {0};
  double limits_zmax = 20.;  ///< faint-end redshift limit used by limits()
  double limits_Mfaint = 0;  ///< faint absolute magnitude limit, computed by
                             ///< limits()
  /// Marginalized PDF for each physical/redshift quantity, keyed by the
  /// index defined in #maptype (e.g. 0 for stellar mass, 9 for the galaxy
  /// redshift); see the onesource(pos, gridz) constructor for the binning
  unordered_map<int, PDF> pdfmap;
  double zsec,     ///< redshift of the secondary chi2 minimum
      zsecChi2,    ///< chi2 of the secondary minimum
      zsecEbv,     ///< E(B-V) of the secondary minimum
      zsecScale,   ///< template scaling of the secondary minimum
      zsecProb,    ///< PDF probability of the secondary minimum
      zsecAge;     ///< age of the secondary minimum
  int zsecMod,     ///< model (template) number of the secondary minimum
      zsecExtlaw;  ///< extinction law index of the secondary minimum

  /// Milky Way E(B-V) for this source (galactic coordinates), used by
  /// correct_classic_mw()/correct_galametz_mw(); -99. if not set
  double mw_ebv = -99.;

  // Minimal constructor of the source
  onesource() {
    spec = "1";      // ident
    zs = INVALID_Z;  // spectroscopic redshift
    cont = 0;        // context
    str_inp = ' ';
    for (int k = 0; k < 3; k++) {
      zmin[k] = INVALID_Z;
      indmin[k] = INVALID_INDEX;
      chimin[k] = HIGH_CHI2;
      imasmin[k] = INVALID_INDEX;
      dmmin[k] = 0.;
    }
    zminIR = INVALID_Z;
    indminIR = INVALID_INDEX;
    chiminIR = HIGH_CHI2;
    imasminIR = INVALID_INDEX;
    dmminIR = 0.;
    LIRmed.fill(INVALID_PHYS);
    nbused = 0;
    pos = 0;
  }

  // Need to initialize the PDF in the constructor after the ":"
  /*! Build a source at a given catalogue position and initialise #pdfmap
   * with the physical-parameter grids and a redshift grid derived from
   * @p gridz
   * @param posC: position (row index) of the source in the catalogue,
   * stored in #pos
   * @param gridz: redshift grid of the SED library, used to set the range
   * and step of the redshift PDFs (a single-element grid is treated as the
   * STAR-only case)
   */
  onesource(const int posC, const vector<double>& gridz) : onesource() {
    pos = posC;  // position in the file

    // 0:["MASS"] / 1:["SFR"] / 2:["SSFR"] / 3:["LDUST"] / 4:["LIR"] / 5:["AGE"]
    // / 6:["COL1"] / 7:["COL2"] / 8:["MREF"]/ 9:["MIN_ZG"] / 10:["MIN_ZQ"] /
    // 11:["BAY_ZG"] / 12:["BAY_ZQ"]
    pdfmap[0] =
        PDF(3., 13, 201);  // log Stellar Mass  [Mo] : 3 to 13 by 0.05 step
    pdfmap[1] = PDF(-10., 5, 151);  // log SFR [Mo/yr] : -10 to 5 by 0.1 step
    pdfmap[2] =
        PDF(-25., -5, 201);  // log SFR/Mass [SFR/Mo] : -25 to -5 by 0.1 step
    pdfmap[3] = PDF(6., 14,
                    161);  // Ldust [L0]  ! from mag_gal (EB-V, Extlaw)  6 to 14
                           // by 0.05 step
    pdfmap[4] = PDF(
        6., 14,
        161);  // Llir [L0]   ! from  sedtolib (*.phys)  6 to 14 by 0.05 step
    pdfmap[5] = PDF(6, 10.5, 91);    // log Age  [yr]  ! 6 -> 10.5 by 0.05 step
    pdfmap[6] = PDF(-10., 10, 201);  // rest-frame color : -10 to 10 by 0.1 step
    pdfmap[7] = PDF(-10., 10, 201);  // rest-frame color : -10 to 10 by 0.1 step
    pdfmap[8] = PDF(-30., -5, 250);  // Mref : -30 to -5 by 0.1 step

    int zgnb;
    double zgmax;
    // case where we only treat STAR libraries
    if (gridz.size() == 1) {
      zgmax = EPS_Z;
      zgnb = 2;
    } else {
      // find the smallest interval in the grid provided by the user
      vector<double> diff(gridz.size(), 1.0);
      adjacent_difference(gridz.begin(), gridz.end(), diff.begin());
      zgmax = gridz.back();
      // the first and last element of adjacent_diff need to be removed from
      // this difference double zgstep = *min_element(diff.begin()+1,
      // diff.end()-1);
      double zgstep = *min_element(diff.begin() + 1, diff.end() - 1);
      zgnb = int(zgmax / zgstep) + 1;
    }
    // create PDF with a linear z grid.
    pdfmap[9] = PDF(0., zgmax, zgnb);
    pdfmap[10] = PDF(0., zgmax, zgnb);
    pdfmap[11] = PDF(0., zgmax, zgnb);
    pdfmap[12] = PDF(0., zgmax, zgnb);
  }

  // erase all entries in onesource
  ~onesource() {
    chibay.clear();
    ab.clear();
    sab.clear();
    abIR.clear();
    sabIR.clear();
    mab.clear();
    msab.clear();
    kap.clear();
    mabs.clear();
    absfilt.clear();
    busnorma.clear();
  }

  //! Set verbosity
  /*!
    @param v bool specifying whether output should be verbose or not.
   */
  inline void set_verbosity(const bool v) { verbose = v; }
  //! Get verbosity
  /*!
    @return bool specifying whether output is verbose or not.
   */
  inline bool get_verbosity() const { return verbose; }

  // Prototype
  /*! Set the observed fluxes and their errors, and the identifying metadata,
   * from a catalogue row
   * @param identifier: source identifier, stored in #spec
   * @param vals: observed flux (or magnitude, converted upstream) per band,
   * stored in #ab
   * @param err_vals: uncertainty per band, stored in #sab
   * @param context: bitmask of bands to use, stored in #cont
   * @param z_spec: spectroscopic redshift, stored in #zs
   * @param additional_input: free-form string carried through to the
   * output, stored in #str_inp
   */
  void readsource(const string& identifier, const vector<double> vals,
                  const vector<double> err_vals, const long context,
                  const double z_spec, const string additional_input);
  /*! Set the absolute magnitude prior range used by nzprior(), for the
   * galaxy and AGN libraries
   * @param magabsB: [galaxy,AGN] bright absolute magnitude bound
   * @param magabsF: [galaxy,AGN] faint absolute magnitude bound
   */
  void setPriors(const array<double, 2> magabsB,
                 const array<double, 2> magabsF);
  /*! Flag, in #busnorma and #busul, which bands are used in the main fit and
   * which are treated as upper limits, based on the context and on the
   * validity of #ab/#sab
   * @param gbcont: global context overriding the per-source context if >= 0
   * @param contforb: context of bands forbidden from being upper limits
   */
  void fltUsed(const long gbcont, const long contforb);
  /*! Flag, in #busfir and #bscfir, which bands are used in the FIR fit and
   * for its flux scaling
   * @param fir_cont: context of bands eligible for the FIR fit
   * @param fir_scale: context of bands used to scale the FIR fit
   * @param allFilters: filter set, used to check the rest-frame wavelength
   * against @p fir_lmin
   * @param fir_lmin: minimum rest-frame wavelength for a band to be eligible
   */
  void fltUsedIR(const long fir_cont, const long fir_scale,
                 vector<flt> allFilters, const double fir_lmin);
  /*! Convert the input magnitudes/fluxes (#ab/#sab) to fluxes in the LePHARE
   * internal convention, handling invalid/negative values
   * @param catmag: magnitude system of the input catalogue ("AB" or "VEGA")
   * @param allFilters: filter set, used for the VEGA-to-AB conversion
   */
  void convertFlux(const string& catmag, const vector<flt> allFilters);
  /*! Rescale the flux errors #sab, e.g. to add a systematic floor
   * @param min_err: minimum relative error per band (or a single value
   * applied to all bands)
   * @param fac_err: multiplicative factor applied to the error per band (or
   * a single value applied to all bands)
   */
  void rescale_flux_errors(const vector<double> min_err,
                           const vector<double> fac_err);

  /*! Fit the source against the given library by chi2 minimisation
   * @param lightLib: light SED library to fit against
   * @param flux: predicted flux of each library template, in each band
   * @param valid: indices, in @p lightLib, of the templates to consider
   * @param funz0: scaling factor. flat prior in flux, or -1 to disable it
   * @param bp: [min,max] index of the bands used for the fit
   * @param restrict: if true, treat a negative predicted flux (outside the
   * SED's rest-frame wavelength coverage) as zero rather than using it as-is
   */
  void fit(SEDlight& lightLib, const vector<vector<double>>& flux,
           const vector<size_t>& valid, const double& funz0,
           const array<int, 2>& bp, const bool restrict);
  /*! Fit the source against a far-infrared library, in the bands flagged by
   * fltUsedIR()
   * @param fulllib: FIR SED library to fit against
   * @param flux: predicted flux of each library template, in each band
   * @param valid: indices, in @p fulllib, of the templates to consider
   * @param fit_frsc: FIR flux-scaling method
   * @param lcdm: cosmology used to rescale distances
   */
  void fitIR(vector<SED*>& fulllib, const vector<vector<double>>& flux,
             const vector<size_t>& valid, const string fit_frsc, cosmo lcdm);
  /*! Compute the N(z) prior weight applied to the chi2 during the fit,
   * following Ilbert et al. (2006)
   * @param luv: rest-frame UV luminosity of the template
   * @param lnir: rest-frame NIR luminosity of the template
   * @param reds: redshift of the template
   * @param bp: [min,max] index of the bands used to define the apparent
   * magnitude the prior is conditioned on
   * @return the prior weight, to be combined with the chi2
   */
  double nzprior(const double luv, const double lnir, const double reds,
                 const array<int, 2> bp);
  /*! Iteratively remove the band contributing most to the chi2 of the
   * best-fit solution, and re-fit, as long as this improves the chi2 by more
   * than @p thresholdChi2
   * @param lightLib: light SED library to fit against
   * @param flux: predicted flux of each library template, in each band
   * @param valid: indices, in @p lightLib, of the templates to consider
   * @param funz0: scaling factor prior, forwarded to fit()
   * @param bp: [min,max] index of the bands used for the fit
   * @param thresholdChi2: minimum chi2 improvement required to discard a
   * band
   * @param restrict: forwarded to fit(), see there
   */
  void rm_discrepant(SEDlight& lightLib, const vector<vector<double>>& flux,
                     const vector<size_t>& valid, const double funz0,
                     const array<int, 2> bp, double thresholdChi2,
                     const bool restrict);
  void deredden_observed_mag(const vector<double>& ext_values);

  /*! Apply a per-template, per-band reddening correction to a predicted
   * flux array
   * @param flux: predicted flux of each library template, in each band
   * @param reddening: per-template, per-band reddening correction to apply
   * @return the reddened flux array
   */
  vector<vector<double>> redden_flux(
      const vector<vector<double>>& flux,
      const vector<vector<double>>& reddening) const;

  /*! Write output in the lephare ascii format
   * @param stout: stream object pointing to the output file
   * @param outkeywords: list of keywords to be output
   */
  void write_out(ofstream& stout, const vector<string>& outkeywords);
  /*! Write the header of the PDF output file(s)
   * @param pdztype: list of PDF types to write a header for (keys of #maptype)
   * @param stpdz: map, keyed by PDF type, of the open output streams
   */
  void write_pdz_header(vector<string> pdztype,
                        unordered_map<string, ofstream>& stpdz);
  /*! Write this source's marginalized PDF(s) to the corresponding output
   * stream(s)
   * @param pdztype: list of PDF types to write (keys of #maptype)
   * @param stpdz: map, keyed by PDF type, of the open output streams
   */
  void write_pdz(vector<string> pdztype,
                 unordered_map<string, ofstream>& stpdz);
  /// Convert the fluxes #ab/#sab to AB magnitudes #mab and errors #msab
  void convertMag();
  /// Save a copy of #ab, #sab and #mab, before corrections, into
  /// #ab_ori/#sab_ori/#mab_ori
  void keepOri();

  //! Update the solution of the fit based on execution flags
  /*!
   * @param zfix: bool that sets whether to set solution to a given redshift,
   * typically a true or spectroscopic redshift
   * @param zintp: bool that sets whether to improve the determination
   * of the minimum on the chi2 curve, by parabolic approximation
   * @param lcdm: `cosmo` object to access the `distMod`
   * and `distDet` functions
   *
   * #zmin and #dmmin are updated in place, for the GAL and QSO solutions.
   * Note that zfix and zintp are not supposed to both be set. In case it
   * happens, zintp is discarded here.
   */
  void interp(const bool zfix, const bool zintp, const cosmo& lcdm);
  /// Compute the frequentist (chi2-based) confidence intervals around the
  /// GAL and QSO chi2-minimum redshifts, filling #zgmin/#zqmin
  void uncertaintiesMin();
  /// Compute the Bayesian confidence intervals (median and credible
  /// intervals) of the redshift and physical parameters from their
  /// marginalized PDFs, filling #zgmed/#zqmed/#massmed/etc.
  void uncertaintiesBay();
  /// Compute the Bayesian confidence intervals of the IR luminosity from
  /// its marginalized PDF (#pdfmap[4]), filling #LIRmed
  void uncertaintiesBayIR();
  /*! Detect a secondary peak in the marginalized galaxy redshift PDF and,
   * if found, store its properties in #zsec and related members
   * @param lightLib: light SED library the source was fit against
   * @param dz_win: minimum redshift separation from the main peak for a
   * secondary peak to be considered
   * @param min_thres: minimum relative probability for a secondary peak to
   * be considered
   */
  void secondpeak(SEDlight& lightLib, const double dz_win,
                  const double min_thres);
  /*! Build the marginalized PDF of the redshift and of the physical
   * parameters (mass, SFR, age, colors...) by summing the chi2-based
   * probability of every valid template, and fill #pdfmap
   * @param lightLib: light SED library the source was fit against
   * @param va: indices, in @p lightLib, of the valid templates
   * @param colAnalysis: whether to also build the rest-frame color PDFs
   * (#pdfmap[6]/[7])
   * @param zfix: whether the redshift was fixed for this fit
   */
  void generatePDF(SEDlight& lightLib, const vector<size_t>& va,
                   const bool colAnalysis, const bool zfix);
  /// Build the marginalized IR luminosity PDF (#pdfmap[4]) from the FIR fit
  /// @param fulllib: FIR SED library the source was fit against
  void generatePDF_IR(vector<SED*>& fulllib);
  /// Compute the mode of the marginalized GAL/QSO redshift PDFs and their
  /// confidence intervals, filling #zgmode/#zqmode
  void mode();
  /*! Interpolate the predicted magnitudes (#magm) of the best-fit template
   * at the redshift #consiz, between the two adjacent grid points of
   * @p fulllib
   * @param fulllib: SED library the source was fit against
   * @param flux: predicted flux of each library template, in each band
   */
  void interp_lib(vector<SED*>& fulllib, const vector<vector<double>>& flux);
  /*! Apply a per-band zero-point offset to the observed magnitudes/fluxes
   * @param a0: magnitude offset to add, one value per band
   */
  void adapt_mag(vector<double> a0);
  /*! Apply the "classic" Milky Way dust correction to the observed flux and
   * error, using a per-band correction coefficient (galaxy-independent)
   * @param Alamb_corr: MW dust correction coefficient, one value per band
   * @param mw_global_ebv: global MW E(B-V) to use for every source if >= 0,
   * overriding #mw_ebv; if negative, #mw_ebv is used instead
   */
  void correct_classic_mw(const vector<double>& Alamb_corr,
                          const double mw_global_ebv);
  /*! Apply the Galametz et al. Milky Way dust correction to the observed
   * flux and error, using the model-dependent reddening of the best-fit
   * (chi2-minimum) galaxy template. No-op if there is no valid galaxy
   * solution or if #mw_ebv is not positive.
   * @param reddening: per-template, per-band reddening correction, indexed
   * as #indmin
   */
  void correct_galametz_mw(const vector<vector<double>>& reddening);
  /*! Allow for stellar component substraction before fitting an IR template
   * When fitting an IR component after the nominal fit, there is an interval
   * in wavelength where both nominal and IR template would contribute.
   * In order to correctly fit the IR template in this interval, one may want
   * to subtract first the nominal stellar component as obtained from the
   * best fit template. This common interval in lambda (in the rest frame)
   * is bounded by the hardcoded value of 250 um at the high end, and by the
   * value set by the FIR_LMIN keyword (in um, default 7 um) at the low end.
   * @param substar : bool value set by the FIR_SUBSTELLAR keyword, defining
   * whether to do this subtraction or not.
   * @param allFilters : the list of filters, needed to discard filters
   * which have \f$\lambda_{mean}/(1+z) > 250\,\mu m\f$.
   */
  void subtract_stellar_component(const bool substar, vector<flt> allFilters);
  /*! Compute the absolute magnitude(s) of the source from its best-fit
   * template(s), filling #mabs/#emabs/#absfilt
   * @param bestFlt: per band-pass, indices of the filters bracketing the
   * rest-frame reference wavelength at the source's redshift
   * @param maxkcolor: maximum rest-frame color allowed for the k-correction
   * @param lcdm: cosmology used to compute the distance modulus
   * @param gridz: redshift grid of the library
   */
  void absmag(const vector<vector<int>>& bestFlt,
              const vector<vector<double>>& maxkcolor, cosmo lcdm,
              const vector<double> gridz);
  /*! Write the best-fit spectra (and, if available, the FIR spectrum) of
   * this source to per-source output files
   * @param fulllib: SED library the source was fit against
   * @param fulllibIR: FIR SED library the source was fit against
   * @param lcdm: cosmology used to compute the spectra
   * @param allFilters: filter set, written alongside the spectra
   * @param outspdir: output directory for the spectra files
   */
  void writeSpec(vector<SED*>& fulllib, vector<SED*>& fulllibIR, cosmo lcdm,
                 const vector<flt>& allFilters, const string outspdir) const;
  /*! Write out in a Id<source id>.chi file the chi2 of all the templates
   * participating to the fit.
   * @param lightLib : the light library of SED objects.
   */
  void writeFullChi(const SEDlight& lightLib);
  /*! Compute the predicted apparent magnitudes of the best-fit template(s)
   * in an arbitrary set of additional filters, filling #magPred
   * @param fulllib: SED library the source was fit against
   * @param lcdm: cosmology used to compute the magnitudes
   * @param allFltAdd: additional filters to predict magnitudes for
   */
  void computePredMag(vector<SED*>& fulllib, cosmo lcdm, vector<flt> allFltAdd);
  /*! Compute the predicted absolute magnitudes of the best-fit template(s)
   * in an arbitrary set of additional filters, filling #absmagPred
   * @param fulllib: SED library the source was fit against
   * @param lcdm: cosmology used to compute the magnitudes
   * @param allFltAdd: additional filters to predict magnitudes for
   */
  void computePredAbsMag(vector<SED*>& fulllib, cosmo lcdm,
                         vector<flt> allFltAdd);
  /*! Compute the emission-line flux and equivalent width of the best-fit
   * galaxy template rescaled to this source, filling #fluxEL_SED and
   * #results_emission_lines
   * @param fulllib: SED library the source was fit against
   * @param lcdm: cosmology used to rescale the flux
   */
  void computeEmFlux(vector<SED*>& fulllib, cosmo lcdm);
  /*! Compute the faint-end absolute-magnitude limit reachable for this
   * source given the survey depth, filling #limits_zmax/#limits_Mfaint
   * @param fulllib: SED library the source was fit against
   * @param limits_zbin: redshift bins in which the limit is computed
   * @param limits_ref: reference filter index
   * @param limits_sel: selection filter indices
   * @param limits_cut: flux limit cut per selection filter
   */
  void limits(vector<SED*>& fulllib, vector<double>& limits_zbin,
              int limits_ref, vector<int>& limits_sel,
              vector<double>& limits_cut);
  /*! Compute the flux-averaged best-fit spectrum (and its uncertainty) of a
   * given solution, between two wavelengths
   * @param sol: solution to use: 0/1 for the primary/secondary galaxy
   * solution, 2 for the FIR solution, 3 for QSO, 4 for STAR
   * @param fulllib: SED library the source was fit against
   * @param lcdm: cosmology used to compute the spectrum
   * @param minl: lower wavelength bound
   * @param maxl: upper wavelength bound
   * @return [flux, flux uncertainty] averaged over the requested range
   */
  pair<vector<double>, vector<double>> best_spec_vec(short sol,
                                                     vector<SED*>& fulllib,
                                                     cosmo lcdm, double minl,
                                                     double maxl) const;

  /// Compute the physical quantities (mass, SFR, luminosities...) of the
  /// best-fit template(s), filling #results
  /// @param fulllib: SED library the source was fit against
  void compute_best_fit_physical_quantities(vector<SED*>& fulllib);
};

#endif
