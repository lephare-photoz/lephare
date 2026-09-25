#ifndef PHOTOZ_LIB_H
#define PHOTOZ_LIB_H

#include <string>  // use string
#include <vector>  // manipulate vector

#include "SED.h"  // to read the libraries
#include "cosmology.h"
#include "flt.h"  // to read the libraries
#include "mag.h"

/*! \brief Class managing photo-z computation.
 *
 * The PhotoZ class is the central executor of photo-z estimation : it manages
 * most configuration parameters, reads synthetic magnitudes inputs and object
 * catalogs, performs the fit, and is in charge of producing results and saving
 * them to files.
 */
class PhotoZ {
 private:
  keymap keys;
  unsigned int rowmin, rowmax;
  int cat_fmt, fl_auto, method, babs;
  size_t nlibext;
  array<int, 2> bp;
  long gbcont, contforb, bdscal;
  double funz0, adzmin, adzmax, auto_thresmin, auto_thresmax, min_thres, dz_win;
  double fir_lmin, fir_cont, fir_scale;
  bool verbose, restrict_rf = false;
  array<double, 2> magabsB, magabsF, zrange, ebvrange;
  bool outchi,  /// Whether or not to output th chi2 values in ascii
      zintp,    /// Whether or not to interpolate the z solution from the grid
                /// result
      zfix,     /// Whether or not to run with a fixed redshift, typically the
                /// true/spectro z
      methz,    /// If true, set z to the MEDIAN solution rather than the BEST
                /// solution, when computing physical parameters.
      mw_galametz = false,  /// Whether or not to apply a Milky Way extinction
                            /// correction to the SEDs with Galametz method
      mw_classic_extinction = false,
      /// apply equal extinction to all SEDs, as opposed
      /// to a different value for each SED
      colAnalysis;  /// if true, measure the PDF to get unceratinties
                    /// associated to the rest-frame colors.

  string cat, typm, catmag, cattyp, zmulti, outf, outsp, outpdz, outpdm;
  vector<double> shifts0, min_err, fac_err, int_pdz, zbmin, zbmax;
  vector<flt> allFiltersAdd;
  vector<vector<int>> goodFlt;
  vector<vector<double>> maxkcol;
  vector<long> magabscont;
  vector<string> colib;
  vector<int> bapp, bappOp, pdz_fabs, emMod;
  cosmo lcdm;
  vector<opa> opaOut;
  double mw_global_ebv;
  string mwExtCurve;

 public:
  vector<vector<double>> flux,  ///< predicted flux of each GAL/QSO library
                                ///< template, in each band, after any MW
                                ///< dust correction
      flux_no_mw,               ///< same as #flux, before MW dust correction
      fluxIR,                   ///< predicted flux of each FIR library
                                ///< template, in each band
      reddening;                ///< per-template, per-band Galametz MW
                                ///< reddening correction (see
                                ///< onesource::redden_flux)
  string mw_ref_mod;     ///< path (relative to $LEPHAREDIR) of the reference
                         ///< stellar SED used to normalise the Milky Way
                         ///< extinction correction (MW_REFERENCE_MODEL keyword)
  vector<double> zLib,   ///< redshift of each GAL/QSO library template
      zLibIR;            ///< redshift of each FIR library template
  vector<SED*> fullLib,  ///< GAL/QSO/STAR library templates (owned)
      fullLibIR;         ///< FIR library templates (owned)
  SEDlight lightLib;     ///< lightweight (memory-reduced) copy of #fullLib used
                         ///< during the fit
  vector<flt> allFilters;      ///< filter set the libraries were built with
  vector<double> gridz;        ///< redshift grid of the GAL/QSO libraries
  vector<string> outkeywords,  ///< requested output column keywords
                               ///< (OUTPUT_PARA / CAT_OUT_PARA)
      pdftype;                 ///< requested output PDF types (PDZ_OUT)
  int imagm;                   ///< number of filters/bands
  string outputHeader,         ///< accumulated header text describing the run
                        ///< configuration, written atop the output catalogue
      outpara;  ///< accumulated documentation of the output columns
  bool one_mw_ebv =
      false;  ///< Whether or not a single MW E(B-V) value is applied to
              ///< all sources, as opposed to a different value for each
              ///< source read from a file
  vector<double>
      mw_classic_extinction_values;  ///< If mw_classic_extinction is true, this
                                     ///< vector contains the extinction values
                                     ///< to apply to each filter for all SEDs.
                                     ///< It is computed once at the beginning
                                     ///< of the code, based on the selected
                                     ///< Milky Way extinction curve and the
                                     ///< filter transmission curves.
  /*! Build a PhotoZ instance: parse every keyword needed to run the fit
   * (input catalogue format, priors, Milky Way extinction options, output
   * format...), then read and merge the GAL/QSO/STAR (and, if configured,
   * FIR) binary magnitude libraries built beforehand by mag_gal, checking
   * their cosmology/redshift-grid consistency
   * @param key_analysed: map of keyword/value pairs
   */
  PhotoZ(keymap& key_analysed);

  // LCOV_EXCL_START
  virtual ~PhotoZ() {
    for (auto& sed : fullLib) delete sed;
    for (auto& sed : fullLibIR) delete sed;
    fullLib.clear();
    fullLibIR.clear();
  }
  // LCOV_EXCL_STOP

  /*! Determine the per-band magnitude zero-point offsets (a0) to apply
   * before fitting, either from the APPLY_SYSSHIFT keyword if it matches
   * the number of filters, or by running the auto-adaptation procedure
   * (AUTO_ADAPT=YES) on @p adaptSources, or 0 for every band otherwise
   * @param adaptSources: sources (typically with a reliable spec-z) used
   * for the auto-adaptation procedure
   * @return the per-band offset, one value per filter
   */
  vector<double> compute_offsets(vector<onesource*> adaptSources);

  /*! Iteratively fit @p adaptSources at their spec-z and derive, for each
   * band, the median magnitude offset between the observed and predicted
   * magnitudes, until convergence or 10 iterations
   * @param adaptSources: sources (typically with a reliable spec-z) used
   * for the auto-adaptation procedure
   * @return the per-band offset, one value per filter
   */
  vector<double> run_autoadapt(vector<onesource*> adaptSources);

  /*! Fit every source in the given list against the GAL/QSO/STAR (and, if
   * configured, FIR) libraries, compute their uncertainties and physical
   * parameters, and write the per-source outputs (ascii catalogue, spectra,
   * PDFs) as configured
   * @param sources: sources to fit
   * @param a0: per-band magnitude zero-point offset to apply before fitting
   */
  void run_photoz(vector<onesource*> sources, const vector<double>& a0);

  /*! Fit a source based on the PhotoZ configuration
   * \param source: the onesource object under consideration
   * \param mag_shifts: magnitude shifts to apply prior to fit (zero points)
   */
  void fit(onesource& source, const vector<double>& mag_shifts);

  /*! Compute the fit uncertainties of a source
   * \param source: the onesource object under consideration
   */
  void fit_uncertainties(onesource& source);

  /*! Compute the physical parameters derived from the fit to a source
   * \param src: the onesource object under consideration
   */
  void physical_parameters(onesource& src);

  /*! Get the spectrum of one of the template of one solution to the fit of
   * a source
   * \param source: the onesource object under consideration
   * \param templateType: which object type (GAL:0, GAL 2nd solution:1, GAL
   * FIR:2, QSO:3, STAR:4)
   * \param minl: lower bound in wavelength for the returned spectrum
   * \param maxl: upper bound in wavelength for the returned
   * spectrum
   *
   * \return a pair of vectors for the wavelengths and spectral energy density.
   */
  pair<vector<double>, vector<double>> best_template(onesource& source,
                                                     int const templateType,
                                                     double const minl,
                                                     double const maxl);

  /*! Write out the spectrum solutions to ascii
   * \param src: the onesource object under consideration
   * \param outputDir: output directory
   * The ascii file will be outputDir/Id<source.spec>.spec
   */
  inline void write_spectrum(const onesource& src,
                             const string outputDir = ".") {
    src.writeSpec(fullLib, fullLibIR, lcdm, allFilters, outputDir);
  }

  /*! Build the header line(s) of the output ascii catalogue from the list
   * of requested output keywords
   * @param outkeywords: output column keywords (OUTPUT_PARA / CAT_OUT_PARA)
   * @return the formatted header, ready to be written to the output file
   */
  string prep_header(vector<string> outkeywords);

  /*! Write the per-source outputs (ascii catalogue line, spectra, PDFs,
   * .chi files) for every source in the list, as configured by the
   * relevant keywords (CAT_OUT, SPEC_OUT, PDZ_OUT, FULL_CHI_OUT...)
   * @param sources: sources to write out
   */
  void write_outputs(vector<onesource*> sources);

  /*! Read a binary SED/magnitude library (and its .doc file) into memory,
   * checking cosmology/redshift-grid consistency with any library already
   * read via check_consistency()
   * @param libFull: SED library to append the read templates to
   * @param ind: running count of templates read so far, updated in place
   * @param nummodpre: per-type (GAL/QSO/STAR) running count of models read
   * so far, used to renumber templates across successive libraries
   * @param libName: base name of the library to read (in
   * $LEPHAREWORK/lib_mag/)
   * @param filtname: filled with the FILTER_FILE used to build this library
   * @param emMod: [min,max] model index range with emission lines, filled
   * from the library's EM_LINES keyword
   * @param babs: filled with the reference-band index for absolute
   * magnitudes (MAG_REF keyword)
   */
  void read_lib(vector<SED*>& libFull, int& ind, int nummodpre[3],
                const string libName, string& filtname, vector<int> emMod,
                int& babs);

  /*! Verify that the cosmology and redshift grid of the library being read
   * are consistent with any library already read (all GAL/QSO libraries
   * used together must share the same cosmology and redshift grid), and
   * that the Milky Way Galametz option matches between the library and the
   * current run
   * @param keys: keywords read from the library's .doc file
   */
  void check_consistency(keymap& keys);

  /*! Parse one catalogue line into a source's identifier, fluxes/magnitudes
   * and errors, and (for the LONG format) its context and spectroscopic
   * redshift
   * @param oneObj: the source to fill in place
   * @param line: one line of the input catalogue (CAT_IN)
   */
  void readsource(onesource* oneObj, const string line);
  // LCOV_EXCL_START
  /*! Build, read and prepare a single source from one catalogue line
   * @param nobj: position (row index) to give the new source
   * @param line: one line of the input catalogue (CAT_IN)
   * @return the newly allocated, ready-to-fit source (ownership passed to
   * the caller)
   */
  onesource* yield(const int nobj, const string line) {
    onesource* oneObj = new onesource(nobj, gridz);
    readsource(oneObj, line);
    prep_data(oneObj);
    return oneObj;
  };
  // LCOV_EXCL_STOP

  /*! Read the input catalogue (CAT_IN) and keep only the sources eligible
   * for the auto-adaptation of zero-points: those with a spec-z within
   * [adzmin, adzmax] and a magnitude in the #fl_auto band within
   * [auto_thresmin, auto_thresmax]. Also applies read_mw_ebv() and
   * read_externalz() to the selected sources.
   * @return the selected sources (heap-allocated; ownership passed to the
   * caller)
   */
  vector<onesource*> read_autoadapt_sources();
  /*! Read every source of the input catalogue (CAT_IN), regardless of
   * spec-z or magnitude, applying read_mw_ebv() and read_externalz()
   * @return all sources (heap-allocated; ownership passed to the caller)
   */
  vector<onesource*> read_photoz_sources();
  /*! Set each source's Milky Way E(B-V) (onesource::mw_ebv), either from a
   * single global value (MW_GLOBAL_EBV, #one_mw_ebv true) or by matching
   * each source's Id against a per-source file (MW_EBV_FILE)
   * @param sources: sources to set the E(B-V) of, in place
   */
  void read_mw_ebv(vector<onesource*> sources);
  /*! Override each source's spectroscopic redshift (onesource::zs) by
   * matching its Id against an external file (EXTERNALZ_FILE); sources not
   * found in the file keep their catalogue redshift. No-op if
   * EXTERNALZ_FILE is "NONE".
   * @param sources: sources to override the redshift of, in place
   */
  void read_externalz(vector<onesource*> sources);
  /*! Prepare every source in the list for fitting (see the single-source
   * overload)
   * @param sources: sources to prepare, in place
   */
  void prep_data(vector<onesource*> sources);
  /*! Prepare a source for fitting: convert magnitudes to fluxes if needed,
   * rescale the flux errors, derive magnitudes, keep a copy of the
   * original values, and flag which bands are used based on the context
   * @param oneObj: the source to prepare, in place
   */
  void prep_data(onesource* oneObj);

  //! Return the indexes over zlib vector on which to run the fit
  /*!
    \param redshift: the selected redshift
    \param ir whether the selection is done on the main template library or on
    the IR one.

    \return Vector of indexes for templates set with `redshift` as redshift.
  */
  vector<size_t> validLib(const double& redshift, const bool& ir = false);

  /*! Check that the source belon to the auto-adapt sample
   * i.e. selected magnitude and redshift range
   *
   * \return a boolean
   */
  bool belong_autoadapt(onesource* src);
};

keymap read_keymap_from_doc(const string libName);

vector<string> readOutKeywords(const string outpara);

void auto_adapt(const vector<onesource*> adaptSources, vector<double>& a0,
                int& converge, int& iteration);

vector<vector<int>> bestFilter(int nbFlt, vector<double> gridz,
                               vector<SED*> fullLib, int method,
                               vector<long> magabscont, vector<int> bapp,
                               vector<int> bappOp, vector<double> zbmin,
                               vector<double> zbmax);

vector<vector<double>> maxkcolor(vector<double> gridz, vector<SED*> fullLib,
                                 vector<vector<int>> bestFlt);

void minimizekcolor(vector<double> gridz, vector<SED*> fulllib,
                    vector<vector<int>>& bestFlt, vector<long> magabscont);

#endif
