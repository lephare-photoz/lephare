#ifndef PHOTOZ_LIB_H
#define PHOTOZ_LIB_H

#include <ctime>   // date
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
  array<int, 2> bp;
  long gbcont, contforb, bdscal;
  double funz0, adzmin, adzmax, auto_thresmin, auto_thresmax, min_thres, dz_win;
  double fir_lmin, fir_cont, fir_scale;
  bool verbose;
  array<double, 2> magabsB, magabsF, zrange, ebvrange;
  bool outchi,  /// Whether or not to output th chi2 values in ascii
      zintp,    /// Whether or not to interpolate the z solution from the grid
                /// result
      zfix,     /// Whether or not to run with a fixed redshift, typically the
                /// true/spectro z
      methz;    /// If true, set z to the MEDIAN solution rather than the BEST
                /// solution, when computing physical parameters.

  string cat, typm, catmag, cattyp, zmulti, outf, outsp, outpdz, outpdm;
  vector<double> shifts0, min_err, fac_err, int_pdz, zbmin, zbmax;
  vector<flt> allFiltersAdd;
  vector<vector<int>> goodFlt;
  vector<vector<double>> maxkcol;
  vector<long> magabscont;
  vector<string> colib;
  vector<int> bapp, bappOp, pdz_fabs, emMod;
  cosmo lcdm;

 public:
  vector<vector<double>> flux, fluxIR;
  vector<double> zLib, zLibIR;
  vector<SED *> fullLib, fullLibIR, lightLib;
  vector<flt> allFilters;
  vector<double> gridz;
  vector<string> outkeywords, pdftype;
  int imagm;
  time_t ti1;
  string outputHeader, outpara;

  PhotoZ(keymap &key_analysed);

  virtual ~PhotoZ() {
    for (auto &sed : fullLib) delete sed;
    for (auto &sed : fullLibIR) delete sed;
    for (auto &sed : lightLib) delete sed;
    fullLib.clear();
    fullLibIR.clear();
    lightLib.clear();
  }

  vector<double> compute_offsets(vector<onesource *>);
  vector<double> run_autoadapt(vector<onesource *>);

  void run_photoz(vector<onesource *> sources, const vector<double> &a0);

  string prep_header(vector<string> outkeywords);

  void write_outputs(vector<onesource *> sources, const time_t &ti1);

  void read_lib(vector<SED *> &fullLib, int &ind, int nummodpre[3],
                const string libName, string &filtname, vector<int> emMod,
                int &babs);

  void check_consistency(keymap &keys);

  void readsource(onesource *oneObj, const string line);
  onesource *yield(const int nobj, const string line) {
    onesource *oneObj = new onesource(nobj, gridz);
    readsource(oneObj, line);
    prep_data(oneObj);
    return oneObj;
  };

  vector<onesource *> read_autoadapt_sources();
  vector<onesource *> read_photoz_sources();
  void prep_data(vector<onesource *> sources);
  void prep_data(onesource *oneObj);

  //! Return the indexes over zlib vector on which to run the fit
  /*!
    \param redshift: the selected redshift
    \param ir whether the selection is done on the main template library or on
    the IR one.

    \return Vector of indexes for templates set with `redshift` as redshift.
  */
  vector<size_t> validLib(const double &redshift, const bool &ir = false);
};

keymap read_keymap_from_doc(const string libName);

vector<string> readOutKeywords(const string outpara);

void auto_adapt(const vector<onesource *> adaptSources, vector<double> &a0,
                int &converge, int &iteration);

vector<vector<int>> bestFilter(int nbFlt, vector<double> gridz,
                               vector<SED *> fullLib, int method,
                               vector<long> magabscont, vector<int> bapp,
                               vector<int> bappOp, vector<double> zbmin,
                               vector<double> zbmax);

vector<vector<double>> maxkcolor(vector<double> gridz, vector<SED *> fullLib,
                                 vector<vector<int>> bestFlt);

void minimizekcolor(vector<double> gridz, vector<SED *> fulllib,
                    vector<vector<int>> &bestFlt, vector<long> magabscont);

#endif
