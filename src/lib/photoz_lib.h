#ifndef PHOTOZ_LIB_H
#define PHOTOZ_LIB_H

#include <ctime>   // date
#include <string>  // use string
#include <tuple>   // use std::make_tuple
#include <vector>  // manipulate vector

#include "SED.h"  // to read the libraries
#include "cosmology.h"
#include "flt.h"  // to read the libraries
#include "mag.h"

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
  bool outchi, zintp, zfix;
  string cat, methz, typm, catmag, cattyp, zmulti, outf, outsp, outpdz, outpdm;
  vector<double> shifts0, shifts1, min_err, fac_err, int_pdz, zbmin, zbmax;
  vector<flt> allFiltersAdd;
  vector<vector<int>> goodFlt;
  vector<vector<double>> maxkcol;
  vector<long> magabscont;
  vector<string> colib;
  vector<int> bapp, bappOp, pdz_fabs, emMod;
  cosmo lcdm;
  vector<size_t> valid;

 public:
  vector<vector<double>> flux, fluxIR;
  vector<SED *> fullLib, fullLibIR, lightLib;
  vector<flt> allFilters;
  vector<double> gridz;
  vector<string> outkeywords, pdftype;
  int imagm;
  time_t ti1;
  string outputHeader, outpara;

  PhotoZ(keymap &key_analysed);

  virtual ~PhotoZ() {
    fullLib.clear();
    fullLibIR.clear();
  }

  std::tuple<vector<double>, vector<double>> run_autoadapt(vector<onesource *>);

  void run_photoz(vector<onesource *> sources, const vector<double> &a0,
                  const vector<double> &a1);

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
};

keymap read_keymap_from_doc(const string libName);

vector<string> readOutKeywords(const string outpara);

void auto_adapt(const vector<onesource *> adaptSources, vector<double> &a0,
                vector<double> &a1, int &converge, int &iteration);

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
