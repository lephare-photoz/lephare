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
#include "opa.h"

using namespace std;

class SED;

static vector<string> phys_par_names = {"AGE",  "LDUST", "LIR",  "MASS", "SFR",
                                        "SSFR", "COL1",  "COL2", "MREF"};
static vector<string> photoz_par_names = {"MIN_ZG", "MIN_ZQ", "BAY_ZG",
                                          "BAY_ZQ"};
// How to associate a name to the corresponding interger
static unordered_map<string, int> maptype = {
    {"MASS", 0},    {"SFR", 1},     {"SSFR", 2},   {"LDUST", 3}, {"LIR", 4},
    {"AGE", 5},     {"COL1", 6},    {"COL2", 7},   {"MREF", 8},  {"MIN_ZG", 9},
    {"MIN_ZQ", 10}, {"BAY_ZG", 11}, {"BAY_ZQ", 12}};

class onesource {
 private:
 public:
  //{"MASS_BEST","SFR_BEST","SSFR_BEST","LDUST_BEST","LUM_TIR_BEST","AGE_BEST","EBV_BEST","EXTLAW_BEST","LUM_NUV_BEST","LUM_R_BEST","LUM_K_BEST",}
  unordered_map<string, double> results = {
      {"MASS_BEST", INVALID_PHYS},    {"SFR_BEST", INVALID_PHYS},
      {"SSFR_BEST", INVALID_PHYS},    {"LDUST_BEST", INVALID_PHYS},
      {"LUM_TIR_BEST", INVALID_PHYS}, {"AGE_BEST", INVALID_PHYS},
      {"EBV_BEST", INVALID_PHYS},     {"EXTLAW_BEST", INVALID_PHYS},
      {"LUM_NUV_BEST", INVALID_PHYS}, {"LUM_R_BEST", INVALID_PHYS},
      {"LUM_K_BEST", INVALID_PHYS},
  };
  unordered_map<string, double> results_emission_lines;

  long cont, new_cont;
  vector<double> ab, sab, mab, msab, magm, magm0, absmagPred, magPred, kap,
      mabs, emabs, ab_ori, mab_ori, abIR, sabIR;
  vector<int> busnorma, busul, busfir, bscfir, absfilt;
  string spec, str_inp;
  int pos, nbused, nbul, nbusIR, indminSec, indminIR;
  double zs, dm, consiz, closest_red;
  array<double, 3> zmin, chimin, dmmin;
  array<int, 3> indmin, imasmin;
  double zminIR, chiminIR, dmminIR, imasminIR;
  array<double, 4>
      priorLib;  // Prior with the range in abs mag gal, abs mag AGN

  vector<double> chibay;
  vector<double> gridzg, gridLdustIR, gridEbv, gridLIR;
  PDF PDFebv;

  vector<double> zgmed, zgmin, zgmode, zqmed, zqmin, zqmode;
  vector<double> massmed, SFRmed, sSFRmed, agemed, Ldustmed, LIRmed, col1med,
      col2med, ebvmed, Mrefmed;

  array<double, 65> fluxEL_SED = {0};
  double limits_zmax = 20.;
  double limits_Mfaint = 0;
  unordered_map<int, PDF> pdfmap;
  double zsec, zsecChi2, zsecEbv, zsecScale, zsecProb, zsecAge;
  int zsecMod, zsecExtlaw;

  // Minimal constructor of the source
  onesource() {
    spec = "1";  // ident
    zs = -99.9;  // spectroscopic redshift
    cont = 0;    // context
    closest_red = 0.;
    str_inp = ' ';
    for (int k = 0; k < 3; k++) {
      zmin[k] = -99.9;
      indmin[k] = -99;
      chimin[k] = 1.e9;
      imasmin[k] = -99;
    }
    zminIR = -99.9;
    indminIR = -99;
    chiminIR = 1.e9;
    imasminIR = -99;
    nbused = 0;
    pos = 0;
  }

  // Need to initialize the PDF in the constructor after the ":"
  onesource(const int posC, const vector<double> &gridz) : onesource() {
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
    mab.clear();
    msab.clear();
    kap.clear();
    mabs.clear();
    absfilt.clear();
    busnorma.clear();
  }

  // Prototype
  void readsource(const string &identifier, const vector<double> vals,
                  const vector<double> err_vals, const long context,
                  const double z_spec, const string additional_input);
  void considered_red(const bool zfix, const string methz);
  void setPriors(const array<double, 2> magabsB,
                 const array<double, 2> magabsF);
  void fltUsed(const long gbcont, const long contforb, const int imagm);
  void fltUsedIR(const long fir_cont, const long fir_scale, const int imagm,
                 vector<flt> allFilters, const double fir_lmin);
  void convertFlux(const string &catmag, const vector<flt> allFilters);
  void rescale_flux_errors(const vector<double> min_err,
                           const vector<double> fac_err);
  vector<size_t> validLib(vector<SED *> &fulllib, const bool &zfix,
                          const double &consideredZ);
  void fit(vector<SED *> &fulllib, const vector<vector<double>> &flux,
           const vector<size_t> &valid, const double &funz0,
           const array<int, 2> &bp);
  void fitIR(vector<SED *> &fulllib, const vector<vector<double>> &flux,
             const int imagm, const string fit_frsc, cosmo lcdm);
  double nzprior(const double luv, const double lnir, const double reds,
                 const array<int, 2> bp);
  void rm_discrepant(vector<SED *> &fulllib, const vector<vector<double>> &flux,
                     const vector<size_t> &valid, const double funz0,
                     const array<int, 2> bp, double thresholdChi2,
                     bool verbose);
  void write_out(vector<SED *> &fulllib, vector<SED *> &fulllibIR,
                 ofstream &stout, vector<string> outkeywords);
  void write_pdz_header(vector<string> pdztype,
                        unordered_map<string, ofstream> &stpdz,
                        const time_t &ti1);
  void write_pdz(vector<string> pdztype,
                 unordered_map<string, ofstream> &stpdz);
  void convertMag();
  void keepOri();
  void interp(const bool zfix, const bool zintp, cosmo lcdm);
  void uncertaintiesMin();
  void uncertaintiesBay();
  void secondpeak(vector<SED *> &fulllib, const double dz_win,
                  const double min_thres);
  void generatePDF(vector<SED *> &fulllib, const vector<size_t> &va,
                   const vector<int> fltColRF, int fltREF, const bool zfix);
  void generatePDF_IR(vector<SED *> &fulllib);
  void mode();
  void interp_lib(vector<SED *> &fulllib, const int imagm, cosmo lcdm);
  void adapt_mag(vector<double> a0, vector<double> a1);
  void substellar(const bool substar, vector<flt> allFilters);
  void absmag(const vector<vector<int>> &bestFlt,
              const vector<vector<double>> &maxkcolor, cosmo lcdm,
              const vector<double> gridz);
  void writeSpec(vector<SED *> &fulllib, vector<SED *> &fulllibIR, cosmo lcdm,
                 vector<opa> opaAll, const vector<flt> &allFilters,
                 const string outspdir);
  void writeFullChi(vector<SED *> &fulllib);
  void computePredMag(vector<SED *> &fulllib, cosmo lcdm, vector<opa> opaAll,
                      vector<flt> allFltAdd);
  void computePredAbsMag(vector<SED *> &fulllib, cosmo lcdm, vector<opa> opaAll,
                         vector<flt> allFltAdd);
  void computeEmFlux(vector<SED *> &fulllib, cosmo lcdm, vector<opa> opaAll);
  void limits(vector<SED *> &fulllib, vector<double> &limits_zbin,
              int limits_ref, vector<int> &limits_sel,
              vector<double> &limits_cut);
  pair<vector<double>, vector<double>> best_spec_vec(short sol,
                                                     vector<SED *> &fulllib,
                                                     cosmo lcdm,
                                                     vector<opa> opaAll,
                                                     double minl, double maxl);

  void compute_best_fit_physical_quantities(vector<SED *> &fulllib);
};

#endif
