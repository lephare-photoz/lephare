#include <math.h>    /* pow */
#include <unistd.h>  //posix interface, prooviding path access method

#include <algorithm>
#include <array>
#include <cmath>
#include <cstring>
#include <ctime>  // date
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

// Le Phare
#include <omp.h>

#include "SED.h"        //our own class to read the keywords
#include "cosmology.h"  // in order to measure the distance modulus
#include "globals.h"    // global variables
#include "keyword.h"    //our own class to read the keywords
#include "mag.h"  // to create the predicted magnitudes/k-corrections along the grid
#include "onesource.h"
#include "photoz_lib.h"

using namespace std;

/*
  Construct the photo-z run.
  1 - Analyse the keywords
  2 - Read the magnitude's library and other documentation
*/
PhotoZ::PhotoZ(keymap &key_analysed) {
  keys = key_analysed;
  imagm = 0;

  // Configuration file
  string config = key_analysed["c"].value;

  outputHeader += "####################################### \n";
  outputHeader += "# PHOTOMETRIC REDSHIFT with OPTIONS   # \n";
  outputHeader += "# Config file            : " + config + '\n';

  /* INPUT CATALOG OPTIONS  */
  // CAT_IN input catalogue
  cat = key_analysed["CAT_IN"].value;

  // INP_TYPE input type, flux or magnitude
  typm = key_analysed["INP_TYPE"].value;

  // CAT_MAG mag type AB/VEGA - AB default
  catmag = ((key_analysed["CAT_MAG"]).split_string("AB", 1))[0];

  // CAT_FMT format of the input file MEME/MMEE - MEME default
  string meme = ((key_analysed["CAT_FMT"]).split_string("MEME", 1))[0];
  cat_fmt = 0;
  if (meme[1] == 'M') cat_fmt = 1;

  // CAT_LINE define the line range to be considered in the input catalogue
  // NOTE: commented lines are NOT considered while reading the catalogue,
  //       so this range should be intended as the number of entries, not rows
  int rowmin_tmp = ((key_analysed["CAT_LINES"]).split_int("0", 2))[0];
  int rowmax_tmp = ((key_analysed["CAT_LINES"]).split_int("20000000000", 2))[1];
  rowmin = rowmin_tmp < 0 ? 0 : rowmin_tmp;
  rowmax = rowmax_tmp < 0 ? 2000000000 : rowmax_tmp;

  // ZPHOTLIB values, multiple librairies are possible,
  // number of expected values unknown in advance -> -1
  colib = (key_analysed["ZPHOTLIB"]).split_string("GAL_LIB", -1);
  int numlib = int(colib.size());

  // PARA_OUT output parameter file - output.para default
  outpara = ((key_analysed["PARA_OUT"]).split_string("output.para", 1))[0];

  // CAT_OUT output  file -  zphot.out default
  outf = ((key_analysed["CAT_OUT"]).split_string("zphot.out", 1))[0];

  // Want to display the template number on the screen
  // VERBOSE output  file -  YES default
  verbose = key_analysed["VERBOSE"].split_bool("YES", 1)[0];

  /* SECONDARY OPTIONS */

  // CAT_TYPE Type of catalogue (short: read only Id mag err; long: add context
  // zs ...) - SHORT default
  cattyp = ((key_analysed["CAT_TYPE"]).split_string("SHORT", 1))[0];

  // ERR_SCALE Minimal uncertainties to be added in quadrature - 0.0 default
  min_err = (key_analysed["ERR_SCALE"]).split_double("0.0", -1);
  int nerr = int(min_err.size());

  // ERR_FACTOR Multiply the flux uncertainties by a given factor - 1.0 default
  fac_err = ((key_analysed["ERR_FACTOR"]).split_double("1.0", -1));
  int nfac = int(fac_err.size());
  if ((nerr > 1) && (nfac > 1) && (nfac != nerr)) {
    cout << "The number of filters in ERR_SCALE and ERR_FACTOR do not "
            "correspond."
         << endl;
    cout << " Error in the keywords. " << endl;
    cout << " Stop " << endl;
  }

  // GLB_CONTEXT Global context to be used for all objects - 0 default (all
  // bands)
  gbcont = ((key_analysed["GLB_CONTEXT"]).split_long("0", 1))[0];

  // FORB_CONTEXT Context to reject some bands for all sources - 0 default
  contforb = ((key_analysed["FORB_CONTEXT"]).split_long("0", 1))[0];

  /*  PRIOR  */

  // Limits the size of the library in redshift and E(B-V)
  zrange[0] = ((keys["Z_RANGE"]).split_double("0", 2))[0];
  zrange[1] = ((keys["Z_RANGE"]).split_double("10000", 2))[1];
  if (zrange[0] < 0 || zrange[1] < 0) {
    zrange[0] = 0;
    zrange[1] = 10000;
  }
  ebvrange[0] = ((keys["EBV_RANGE"]).split_double("0", 2))[0];
  ebvrange[1] = ((keys["EBV_RANGE"]).split_double("10000", 2))[1];
  if (ebvrange[0] < 0 || ebvrange[1] < 0) {
    ebvrange[0] = 0;
    ebvrange[1] = 10000;
  }

  // MAG_ABS allowed range in absolute magnitude for galaxies - 0 default
  magabsB[0] = ((key_analysed["MAG_ABS"]).split_double("0.", 2))[0];
  magabsF[0] = ((key_analysed["MAG_ABS"]).split_double("0.", 2))[1];
  // be sure that magabsB is the minimum value
  if (magabsB[0] > magabsF[0]) {
    double inter = magabsF[0];
    magabsF[0] = magabsB[0];
    magabsB[0] = inter;
  }

  // MAG_ABS_QSO allowed range in absolute magnitude for QSO library - 0 default
  // (all bands)
  magabsB[1] = ((key_analysed["MAG_ABS_QSO"]).split_double("0.", 2))[0];
  magabsF[1] = ((key_analysed["MAG_ABS_QSO"]).split_double("0.", 2))[1];
  // be sure that magabsB is the minimum value
  if (magabsB[1] > magabsF[1]) {
    double inter = magabsF[1];
    magabsF[1] = magabsB[1];
    magabsB[1] = inter;
  }

  // MAG_REF allowed range in absolute magnitude: define the band - 1 default
  babs = ((key_analysed["MAG_REF"]).split_int("1", 1))[0];
  // Shift of 1 because of the convention (1 to start in the para, 0 in the
  // array)
  babs = babs - 1;
  if (babs < 0) {
    magabsB[0] = 99.;
    magabsF[0] = 99.;
    magabsB[1] = 99.;
    magabsF[1] = 99.;
  }

  // NZ_PRIOR prior on N(z) based on z-VVDS: I mag (the second number is in the
  // case I band is not defined) - -1 by default number of expected values could
  // be one or two
  vector<int> test;
  test = (key_analysed["NZ_PRIOR"]).split_int("-1", -1);
  int testNb = int(test.size());
  if (testNb == 2) {
    bp[0] = ((key_analysed["NZ_PRIOR"]).split_int("-1", 2))[0];
    bp[1] = ((key_analysed["NZ_PRIOR"]).split_int("-1", 2))[1];
    // Shift of 1 because of the convention (1 to start in the para, 0 in the
    // array)
    bp[0] = bp[0] - 1;
    bp[1] = bp[1] - 1;
  } else {
    bp[0] = ((key_analysed["NZ_PRIOR"]).split_int("-1", 1))[0];
    bp[0] = bp[0] - 1;
    bp[1] = bp[0];
  }
  /* Fix redshift */
  // ZFIX  fixed the redshift based on the inoput catalogue - NO default
  zfix = key_analysed["ZFIX"].split_bool("NO", 1)[0];

  // Z_INTERP redshift interpolation between library Z-STEP
  zintp = key_analysed["Z_INTERP"].split_bool("NO", 1)[0];

  // Z_METHOD compute the absolute magnitude at a given redshift solution ML or
  // BEST - BEST by default
  methz = ((key_analysed["Z_METHOD"]).split_string("BEST", 1))[0];

  /* search for secondary solution  */

  // DZ_WIN minimal delta z window to search - 0.25 by default
  dz_win = ((key_analysed["DZ_WIN"]).split_double("0.25", 1))[0];

  // MIN_THRES threshold to trigger the detection - 0.1 by default
  min_thres = ((key_analysed["MIN_THRES"]).split_double("0.1", 1))[0];

  // PROB_INTZ Integrated PDFz over Z ranges - 0. by default
  int_pdz = (key_analysed["PROB_INTZ"]).split_double("0.", -1);
  int npdz = int(int_pdz.size());

  /* Output */

  // SPEC_OUT Output individual spectra - NO default
  outsp = ((key_analysed["SPEC_OUT"]).split_string("NO", 1))[0];
  // CHI2_OUT Output the full chi2 library - NO default
  outchi = ((key_analysed["CHI2_OUT"]).split_bool("NO", 1))[0];

  /* FIR Libraries */

  // FIR_LIB name of the library in IR
  vector<string> libext = (key_analysed["FIR_LIB"]).split_string("NONE", -1);
  int nlibext = int(libext.size());
  // if the library is NONE, put the number of library at -1
  if (nlibext == 1 && libext[0] == "NONE") nlibext = -1;

  // FIR_LMIN Lambda min given in micron  (um)
  fir_lmin = ((key_analysed["FIR_LMIN"]).split_double("7.0", 1))[0];

  // FIR_CONT context for the far-IR
  fir_cont = ((key_analysed["FIR_CONT"]).split_long("-1", 1))[0];

  // FIR_SCALE context of the bands used for the rescaling in IR
  fir_scale = ((key_analysed["FIR_SCALE"]).split_int("-1", 1))[0];

  // FIR_FREESCALE possible free rscaling in IR, when several bands. Otherwise,
  // model imposed by its LIR
  string fir_frsc = (key_analysed["FIR_FREESCALE"]).split_string("NO", 1)[0];

  // FIR_SUBSTELLAR remove the stellar component
  bool substar = key_analysed["FIR_SUBSTELLAR"].split_bool("NO", 1)[0];

  // MABS_METHOD method to compute the absolute magnitudes
  method = ((key_analysed["MABS_METHOD"]).split_int("0", -1))[0];

  // MABS_CONTEXT context for method 1
  magabscont = (key_analysed["MABS_CONTEXT"]).split_long("0", -1);
  int nmagabscont = magabscont.size();

  // MABS_REF reference filter in case of mag abs method 2
  bapp = (key_analysed["MABS_REF"]).split_int("1", -1);
  int nbapp = int(bapp.size());
  // Need to substract one because the convention in the .para file start at 1,
  // but 0 in the code
  for (int k = 0; k < nbapp; k++) bapp[k]--;

  // MABS_ZBIN give the redshift bins corresponding to MABS_FILT
  // MABS_FILT choose filters per redshift bin (MABS_ZBIN) if method 4
  bappOp = (key_analysed["MABS_FILT"]).split_int("1", -1);
  int nbBinZ = int(bappOp.size());
  zbmin = (key_analysed["MABS_ZBIN"]).split_double("0", nbBinZ + 1);
  zbmin.erase(zbmin.end() - 1);
  zbmax = (key_analysed["MABS_ZBIN"]).split_double("6", nbBinZ + 1);
  zbmax.erase(zbmax.begin(), zbmax.begin() + 1);

  /*  AUTO-ADAPT */
  // APPLY_SYSSHIFT shift to be applied in each filter. Create a fake shifts1
  // vector
  shifts0 = (key_analysed["APPLY_SYSSHIFT"]).split_double("0.", -1);
  for (size_t i = 0; i < shifts0.size(); i++) shifts1.push_back(0.);

  /* PDZ OUTPUT */
  // PDZ_OUT pdz output file
  outpdz = (key_analysed["PDZ_OUT"]).split_string(nonestring, 1)[0];
  // PDZ_TYPE type of PZD (BAY_ZG,BAY_ZQ,MIN_ZG,MIN_ZQ,MASS,SFR,SSFR,AGE)
  pdftype = ((key_analysed["PDZ_TYPE"]).split_string("BAY_ZG", -1));
  // PDZ_MABS_FILT filter in which we want the mag abs in each pdz
  pdz_fabs = (key_analysed["PDZ_MABS_FILT"]).split_int("0", -1);
  // // PDM_OUT pdf(stellar mass) output file
  // outpdm = ((key_analysed["PDM_OUT"]).split_string(nonestring,1))[0];

  // ADD_EMLINES
  // minimum and maximum models to add emission lines
  emMod.push_back(((key_analysed["ADD_EMLINES"]).split_int("-9999", 2))[0]);
  emMod.push_back(((key_analysed["ADD_EMLINES"]).split_int("-9999", 2))[1]);

  /*
    INFO PARAMETERS ON SCREEN AND DOC
  */

  outputHeader += "# CAT_IN                 : " + cat + '\n';
  outputHeader += "# CAT_OUT                : " + outf + '\n';
  outputHeader += "# CAT_LINES              : " + to_string(rowmin) + ' ' +
                  to_string(rowmax) + '\n';
  outputHeader += "# PARA_OUT               : " + outpara + '\n';
  outputHeader += "# INP_TYPE               : " + typm + '\n';
  outputHeader += "# CAT_FMT[0:MEME 1:MMEE] : " + to_string(cat_fmt) + '\n';
  outputHeader += "# CAT_MAG                : " + catmag + '\n';
  outputHeader += "# ZPHOTLIB               : ";
  for (int k = 0; k < numlib; k++) {
    outputHeader += colib[k] + ' ';
  };
  outputHeader += '\n';
  outputHeader += "# FIR_LIB                : ";
  for (int k = 0; k < nlibext; k++) {
    outputHeader += libext[k] + ' ';
  };
  outputHeader += '\n';
  outputHeader += "# FIR_LMIN               : " + to_string(fir_lmin) + '\n';
  outputHeader += "# FIR_CONT               : " + to_string(fir_cont) + '\n';
  outputHeader += "# FIR_SCALE              : " + to_string(fir_scale) + '\n';
  outputHeader += "# FIR_FREESCALE          : " + fir_frsc + '\n';
  outputHeader += "# FIR_SUBSTELLAR         : " + bool2string(substar) + '\n';
  outputHeader += "# ERR_SCALE              : ";
  for (int k = 0; k < nerr; k++) {
    outputHeader += to_string(min_err[k]) + ' ';
  };
  outputHeader += '\n';
  outputHeader += "# ERR_FACTOR             : ";
  for (int k = 0; k < nfac; k++) {
    outputHeader += to_string(fac_err[k]) + ' ';
  };
  outputHeader += '\n';
  outputHeader += "# GLB_CONTEXT            : " + to_string(gbcont) + '\n';
  outputHeader += "# FORB_CONTEXT           : " + to_string(contforb) + '\n';
  outputHeader += "# DZ_WIN                 : " + to_string(dz_win) + '\n';
  outputHeader += "# MIN_THRES              : " + to_string(min_thres) + '\n';
  outputHeader += "# MAG_ABS                : " + to_string(magabsB[0]) + ' ' +
                  to_string(magabsF[0]) + '\n';
  outputHeader += "# MAG_ABS_AGN            : " + to_string(magabsB[1]) + ' ' +
                  to_string(magabsF[1]) + '\n';
  outputHeader += "# MAG_REF                : " + to_string(babs) + '\n';
  outputHeader += "# NZ_PRIOR               : " + to_string(bp[0]) + ' ' +
                  to_string(bp[1]) + '\n';
  outputHeader += "# Z_INTERP               : " + bool2string(zintp) + '\n';
  outputHeader += "# Z_METHOD               : " + methz + '\n';
  outputHeader += "# PROB_INTZ              : ";
  for (int k = 0; k < npdz; k++) {
    outputHeader += to_string(int_pdz[k]) + ' ';
  };
  outputHeader += '\n';

  outputHeader += "# MABS_METHOD            : " + to_string(method) + '\n';
  outputHeader += "# MABS_CONTEXT           : ";
  for (int k = 0; k < nmagabscont; k++) {
    outputHeader += to_string(magabscont[k]) + ' ';
  };
  outputHeader += '\n';
  outputHeader += "# MABS_REF               : ";
  for (int k = 0; k < nbapp; k++) {
    outputHeader += to_string(bapp[k]) + ' ';
  };
  outputHeader += '\n';
  // AUTO-ADAPT
  outputHeader += "# AUTO_ADAPT             : " +
                  bool2string(keys["AUTO_ADAPT"].split_bool("NO", 1)[0]) + '\n';

  // ADAPT_BAND selection in one band
  fl_auto = ((key_analysed["ADAPT_BAND"]).split_int("1", 1))[0];
  // Need to substract one because the convention in the .para file start at 1,
  // but 0 in the code
  fl_auto--;
  outputHeader += "# ADAPT_BAND             : " + to_string(fl_auto) + '\n';
  // ADAPT_LIM limit for the selection in auto-adapt
  auto_thresmin = ((key_analysed["ADAPT_LIM"]).split_double("15", 2))[0];
  auto_thresmax = ((key_analysed["ADAPT_LIM"]).split_double("35", 2))[1];
  outputHeader += "# ADAPT_LIM              : " + to_string(auto_thresmin) +
                  ' ' + to_string(auto_thresmax) + '\n';
  // ADAPT_ZBIN minimum and maximum redshift for the adaptation
  adzmin = ((key_analysed["ADAPT_ZBIN"]).split_double("0.001", 2))[0];
  adzmax = ((key_analysed["ADAPT_ZBIN"]).split_double("6", 2))[1];
  outputHeader += "# ADAPT_ZBIN             : " + to_string(adzmin) + ' ' +
                  to_string(adzmax) + '\n';

  outputHeader += "# ZFIX                   : " + bool2string(zfix) + '\n';
  outputHeader += "# SPEC_OUT               : " + outsp + '\n';
  outputHeader += "# CHI_OUT                : " + bool2string(outchi) + '\n';
  outputHeader += "# PDZ_OUT                : " + outpdz + '\n';
  outputHeader += "####################################### \n";

  // Write the header on the screen
  cout << outputHeader;

  /*
    PREPARE THE WORK
  */

  /* Reading input libraries */
  // typlib : 1 for GAL  / 2 for QSO / 3 for STAR
  cout << "Reading input librairies ..." << endl;
  string filtName, filtNameFIR, line;
  cout << "Read lib " << endl;
  // vector containing all SEDs
  // loop over each type of SED (star/gal/AGN)
  int ind = 0, nummodpre[3] = {0, 0, 0};

  // Open the galaxies last to keep the grid z of their library
  for (int j = numlib - 1; j >= 0; j--) {
    read_lib(fullLib, ind, nummodpre, colib[j], filtName, emMod, babs);
  }
  cout << "Read lib out " << endl;

  // reading FIR SED libraries
  // loop over each type of SED (star/gal/AGN)
  int indFIR = 0, nummodpreFIR[3];

  // Infrared library
  for (int j = nlibext - 1; j >= 0; j--) {
    read_lib(fullLibIR, indFIR, nummodpreFIR, libext[j], filtNameFIR, emMod,
             babs);
  }
  cout << "Read filt " << endl;

  // we need to define a zgrid singleton in case zphota ends up
  // running with only STAR templates.
  if (gridz.empty()) gridz = {0.};

  /* Reading filters */
  allFilters = read_doc_filters(filtName);
  imagm = allFilters.size();

  /* Reading additional filters in which to compute the magnitude*/
  string filtNameAdd =
      (key_analysed["ADDITIONAL_MAG"]).split_string("none", 1)[0];
  if (filtNameAdd != "none" && filtNameAdd != "NONE") {
    allFiltersAdd = read_doc_filters(filtNameAdd);
  }

  /* Create a 2D array with the predicted flux.
  Done to improve the performance in the fit*/
  flux.resize(fullLib.size(), vector<double>(imagm, 0.));
  fluxIR.resize(fullLibIR.size(), vector<double>(imagm, 0.));
// Convert the magnitude library in flux
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  // Initialize the chi2
  for (size_t i = 0; i < fullLib.size(); i++) {
    // Loop over the filters
    for (size_t k = 0; k < allFilters.size(); k++) {
      flux[i][k] = pow(10., -0.4 * (fullLib[i]->mag[k] + 48.6));
    }
  }
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (size_t i = 0; i < fullLibIR.size(); i++) {
    for (size_t k = 0; k < allFilters.size(); k++) {
      fluxIR[i][k] = pow(10., -0.4 * (fullLibIR[i]->mag[k] + 48.6));
    }
  }
}

keymap read_keymap_from_doc(const string libName) {
  // List of the keywords to be found in the config file/command line
  string list_keywords[] = {
      "LIB_TYPE",   "NUMBER_ROWS", "FILTER_FILE", "FILTERS", "EM_LINES",
      "LIB_NAME",   "NUMBER_SED",  "ZGRID_TYPE",  "Z_STEP",  "COSMOLOGY",
      "EXTINC_LAW", "EB_V",        "MOD_EXTINC",  "Z_FORM"};
  // Number of keywords
  int nb_doc_key = (int)(sizeof(list_keywords) / sizeof(list_keywords[0]));
  cout << "Number of keywords to be read in the doc: " << nb_doc_key << endl;

  // Find the reference keywords within the command line + config file -> store
  // them into a vector
  string docOutFile = lepharework + "/lib_mag/" + libName + ".doc";
  // Check the exitence of the file. Stop if doesn't exist
  if (access(docOutFile.c_str(), F_OK) == -1) {
    throw invalid_argument("The doc file " + docOutFile +
                           " does not exist. Stop.");
  }
  // need to create an array of char pointer, exactly as argv, in order to use
  // analyse_keywords
  char *writable = new char[docOutFile.size() + 1];
  copy(docOutFile.begin(), docOutFile.end(), writable);
  writable[docOutFile.size()] = '\0';
  char *argv[] = {(char *)"bid", (char *)"-c", writable, NULL};
  keymap key_analysed = analyse_keywords(2, argv, list_keywords, nb_doc_key);
  delete[] writable;
  return key_analysed;
}

void PhotoZ::check_consistency(keymap &keys) {
  string valc = keys["LIB_TYPE"].split_string("GALAXY", 1)[0];
  // if we are not looking at a STAR library, we need
  // to assert that cosmology and z grids have been defined identically
  // when building the input libraries with mag_gal
  if (valc[0] != 'S' && valc[0] != 's') {
    vector<double> cosmo_keys = keys["COSMOLOGY"].split_double("70", 3);
    double h0 = cosmo_keys[0];
    double om0 = cosmo_keys[1];
    double l0 = cosmo_keys[2];
    int gridType = keys["ZGRID_TYPE"].split_int("0", 1)[0];
    vector<double> zstep_keys = keys["Z_STEP"].split_double("0.1", 3);
    double zstep = zstep_keys[0];
    double zmin = zstep_keys[1];
    double zmax = zstep_keys[2];
    if (gridz.size() == 0) {
      lcdm = cosmo(h0, om0, l0);
      gridz = zgrid(gridType, zstep, zmin, zmax);
    } else {
      cosmo lcdm2(h0, om0, l0);
      vector<double> gridz2 = zgrid(gridType, zstep, zmin, zmax);
      if (lcdm2 != lcdm || gridz2 != gridz) {
        throw runtime_error(
            "Cosmology and redshift grids parameters are not "
            "consistent in the input libraries!");
      }
    }
  }
}

/*
  Read the magnitude library
*/
void PhotoZ::read_lib(vector<SED *> &libFull, int &ind, int nummodpre[3],
                      const string libName, string &filtname, vector<int> emMod,
                      int &babs) {
  keyword oneel;
  vector<keyword> key_doc;
  string line;

  // get the keywords from the saved doc files in order to ensure library
  // consistency
  keymap key_analysed = read_keymap_from_doc(libName);
  // build cosmology and gridz, and verify that they are identical across GAL
  // and QSO libraries
  check_consistency(key_analysed);

  // LIB_TYPE type of library GAL/QSO/STAR
  string valc = ((key_analysed["LIB_TYPE"]).split_string("GALAXY", 1))[0];
  // FILTER_FILE
  filtname = ((key_analysed["FILTER_FILE"]).split_string("filters.dat", 1))[0];
  // EM_LINES
  string emlines = ((key_analysed["EM_LINES"]).split_string("NO", 1))[0];

  // EXTINC_LAW"
  vector<string> extlaw =
      (key_analysed["EXTINC_LAW"]).split_string("calzetti.dat", -1);
  int nextlaw = int(extlaw.size());
  // EB_V
  vector<double> ebv = (key_analysed["EB_V"]).split_double("0", -1);
  // MOD_EXTINC
  vector<int> modext =
      (key_analysed["MOD_EXTINC"]).split_int("0,0", nextlaw * 2);

  // Z_FORM
  vector<double> zf = (key_analysed["Z_FORM"]).split_double("10000", -1);

  /*  LIBRARY */

  // read the library in binary file
  string binLibFile = lepharework + "/lib_mag/" + libName + ".bin";
  cout << "Reading library: " << binLibFile << endl;
  ifstream slibIn(binLibFile.c_str(), ios::binary);
  if (!slibIn) {
    throw invalid_argument("Library doesn't exist ");
  }

  // find the end of the stream for the binary input file
  slibIn.seekg(0, slibIn.end);
  long length = slibIn.tellg();
  slibIn.seekg(0, slibIn.beg);
  // read each step in the library in the binary file
  // and store it as a SED
  int ind0 = 0;
  vector<double> mag_z0;
  while (slibIn.tellg() < length) {
    SED *oneSED;

    if (valc[0] == 'G' || valc[0] == 'g') {
      oneSED = new GalSED("bid.dat", 0);
      // If emission lines used, only for galaxies
      if (emlines[0] == 'P' || emlines[0] == 'E') oneSED->has_emlines = true;
    } else if (valc[0] == 'Q' || valc[0] == 'q') {
      oneSED = new QSOSED("bid.dat", 0);
    } else if (valc[0] == 'S' || valc[0] == 's') {
      oneSED = new StarSED("bid.dat", 0);
    } else {
      throw invalid_argument("There is no such SED type defined: " + valc);
    }
    // read each SED in the library binary file
    oneSED->readMagBin(slibIn);

    // Keep in the library only redshift or EBV in of the prior range (except 0
    // in all cases).
    if ((oneSED->red >= zrange[0] && oneSED->red <= zrange[1]) ||
        oneSED->red < 1.e-20) {
      if ((oneSED->ebv >= ebvrange[0] && oneSED->ebv <= ebvrange[1]) ||
          oneSED->ebv < 1.e-20) {
        // Case in which we add the emission lines (only for galaxies), if
        // option turned on
        if ((valc[0] == 'G' || valc[0] == 'g') &&
            (emlines[0] == 'P' || emlines[0] == 'E')) {
          // Check if we are well in the considered model range to add lines
          if (emMod[0] >= 0 && emMod[1] >= 0 && (oneSED->nummod) >= emMod[0] &&
              (oneSED->nummod) <= emMod[1]) {
            oneSED->sumEmLines();
          } else {
            // Loop over the lines: put them at 0 since out of the considered
            // model range
            if (oneSED->red == 0) {
              for (int k = 0; k < 65; k++) {
                oneSED->fac_line[k].val = 0.;
              }
            }
          }
        }
        // Index la SED
        (oneSED->index) = ind;
        // keep the index and mag at z=0
        if (oneSED->red == 0) {
          ind0 = oneSED->index;
          // Keep the magnitude in each band
          mag_z0.clear();
          for (int k = 0; k < int((oneSED->mag).size()); k++)
            mag_z0.push_back(oneSED->mag[k]);
        }
        // Store the index at z=0 in the SED
        oneSED->index_z0 = ind0;
        // Store the magnitude in the band used for the prior
        oneSED->mag0 = mag_z0[babs];
        // Recompute the kcorr because of the emission lines not included in
        // mag_gal
        oneSED->kcorrec(mag_z0);
        // Add the model number of the previous library, for each library
        oneSED->nummod += nummodpre[oneSED->nlib];
        // Add a new SED to the output
        libFull.push_back(oneSED);
        ind++;
      }
    }
  }
  // keep the last model number to start the numbering with it at the next
  // library
  nummodpre[libFull.back()->nlib] = (libFull.back())->nummod;
  slibIn.close();
  cout << " Done with the library reading with " << libFull.size()
       << " SED read. " << endl;
}

/*
   read the keywords of the ouput parameter file
*/
vector<string> readOutKeywords(const string outpara) {
  ifstream stoutpara;
  string lit;
  vector<string> vecKey;
  int nbLines = 0;
  bool inString = false;

  // check PARA_OUT config parameter
  if (outpara == nonestring) {
    cout << "WARNING: PARA_OUT not set or set to NONE, output will contain all "
            "the parameters resulting from the fit."
         << endl;
    return readOutKeywords(lepharedir + "/alloutputkeys.txt");
  } else {
    // open the output parameter file into a stream
    stoutpara.open(outpara.c_str());
    if (!stoutpara) {
      throw invalid_argument(
          "Incorrect setting for PARA_OUT: can't open output parameter file " +
          outpara);
    }

    // using a set to find and discard duplicates
    set<string> tmp;

    // Take the stream line by line
    while (getline(stoutpara, lit)) {
      // If the first character of the line is not #
      if (check_first_char(lit)) {
        // put the line into the stream ss again
        stringstream ss(lit);

        // fill the lambda/trans values of the SED
        ss >> fakeString;

        if (fakeString == "STRING_INPUT") {
          inString = true;
          continue;
        }
        if (fakeString == "Z_QSO") {
          fakeString = "ZQ_BEST";
        }
        //_ML is legacy, replaced by _MED
        string deprecated = "_ML";
        size_t start = fakeString.find(deprecated);
        if (start != string::npos) {
          fakeString.replace(start, deprecated.length(), "_MED");
        }
        // If insertion to the set occurs, then
        // it is not a duplicate; res.second is bool true if
        // insertion occurred
        auto res = tmp.insert(fakeString);
        if (res.second) {
          vecKey.push_back(fakeString);
          nbLines++;
        }
      }
    }
  }

  if (nbLines == 0) cout << "No keyword for the output parameter file " << endl;

  // Put "STRING_INPUT" at the end
  if (inString) vecKey.push_back("STRING_INPUT");

  // Close the stream
  stoutpara.close();

  return vecKey;
}

/*
  create the format to be included in the header
*/
string PhotoZ::prep_header(vector<string> outkeywords) {
  const int imagm = allFilters.size();
  const int imagmAdd = int(allFiltersAdd.size());
  int numCol = 1;
  string form = "# Format: \n#";

  // Loop over each keyword
  for (const auto &outkey : outkeywords) {
    for (const string k : {"IDENT", "CONTEXT", "NBAND_USED", "NBAND_ULIM",
                           "ZSPEC", "STRING_INPUT"}) {
      if (outkey == k) {
        form += " " + k + " " + to_string(numCol);
        numCol++;
      }
    }

    for (const string k : {"Z", "ZQ"}) {
      for (const string k1 : {"BEST", "MED", "MODE"}) {
        for (const string k2 : {"", "68_LOW", "68_HIGH", "90_LOW", "90_HIGH",
                                "99_LOW", "99_HIGH"}) {
          if (outkey == k + "_" + k1 + k2) {
            form += " " + outkey + " " + to_string(numCol);
            numCol++;
          }
        }
      }
    }

    for (const string k :
         {"MOD", "EXTLAW", "EBV", "MASS", "LDUST", "LUM_TIR", "LUM_NUV",
          "LUM_R", "LUM_K", "PDZ", "CHI", "SCALE", "AGE", "SFR", "SSFR"}) {
      if (outkey == k + "_BEST") {
        form += " " + outkey + " " + to_string(numCol);
        numCol++;
      }
    }

    for (const string k : {"MAG_OBS", "ERR_MAG_OBS", "MAG_MOD", "MABS_FILT",
                           "K_COR", "MAG_ABS", "EMAG_ABS"}) {
      if (outkey == k + "()") {
        form += " " + outkey + to_string(numCol);
        numCol += imagm;
        form += "-" + to_string(numCol - 1);
      }
    }

    for (const string k : {
             "MAG_PRED",
             "ABSMAG_PRED",
         }) {
      if (outkey == k + "()" && imagmAdd != 0) {
        form += " " + outkey + to_string(numCol);
        numCol += imagmAdd;
        form += "-" + to_string(numCol - 1);
      }
    }

    if (outkey == "EM_FLUX()") {
      form += " EM_FLUX() " + to_string(numCol);
      numCol += 65;
      form += "-" + to_string(numCol - 1);
    }

    for (const string k : {"MOD_QSO", "MOD_STAR", "CHI_QSO", "CHI_STAR",
                           "LIMITS_ZMAX", "LIMITS_MFAINT"}) {
      if (outkey == k) {
        form += " " + k + " " + to_string(numCol);
        numCol++;
      }
    }

    for (const string k :
         {"MOD", "EXTLAW", "Z", "EBV", "PDZ", "CHI", "SCALE", "AGE"}) {
      if (outkey == k + "_SEC") {
        form += " " + outkey + " " + to_string(numCol);
        numCol++;
      }
    }

    for (const string k1 : {"AGE", "LDUST", "LUM_TIR", "MASS", "SFR", "SSFR",
                            "COL1", "COL2", "MREF"}) {
      for (const string k2 : {"MED", "INF", "SUP"}) {
        if (outkey == k1 + "_" + k2) {
          form += " " + outkey + " " + to_string(numCol);
          numCol++;
        }
      }
    }

    for (const string k1 : {"EM_FLUX", "EM_EW"}) {
      for (const string k2 :
           {"LYA", "OII", "HB", "OIIIA", "OIIIB", "HA", "SIIIA", "SIIIB"}) {
        if (outkey == k1 + "_" + k2) {
          form += " " + outkey + " " + to_string(numCol);
          numCol++;
        }
      }
    }
  }
  form += '\n';

  ////////////////////////////////////////////////////////////////////////
  // Format for topcat
  form += "# Format topcat: \n#";
  for (const auto &outkey : outkeywords) {
    for (const string k : {"IDENT", "CONTEXT", "NBAND_USED", "NBAND_ULIM",
                           "ZSPEC", "STRING_INPUT"}) {
      if (outkey == k) {
        form += " " + k + " ";
      }
    }

    for (const string k : {"Z", "ZQ"}) {  // ML is legacy, replaced by MED
      for (const string k1 : {"BEST", "ML", "MED", "MODE"}) {
        for (const string k2 : {"", "68_LOW", "68_HIGH", "90_LOW", "90_HIGH",
                                "99_LOW", "99_HIGH"}) {
          if (outkey == k + "_" + k1 + k2) {
            form += " " + outkey + " ";
          }
        }
      }
    }
    // legacy
    if (outkey == "Z_QSO") {
      form += " ZQ_BEST ";
    }

    for (const string k :
         {"MOD", "EXTLAW", "EBV", "MASS", "LDUST", "LUM_TIR", "LUM_NUV",
          "LUM_R", "LUM_K", "PDZ", "CHI", "SCALE", "AGE", "SFR", "SSFR"}) {
      if (outkey == k + "_BEST") {
        form += " " + outkey + " ";
      }
    }

    for (const string k : {"MAG_OBS", "ERR_MAG_OBS", "MAG_MOD", "MABS_FILT",
                           "K_COR", "MAG_ABS", "EMAG_ABS"}) {
      if (outkey == k + "()") {
        for (int l = 0; l < imagm; l++) form += " " + k + to_string(l) + " ";
      }
    }

    for (const string k : {
             "MAG_PRED",
             "ABSMAG_PRED",
         }) {
      if (outkey == k + "()") {
        for (int l = 0; l < imagmAdd; l++) form += " " + k + to_string(l) + " ";
      }
    }

    if (outkey == "EM_FLUX()") {
      for (int l = 0; l < 65; l++) form += " EM_FLUX" + to_string(l) + " ";
    }

    for (const string k : {"MOD_QSO", "MOD_STAR", "CHI_QSO", "CHI_STAR",
                           "LIMITS_ZMAX", "LIMITS_MFAINT"}) {
      if (outkey == k) {
        form += " " + k + " ";
      }
    }

    for (const string k :
         {"MOD", "EXTLAW", "Z", "EBV", "PDZ", "CHI", "SCALE", "AGE"}) {
      if (outkey == k + "_SEC") {
        form += " " + outkey + " ";
      }
    }

    for (const string k1 : {"AGE", "LDUST", "LUM_TIR", "MASS", "SFR", "SSFR",
                            "COL1", "COL2", "MREF"}) {
      for (const string k2 : {"MED", "INF", "SUP"}) {
        if (outkey == k1 + "_" + k2) {
          form += " " + outkey + " ";
        }
      }
    }

    for (const string k1 : {"EM_FLUX", "EM_EW"}) {
      for (const string k2 :
           {"LYA", "OII", "HB", "OIIIA", "OIIIB", "HA", "SIIIA", "SIIIB"}) {
        if (outkey == k1 + "_" + k2) {
          form += " " + outkey + " ";
        }
      }
    }

  }  // for loop on outkeywords
  form += '\n';

  return form;
}

/*
 READ THE SOURCE IN THE INPUT CATALOGUE, CHOICE BETWEEN VARIOUS FORMATS
*/
void PhotoZ::readsource(onesource *src, const string line) {
  double dab, dsab;

  // put the line into the stream ss
  stringstream ss(line);
  // Read the identifiant
  ss >> src->spec;

  // Read the flux/mag and associated errors
  if (cat_fmt == 0) {  // Case mag err mag err MEME

    // Loop over each filter and store the flux+err
    for (int k = 0; k < imagm; k++) {
      ss >> dab >> dsab;
      src->ab.push_back(dab);
      src->sab.push_back(dsab);
    }

  } else {  // Case mag mag err err MMEE

    // Loop over each filter and store the flux
    for (int k = 0; k < imagm; k++) {
      ss >> dab;
      src->ab.push_back(dab);
    }

    // Loop over each filter and store the errors
    for (int k = 0; k < imagm; k++) {
      ss >> dsab;
      src->sab.push_back(dsab);
    }
  }

  // Read the context/spec-z/final string in case of the LONG format
  if (cattyp[0] == 'L' || cattyp[0] == 'l') {
    ss >> src->cont;            // context
    ss >> src->zs;              // spectro-z
    getline(ss, src->str_inp);  // store the end of the file into a string
  }

  return;
}

/*
  read the sources which are used for the adaotation of the zero-points
 */
vector<onesource *> PhotoZ::read_autoadapt_sources() {
  string line;
  // Vector of the objects with a spec-z
  vector<onesource *> adaptSources;
  ifstream sin(cat.c_str());
  // Read all the sources for auto-adapt, store them
  int nobj = 0;
  while (getline(sin, line)) {
    // If the first character of the line is not #
    if (check_first_char(line)) {
      // Construct one objet
      onesource *oneObj = yield(nobj, line);

      // Keep only sources with a spectroscopic redshift
      if (oneObj->zs > adzmin && oneObj->zs < adzmax) {
        // Correct the observed magnitudes and fluxes with the coefficients
        // given in APPLY_SHIFTS
        if (shifts0.size() == (size_t)imagm)
          oneObj->adapt_mag(shifts0, shifts1);
        // Keep all the objects within a vector if in the right mag range
        double magSel;
        if (oneObj->ab[fl_auto] > 0)
          magSel = oneObj->mab[fl_auto];
        else
          magSel = HIGH_MAG;
        if (magSel > auto_thresmin && magSel < auto_thresmax) {
          nobj++;
          adaptSources.push_back(oneObj);
        } else {
          delete oneObj;
        }
      } else {
        delete oneObj;
      }
    }
  }
  return adaptSources;
}

/*
  Run the fit in order to get an adaptation of the zero-points
  Median of the difference between the modeled magnitudes and the observed ones
*/
std::tuple<vector<double>, vector<double>> PhotoZ::run_autoadapt(
    vector<onesource *> adaptSources) {
  double funz0 = lcdm.distMod(gridz[1] / 20.);
  vector<double> a0, a1;
  a0.assign(imagm, 0.);
  a1.assign(imagm, 0.);
  // Use the spec-z for the adpation
  const bool zfix = true;
  if (verbose)
    cout << "\n Number of sources read for auto adapt " << adaptSources.size()
         << endl;
  // Adaptation only if some sources selected for auto-adapt
  if (adaptSources.size() > 0) {
    int iteration = 0;
    int converge = 0;
    // While the convergence is not reached and we have less than 10 iterations
    while (converge == 0 && iteration < 10) {
      // Loop over the sources
      for (auto &oneObj : adaptSources) {
        // Correct the observed magnitudes and fluxes with the coefficients
        // found by auto-adapt
        oneObj->adapt_mag(a0, a1);
        // set the prior on the redshift, abs mag, ebv, etc on the object
        oneObj->setPriors(magabsB, magabsF);
        // Fixed the redshift considered
        oneObj->considered_red(zfix, methz);
        oneObj->closest_red = gridz[indexz(oneObj->zs, gridz)];
        // Select the valid index of the library in case of ZFIX=YES to save
        // computational time
        valid = oneObj->validLib(fullLib, zfix, oneObj->zs);
        // Fit the source at the spec-z value
        oneObj->fit(fullLib, flux, valid, funz0, bp);
        // Interpolation of the predicted magnitudes, scaling at zs, checking
        // first that the fit was sucessfull
        if (oneObj->indmin[0] >= 0) {
          oneObj->zmin[0] = oneObj->zs;
          oneObj->interp_lib(fullLib, imagm, lcdm);
        }
        if (verbose)
          cout << " Fit source for adapt " << oneObj->spec << "  \r " << flush;
      }
      // run auto-adapt
      auto_adapt(adaptSources, a0, a1, converge, iteration);
      // Display the results
      if (verbose) {
        cout << "Offsets:  ";
        for (int k = 0; k < imagm - 1; k++) cout << a0[k] << ",";
        cout << a0[imagm - 1] << endl;
      }
    }
  }
  if (verbose) cout << endl << "Done with auto-adapt " << endl;
  return std::make_tuple(a0, a1);
}

/*
  function to compare observed magnitudes and predicted ones
*/
void auto_adapt(const vector<onesource *> adaptSources, vector<double> &a0,
                vector<double> &a1, int &converge, int &iteration) {
  vector<double> diff, a0pre, a1pre;
  double inter;

  // Keep the values of a0 and a1 to check the convergence
  a0pre.swap(a0);
  a1pre.swap(a1);
  a0.clear();
  a1.clear();
  // Number of filters
  int imagm = (adaptSources[0]->ab).size();

  // Store the median of the difference in each band between the predicted and
  // the observed magnitudes
  for (int k = 0; k < imagm; k++) {
    diff.clear();
    // define a vector difference between the observed and predicted mag
    for (auto &oneObj : adaptSources) {
      // Only in the case of a positive flux and a fit successfully performed
      if (oneObj->mab_ori[k] > 0 && oneObj->indmin[0] > 0 &&
          oneObj->busnorma[k] == 1) {
        // Same convention as before. Need to be added to predicted magnitude,
        // or substracted to observed magnitudes
        inter = oneObj->mab_ori[k] - oneObj->magm[k];
        // If difference smaller than 1
        if (abs(inter) < 1) diff.push_back(inter);
      }
    }
    // If less than 10 elements, don't do it
    if (diff.size() > 10) {
      // Sort the vector diff in increasing order
      sort(diff.begin(), diff.end(),
           [](double a, double b) { return (a > b); });
      // Store the median in a0
      a0.push_back(diff[int(diff.size() / 2.)]);
      a1.push_back(0.);
    } else {
      a0.push_back(0.);
      a1.push_back(0.);
    }
    // cout << "Adaptation a0 : " << k << " " << a0[k] << " " << endl;
  }

  // Compare the precedent value of a0 with the new ones to check for
  // convergence
  double somme = 0.;
  for (int k = 0; k < imagm; k++) somme += abs(a0[k] - a0pre[k]);
  if (somme / double(imagm) < 0.02)
    converge = 1;  // average diff over the filter lower than 0.02
  // If one band varies by more than 0.03 mag, go on with iteration
  for (int k = 0; k < imagm; k++)
    if (abs(a0[k] - a0pre[k]) > 0.03) converge = 0;

  // Stop before iteratio 10
  if (iteration == 10) converge = 1;
  iteration++;

  cout << " Done with iteration " << iteration << " and converge " << converge
       << endl;

  return;
}

/*
   Determine the best filter to be used as a function of redshift
*/
vector<vector<int>> bestFilter(int nbFlt, vector<double> gridz,
                               vector<SED *> fulllib, int method,
                               vector<long> magabscont, vector<int> bapp,
                               vector<int> bappOp, vector<double> zbmin,
                               vector<double> zbmax) {
  vector<vector<int>> bestFlt;
  vector<int> bestFlt_z;
  int j = 0;
  // Loop over the grid in redshift
  for (size_t i = 0; i < gridz.size(); i++) {
    bestFlt_z.clear();

    // Define the best filter depending on the method
    switch (method) {
        // Use the same filter for the apparent magnitude as for the  absolute
        // one
      case 0:
        for (int k = 0; k < nbFlt; ++k) {
          bestFlt_z.push_back(k);
        }
        break;

        // method 2 : Fix bapp as observed-frame
      case 2:
        if (bapp.size() == (size_t)nbFlt) {
          for (int k = 0; k < nbFlt; ++k) {
            bestFlt_z.push_back(bapp[k]);
          }
        } else {
          for (int k = 0; k < nbFlt; ++k) {
            bestFlt_z.push_back(bapp[0]);
          }
        }
        break;

        // method 4 : Fix bapp as observed-frame depending on the redshift
      case 4:
        j = 0;
        // change the redshift bin as long as we are not in the right redshift
        // bin
        while (!(gridz[i] >= zbmin[j] && gridz[i] < zbmax[j]) &&
               (size_t)j < zbmax.size() - 1)
          j++;
        for (int k = 0; k < nbFlt; ++k) {
          bestFlt_z.push_back(bappOp[j]);
        }
        break;

        // method 3, 1 and default
      default:
        for (int k = 0; k < nbFlt; ++k) {
          bestFlt_z.push_back(-1);
        }
    }

    // store the best filter at this redshift
    bestFlt.push_back(bestFlt_z);

  }  // end loop over the redshift grid

  // Method 1 to change observe frame as a function of redshift
  if (method == 1) {
    minimizekcolor(gridz, fulllib, bestFlt, magabscont);
  }

  return bestFlt;
}

/*
  Derive the k-color term (rest-frame color - k-correction) maximum, in order to
  conserve the error bars for the absolute magnitudes
*/
vector<vector<double>> maxkcolor(vector<double> gridz, vector<SED *> fulllib,
                                 vector<vector<int>> bestFlt) {
  vector<vector<double>> extremeDiff, extremeDiffMin, extremeDiffMax;
  vector<double> extremeValMin, extremeValMax, extremeVal, di;
  int indz;
  double kcolor;

  // initialize the minimum and maximum values
  for (size_t i = 0; i < bestFlt.size(); i++) {       // Loop over all redshifts
    for (size_t k = 0; k < bestFlt[0].size(); ++k) {  // loop over all filters
      extremeValMin.push_back(1.e4);
      extremeValMax.push_back(-1.e4);
    }
    extremeDiffMin.push_back(extremeValMin);
    extremeDiffMax.push_back(extremeValMax);
  }

  // Find the extreme value of k-color
  // Loop over all SEDs from the library
  // parrallellize the loop
#ifdef _OPENMP
#pragma omp parallel private(indz, kcolor)
  {
#pragma omp for
#endif
    for (vector<SED *>::iterator it = fulllib.begin(); it < fulllib.end();
         ++it) {
      // computed only for galaxies
      if ((*it)->nlib == 0) {
        // loop over the rest-frame filters to be considered
        for (size_t k = 0; k < bestFlt[0].size(); ++k) {
          // obtain the index of the redshift grid
          indz = indexz((*it)->red, gridz);
          // derive the k-color term at this redshift: m(ref)-m(obs)-kcorr(ref)
          if (bestFlt[indz][k] >= 0) {
            kcolor = ((fulllib[(*it)->index_z0])->mag[k]) -
                     ((fulllib[(*it)->index_z0])->mag[bestFlt[indz][k]]) -
                     (*it)->kcorr[bestFlt[indz][k]];
            // kcolor= (*it)->mag[k] - (*it)->mag[bestFlt[indz][k]] -
            // (*it)->kcorr[k]; Store it if the value is maximum or minimum
            if (extremeDiffMin[indz][k] > kcolor)
              extremeDiffMin[indz][k] = kcolor;
            if (extremeDiffMax[indz][k] < kcolor)
              extremeDiffMax[indz][k] = kcolor;
          }
        }
      }
    }
#ifdef _OPENMP
  }
#endif

  // Do the difference of the extreme values
  for (size_t i = 0; i < bestFlt.size(); i++) {  // Loop over all redshift
    di.clear();
    for (size_t k = 0; k < bestFlt[0].size(); ++k) {  // loop over all filters
      di.push_back(extremeDiffMax[i][k] - extremeDiffMin[i][k]);
    }
    extremeDiff.push_back(di);
  }

  return extremeDiff;
}

/*
  Define the filters to pick (depending on the redshift) to minimize the
  k-correction + rest-frame color term in the absolute magnitude computation
*/
void minimizekcolor(vector<double> gridz, vector<SED *> fulllib,
                    vector<vector<int>> &bestFlt, vector<long> magabscont) {
  vector<vector<double>> extremeDiffMin, extremeDiffMax;
  vector<double> extremeValMin, extremeValMax, extremeVal, di;
  vector<int> possiFlt;
  int indz;
  double kcolor, mini;

  // bestFlt has been initialize at -1 before the call, with a dimention
  // redshift grid/nb filt
  int imagm = bestFlt[0].size();
  int nbz = bestFlt.size();

  // initialize the maximum and minimum values
  for (int i = 0; i < nbz; i++) {      // Loop over all redshifts
    for (int k = 0; k < imagm; k++) {  // loop over all filters
      extremeValMin.push_back(1.e4);
      extremeValMax.push_back(-1.e4);
    }
    extremeDiffMin.push_back(extremeValMin);
    extremeDiffMax.push_back(extremeValMax);
  }

  // List of possible filters for the observed magnitudes, given the context
  for (int k = 0; k < imagm; ++k) {
    if (magabscont[0] > 0) {
      // Only if the context allows this filter
      if (magabscont.size() == (size_t)imagm) {
        if (CHECK_CONTEXT_BIT(magabscont[k], k)) possiFlt.push_back(k);
      } else {
        if (CHECK_CONTEXT_BIT(magabscont[0], k)) possiFlt.push_back(k);
      }
    } else {
      possiFlt.push_back(k);
    }
  }
  int nbPossi = possiFlt.size();

  // loop over the reference filters to be considered
  for (int r = 0; r < imagm; ++r) {
    // reinitialize
    for (int i = 0; i < nbz; i++) {      // Loop over all redshifts
      for (int k = 0; k < imagm; ++k) {  // loop over all obs filters
        extremeDiffMin[i][k] = 1.e4;
        extremeDiffMax[i][k] = -1.e4;
      }
    }

    // Find the extreme value of k-color
    // Loop over all SEDs from the library
    // parrallellize the loop
#ifdef _OPENMP
#pragma omp parallel private(indz)
    {
#pragma omp for
#endif
      for (vector<SED *>::iterator it = fulllib.begin(); it < fulllib.end();
           ++it) {
        // computed only for galaxies
        if ((*it)->nlib == 0) {
          // obtain the index of the redshift grid
          indz = indexz((*it)->red, gridz);

          // loop over the possible filters for the apparent magnitudes
          for (int l = 0; l < nbPossi; ++l) {
            // derive the k-color term at this redshift:
            // (m(ref)-m(obs))(z=0)-kcorr(obs)
            kcolor = ((fulllib[(*it)->index_z0])->mag[r]) -
                     ((fulllib[(*it)->index_z0])->mag[possiFlt[l]]) -
                     (*it)->kcorr[possiFlt[l]];

            // Store it if the value is maximum or minimum
            if (extremeDiffMin[indz][possiFlt[l]] > kcolor)
              extremeDiffMin[indz][possiFlt[l]] = kcolor;
            if (extremeDiffMax[indz][possiFlt[l]] < kcolor)
              extremeDiffMax[indz][possiFlt[l]] = kcolor;
          }
        }
      }
#ifdef _OPENMP
    }
#endif

    // Do the difference between the extreme values
    for (int i = 0; i < nbz; i++) {  // Loop over all redshifts
      mini = 1.e9;
      for (int k = 0; k < nbPossi; ++k) {  // loop over all filters
        // keep the filter which is able to minimize the k-color, save it into
        // bestFlt
        if ((extremeDiffMax[i][possiFlt[k]] - extremeDiffMin[i][possiFlt[k]]) <
            mini) {
          bestFlt[i][r] = possiFlt[k];
          mini =
              (extremeDiffMax[i][possiFlt[k]] - extremeDiffMin[i][possiFlt[k]]);
        }
      }
    }
  }

  return;
}

vector<onesource *> PhotoZ::read_photoz_sources() {
  vector<onesource *> photoz_sources;
  // open the external file with zspec
  ifstream szex;
  string externalzfile = ((keys["EXTERNALZ_FILE"]).split_string("NONE", 1))[0];
  if (externalzfile.substr(0, 4) != "NONE") {
    szex.open(externalzfile.c_str());
    if (!szex) {
      cout << "External spec-z option, but no file " << externalzfile << endl;
      exit(0);
    }
    string linezex;
    // Ignore the comments
    int nbcomments = 0;
    while (!(check_first_char(linezex))) {
      getline(szex, linezex);
      nbcomments++;
    }
    // back to the beginning of the file
    szex.seekg(0, ios::beg);
    ;
    // Go directly to the right lines, skip commented lines
    for (int k = 1; k < nbcomments; k++) {
      getline(szex, linezex);
      cout << "skip comments " << '\n';
    }  // go to the right starting row of the file
    // Go directly to the right lines, skipping lines if CAT_LINES
    for (unsigned int k = 1; k < rowmin; k++) {
      getline(szex, linezex);
      cout << "done skip " << k << " " << rowmin << '\n';
    }  // go to the right starting row of the file
  }

  // Take the stream line by line
  unsigned int nobj = 0;
  string line;

  // open the input file
  ifstream sin(cat.c_str());
  if (!sin) {
    cout << "No input file " << cat << endl;
    exit(0);
  }
  // Back to the beginning of the file
  sin.clear();
  sin.seekg(0, sin.beg);
  while (getline(sin, line)) {
    // If the first character of the line is not #
    if (check_first_char(line)) {
      nobj++;
      if ((nobj < rowmin && rowmin > 0) || (nobj > rowmax && rowmax > 0)) {
        continue;
      }  // CAT_LINES option

      // Generate one objet
      onesource *oneObj = yield(nobj, line);

      // Use zspec from external file
      // open the external file with zspec
      if (externalzfile.substr(0, 4) != "NONE") {
        string idzex, linezex;
        getline(szex, linezex);
        stringstream sszex(linezex);
        sszex >> idzex;
        if (idzex != oneObj->spec)
          cout << endl
               << "ERROR: mismatch in the external file " << idzex << " "
               << oneObj->spec << endl;
        sszex >> oneObj->zs;
      }

      // Add the source
      photoz_sources.push_back(oneObj);
    }
  }
  return photoz_sources;
}

/*
  Additional layer to prepare the data for the run (context, flux and asscoiated
  uncertainties)
*/
void PhotoZ::prep_data(onesource *oneObj) {
  // Convert the magnitude in fluxes if needed
  if (typm[0] == 'M') oneObj->convertFlux(catmag, allFilters);
  // Rescale the flux errors if needed
  oneObj->rescale_flux_errors(min_err, fac_err);
  // turn fluxes and possibly rescaled flux errors into magnitudes
  oneObj->convertMag();
  // Keep original magnitudes
  oneObj->keepOri();
  // Define the filters used for the fit based on the context
  oneObj->fltUsed(gbcont, contforb, imagm);
  return;
}

void PhotoZ::prep_data(vector<onesource *> sources) {
  // Loop over all sources
  for (auto &oneObj : sources) {
    prep_data(oneObj);
  }
  return;
}

/*
  Central part of the code to fit the templates and measure the photo-z
*/
void PhotoZ::run_photoz(vector<onesource *> sources, const vector<double> &a0,
                        const vector<double> &a1) {
  // Open the output file
  // AUTO_ADAPT adaptation of the zero-points
  bool autoadapt = keys["AUTO_ADAPT"].split_bool("NO", 1)[0];
  // RM_DISCREPANT_BD
  // Threshold in chi2 to consider. Remove <3 bands, stop when below this chi2
  // threshold
  double thresholdChi2 =
      ((keys["RM_DISCREPANT_BD"]).split_double("1.e9", 2))[0];

  // Parabolic interpolation of the redshift
  bool zintp = keys["Z_INTERP"].split_bool("NO", 1)[0];
  // EXTERNALZ_FILE file with the redshifts to be fixed - NONE default
  string externalzfile = ((keys["EXTERNALZ_FILE"]).split_string("NONE", 1))[0];
  // DZ_WIN minimal delta z window to search - 0.25 by default
  double dz_win = ((keys["DZ_WIN"]).split_double("0.25", 1))[0];

  // RF_COLORS compute 2 rest-frame colors and associated errorbars
  vector<int> fltColRF = (keys["RF_COLORS"]).split_int("-1", 4);
  // Need to substract one because the convention in the .para file start at 1,
  // but 0 in the code
  if (fltColRF.size() > 1) {
    for (int k = 0; k < 4; k++) fltColRF[k]--;
  }
  // M_REF compute the absolute magnitudes and associated errorbars
  int fltREF = ((keys["M_REF"]).split_int("0", -1))[0];
  // Need to substract one because the convention in the .para file start at 1,
  // but 0 in the code
  if (fltREF >= 0) {
    fltREF--;
  }
  // LIMITS_ZBIN Compute the z_max and M_faint in several bins of redshift. Give
  // the z bin.
  vector<double> limits_zbin =
      (keys["LIMITS_ZBIN"]).split_double("0.0,90.", -1);
  int nzbin = int(limits_zbin.size()) - 1;
  // LIMITS_MAPP_REF Compute the z_max and M_faint in several bins of redshift.
  // Give the reference band.
  int limits_ref = ((keys["LIMITS_MAPP_REF"]).split_int("1", 1))[0];
  // LIMITS_MAPP_SEL Compute the z_max and M_faint in several bins of redshift.
  // Give the selection band in each bin.
  vector<int> limits_sel = (keys["LIMITS_MAPP_SEL"]).split_int("1", nzbin);
  // LIMITS_MAPP_CUT Compute the z_max and M_faint in several bins of redshift.
  // Give the cut in magnitude in each bin.
  vector<double> limits_cut =
      (keys["LIMITS_MAPP_CUT"]).split_double("90.", nzbin);

  // FIR_LIB name of the library in IR
  vector<string> libext = (keys["FIR_LIB"]).split_string("NONE", -1);
  int nlibext = int(libext.size());
  // if the library is NONE, put the number of library at -1
  if (nlibext == 1 && libext[0] == "NONE") nlibext = -1;
  // FIR_LMIN Lambda min given in micron  (um)
  double fir_lmin = ((keys["FIR_LMIN"]).split_double("7.0", 1))[0];
  // FIR_CONT context for the far-IR
  long fir_cont = ((keys["FIR_CONT"]).split_long("-1", 1))[0];
  // FIR_SCALE context of the bands used for the rescaling in IR
  int fir_scale = ((keys["FIR_SCALE"]).split_int("-1", 1))[0];
  // FIR_FREESCALE possible free rscaling in IR, when several bands. Otherwise,
  // model imposed by its LIR
  string fir_frsc = ((keys["FIR_FREESCALE"]).split_string("NO", 1))[0];
  // FIR_SUBSTELLAR remove the stellar component
  bool substar = keys["FIR_SUBSTELLAR"].split_bool("NO", 1)[0];
  // MIN_THRES threshold to trigger the detection - 0.1 by default
  double min_thres = ((keys["MIN_THRES"]).split_double("0.1", 1))[0];

  /* Define what are the filters to be used for the absolute magnitude depending
   * on the method adopted */
  // MABS_METHOD method to compute the absolute magnitudes
  int method = ((keys["MABS_METHOD"]).split_int("0", -1))[0];
  // MABS_REF reference filter in case of mag abs method 2
  vector<int> bapp = (keys["MABS_REF"]).split_int("1", -1);
  int nbapp = int(bapp.size());

  // MABS_FILT choose filters per redshift bin (MABS_ZBIN) if method 4
  vector<int> bappOp = (keys["MABS_FILT"]).split_int("1", -1);
  int nbBinZ = int(bappOp.size());

  // Need to substract one because the convention in the .para file start at 1,
  // but 0 in the code
  for (int k = 0; k < nbBinZ; k++) bappOp[k]--;
  // Need to substract one because the convention in the .para file start at 1,
  // but 0 in the code
  for (int k = 0; k < nbapp; k++) bapp[k]--;

  // MABS_ZBIN give the redshift bins corresponding to MABS_FILT
  vector<double> zbmin = (keys["MABS_ZBIN"]).split_double("0", nbBinZ + 1);
  zbmin.erase(zbmin.end() - 1);
  vector<double> zbmax = (keys["MABS_ZBIN"]).split_double("6", nbBinZ + 1);
  zbmax.erase(zbmax.begin(), zbmax.begin() + 1);
  // MABS_CONTEXT context for method 1
  vector<long> magabscont = (keys["MABS_CONTEXT"]).split_long("0", -1);
  vector<vector<int>> goodFlt = bestFilter(
      imagm, gridz, fullLib, method, magabscont, bapp, bappOp, zbmin, zbmax);
  vector<vector<double>> maxkcol = maxkcolor(gridz, fullLib, goodFlt);

  double funz0 = lcdm.distMod(gridz[1] / 20.);
  vector<opa> opaOut = Mag::read_opa();

  // Specify the offsets in the header
  string offsets;
  if (autoadapt) {
    for (int k = 0; k < imagm; k++) offsets = offsets + to_string(a0[k]) + ",";
    offsets = "# Offsets from auto-adapt: " + offsets + '\n';
    outputHeader += offsets;
  }
  if (shifts0.size() == (size_t)imagm) {
    offsets = "";
    for (int k = 0; k < imagm; k++)
      offsets = offsets + to_string(shifts0[k]) + ",";
    offsets = "# Offsets applied directly from keyword: " + offsets + '\n';
    outputHeader += offsets;
  }

  unsigned int nobj = 0;
  for (auto &oneObj : sources) {
    if (verbose)
      cout << "Fit source " << nobj << " with Id " << oneObj->spec << " \r "
           << flush;
    nobj++;
    // Correct the observed magnitudes and fluxes with the coefficients given in
    // APPLY_SHIFTS
    if (shifts0.size() == (size_t)imagm) oneObj->adapt_mag(shifts0, shifts1);
    // Correct the observed magnitudes and fluxes with the coefficients found by
    // auto-adapt
    if (autoadapt) oneObj->adapt_mag(a0, a1);
    oneObj->closest_red = gridz[indexz(oneObj->zs, gridz)];
    // set the prior on the redshift, abs mag, ebv, etc on the object
    oneObj->setPriors(magabsB, magabsF);
    // Select the valid index of the library in case of ZFIX=YES to save
    // computational time
    valid = oneObj->validLib(fullLib, zfix, oneObj->zs);
    // Core of the program: compute the chi2
    oneObj->fit(fullLib, flux, valid, funz0, bp);
    // Try to remove some bands to improve the chi2, only as long as the chi2 is
    // above a threshold
    oneObj->rm_discrepant(fullLib, flux, valid, funz0, bp, thresholdChi2,
                          verbose);
    // Generate the marginalized PDF (z+physical parameters) from the chi2
    // stored in each SED
    oneObj->generatePDF(fullLib, valid, fltColRF, fltREF, zfix);
    // Interpolation of Z_BEST and ZQ_BEST (zmin) via Chi2 curves, put z-spec if
    // ZFIX YES  (only gal for the moment)
    oneObj->interp(zfix, zintp, lcdm);
    // Uncertainties from the minimum chi2 + delta chi2
    oneObj->uncertaintiesMin();
    // Uncertainties from the bayesian method
    oneObj->uncertaintiesBay();
    // find a second peak in the PDZ
    oneObj->secondpeak(fullLib, dz_win, min_thres);
    // find the mode of the marginalized PDF and associated uncertainties
    oneObj->mode();
    // Fixed the redshift considered for abs mag, etc depending on the option
    oneObj->considered_red(zfix, methz);
    // If use the median of the PDF for the abs mag, etc, need to redo the fit
    // Need to redo the fit to get the right scaling. It would change ZMIN, etc
    if ((methz[0] == 'M' || methz[0] == 'm') && (!zfix)) {
      oneObj->chimin[0] = 1.e9;
      oneObj->closest_red = gridz[indexz(oneObj->zgmed[0], gridz)];
      oneObj->fit(fullLib, flux, valid, funz0, bp);
    }
    // Interpolation at the new redshift  (only gal for the moment)
    oneObj->interp_lib(fullLib, imagm, lcdm);
    // Compute absolute magnitudes
    oneObj->absmag(goodFlt, maxkcol, lcdm, gridz);
    // Compute zmax and M_faint
    oneObj->limits(fullLib, limits_zbin, limits_ref, limits_sel, limits_cut);
    // Compute predicted magnitude in new filters
    if (allFiltersAdd.size() > 0) {
      oneObj->computePredMag(fullLib, lcdm, opaOut, allFiltersAdd);
      oneObj->computePredAbsMag(fullLib, lcdm, opaOut, allFiltersAdd);
    }
    // Compute flux of emission lines
    oneObj->computeEmFlux(fullLib, lcdm, opaOut);
    if (nlibext > 0) {
      // FIR FIT
      // Define the filters used for the FIR fit based on the FIR context
      oneObj->fltUsedIR(fir_cont, fir_scale, imagm, allFilters, fir_lmin);
      // Substract the stellar component to the FIR observed flux
      oneObj->substellar(substar, allFilters);
      oneObj->closest_red = gridz[indexz(oneObj->consiz, gridz)];
      // Fit the SED on FIR data, with the redshift fixed at zmin or zmed
      oneObj->fitIR(fullLibIR, fluxIR, imagm, fir_frsc, lcdm);
      // Compute the IR luminosities
      oneObj->generatePDF_IR(fullLibIR);
    }
    // compute physical quantities for the best fit GAL solution
    oneObj->compute_best_fit_physical_quantities(fullLib);
  }
  return;
}

void PhotoZ::write_outputs(vector<onesource *> sources, const time_t &ti1) {
  // CAT_OUT output  file -  zphot.out default
  string outf = ((keys["CAT_OUT"]).split_string("zphot.out", 1))[0];
  ofstream stout;
  stout.open(outf.c_str());
  // Start the header
  stout << "# Creation date: " << asctime(localtime(&ti1));
  /* Read the output parameter file */
  vector<string> outkeywords = readOutKeywords(outpara);
  /* Add the format to the header */
  outputHeader += prep_header(outkeywords);
  // Write the header
  stout << outputHeader;

  // If the pdf(z) is requested in output, open the several streams
  unordered_map<string, ofstream> pdf_streams;
  if (outpdz.compare(nonestring) != 0) {
    for (const auto &type : pdftype) {
      string output = outpdz + "_" + type + ".prob";
      pdf_streams[type].open(output.c_str());
    }
  }

  vector<opa> opaOut = Mag::read_opa();

  static bool first_obj = true;
  for (auto &oneObj : sources) {
    // write the object in output
    oneObj->write_out(fullLib, fullLibIR, stout, outkeywords);
    // Write an ascii file with the best fit template
    if (outsp.compare("NO") != 0)
      oneObj->writeSpec(fullLib, fullLibIR, lcdm, opaOut, allFilters, outsp);
    // write the full chi2 if asked (could take a lot of space)
    if (outchi) oneObj->writeFullChi(fullLib);

    // write the PDF
    if ((outpdz.compare(nonestring) != 0) && first_obj)
      oneObj->write_pdz_header(pdftype, pdf_streams, ti1);
    if (outpdz.compare(nonestring) != 0)
      oneObj->write_pdz(pdftype, pdf_streams);
    first_obj = false;
  }

  if (outpdz.compare(nonestring) != 0)
    for (const auto &type : pdftype) {
      pdf_streams[type].close();
    }
  stout.close();

  for (auto it = fullLib.begin(); it != fullLib.end(); it++) {
    delete *it;
  }
  fullLib.clear();
  for (auto it = fullLibIR.begin(); it != fullLibIR.end(); it++) {
    delete *it;
  }
  fullLibIR.clear();

  return;
}
