/*
    01/04/2015
    Main code to estimate the photometric redshifts
    zphota   -c zphot.para
    zphot.para  : config file
*/

#include <getopt.h>  // get option line commands getopt_long

#include <cmath>     // for the log
#include <cstring>   // Use string in the old c format
#include <ctime>     // date
#include <fstream>   // print output file
#include <iomanip>   // std::set precision
#include <iostream>  // print standard file
#include <sstream>
#include <string>  // use string
#include <vector>  // manipulate vector

// Le Phare
#include "SED.h"        // to read the libraries
#include "cosmology.h"  // in order to measure the distance modulus
#include "flt.h"        // filter class
#include "globals.h"    // global variables
#include "keyword.h"    // to read the keyword in the command line/para file
#include "onesource.h"  // class including all the attributes and the operations on one object
#include "photoz_lib.h"  // class including helper functions

using namespace std;

// START
int main(int argc, char *argv[]) {
  // to estimate the times needed to run the code
  clock_t t1, t2;
  t1 = clock();
  cout << endl
       << "Starting times in sec: " << (((float)t1) / CLOCKS_PER_SEC) << endl;
  time_t ti1 = time(nullptr);
  cout << "local:     " << asctime(localtime(&ti1)) << '\n';

  /*
  ENVIRONMENT VARIABLES LEPHAREDIR and LEPHAREWORK
  */
  get_lephare_env();

  /*
  ANALYSE KEYWORDS
  */

  // List of the keywords to be found in the config file/command line
  string list_keywords[] = {
      "CAT_IN",          "INP_TYPE",       "CAT_MAG",
      "CAT_FMT",         "CAT_LINES",      "ZPHOTLIB",
      "PARA_OUT",        "CAT_OUT",        "VERBOSE",
      "CAT_TYPE",        "ERR_SCALE",      "ERR_FACTOR",
      "BD_SCALE",        "GLB_CONTEXT",    "FORB_CONTEXT",
      "MASS_SCALE",      "MAG_ABS",        "MAG_ABS_AGN",
      "MAG_REF",         "NZ_PRIOR",       "ZFIX",
      "Z_INTERP",        "EXTERNALZ_FILE", "RM_DISCREPANT_BD",
      "Z_RANGE",         "EBV_RANGE",      "DZ_WIN",
      "MIN_THRES",       "PROB_INTZ",      "SPEC_OUT",
      "FIR_LIB",         "FIR_LMIN",       "FIR_CONT",
      "FIR_SCALE",       "FIR_FREESCALE",  "FIR_SUBSTELLAR",
      "MABS_METHOD",     "MABS_CONTEXT",   "MABS_REF",
      "MABS_FILT",       "MABS_ZBIN",      "RF_COLORS",
      "M_REF",           "Z_METHOD",       "APPLY_SYSSHIFT",
      "AUTO_ADAPT",      "ADAPT_BAND",     "ADAPT_LIM",
      "ADAPT_CONTEXT",   "ADAPT_ZBIN",     "PDZ_OUT",
      "PDZ_TYPE",        "PDZ_MABS_FILT",  "ADD_EMLINES",
      "ADDITIONAL_MAG",  "LIMITS_ZBIN",    "LIMITS_MAPP_REF",
      "LIMITS_MAPP_SEL", "LIMITS_MAPP_CUT"};
  // Number of keywords
  int nb_ref_key = (int)(sizeof(list_keywords) / sizeof(list_keywords[0]));

  // Find the reference keywords within the command line + config file -> store
  // them into a vector
  keymap key_analysed = analyse_keywords(argc, argv, list_keywords, nb_ref_key);
  bool autoadapt = key_analysed["AUTO_ADAPT"].split_bool("NO", 1)[0];

  // Construct the photo-z run with an analyze of all the keywords
  PhotoZ photoz = PhotoZ(key_analysed);

  /*
   AUTO ADAPT
  */

  // shifts derived by auto-adapt
  vector<double> a0, a1;
  if (autoadapt) {
    // Read the sources to be used for auto-adapt
    vector<onesource *> adaptSrcs = photoz.read_autoadapt_sources();
    // Compute the offsets
    std::tie(a0, a1) = photoz.run_autoadapt(adaptSrcs);
    cout << " Done " << endl;
    // Clean the vector
    for (auto &src : adaptSrcs) delete src;
    adaptSrcs.shrink_to_fit();
  } else {
    for (int k = 0; k < photoz.imagm; k++) {
      a0.push_back(0.);
      a1.push_back(0.);
    }
  }

  /*
   PHOTOZ ESTIMATION
  */

  // Read the sources for which we want photo-z
  vector<onesource *> fitSrcs = photoz.read_photoz_sources();
  photoz.run_photoz(fitSrcs, a0, a1);
  photoz.write_outputs(fitSrcs, ti1);

  // display the time needed to run the code
  t2 = clock();
  float diff = (((float)t2 - (float)t1) / CLOCKS_PER_SEC);
  cout << endl
       << "Time to run the code in sec, sum over all threads: " << diff << endl;
  time_t ti2 = time(nullptr);
  cout << "Local time at beginning:     " << asctime(localtime(&ti1));
  cout << "Local time at the end:     " << asctime(localtime(&ti2));

  for (auto &src : fitSrcs) delete src;
  fitSrcs.shrink_to_fit();

  return 0;
}
