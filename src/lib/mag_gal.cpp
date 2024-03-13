/*
    10/12/2014
    Construct the different models of evolution used to measure the photo-z

    2 options : mag_gal
                        -c zphot.para    : config file
                        -t (G/g or Q/q)  : GAL/QSO library
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
#include "SED.h"
#include "globals.h"  // global variables
#include "keyword.h"  // to read the keyword in the command line/para file
#include "mag.h"      // SED class

using namespace std;

// START
int main(int argc, char *argv[]) {
  int num_threads;

#ifdef _OPENMP
#pragma omp parallel
  num_threads = omp_get_num_threads();
  cout << "Running with OpenMP on " << num_threads << " threads!" << endl;
#endif

  // to estimate the times needed to run the code
  // since clock() returns the sum of all threads, we need to use the openmp
  // timer to measure the time
#ifdef _OPENMP
  double start = omp_get_wtime();
#else
  clock_t t1 = clock();
#endif

  /*
  ENVIRONMENT VARIABLES LEPHAREDIR and LEPHAREWORK
  */
  get_lephare_env();

  /*
  ANALYSE KEYWORDS
  */

  // List of the keywords to be found in the config file/command line
  string list_keywords[] = {
      "t",           "COSMOLOGY",     "FILTER_FILE",  "MAGTYPE",
      "EXTINC_LAW",  "EB_V",          "MOD_EXTINC",   "ZGRID_TYPE",
      "Z_STEP",      "GAL_LIB_IN",    "QSO_LIB_IN",   "STAR_LIB_IN",
      "GAL_LIB_OUT", "QSO_LIB_OUT",   "STAR_LIB_OUT", "LIB_ASCII",
      "EM_LINES",    "EM_DISPERSION", "ADD_DUSTEM",   "VERBOSE"};
  // Number of keywords
  int nb_ref_key = (int)(sizeof(list_keywords) / sizeof(list_keywords[0]));

  // Find the reference keywords within the command line + config file -> store
  // them into a vector
  keymap key_analysed = analyse_keywords(argc, argv, list_keywords, nb_ref_key);

  // type of source which is read (Galaxy G, QSO Q, Star S)
  string typ = key_analysed["t"].value;

  // keyword to add the LDUST component to the stellar component (e.g. in BC03)
  string addDust = key_analysed["ADD_DUSTEM"].value;

  // Define a pointer of the basis class "Mag" which encompasses all the
  // elements to create the library
  Mag *Magnitude;

  // GALAXY CASE
  if (typ[0] == 'G' || typ[0] == 'g') {
    Magnitude = new GalMag(key_analysed);

    // QSO CASE
  } else if (typ[0] == 'Q' || typ[0] == 'q') {
    Magnitude = new QSOMag(key_analysed);

    // STAR CASE
  } else if (typ[0] == 'S' || typ[0] == 's') {
    Magnitude = new StarMag(key_analysed);

  } else {
    cout << "The type is not correctly indicated with -t " << endl;
    exit(1);
  }

  /*
  OPEN FILES
  */

  Magnitude->open_files();

  /*
  INFO PARAMETERS ON SCREEN AND DOC
  */

  Magnitude->print_info();

  /*
  READ DUST EXTINCTION LAWS
  */

  Magnitude->read_ext();

  /*
  DEFINE THE REDSHIFT GRID
  */

  Magnitude->def_zgrid();

  /*
  READ B12 TEMPLATES TO ADD DUST EMISSION TO BC03
  */

  if (addDust[0] == 'Y' || addDust[0] == 'y') Magnitude->read_B12();

  /*
  READ SED, APPLY EXTINCTION AND IGM OPACITY
  */

  Magnitude->read_SED();

  /*
  WRITE DOCUMENTATION
  */
  Magnitude->write_doc();

  /*
  DONE
  */

  delete Magnitude;

  // display the time needed to run the code
#ifdef _OPENMP
  double end = omp_get_wtime();
  double diff = end - start;
#else
  clock_t t2 = clock();
  float diff(((float)t2 - (float)t1) / CLOCKS_PER_SEC);
#endif

  cout << endl << "Times to run the code in sec: " << diff << endl;

  return 0;
}
