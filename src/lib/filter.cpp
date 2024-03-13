

/*
    10/11/2014

    Create a unique filter file (FILTER_FILE) with
    all the single files combined in the list FILTER_LIST

    Derive useful quantities as the AB corrections
*/

#include <getopt.h>  // get option line commands getopt_long
#include <stdlib.h>  // absolute value, exit

#include <cstring>   // Use string in the old c format
#include <fstream>   // print output file
#include <iostream>  // print standard file
#include <sstream>
#include <string>  // use string
#include <vector>  // manipulate vector

// Le Phare
#include "SED.h"
#include "flt.h"      // define the filters
#include "globals.h"  // global variables
#include "keyword.h"  // to read the keyword in the command line/para file

using namespace std;

// START
int main(int argc, char *argv[]) {
  // Declare variables
  vector<flt> vecFlt;

  // to estimate the times needed to run the code
  clock_t t1, t2;
  t1 = clock();

  /*
  ENVIRONMENT VARIABLES LEPHAREDIR and LEPHAREWORK
  */
  get_lephare_env();

  /*
  ANALYSE KEYWORDS
  */

  // List of the keywords to be found in the config file/command line
  string list_keywords[] = {"FILTER_REP", "FILTER_LIST", "TRANS_TYPE",
                            "FILTER_CALIB", "FILTER_FILE"};
  // Number of keywords
  int nb_ref_key = (int)(sizeof(list_keywords) / sizeof(list_keywords[0]));

  // Find the reference keywords within the command line + config file -> store
  // them into a vector
  keymap key_analysed = analyse_keywords(argc, argv, list_keywords, nb_ref_key);

  // For each keyword from the reference list, fill the corresponding variable

  // Configuration file
  string config = key_analysed["c"].value;
  // Repository in which the filters are stored
  string fltRep =
      ((key_analysed["FILTER_REP"]).split_string(lepharedir + "/filt/", 1))[0];
  // filter list with the name of the filter file to be oppened
  vector<string> fltFiles =
      (key_analysed["FILTER_LIST"]).split_string("flt.pb", -99);
  // Use the number of filters as the reference number of keyword (the following
  // list of values should match 1 or this number)
  int ref_size = fltFiles.size();
  // Transmission in energy or photons
  vector<int> transtyp = (key_analysed["TRANS_TYPE"]).split<int>("0", ref_size);
  // calibration depending on the instrument
  vector<int> calibtyp =
      (key_analysed["FILTER_CALIB"]).split<int>("0", ref_size);
  // Output file
  string outputName =
      (key_analysed["FILTER_FILE"]).split<string>("filters", 1)[0];
  string filtfile = lepharework + "/filt/" + outputName + ".dat";
  string filtdoc = lepharework + "/filt/" + outputName + ".doc";

  /*
  INFO PARAMETERS ON SCREEN
  */

  cout << "#######################################" << endl;
  cout << "# Build the filter file with the following options: " << endl;
  cout << "# Config file: " << config << endl;
  cout << "# FILTER_REP: " << fltRep << endl;
  cout << "# FILTER_LIST: ";
  for (int k = 0; k < ref_size; k++) {
    cout << fltFiles[k] << " ";
  }
  cout << endl << "# TRANS_TYPE: ";
  for (int k = 0; k < ref_size; k++) {
    cout << transtyp[k] << " ";
  }
  cout << endl << "# FILTER_CALIB: ";
  for (int k = 0; k < ref_size; k++) {
    cout << calibtyp[k] << " ";
  }
  cout << endl << "# FILTER_FILE: " << filtfile << endl;
  cout << "# FILTER_FILE.doc: " << filtdoc << endl;
  cout << "#######################################" << endl;

  /*
  READ THE FILTER DATA FILE
  */

  // Loop over each filter from the FILTER_LIST
  for (int k = 0; k < ref_size; k++) {
    // Open this filter file
    string fltFile = fltRep + "/" + fltFiles[k];
    // Generate one object "flt"
    flt oneFilt(k, fltFile, transtyp[k], calibtyp[k]);
    // store all filters in a vector
    vecFlt.push_back(oneFilt);
  }

  /*
  WRITE THE FILTER FILE
  */

  // function to write the filters
  write_output_filter(filtfile, filtdoc, vecFlt);

  /*
  DONE
  */

  // display the time needed to run the code
  t2 = clock();
  float diff(((float)t2 - (float)t1) / CLOCKS_PER_SEC);
  cout << "Times to run the code in sec: " << diff << endl;

  return 0;
}
