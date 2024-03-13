/*
    26/11/2014
    Create a unique binary file with all the SED stored in the right format

*/

#include <getopt.h>  // get option line commands getopt_long

#include <cstring>   // Use string in the old c format
#include <ctime>     // date
#include <fstream>   // print output file
#include <iostream>  // print standard file
#include <sstream>
#include <string>  // use string
#include <vector>  // manipulate vector

// Le Phare
#include "SED.h"      // SED objects
#include "SEDLib.h"   // SED libraries
#include "globals.h"  // global variables
#include "keyword.h"  // to read the keyword in the command line/para file

using namespace std;

vector<SED> readSED(string sedFile, string sedFormat, int nummod, string type,
                    vector<double> ageSel);

template <class T>
void runner(SEDLib<T> &sedlib) {
  // INFO PARAMETERS ON SCREEN AND DOC
  sedlib.print_info();
  // READ THE MODEL LIST
  sedlib.read_model_list();
  // COMPUTE SED PROPERTIES AND WRITE THE BINARY
  sedlib.write_SED_lib();
  // FINALISE DOC
  time_t result = time(nullptr);
  sedlib.print_time_tofile(result);
}

// START
int main(int argc, char *argv[]) {
  // to estimate the times needed to run the code
  clock_t t1, t2;
  t1 = clock();

  /*
  ANALYSE KEYWORDS
  */

  // List of the keywords to be found in the config file/command line
  string list_keywords[] = {"t",        "GAL_SED",    "GAL_FSCALE",
                            "GAL_LIB",  "SEL_AGE",    "AGE_RANGE",
                            "QSO_SED",  "QSO_FSCALE", "QSO_LIB",
                            "STAR_SED", "STAR_LIB",   "STAR_FSCALE"};
  // Number of keywords
  int nb_ref_key = (int)(sizeof(list_keywords) / sizeof(list_keywords[0]));

  // Find the reference keywords within the command line + config file -> store
  // them into a vector
  keymap key_analysed = analyse_keywords(argc, argv, list_keywords, nb_ref_key);

  // Configuration file
  string config = key_analysed["c"].value;
  // type of source which is read (Galaxy G, QSO Q, Star S)
  string typ = key_analysed["t"].value;

  // GALAXY CASE
  if (typ[0] == 'G' || typ[0] == 'g') {
    GalSEDLib SEDLibrary(key_analysed, config, typ);
    runner(SEDLibrary);

    // QSO CASE
  } else if (typ[0] == 'Q' || typ[0] == 'q') {
    QSOSEDLib SEDLibrary(key_analysed, config, typ);
    runner(SEDLibrary);

    // STAR CASE
  } else if (typ[0] == 'S' || typ[0] == 's') {
    StarSEDLib SEDLibrary(key_analysed, config, typ);
    runner(SEDLibrary);

  } else {
    cout << "The type (-t) is not correctly indicated " << endl;
    exit(1);
  }

  /*
  DONE
  */

  // display the time needed to run the code
  t2 = clock();
  float diff(((float)t2 - (float)t1) / CLOCKS_PER_SEC);
  cout << endl << "Times to run the code in sec: " << diff << endl;

  return 0;
}
