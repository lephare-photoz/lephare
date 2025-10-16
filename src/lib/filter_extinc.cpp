/*
    Compute the mean extinction for a given set of filters
*/

#include <getopt.h>  // get option line commands getopt_long
#include <math.h>    /* pow */
#include <stdlib.h>  // absolute value, exit

#include <cstring>   // Use string in the old c format
#include <fstream>   // print output file
#include <iomanip>   // std::set precision, et al.
#include <iostream>  // print standard file
#include <sstream>
#include <string>  // use string
#include <vector>  // manipulate vector

// Le Phare
#include "SED.h"
#include "ext.h"
#include "flt.h"      // define the filters
#include "globals.h"  // global variables
#include "keyword.h"  // to read the keyword in the command line/para file
#include "mag.h"

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
  string list_keywords[] = {"FILTER_FILE", "EXT_CURVE", "GAL_CURVE", "OUTPUT"};
  // Number of keywords
  int nb_ref_key = (int)(sizeof(list_keywords) / sizeof(list_keywords[0]));

  // Find the reference keywords within the command line + config file -> store
  // them into a vector
  keymap key_analysed = analyse_keywords(argc, argv, list_keywords, nb_ref_key);

  // For each keyword from the reference list, fill the corresponding variable

  // Filter file
  string filters =
      (key_analysed["FILTER_FILE"]).split<string>("filters.dat", 1)[0];
  // Atmospheric extinction curve
  string atmec =
      (key_analysed["EXT_CURVE"]).split<string>("extinc_etc.dat", 1)[0];
  // Galactic extinction curve
  string galec = (key_analysed["GAL_CURVE"]).split<string>("CARDELLI", 1)[0];
  // Output file
  string output =
      (key_analysed["OUTPUT"]).split<string>("filter_extinc.dat", 1)[0];

  // Atmospheric extinction curve file path
  string atmec_file = lepharedir + "/ext/" + atmec;
  // Galactic extinction curve file path
  string galec_file = lepharedir + "/ext/" + galec;
  // Filter file path
  string filtfile = lepharework + "/filt/" + filters;

  /*
  INFO PARAMETERS ON SCREEN
  */

  cout << "#######################################" << endl;
  cout << "# Computing ATMOSPHERIC AND GALACTIC EXTINCTION " << endl;
  cout << "# with the following options: " << endl;
  cout << "# Filters: " << filters << endl;
  cout << "# Atmospheric extinction curve: " << atmec << endl;
  cout << "# Galactic extinction curve: " << galec << endl;
  cout << "# Output file: " << output << endl;
  cout << "#######################################" << endl;
  cout << " Filters Ext(mag/airmass) Albda/Av Albda/E(B-V) " << endl;

  // read atmospheric extinction data
  ext Atmospheric_Ext(atmec, 0);
  Atmospheric_Ext.read(atmec_file);

  // read galactic extinction data
  ext Galactic_Ext(galec, 1);
  if (galec.compare("CARDELLI") != 0) Galactic_Ext.read(galec_file);

  // The input filter file
  string fltFile = lepharework + "/filt/" + filters;
  ifstream sfiltIn;
  sfiltIn.open(fltFile.c_str());
  // Check if file is opened
  if (!sfiltIn) {
    throw invalid_argument("Can't open file " + fltFile);
  }
  vector<flt> allFlt;
  allFlt = read_flt(sfiltIn);

  // check Rv to be applied
  double rv = 3.1;  // Use by default Cardelli
  if (galec.compare("SMC_prevot") == 0) rv = 2.72;
  if ((galec.compare("SB_calzetti") == 0) || (galec.compare("calzetti") == 0))
    rv = 4.05;
  cout << " assuming Rv=" << rv << " for this Extinction law " << galec << endl;

  // Compute galactic extinction
  vector<double> albd, albdav, aint;
  double inter;
  // Loop over all filters
  for (int k = 0; k < int(allFlt.size()); k++) {
    // Compute atmospheric extinction
    if (atmec.compare("NONE") != 0)
      aint.push_back(compute_filter_extinction(allFlt[k], Atmospheric_Ext));
    else
      aint.push_back(99.);

    // If not Cardelli
    if (galec.compare("CARDELLI") != 0) {
      // galactic curves given in k(lbda) (=A(lbda)/E(B-V))
      //  -> A(lbda)/Av = A(lbda)/E(B-V) / Rv)
      //  Rv=3.1 except for Calzetti law (4.05) and SMC Prevot (2.72)
      inter = compute_filter_extinction(allFlt[k], Galactic_Ext);
      albd.push_back(inter);
      albdav.push_back(inter / rv);

    } else {
      // If cardelli (hardcoded)
      // output A(lbd)/Av->A(lbd)/(Rv*E(B-V))->A(lbd)/E(B-V)=Rv*A(lbd)/Av
      inter = cardelli_ext(allFlt[k]);
      albdav.push_back(inter);
      albd.push_back(inter * 3.1);
    }
  }

  /*
 WRITE THE FILTER FILE
 */

  // create a stream with the two output files
  ofstream stextoutput;
  stextoutput.open(output.c_str());
  // Check if file is opened
  if (!stextoutput) {
    cerr << "Can't open file " << output << endl;
    exit(1);
  }
  stextoutput << "#######################################" << endl;
  stextoutput << "# Computing ATMOSPHERIC AND GALACTIC EXTINCTION " << endl;
  stextoutput << "# with the following options: " << endl;
  stextoutput << "# Filters: " << filters << endl;
  stextoutput << "# Atmospheric extinction curve: " << atmec << endl;
  stextoutput << "# Galactic extinction curve: " << galec << endl;
  stextoutput << "# Output file: " << output << endl;
  stextoutput << "#######################################" << endl;
  stextoutput << " Filters Ext(mag/airmass) Albda/Av Albda/E(B-V) " << endl;
  for (int k = 0; k < int(allFlt.size()); k++) {
    string shortName =
        (allFlt[k].name).substr(((allFlt[k].name).rfind("/") + 1));
    cout << setw(30) << shortName << " " << setw(12) << aint[k] << " "
         << setw(12) << albdav[k] << " " << setw(12) << albd[k] << endl;
  }

  /*
  DONE
  */

  // display the time needed to run the code
  t2 = clock();
  float diff(((float)t2 - (float)t1) / CLOCKS_PER_SEC);
  cout << "Times to run the code in sec: " << diff << endl;

  return 0;
}
