
/*

    Compute the mean extinction for a given set of filters

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

// Declare prototypes
void get_lephare_env();
vector<keyword> analyse_keywords(int argc, char *argv[],
                                 const string *list_keywords, int nb_ref_key);
void write_output_filter(ofstream &stfiltfile, ofstream &stfiltdoc,
                         vector<flt> vecFlt);

int SED::nrec = 0;  // Even if SED used in the class FLT, still need to
                    // instantiate this variable (used in sedtolib)

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
  vector<keyword> key_analysed =
      analyse_keywords(argc, argv, list_keywords, nb_ref_key);

  // For each keyword from the reference list, fill the corresponding variable

  // Filter name
  string filter = key_analysed[0].value;
  // Atmospheric extinction curve
  string atmec = key_analysed[1].value;
  // Galactic extinction curve
  string galec = key_analysed[2].value;
  // Output file
  string output = key_analysed[3].value;

  // Atmospheric extinction curve file path
  string atmec_file = lepharedir + "/ext/" + atmec + ".dat";
  // Galactic extinction curve file path
  string galec_file = lepharedir + "/ext/" + galec + ".dat";
  // Filter file path
  string filtfile = lepharework + "/filt/" + filter + ".dat";

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
  ext Atmospheric_Ext(atmec);
  Atmospheric_Ext.read(atmec_file);

  // read galactic extinction data
  ext Galactic_Ext(galec);
  Galactic_Ext.read(galec_file);

  // Filter object
  flt Filt(filter);
  // read the filter
  Filt.read(filtfile);

  // check Rv to be applied

  double rv = 3.1;
  if (galec.compare("CARDELLI") != 0)
    if (galec.compare("SMC_prevot") == 0)
      rv = 2.72;
    else if ((galec.compare("SB_calz") == 0) || (galec.compare("calzet") == 0))
      rv = 4.05;

  cout << " assuming Rv=" << rv << " for this Extinction law " << endl;

  // Compute atmospheric extinction
  double aint = comp_atm(&Filt, &Atmospheric_Ext);

  // Compute galactic extinction
  double albd, albdav;
  if (galec.compare("CARDELLI") != 0) {
    albd = comp_atm(&Filt, &Atmospheric_Ext);
    if (galec.compare("SMC_prevot") == 0) {
      albdav = albd / 2.72;
      rv = 2.72;
    } else if ((galec.compare("SB_calz") == 0) ||
               (galec.compare("calzet") == 0)) {
      albdav = albd / 4.05;
      rv = 4.05;
    } else {
      albdav = albd / 3.1;
      rv = 3.1;
    }
  } else {
    albdav = comp_gal(&Filt);
    albd = albdav * 3.1;
    rv = 3.1;
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
  stextoutput << "# Filters: " << filter << endl;
  stextoutput << "# Atmospheric extinction curve: " << atmec << endl;
  stextoutput << "# Galactic extinction curve: " << galec << endl;
  stextoutput << "# Output file: " << output << endl;
  stextoutput << "#######################################" << endl;
  stextoutput << " Filters Ext(mag/airmass) Albda/Av Albda/E(B-V) " << endl;

  /*
  DONE
  */

  // display the time needed to run the code
  t2 = clock();
  float diff(((float)t2 - (float)t1) / CLOCKS_PER_SEC);
  cout << "Times to run the code in sec: " << diff << endl;

  return 0;
}

// compute atmospheric extinction
double comp_atm(flt *oneFlt, ext *oneExt) { return comp_ext(oneFlt, &oneExt); }

// compute galactic extinction based on Cardelli et al., 1989, ApJ 345
double comp_gal(flt *oneFlt) {
  double lmin = oneFlt->lmin();
  double lmax = oneFlt->lmax();
  ext oneExt("DUMMY", 2);

  // computes the galactic extinction
  double dlbd = (lmax - lmin) / 400.;

  double lextg, extg;

  for (size_t i = 0; i < 402; i++) {
    lextg = lmin + (i - 1) * dlbd;
    extg = ext_gal(lextg);
    oneExt.add(lextg, extg, 2);
  }

  return comp_ext(oneFlt, &oneExt);
}

// compute extinction
double comp_ext(flt *oneFlt, ext *oneExt) {
  // work with the original lamb_flux
  vector<oneElLambda> lamb_all = oneFlt->lamb_trans;

  // Concatenate two vectors composed of "oneElLambda" including this spectra
  // and the one to be added
  lamb_all.insert(lamb_all.end(), oneExt->lamb_ext.begin(),
                  oneExt->lamb_ext.end());

  // Sort the vector in increasing lambda
  sort(lamb_all.begin(), lamb_all.end());

  // Resample the filter into a common lambda range
  vector<oneElLambda> new_lamb_flt;
  resample(lamb_all, new_lamb_flt, 0, 0, 1.e50);
  // Resample the extinction into a common lambda range
  vector<oneElLambda> new_lamb_ext;
  resample(lamb_all, new_lamb_ext, 2, 0, 1.e50);

  double fint = 0;
  double aint = 0;

  // integrate the extinction curve through the filter
  for (size_t i = 0; i < new_lamb_flt.size() - 1; i++) {
    fint += (new_lamb_flt[i].lamb + new_lamb_flt[i + 1].lamb) / 2. *
            (new_lamb_ext[i + 1].val - new_lamb_ext[i].val);
    aint += (new_lamb_flt[i].lamb + new_lamb_flt[i + 1].lamb) *
            (new_lamb_flt[i + 1].val - new_lamb_flt[i].val) *
            (new_lamb_flt[i].val + new_lamb_flt[i + 1].val) / 4.;
  }

  // clean
  lamb_all.clear();
  new_lamb_flt.clear();
  new_lamb_ext.clear();

  return (aint /= fint);
}

//  compute albd/av at a given lambda (A)
double ext_gal(double lb) {
  double rv = 3.1;
  double x = 10000. / lb;
  double y = x - 1.82;
  double f1 = 0;
  double f2 = 0;

  if (x <= 1.1) {
    f1 = 0.574 * pow(x, 1.61);
    f2 = -0.527 * pow(x, 1.61);
    elseif(x > 1.1 && x < 3.3) {
      f1 = 1 + 0.17699 * y - 0.50447 * y * y - 0.02427 * y * y * y +
           0.72085 * y * y * y * y;
      f1 = f1 + 0.01979 * pow(y, 5) - 0.77530 * pow(y, 6) + 0.32999 * pow(y, 7);
      f2 = 1.41338 * y + 2.28305 * y * y + 1.07233 * y * y * y;
      f2 = f2 - 5.38434 * pow(y, 4) - 0.62251 * pow(y, 5) +
           5.30260 * pow(y, 6) - 2.09002 * pow(y, 7);
    }
    elseif(x >= 3.3 && x < 5.9) {
      f1 = 1.752 - 0.316 * x - 0.104 / ((x - 0.467) * (x - 0.467) + 0.341);
      f2 = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) * (x - 4.62) + 0.262);
    }
    elseif(x >= 5.9 && x < 8) {
      fa = -0.04473 * (x - 5.9) * (x - 5.9) -
           0.009779 * (x - 5.9) * (x - 5.9) * (x - 5.9);
      fb = 0.2130 * (x - 5.9) * (x - 5.9) +
           0.1207 * (x - 5.9) * (x - 5.9) * (x - 5.9);
      f1 = 1.752 - 0.316 * x - 0.104 / ((x - 0.467) * (x - 0.467) + 0.341) + fa;
      f2 = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) * (x - 4.62) + 0.262) + fb;
    }
    elseif(x >= 8) {
      f1 = -1.073 - 0.628 * (x - 8) + 0.137 * (x - 8) * (x - 8) -
           0.070 * (x - 8) * (x - 8) * (x - 8);
      f2 = 13.670 + 4.257 * (x - 8) - 0.420 * (x - 8) * (x - 8) +
           0.374 * (x - 8) * (x - 8) * (x - 8);
    }
  }

  return (f1 + f2 / rv);
}

/*
  resample the vector

  INPUT
   lamb_all= all elements concatenated (filter+SED)
   origine=  indicate how the vector needs to be filled (0 for filter, 1 for
  SED, 2 for attenuation laws, 3 opacity) lmin= lambda min to be considered
   lmax= lambda max to be considered

  OUTPUT
   lamb_interp= vector with a value in each lambda element, obtained with linear
  interpolation
*/
void resample(vector<oneElLambda> &lamb_all, vector<oneElLambda> &lamb_interp,
              const int origine, const double lmin, const double lmax) {
  // Initialize the previous and next element used for the interpolation
  oneElLambda prevEl(-999, -999, -999);
  oneElLambda nextEl(-999, -999, -999);

  // Loop over the full vector with the two vectors concatenated
  for (vector<oneElLambda>::iterator it = lamb_all.begin(); it < lamb_all.end();
       ++it) {
    // If well within the considered lambda range
    if ((it->lamb) >= lmin && (it->lamb) <= lmax) {
      // Add the new element
      lamb_interp.push_back(*it);
    }
  }

  // Loop over the full vector with the two vectors concatenated
  for (vector<oneElLambda>::iterator it = lamb_interp.begin();
       it < lamb_interp.end(); ++it) {
    // If the considered origin is the one we want to fill -> don't do anything
    if (it->ori == origine) {
      // Store this element for the next interpolation
      prevEl = *it;

      // If the considered origin is not the one we want to fill, replace it
      // with a linear interpolation of the others
    } else {
      // Check if the lower value is already defined, which is required for the
      // interpolation
      if (prevEl.lamb > 0) {
        // Increase the iterator to reach the next item with the right origin
        for (vector<oneElLambda>::iterator kt = it; kt < lamb_interp.end();
             ++kt) {
          // Store it as the upper value for the interpolation
          if (kt->ori == origine) {
            nextEl = *kt;
            break;
          }

          // Case when the final position in the vector is examined without
          // finding any item with the right origin Not possible to find an
          // upper value for the interpolation
          if (kt == lamb_interp.end()) nextEl.lamb = -999;
        }
      }

      // Change the vector as the result of the interpolation
      // Mofify the origin variable
      it->interp(prevEl, nextEl);
      it->ori = origine;
    }
  }
}
