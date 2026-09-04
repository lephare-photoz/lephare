/*

  15/12/14
  Implementation of the functions of the ext class

*/

#include "ext.h"

#include <algorithm>
#include <cmath>
#include <fstream>   // print output file
#include <iostream>  // print standard file
#include <sstream>
#include <string>
#include <vector>

#include "SED.h"
#include "globals.h"
#include "oneElLambda.h"

using namespace std;

/*
 READ THE EXTINCTION LAW
*/
void ext::read(string extFile) {
  ifstream sext;
  string lit;

  // open the ascii file with the extinction law into a stream
  sext.open(extFile.c_str());
  // Check if file is opened
  if (!sext) {
    throw invalid_argument("Can't open file with the attenuation curve in " +
                           extFile);
  }

  // Take the stream line by line
  while (getline(sext, lit)) {
    // If the first character of the line is not #
    if (check_first_char(lit)) {
      // put the line into the stream ss again
      stringstream ss(lit);

      // fill the lambda/trans values of the SED
      double l, v;
      ss >> l;
      ss >> v;
      lamb_ext.emplace_back(l, v);
    }
  }

  // Close the stream
  sext.close();

  // Sort the vector according to ascending lambda
  sort(lamb_ext.begin(), lamb_ext.end());
  // minimum and maximum lambda of the extinction curves
  lmin = (lamb_ext[0]).lamb;
  lmax = (lamb_ext[lamb_ext.size() - 1]).lamb;

  return;
}

/*
Add a value to the extinction law
*/
void ext::add_element(double lam, double val) {
  // fill the lambda/trans values of the SED
  lamb_ext.emplace_back(lam, val);

  // minimum and maximum lambda of the extinction curves
  lmin = (lamb_ext[0]).lamb;
  lmax = (lamb_ext[lamb_ext.size() - 1]).lamb;

  return;
}

double compute_filter_extinction(const flt& oneFlt, const ext& oneExt) {
  // extinction curves required to be defined over the whole filter range
  if (oneExt.lamb_ext.front().lamb > oneFlt.lmin() ||
      oneExt.lamb_ext.back().lamb < oneFlt.lmax()) {
    return INVALID_VAL;
  }

  SED flatSED("flatSED", 0, "GAL");

  // Clear anything just in case
  flatSED.lamb_flux.clear();

  // Make the flat SED with a constant flux of 1.0 across the wavelength range
  // of interest
  double lmin = 0.0;
  double lmax = 200000.0;
  int N = 20000;

  for (int i = 0; i < N; ++i) {
    double l = lmin + (lmax - lmin) * i / (N - 1);
    flatSED.lamb_flux.emplace_back(l, 1.0);
  }

  // Ensure sorted order (probably already fine, but consistent with read())
  std::sort(flatSED.lamb_flux.begin(), flatSED.lamb_flux.end());

  // Compute integrals prior to extinction
  auto result_before_extinction = flatSED.integrateSED(oneFlt);

  // Apply extinction to the SED
  // we use an example pD/EBV value assuming linear relations as in Galametz
  // Appendix A
  double reference_ebv = 0.1;
  flatSED.apply_extinction(reference_ebv, oneExt);

  // Compute integrals afterextinction
  auto result_after_extinction = flatSED.integrateSED(oneFlt);

  return -2.5 *
         LOG10D(result_after_extinction[3] / result_before_extinction[3]) /
         reference_ebv;
  ;
}

// compute galactic extinction in the filter based on Cardelli et al., 1989, ApJ
// 345
double cardelli_ext(flt& oneFlt) {
  // Define the limits of this filter
  double lmin = oneFlt.lmin();
  double lmax = oneFlt.lmax();
  ext oneExt("CARDELLI", 2);

  double lextg, extg;

  // computes the galactic extinction
  double dlbd = (lmax - lmin) / 400.;
  // i needs to be an int else double(i-1) is not correctly cast
  for (int i = 0; i < 401; i++) {
    lextg = lmin + double(i) * dlbd;
    extg = cardelli_law(lextg);
    oneExt.add_element(lextg, extg);
  }

  return compute_filter_extinction(oneFlt, oneExt);
}

//  compute albd/av at a given lambda (A) for the Cardelli law
// value straight from The Astrophysical Journal, 345:245-256,1989
// https://articles.adsabs.harvard.edu/pdf/1989ApJ...345..245C
double cardelli_law(double lb) {
  double rv = 3.1;
  double x = 10000. / lb;
  double y = x - 1.82;
  double f1, f2, fa, fb;

  if (x <= 1.1) {
    f1 = 0.574 * pow(x, 1.61);
    f2 = -0.527 * pow(x, 1.61);
  } else if (x > 1.1 && x < 3.3) {
    f1 = 1 + 0.17699 * y - 0.50447 * y * y - 0.02427 * y * y * y +
         0.72085 * y * y * y * y;
    f1 = f1 + 0.01979 * pow(y, 5) - 0.77530 * pow(y, 6) + 0.32999 * pow(y, 7);
    f2 = 1.41338 * y + 2.28305 * y * y + 1.07233 * y * y * y;
    f2 = f2 - 5.38434 * pow(y, 4) - 0.62251 * pow(y, 5) + 5.30260 * pow(y, 6) -
         2.09002 * pow(y, 7);
  } else if (x >= 3.3 && x < 5.9) {
    f1 = 1.752 - 0.316 * x - 0.104 / ((x - 0.467) * (x - 0.467) + 0.341);
    f2 = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) * (x - 4.62) + 0.262);
  } else if (x >= 5.9 && x < 8) {
    fa = -0.04473 * (x - 5.9) * (x - 5.9) -
         0.009779 * (x - 5.9) * (x - 5.9) * (x - 5.9);
    fb = 0.2130 * (x - 5.9) * (x - 5.9) +
         0.1207 * (x - 5.9) * (x - 5.9) * (x - 5.9);
    f1 = 1.752 - 0.316 * x - 0.104 / ((x - 0.467) * (x - 0.467) + 0.341) + fa;
    f2 = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) * (x - 4.62) + 0.262) + fb;
  } else {
    f1 = -1.073 - 0.628 * (x - 8) + 0.137 * (x - 8) * (x - 8) -
         0.070 * (x - 8) * (x - 8) * (x - 8);
    f2 = 13.670 + 4.257 * (x - 8) - 0.420 * (x - 8) * (x - 8) +
         0.374 * (x - 8) * (x - 8) * (x - 8);
  }

  return (f1 + f2 / rv);
}
