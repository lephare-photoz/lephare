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
      lamb_ext.emplace_back(l, v, 2);
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
void ext::add_element(double lam, double val, double ori) {
  // fill the lambda/trans values of the SED
  lamb_ext.emplace_back(lam, val, ori);

  // minimum and maximum lambda of the extinction curves
  lmin = (lamb_ext[0]).lamb;
  lmax = (lamb_ext[lamb_ext.size() - 1]).lamb;

  return;
}

double compute_filter_extinction(const flt &oneFlt, const ext &oneExt) {
  // work with the original lamb_flux
  vector<oneElLambda> lamb_all = oneFlt.lamb_trans;

  // Concatenate two vectors composed of "oneElLambda" including this spectra
  // and the one to be added
  lamb_all.insert(lamb_all.end(), oneExt.lamb_ext.begin(),
                  oneExt.lamb_ext.end());

  // Sort the vector in increasing lambda
  sort(lamb_all.begin(), lamb_all.end());

  // Resample the filter into a common lambda range
  vector<oneElLambda> new_lamb_flt;
  SED::resample(lamb_all, new_lamb_flt, 0, 0, 1.e50);
  // Resample the extinction into a common lambda range
  vector<oneElLambda> new_lamb_ext;
  SED::resample(lamb_all, new_lamb_ext, 2, 0, 1.e50);

  double fint = 0;
  double aint = 0;

  // integrate the extinction curve through the filter
  for (size_t i = 0; i < new_lamb_flt.size() - 1; i++) {
    // Integral of the transmission by the filter
    oneElLambda flt1 = new_lamb_flt[i];
    oneElLambda flt2 = new_lamb_flt[i + 1];
    oneElLambda ext1 = new_lamb_ext[i];
    oneElLambda ext2 = new_lamb_ext[i + 1];
    if (flt1.ori >= 0 && flt2.ori >= 0 && ext1.ori >= 0 && ext2.ori >= 0) {
      double delta = flt2.lamb - flt1.lamb;
      double mid_flt = (flt1.val + flt2.val) / 2.;
      double mid_ext = (ext1.val + ext2.val) / 2.;
      fint += mid_flt * delta;
      // Integral of the transmission by the filter x extinction
      aint += mid_flt * mid_ext * delta;
    }
  }

  // clean
  lamb_all.clear();
  new_lamb_flt.clear();
  new_lamb_ext.clear();

  return (aint /= fint);
}

// Function of the basis class which read all the filters
vector<flt> read_flt(ifstream &sfiltIn) {
  vector<flt> allFlt;
  string bid;
  int imag;
  // read the number of filter
  sfiltIn >> bid >> imag;

  // Loop over each filter
  for (int k = 0; k < imag; k++) {
    // Generate one object "flt" and read it
    flt oneFilt(k, sfiltIn, 0, 0);
    // Compute fcorr, useful for FIR filters
    oneFilt.fcorrec();
    // store all filters in a vector
    allFlt.push_back(oneFilt);
  }

  return allFlt;
}

// compute galactic extinction in the filter based on Cardelli et al., 1989, ApJ
// 345
double cardelli_ext(flt &oneFlt) {
  // Define the limits of this filter
  double lmin = oneFlt.lmin();
  double lmax = oneFlt.lmax();
  ext oneExt("CARDELLI", 2);

  double lextg, extg;

  // computes the galactic extinction
  double dlbd = (lmax - lmin) / 400.;
  for (size_t i = 0; i < 402; i++) {
    lextg = lmin + double(i - 1) * dlbd;
    extg = cardelli_law(lextg);
    oneExt.add_element(lextg, extg, 2);
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
