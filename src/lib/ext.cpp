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
    oneElLambda pos = new_lamb_flt[i];
    oneElLambda next = new_lamb_flt[i + 1];
    double delta = next.lamb - pos.lamb;
    double mid_flt = (pos.val + next.val) / 2.;
    double mid_ext = (new_lamb_ext[i].val + new_lamb_ext[i + 1].val) / 2.;
    fint += mid_flt * delta;
    // Integral of the transmission by the filter x extinction
    aint += mid_flt * mid_ext * delta;
  }

  // clean
  lamb_all.clear();
  new_lamb_flt.clear();
  new_lamb_ext.clear();

  return (aint /= fint);
}
