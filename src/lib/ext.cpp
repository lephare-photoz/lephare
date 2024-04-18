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
    throw invalid_argument("Can't open file " + extFile);
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
