/*

  17/12/14
  Implementation of the functions of the opacity class

*/

#include "opa.h"

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
 READ THE IGM OPACITY
*/
void opa::read() {
  ifstream sopa;
  string lit;

  // open the ascii file with the extinction law into a stream
  string name = lepharedir + "/opa/" + opaFile;
  sopa.open(name.c_str());

  // Take the stream line by line
  while (getline(sopa, lit)) {
    // If the first character of the line is not #
    if (check_first_char(lit)) {
      // put the line into the stream ss again
      stringstream ss(lit);

      // fill the lambda/trans values of the SED
      double l, v;
      ss >> l;
      ss >> v;
      // Skip 2 lines, with a sampling every 3A
      // getline( sopa, lit);
      // getline( sopa, lit);

      lamb_opa.emplace_back(l, v, 3);
    }
  }

  // Close the stream
  sopa.close();

  // Sort the vector according to ascending lambda
  sort(lamb_opa.begin(), lamb_opa.end());

  // minimum and maximum lambda of the IGM curves
  lmin = (lamb_opa[0]).lamb;
  lmax = (lamb_opa[lamb_opa.size() - 1]).lamb;

  return;
}
