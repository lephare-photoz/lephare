/*
 17/12/2014
 Class to store the opacity
*/

// avoid multiple def of the same class
#ifndef OPA_H  // check that this keyword has been set already
#define OPA_H  // define the keyword to be checked

#include <string>
#include <vector>

#include "oneElLambda.h"

using std::string;
using std::vector;

// Class filter
class opa {
 private:
  string opaFile;

 public:
  vector<oneElLambda> lamb_opa;
  double lmin, lmax;
  double red;

  // minimal constructor of the ext class, with its name and the id of the model
  opa(const double redC, const string opaFileC) {
    opaFile = opaFileC;  // name of the file
    red = redC;          // corresponding redshift
  }

  // Prototype
  void read();
};

#endif
