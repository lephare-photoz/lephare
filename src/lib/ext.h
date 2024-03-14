/*
 15/12/2014
 Class to store one extinction law
*/

// avoid multiple def of the same class
#ifndef EXT_H  // check that this keyword has been set already
#define EXT_H  // define the keyword to be checked

#include <string>
#include <vector>

#include "oneElLambda.h"

using std::string;
using std::vector;

/*! \brief Class Extinction to store the lambda/value vector from the extinction
 * law read
 *
 * The Mag class is in charge of creating and filling the ext instances for each
 * extinction laws to be read from files provided in the config.
 */
class ext {
 public:
  vector<oneElLambda> lamb_ext;  ///< vector of struct oneElLambda
  string name;                   ///< name of the extinction file
  double lmin;                   ///< min lambda value read from the file
  double lmax;                   ///< max lambda value read from the file
  int numext;                    ///< id of this extinction law

  /// minimal constructor of the ext class, with its name and the id of the
  /// extinction law
  ext(const string nameC, int numextC = 0) {
    name = nameC;
    numext = numextC;
  }

  /// read an extinction file, sort the resulting vector according to ascending
  /// lambda, and update class members
  void read(string extFile);

  /// add a single element
  void add_element(double lam, double val, double ori);
};

#endif
