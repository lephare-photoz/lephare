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

/*! \brief Extragalactic Opacity
 *
 * LePHARE can correct for extragalactic opacity based on models that are
 * stored in <a
 * href="https://github.com/lephare-photoz/lephare-data/opa">lephare-data/opa</a>.
 */
class opa {
 private:
  string opaFile;

 public:
  vector<oneElLambda> lamb_opa;  ///< tabulated (lambda, transmission) curve
                                 ///< of the IGM opacity, read from #opaFile
  double lmin,                   ///< minimum wavelength of #lamb_opa
      lmax;                      ///< maximum wavelength of #lamb_opa
  double red;                    ///< redshift this opacity curve applies to

  /// minimal constructor of the opa class, with its redshift and the name
  /// of the opacity model file
  /// @param redC: redshift, stored in #red
  /// @param opaFileC: name of the opacity file (relative to
  /// $LEPHAREDIR/opa/), stored in #opaFile
  opa(const double redC, const string opaFileC) {
    opaFile = opaFileC;  // name of the file
    red = redC;          // corresponding redshift
  }

  /// Read the (lambda, transmission) opacity curve from #opaFile into
  /// #lamb_opa, and set #lmin/#lmax accordingly
  void read();
};

#endif
