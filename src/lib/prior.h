/*

  24/11/2023
  Class to deploy a prior which modifies chi2 as a function of z and model
  properties.

*/

// avoid multiple def of the same class
#ifndef PRIOR_H  // check that this keyword has been set already
#define PRIOR_H  // define the keyword to be checked

#include <array>
#include <functional>
#include <stdexcept>
#include <vector>

using namespace std;

// Use forward declaration and overwrite these empyy declarations
class SED;
class onesource;

class prior
// The prior class should be initialised before the onesource fit method is
// called. We apply a number of generic prior forms typically based on redshift.
// We must initialise the prior before the fit to enable over writing in Python.
{
 private:
  // Variables
  int libtype, luv, lnir;
  double mag0;
  SED* lib;

 public:
  // Variables
  int apply_nz, apply_weights;
  double chi2, reds;
  vector<double> weights;
  // const onesource& source;
  array<int, 2> bp;
  // Constructor
  prior() {
    apply_nz = 1;  // Apply n of z if it is specified unless we switch it off
    apply_weights = 0;
    luv = 0;
    lnir = 0;
    chi2 = 1.e9;
    reds = 0.;
    bp = {0, 0};
  }
  // Methods
  double update_chi2(onesource* source, double chi2, SED* lib, int il,
                     double dm, double funz0, const array<int, 2>& bp,
                     bool mabsGALprior, bool mabsAGNprior);
  double absmag_prior(onesource* source, double chi2, double reds, int libtype,
                      double dm, double mag0, double funz0);
  double nz_prior(onesource* source, double chi2, const double luv,
                  const double lnir, const double reds, const array<int, 2> bp);

  // Python overwriteable prior which takes the fullLib of SEDs and the source
  // and returns a vector of weights equal in length to the fullLib.
  // In general it is safer to simply overwrite the weights directly.

  // General function allowing pybind overwriting
  std::function<vector<double>(vector<SED*>, onesource*)> weights_function =
      [](vector<SED*> fullLib, onesource* source) {
        // Default weight behavior all set to 1.
        vector<double> onesVector(fullLib.size(), 1.0f);
        return onesVector;
      };
  // Default behavior of the method
  int set_weights(vector<SED*> fullLib, onesource* source) {
    // Call the current chi2_function
    weights = weights_function(fullLib, source);
    if (fullLib.size() != weights.size()) {
      throw std::length_error("weights and fullLib are not the same length.");
    }
    return 0;
  };
  // Allow setting a new weights_function from Python
  void set_weights_function(
      std::function<vector<double>(vector<SED*>, onesource*)> func) {
    weights_function = func;
  };
};

#endif