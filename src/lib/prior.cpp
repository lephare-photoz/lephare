/*

  24/11/2023
  Implementation of the functions of the prior class.

  This should allow the current absolute magnitude and N(z) priors to be
  abstracted outside of the onesource class.

  It should be designed in such a way as to allow overwriting via Python or c++
  with a general prior.

  Multiple priors should be possible. In the first instance we implement
  absolute mag and N(z) priors.

*/

#include "prior.h"

#include <functional>
#include <iostream>
#include <vector>

#include "SED.h"
#include "globals.h"
#include "onesource.h"

// #include <pybind11/pybind11.h>
// namespace py = pybind11;

using namespace std;

// The update_chi2 method takes the chi2 and modifies it as a function of the
// SED properties. It must not take source properties which do not change during
// fit.
double prior::update_chi2(onesource* source, double chi2, SED* lib, int il,
                          double dm, double funz0, const array<int, 2>& bp,
                          bool mabsGALprior, bool mabsAGNprior) {
  // Get all required values from the SED
  reds = (lib)->red;
  libtype = (lib)->nlib;
  luv = (lib)->luv;
  lnir = (lib)->lnir;
  mag0 = (lib)->mag0;

  // Apply all priors
  // nz prior
  if (bp[0] >= 0 && libtype == 0 && apply_nz) {
    chi2 = nz_prior(source, chi2, luv, lnir, reds, bp);
  }
  // Abs Magnitude from model @ z=0 for rejection for galaxies and AGN
  if ((mabsGALprior && libtype == 0) || (mabsAGNprior && libtype == 1)) {
    chi2 = absmag_prior(source, chi2, reds, libtype, dm, mag0, funz0);
  }
  // General weight
  if (apply_weights == 1) {
    chi2 = chi2 - 2 * log(weights[il]);
  }
  // Return final updated chi2 value
  return chi2;
}

// The absolute magnitude prior
double prior::absmag_prior(onesource* source, double chi2, double reds,
                           int libtype, double dm, double mag0, double funz0) {
  // predicted magnitudes within the babs filter renormalized by the scaling
  double abs_mag;
  abs_mag = mag0 - 2.5 * log10(dm);
  // specific case when the redshift is 0
  if (reds < 1.e-10) abs_mag = abs_mag - funz0;
  // Galaxy rejection
  if ((abs_mag <= source->priorLib[0] || abs_mag >= source->priorLib[1]) &&
      libtype == 0)
    chi2 = 1.e9;
  // AGN rejection
  if ((abs_mag <= source->priorLib[2] || abs_mag >= source->priorLib[3]) &&
      libtype == 1)
    chi2 = 1.e9;

  return chi2;
}

/*
  Use prior info on N(z) as a function of magnitude (Iband) and of the type.

  For now we simply try to make this a method of the prior class which somehow
  gets access to all the variables in the onesource class
*/
double prior::nz_prior(onesource* source, double chi2, const double luv,
                       const double lnir, const double reds,
                       const array<int, 2> bp) {
  double val = 1.;
  //   return val;
  //
  // }
  int mod;
  double zot = 0., kt = 0., alpt0 = 0.;

  // Define the NUV-R rest-frame color corrected for dust-extinction
  double nuvk = -2.5 * (luv - lnir);

  // Apparent magnitude in i band, if possible
  // Apparent magnitude in i band, if possible
  double mi = -99.;
  if (source->busnorma[bp[0]] == 1) {
    mi = flux2mag(source->ab[bp[0]]);
  } else if (source->busnorma[bp[1]] == 1) {
    mi = flux2mag(source->ab[bp[1]]);
  } else
    return val;

  // Bright case I<20, no solution above 1, don't compute the prior
  if (mi < 20.) {
    if (reds > 1) val = 0;
    return val;
  }
  // Bright case I<22, no solution above 2, but compute the prior
  if (mi < 22. && reds > 2) {
    val = 0;
    return val;
  }

  // Set up the parameters to define the redshift distribution
  if (nuvk > 4.25) {
    // Case E/S0
    // Color UV-K of PHOTO_230506/El_cww.sed.resample.new.resample15.inter
    mod = 0;
    zot = 0.45181;
    kt = 0.13677;
    alpt0 = 3.33078;
  } else if (nuvk > 3.19 && nuvk < 4.25) {
    // Case Sbc
    // Color UV-K of PHOTO_230506/Sbc_cww.sed.resample.new.resample8.inter
    mod = 1;
    zot = 0.16560;
    kt = 0.12983;
    alpt0 = 1.42815;
  } else if (nuvk > 1.9 && nuvk < 3.19) {
    // Case Scd. Color UV-K of
    // PHOTO_230506/Scd_cww.sed.resample.new.resample7.inter -19.4878 + 21.1501
    mod = 2;
    zot = 0.21072;
    kt = 0.14008;
    alpt0 = 1.58310;
  } else {
    // Case Irr
    mod = 3;
    zot = 0.20418;
    kt = 0.13773;
    alpt0 = 1.34500;
  }

  // P(z|T,m0)
  double zmax = zot + kt * (mi - 20);
  double pz = pow(reds, alpt0) * exp(-pow((reds / zmax), alpt0));

  // P(T|m0)
  double ktf[4] = {0.47165, 0.30663, 0.12715, -0.34437};
  double ft[4] = {0.43199, 0.07995, 0.31162, 0.21220};

  // Ratio for each type
  double rappSum =
      ft[0] * exp(-ktf[0] * (mi - 20.0)) + ft[1] * exp(-ktf[1] * (mi - 20.0));
  rappSum +=
      ft[2] * exp(-ktf[2] * (mi - 20.0)) + ft[3] * exp(-ktf[3] * (mi - 20.0));
  double rapp = ft[mod] * exp(-ktf[mod] * (mi - 20.0));

  // Normalisation of the probability function
  // pcal=exp(gammln(1.d0/alpt0(mod)+1.d0))
  double pcal;
  if (mod == 0) pcal = 0.89744;
  if (mod == 1) pcal = 0.90868;
  if (mod == 2) pcal = 0.89747;
  if (mod == 3) pcal = 0.91760;
  pcal = pow(zmax, alpt0 + 1) / alpt0 * pcal;

  // Final value
  val = pz / pcal * rapp / rappSum;

  return chi2 - 2 * log(val);
}
