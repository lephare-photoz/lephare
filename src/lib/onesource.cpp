/*

  05/05/2015
  Implementation of the functions of the onesource class

*/

#include "onesource.h"

#include <omp.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>   // print output file
#include <iomanip>   // std::set precision, et al.
#include <iostream>  // print standard file
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include "PDF.h"
#include "SED.h"
#include "globals.h"
#include "opa.h"

using namespace std;

void onesource::readsource(const string &identifier, const vector<double> vals,
                           const vector<double> err_vals,
                           const long context = 0, const double z_spec = -99.9,
                           const string additional_input = " ") {
  spec = identifier;
  if (ab.size() != sab.size()) {
    throw invalid_argument("vals and err_vals do not have the same dimension");
  }
  ab = vals;
  sab = err_vals;
  zs = z_spec;
  str_inp = additional_input;
}

/*
 Fixed the redshift considered for abs mag, etc depending on the option
*/
void onesource::considered_red(const bool zfix, const string methz) {
  // if zfix is not fixed
  if (!zfix) {
    if (methz[0] == 'M' || methz[0] == 'm') {
      // if considered redshift is the median of the PDF
      consiz = zgmed[0];
    } else {
      // if considered redshift is the minimum of the chi2
      consiz = zmin[0];
    }
  } else {
    // if considered redshift is the spectroscopic redshift
    consiz = zs;
  }
}

/*
Set the priors on the redshift min-max, E(B-V) min-max, mass min-max, abs mag
min-max+filter
*/
void onesource::setPriors(const array<double, 2> magabsB,
                          const array<double, 2> magabsF) {
  // Prior on the abs mag range for galaxies
  priorLib[0] = magabsB[0];
  priorLib[1] = magabsF[0];

  // Prior on the abs mag range for AGN
  priorLib[2] = magabsB[1];
  priorLib[3] = magabsF[1];
}

/*
 DEFINE THE FILTERS WHICH SHOULD BE USED IN THE FIT
*/
void onesource::fltUsed(const long gbcont, const long contforb,
                        const int imagm) {
  vector<int> bused;

  // Replace the context by the global context if defined
  if (gbcont >= 0) cont = gbcont;

  int nf = 0;
  nbul = 0;
  // Loop over each filter
  for (int k = 0; k < imagm; k++) {
    // Define if the band should be used based on the context
    if (cont <= 0) {
      // if context=0 or badly define, use all bands
      bused.push_back(1);
    } else {
      // define if the band should be used or not, based on the context
      bused.push_back(CHECK_CONTEXT_BIT(cont, k));
    }

    // reject band if forbiden context defined
    if (contforb > 0 && CHECK_CONTEXT_BIT(contforb, k) == 1) bused[k] = 0;

    // reject band if error and flux lower than 0
    if (sab[k] < 0 && ab[k] < 0) bused[k] = 0;

    // reject band if error=0. Not possible to include in the chi2
    if (sab[k] == 0) {
      bused[k] = 0;
      ab[k] = -99.9;
      sab[k] = -99.9;
    }

    // Initialise the mask used for normalisation at the same value than bused
    busnorma.push_back(bused[k]);
    // reject the band for the normalisation if upper-limit sab<0
    if (sab[k] < 0) busnorma[k] = 0;

    // Initialise the mask used to test upper-limits. No upper-limits first
    busul.push_back(0);
    // reject the band for the normalisation if upper-limit sab<0
    if (bused[k] == 1 && sab[k] < 0) busul[k] = 1;

    // Count the band which should be use in the fit and save the new context
    if (bused[k] == 1) {
      nbused++;
      new_cont += pow(2, k);
    }

    // Count the bands used in the normalisation
    if (busnorma[k] == 1) nf++;
    // Count the bands used as upper-limits
    if (busul[k] == 1) nbul++;
  }
  if (nf == 0) cout << "WARNING: No scaling --> No z " << spec << endl;
}

/*
 DEFINE THE FILTERS WHICH SHOULD BE USED IN THE FIR FIT
*/
void onesource::fltUsedIR(const long fir_cont, const long fir_scale,
                          const int imagm, vector<flt> allFilters,
                          const double fir_lmin) {
  // Loop over each filter
  for (int k = 0; k < imagm; k++) {
    // Define if the band should be used based on the context
    if (fir_cont <= 0) {
      // if context=0 or badly define, use all bands
      busfir.push_back(1);
    } else {
      // define if the band should be used or not, based on the context
      busfir.push_back(CHECK_CONTEXT_BIT(fir_cont, k));
    }

    // reject band if error and flux lower than 0
    if (sab[k] < 0 && ab[k] < 0) busfir[k] = 0;

    // reject the band below the define FIR minimum wavelength
    if (allFilters[k].lmean / (1 + consiz) < fir_lmin) busfir[k] = 0;

    // Bands used for the scaling
    bscfir.push_back(busfir[k]);  // same for scaling as the fit
    if (busfir[k] == 1 && fir_scale > 0)
      bscfir[k] = (CHECK_CONTEXT_BIT(fir_scale, k));
  }
}

/*
 CONVERT THE MAG INTO FLUX IF NECESSARY
*/
void onesource::convertFlux(const string &catmag,
                            const vector<flt> allFilters) {
  size_t imagm = allFilters.size();

  // If both magnitude and error are negative, flux/err are invalid
  for (size_t k = 0; k < allFilters.size(); k++) {
    if (ab[k] < 0 && sab[k] < 0) {
      ab[k] = INVALID_FLUX;
      sab[k] = INVALID_FLUX;
    }
  }

  // #pragma GCC ivdep forces the compiler to not check for aliasing, as it
  // cannot determine that there is none here
  if (catmag[0] == 'V') {  // Vega magnitudes
#pragma GCC ivdep
    for (size_t k = 0; k < allFilters.size(); k++) {
      if (ab[k] != INVALID_FLUX) ab[k] = mag2flux(ab[k], allFilters[k].ab);
    }
  } else {
#pragma GCC ivdep
    for (size_t k = 0; k < imagm; k++) {
      if (ab[k] != INVALID_FLUX) ab[k] = mag2flux(ab[k]);  // AB magnitudes
    }
  }
  // convert error from mag to flux
#pragma GCC ivdep
  for (size_t k = 0; k < imagm; k++) {
    if (ab[k] != INVALID_FLUX) sab[k] = ab[k] * sab[k] / 1.086;
  }
  return;
}

/*
 CONVERT THE FLUX INTO MAG
*/
void onesource::convertMag() {
  // Loop over each filter
  for (size_t k = 0; k < ab.size(); k++) {
    // Case with everything positive
    if (ab[k] > 0 && sab[k] >= 0) {
      msab.push_back(1.086 * sab[k] / ab[k]);  // Error on the mag
      mab.push_back(flux2mag(ab[k]));          // AB magnitudes

      // Case with an upper-limit
    } else if (ab[k] > 0 && sab[k] < 0) {
      msab.push_back(-1);              // Error on the mag
      mab.push_back(flux2mag(ab[k]));  // AB magnitudes

      // case with negative flux
    } else {
      msab.push_back(HIGH_MAG);  // Error on the mag
      mab.push_back(HIGH_MAG);   // AB magnitudes
    }
  }

  return;
}

/*
 KEEP ORIGINAL MAGNITUDE
*/
void onesource::keepOri() {
  ab_ori = ab;
  mab_ori = mab;
  return;
}

/*
 MODIFY THE OBSERVED MAGNITUDES WITH THE VALUE A0 and A1
*/
void onesource::adapt_mag(vector<double> a0, vector<double> a1) {
  double corr;
  // Loop over each filter
  for (size_t k = 0; k < ab.size(); k++) {
    // Define the correction to be applied to the observed magnitudes (same
    // convention as before)
    corr = a0[k] + a1[k] * 0.;
    // Observed magnitudes and flux with the correction
    if (ab_ori[k] > 0 || sab[k] > 0) ab[k] = ab_ori[k] * pow(10., 0.4 * corr);
    if (ab[k] > 0)
      mab[k] = flux2mag(ab[k]);
    else
      mab[k] = HIGH_MAG;
  }
  return;
}

/*
 RETURN THE INDEX OF THE LIBRARY TO BE CONSIDERED (MAINLY ZFIX CASE)
*/
vector<size_t> onesource::validLib(vector<SED *> &thelib, const bool &zfix,
                                   const double &consideredZ) {
  vector<size_t> val;
  // Condition with the redshift set ZFIX YES
  if (zfix) {
    for (size_t i = 0; i < thelib.size(); i++) {
      if ((thelib[i])->red == closest_red) val.push_back(i);
    }
  } else {
    // If not fixed redshift, use everything
    for (size_t i = 0; i < thelib.size(); i++) val.push_back(i);
  }

  return val;
}

/*
 SUBSTRACT THE STELLAR COMPONENT TO THE FIR OBSERVED FLUXES
*/
void onesource::substellar(const bool substar, vector<flt> allFilters) {
  double fluxMod;

  // Initalize with the observed fluxes and uncertainties
  for (size_t k = 0; k < ab.size(); k++) {
    abIR.push_back(ab[k]);
    sabIR.push_back(sab[k]);
  }

  // If the option to remove the stellar component is yes
  if (substar) {
    // Loop over each filter
    for (size_t k = 0; k < ab.size(); k++) {
      // Substract the predicted flux to the observed ones, only for
      // lambda_rf>25000 Check that the filter can be used for that
      if (busfir[k] == 1 && magm[k] > 0 && ab[k] > 0 &&
          allFilters[k].lmean / (1 + consiz) < 25.) {
        fluxMod = mag2flux(magm[k]);

        // Check that the observed flux is larger than the predicted one
        if (((ab[k] - fluxMod) > 0) && sab[k] > 0) {
          abIR[k] = ab[k] - fluxMod;
          sabIR[k] = sqrt(pow(sab[k], 2.) + pow(fluxMod, 2.));
        } else if (((ab[k] - fluxMod) <= 0) && sab[k] > 0) {
          abIR[k] = sab[k];
          sabIR[k] = sqrt(pow(sab[k], 2.) + pow(fluxMod, 2.));
        } else {
          busfir[k] = 0;
          bscfir[k] = 0;
        }
      }
    }
  }

  return;
}

/*
 RESCALE THE UNCERTAINTIES IN FLUX
*/
void onesource::rescale_flux_errors(const vector<double> min_err,
                                    const vector<double> fac_err) {
  size_t imagm = ab.size();

  // check if the number of filters in min_err corresponds to the expected one
  if ((size_t)imagm == min_err.size() || min_err.size() == 1) {
    for (size_t k = 0; k < imagm; k++) {
      if (sab[k] >= 0 && ab[k] != 0.) {
        sab[k] = 1.086 * sab[k] / ab[k];
        if (min_err.size() > 1) {
          sab[k] = sqrt(sab[k] * sab[k] + min_err[k] * min_err[k]);
        } else {
          sab[k] = sqrt(sab[k] * sab[k] + min_err[0] * min_err[0]);
        }
        // Put abs(ab[k]) if negative flux
        sab[k] = abs(ab[k]) * sab[k] / 1.086;
      }
    }
  } else {
    cout << " Can not add error in quadrature: dimension of min_err " << imagm
         << " " << min_err.size() << endl;
  }

  // multiply all the errors in flux by a common factor
  // check if the number of filters in fac_err corresponds to the expected one
  if ((size_t)imagm == fac_err.size() || fac_err.size() == 1) {
    for (size_t k = 0; k < imagm; k++) {
      if (fac_err.size() > 1) {
        if (sab[k] > 0) sab[k] = sab[k] * fac_err[k];
      } else {
        if (sab[k] > 0) sab[k] = sab[k] * fac_err[0];
      }
    }
  } else {
    cout << " Can not add error in quadrature: dimension of fac_err " << imagm
         << " " << fac_err.size() << endl;
  }

  return;
}

/*
 Principal function of the code with the fitting procedure
 */
void onesource::fit(vector<SED *> &fulllib, const vector<vector<double>> &flux,
                    const vector<size_t> &va, const double &funz0,
                    const array<int, 2> &bp) {
  int number_threads = 1, thread_id = 0;
  size_t imagm = ab.size();

  // Do a local minimisation per thread (store chi2 and index)
  // Catch first the number of threads
#ifdef _OPENMP
  number_threads = omp_get_max_threads();
#endif
  vector<vector<double>> locChi2(3, vector<double>(number_threads, 1.e9));
  vector<vector<int>> locInd(3, vector<int>(number_threads, -1));

  // Compute some quantities linked to ab and sab to save computational time in
  // the fit. invsab = busnorma / sab busnorma here ensures that all the
  // precomputed values will be 0 if busnorma=0 (aka not to be used in fit)
  vector<double> s2n(imagm, 0.0), invsab(imagm, 0.0), invsabSq(imagm, 0.0),
      abinvsabSq(imagm, 0.0);
  for (size_t i = 0; i < imagm; i++) {
    invsab[i] = busnorma[i] * 1.0 / sab[i];
  }
  for (size_t i = 0; i < imagm; i++) {
    invsabSq[i] = invsab[i] * invsab[i];
    s2n[i] = ab[i] * invsab[i];
  }
  for (size_t i = 0; i < imagm; i++) {
    abinvsabSq[i] = ab[i] * invsabSq[i];
  }

  // Already check the priors
  bool mabsGALprior = (priorLib[0] < 0 && priorLib[1] < 0);
  bool mabsAGNprior = (priorLib[2] < 0 && priorLib[3] < 0);

  // parrallellize over each SED
  vector<double> chi2loc(fulllib.size(), 1.e9), dmloc(fulllib.size(), -999.);
  size_t il;
#ifdef _OPENMP
  // double start = omp_get_wtime();
#pragma omp parallel shared(fulllib)                                    \
    firstprivate(s2n, invsab, invsabSq, abinvsabSq, imagm, nbul, busul, \
                     priorLib, number_threads, thread_id) private(il)
  {
    // Catch the name of the local thread in the parrallelisation
    thread_id = omp_get_thread_num();
#endif

    // Compute normalisation and chi2
#pragma omp for schedule(static, 10000)
    for (size_t i = 0; i < va.size(); i++) {
      // index to be considered because of ZFIX=YES
      il = va[i];
      // Initialize
      chi2loc[il] = 1.e9;
      dmloc[il] = -999.;

      // Measurement of scaling factor dm only with (fobs>flim), dchi2/ddm = 0
      double avmago = 0., avmagt = 0.;
      for (size_t k = 0; k < imagm; k++) {
        double fluxin = flux[il][k];
        avmago += fluxin * abinvsabSq[k];
        avmagt += fluxin * fluxin * invsabSq[k];
      }
      // Normalisation
      if (avmagt > 0) dmloc[il] = avmago / avmagt;

      // Measurement of chi^2
      chi2loc[il] = 0;
      for (size_t k = 0; k < imagm; k++) {
        double inter = s2n[k] - dmloc[il] * flux[il][k] * invsab[k];
        chi2loc[il] += inter * inter;
      }

      // Upper-limits. Check first if some bands have upper-limits, before
      // applying the condition
      if (nbul > 0) {
        for (size_t k = 0; k < imagm; k++) {
          if ((dmloc[il] * busul[k] * flux[il][k]) > ab[k] && busnorma[k] == 1)
            chi2loc[il] = 1.e9;
        }
      }

      // Model rejection based on prior
      // Abs mag rejection
      double reds;
      int libtype;
      reds = (fulllib[il])->red;
      libtype = (fulllib[il])->nlib;
      // Abs Magnitude from model @ z=0 for rejection for galaxies and AGN
      if ((mabsGALprior && libtype == 0) || (mabsAGNprior && libtype == 1)) {
        // predicted magnitudes within the babs filter renormalized by the
        // scaling
        double abs_mag;
        abs_mag = (fulllib[il])->mag0 - 2.5 * log10(dmloc[il]);
        // specific case when the redshift is 0
        if (reds < 1.e-10) abs_mag = abs_mag - funz0;
        // Galaxy rejection
        if ((abs_mag <= priorLib[0] || abs_mag >= priorLib[1]) && libtype == 0)
          chi2loc[il] = 1.e9;
        // AGN rejection
        if ((abs_mag <= priorLib[2] || abs_mag >= priorLib[3]) && libtype == 1)
          chi2loc[il] = 1.e9;
      }
      // Prior N(z)
      if (bp[0] >= 0 && libtype == 0) {
        double pweight =
            nzprior((fulllib[il])->luv, (fulllib[il])->lnir, reds, bp);
        chi2loc[il] = chi2loc[il] - 2. * log(pweight);
      }

      // keep the minimum chi2
      if (chi2loc[il] < 0.99e9) {
        int nlibloc = (fulllib[il])->nlib;
        // If local minimum inside the thread, store the new minimum for the
        // thread done per type s (0 gal, 1 AGN, 2 stars)
        if (locChi2[nlibloc][thread_id] > chi2loc[il]) {
          locChi2[nlibloc][thread_id] = chi2loc[il];
          locInd[nlibloc][thread_id] = (fulllib[il])->index;
        }
      }

      // Write the chi2
      fulllib[il]->chi2 = chi2loc[il];
      fulllib[il]->dm = dmloc[il];
    }

#ifdef _OPENMP
  }
  // double end = omp_get_wtime();
  // cout << endl << "Times to run the OpenMP fit in sec: " <<  end-start
  // <<endl;
#endif

  // Find the minimum chi2 for each redshift step by collapsing the threads
  // for the 3 types (0 gal, 1 AGN, 2 stars)
  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < number_threads; j++) {
      // Minimum over the full redshift range for the galaxies
      if (chimin[k] > locChi2[k][j]) {
        chimin[k] = locChi2[k][j];
        indmin[k] = locInd[k][j];
      }
    }
  }

  // Store the minimum redshift, dm, model for each type
  for (int k = 0; k < 3; k++) {
    if (indmin[k] >= 0) {
      zmin[k] = fulllib[indmin[k]]->red;
      dmmin[k] = fulllib[indmin[k]]->dm;
      imasmin[k] = fulllib[indmin[k]]->nummod;
    }
  }

  return;
}

void onesource::compute_best_fit_physical_quantities(vector<SED *> &fulllib) {
  // Best fit values for GAL physical parameters
  if (indmin[0] >= 0 && dmmin[0] > 0) {
    auto best_gal_sed = fulllib[indmin[0]];
    double tmp1 = dmmin[0];
    double tmp2 = tmp1 / Lsol;
    // Be sure that the mass is defined
    if (best_gal_sed->mass > 0) {
      results["MASS_BEST"] = LOG10D(best_gal_sed->mass * tmp1);
      results["SFR_BEST"] = LOG10D(best_gal_sed->sfr * tmp1);
      results["SSFR_BEST"] = LOG10D(best_gal_sed->sfr / best_gal_sed->mass);
      results["LDUST_BEST"] = best_gal_sed->ltir + LOG10D(tmp1);
    }
    double age = best_gal_sed->age;
    results["AGE_BEST"] = (age != INVALID_PHYS) ? LOG10D(age) : INVALID_PHYS;
    results["EBV_BEST"] = best_gal_sed->ebv;
    results["EXTLAW_BEST"] = best_gal_sed->extlawId;
    results["LUM_NUV_BEST"] =
        best_gal_sed->luv + LOG10D(3.e18 * 400 / pow(2300, 2) * tmp2);
    results["LUM_R_BEST"] =
        best_gal_sed->lopt + LOG10D(3.e18 * 1000 / pow(6000, 2) * tmp2);
    results["LUM_K_BEST"] =
        best_gal_sed->lnir + LOG10D(3.e18 * 2000 / pow(22000, 2) * tmp2);
  }
  return;
}

/*
  Use prior info on N(z) as a function of magnitude (Iband) and of the type.
*/
double onesource::nzprior(const double luv, const double lnir,
                          const double reds, const array<int, 2> bp) {
  double val = 1.;
  int mod;
  double zot = 0., kt = 0., alpt0 = 0.;

  // Define the NUV-R rest-frame color corrected for dust-extinction
  double nuvk = -2.5 * (luv - lnir);

  // Apparent magnitude in i band, if possible
  double mi = -99.;
  if (busnorma[bp[0]] == 1) {
    mi = flux2mag(ab[bp[0]]);
  } else if (busnorma[bp[1]] == 1) {
    mi = flux2mag(ab[bp[1]]);
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
    // Case Scd
    // Color UV-K of PHOTO_230506/Scd_cww.sed.resample.new.resample7.inter
    // -19.4878 + 21.1501
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
      ft[0] * exp(-ktf[0] * (mi - 20.0)) + ft[1] * exp(-ktf[1] * (mi - 20.0)) +
      ft[2] * exp(-ktf[2] * (mi - 20.0)) + ft[3] * exp(-ktf[3] * (mi - 20.0));
  double rapp = ft[mod] * exp(-ktf[mod] * (mi - 20.0));

  // Normalisation of the probability function
  // pcal=exp(gammln(1.0/alpt0(mod)+1.0))
  double pcal;
  if (mod == 0) pcal = 0.89744;
  if (mod == 1) pcal = 0.90868;
  if (mod == 2) pcal = 0.89747;
  if (mod == 3) pcal = 0.91760;
  pcal = pow(zmax, alpt0 + 1) / alpt0 * pcal;

  // Final value
  val = pz / pcal * rapp / rappSum;

  return val;
}

/*
 Check for discrepant band
 Test if some bands can be removed to improve the chi2, as long as the chi2
 remains above a threshold Done over all libraries
*/
void onesource::rm_discrepant(vector<SED *> &fulllib,
                              const vector<vector<double>> &flux,
                              const vector<size_t> &va, const double funz0,
                              const array<int, 2> bp, double thresholdChi2,
                              bool verbose) {
  size_t imagm = busnorma.size();
  double newmin, improvedChi2;
  // Start with the best chi2 among the libraries
  improvedChi2 = min({chimin[0], chimin[1], chimin[2]});
  double oldchi2 = improvedChi2;

  // Number of bands removed
  int nbRemoved = 0;
  // As long as the chi2 is above the treshold, and that the code didn't remove
  // more than 3 bands
  while (improvedChi2 > thresholdChi2 && nbRemoved < 2) {
    int flDis = -99;
    // Remove every bands one by one
    for (size_t k = 0; k < imagm; k++) {
      if (busnorma[k] == 1) {
        // Turn the filt off and redo the fit
        busnorma[k] = 0;
        this->fit(fulllib, flux, va, funz0, bp);
        newmin = min({chimin[0], chimin[1], chimin[2]});
        // if the chi2 has been improved, keep the best chi2, and keep the
        // disable band index
        if (improvedChi2 > newmin) {
          improvedChi2 = newmin;
          flDis = k;
        }
        busnorma[k] = 1;  // Turn the filt back on
      }
    }
    if (flDis >= 0) {
      // redo the fit without the disable band
      busnorma[flDis] = 0;
      this->fit(fulllib, flux, va, funz0, bp);
      newmin = min({chimin[0], chimin[1], chimin[2]});
      // One band has been removed in the fit
      nbused--;
      if (verbose)
        cout << "Source " << spec << " // Band " << flDis
             << " removed to improve the chi2, with old and new chi2 "
             << oldchi2 << " " << newmin << endl;
      nbRemoved++;
    } else
      break;
  }
  return;
}

/*
 Compute the chi2 for the IR library
*/
void onesource::fitIR(vector<SED *> &fulllibIR,
                      const vector<vector<double>> &fluxIR, const int imagm,
                      const string fit_frsc, cosmo lcdm) {
  int number_threads = 1;
// Do a local minimisation per thread (store chi2 and index)
// Catch first the number of threads
#ifdef _OPENMP
  number_threads = omp_get_max_threads();
#endif
  vector<double> locChi2(number_threads, 1.e9);
  vector<int> locInd(number_threads, -1);

  // Compute some quantities linked to ab and sab to save computational time in
  // the fit. No need to do it every step
  vector<double> s2n, invsab, invsabSq, abinvsabSq;
  for (int k = 0; k < imagm; k++) {
    s2n.push_back(abIR[k] / sabIR[k]);
    invsab.push_back(1. / sabIR[k]);
    invsabSq.push_back(1. / sabIR[k] / sabIR[k]);
    abinvsabSq.push_back(abIR[k] / sabIR[k] / sabIR[k]);
  }

// parrallellize over each SED
#ifdef _OPENMP
// double start = omp_get_wtime();
#pragma omp parallel firstprivate(lcdm) shared(fulllibIR)
  {
#endif

#pragma omp for schedule(static, 10000)
    // Loop over all SEDs from the library
    for (size_t i = 0; i < fulllibIR.size(); i++) {
      if (fulllibIR[i]->red == closest_red) {
        // Difference between the distance modulus at the redshift of the
        // library and the finel one considered
        double dmcor;
        dmcor = pow(10., (0.4 * (lcdm.distMod(fulllibIR[i]->red) -
                                 lcdm.distMod(consiz))));

        // normalization
        // Measurement of scaling factor dm only with (fobs>flim), dchi2/ddm = 0
        double avmago = 0.0, avmagt = 0.0, dmloc = -99.0;
        int nbusIR = 0;
        nbusIR = accumulate(bscfir.begin(), bscfir.end(), 0);
        for (int k = 0; k < imagm; k++) {
          double fluxin = fluxIR[i][k];
          avmago += fluxin * dmcor * bscfir[k] * abinvsabSq[k];
          avmagt += fluxin * fluxin * dmcor * dmcor * bscfir[k] * invsabSq[k];
        }
        // Check that the normalisation should be computed (if the scaling
        // should be free)
        if (nbusIR < 1) {
          cout << "WARNING: No scaling in IR " << spec << " " << nbusIR << " "
               << endl;
        } else {
          dmloc = avmago / avmagt;
        }
        // If one band or no free scale, use directly the scaling from the
        // template to find the best fit template
        double dmEff = dmloc;
        if (nbusIR <= 1 || fit_frsc[0] == 'n' || fit_frsc[0] == 'N')
          dmEff = 1.0;

        double chi2loc = 0;
        // Measurement of chi^2
        for (int k = 0; k < imagm; k++) {
          double inter = (s2n[k] - (dmEff * fluxIR[i][k] * dmcor * invsab[k]));
          chi2loc += busfir[k] * inter * inter;
        }
        // Associate results to the SED
        fulllibIR[i]->dm = dmEff;
        fulllibIR[i]->chi2 = chi2loc;

      } else {
        fulllibIR[i]->chi2 = 1.e9;
      }
    }

// Catch the name of the local thread in the parallelisation
#ifdef _OPENMP
    int thread_id = omp_get_thread_num();
#pragma omp for
#endif
    // Loop over all SEDs from the library
    for (size_t i = 0; i < fulllibIR.size(); i++) {
      // If local minimum inside the thread, store the new minimum for the
      // thread done per type s (0 gal, 1 AGN, 2 stars)
      if (locChi2[thread_id] > fulllibIR[i]->chi2) {
        locChi2[thread_id] = fulllibIR[i]->chi2;
        locInd[thread_id] = fulllibIR[i]->index;
      }
    }

#ifdef _OPENMP
  }
// double end = omp_get_wtime();
// cout << endl << "Times to run the OpenMP code in sec: " <<  end-start <<endl;
#endif

  // Find the minimum chi2 and corresponding index
  for (int j = 0; j < number_threads; j++) {
    if (chiminIR > locChi2[j]) {
      chiminIR = locChi2[j];
      indminIR = locInd[j];
    }
  }
  zminIR = fulllibIR[indminIR]->red;
  dmminIR = fulllibIR[indminIR]->dm;
  imasminIR = fulllibIR[indminIR]->nummod;

  if (indminIR >= 0)
    results["LUM_TIR_BEST"] = fulllibIR[indminIR]->ltir + log10(dmminIR);

  return;
}

/*
 Generate PDF marginalized over several parameters
 based on the chi2 stored in the SED class
*/
void onesource::generatePDF(vector<SED *> &fulllib, const vector<size_t> &va,
                            const vector<int> fltColRF, int fltREF,
                            const bool zfix) {
  // dimension of the redshift grid for the PDF
  int dimzg = pdfmap[9].size();
  // catch first the number of threads
  int number_threads = 1, thread_id = 0;
#ifdef _OPENMP
  number_threads = omp_get_max_threads();
#endif
  // Do a local minimisation per thread (store chi2 and index) dim 1: type, dim
  // 2: thread, 3: index of the redshift grid
  vector<vector<vector<double>>> locChi2(
      3, vector<vector<double>>(dimzg, vector<double>(number_threads, 1.e9)));
  vector<vector<vector<int>>> locInd(
      3, vector<vector<int>>(dimzg, vector<int>(number_threads, -1)));

  // to create the marginalized PDF
  int pos = 0;
  double col1 = -999., col2 = -999.;
  SED *rfSED;
  double prob;
  // need to convert into 1 dimension array for openMP reduction
  // 0:["MASS"] / 1:["SFR"] / 2:["SSFR"] / 3:["LDUST"] / 4:["LIR"] / 5:["AGE"] /
  // 6:["COL1"] / 7:["COL2"] / 8:["MREF"]/ 9:["MIN_ZG"] / 10:["MIN_ZQ"] /
  // 11:["BAY_ZG"] / 12:["BAY_ZQ"]
  double PDFzloc[pdfmap[11].size()];
  for (size_t i = 0; i < pdfmap[11].size(); ++i) PDFzloc[i] = 0.0;

  double PDFzqloc[pdfmap[12].size()];
  for (size_t i = 0; i < pdfmap[12].size(); ++i) PDFzqloc[i] = 0.0;

  double PDFmassloc[pdfmap[0].size()];
  for (size_t i = 0; i < pdfmap[0].size(); ++i) PDFmassloc[i] = 0.0;

  double PDFSFRloc[pdfmap[1].size()];
  for (size_t i = 0; i < pdfmap[1].size(); ++i) PDFSFRloc[i] = 0.0;

  double PDFsSFRloc[pdfmap[2].size()];
  for (size_t i = 0; i < pdfmap[2].size(); ++i) PDFsSFRloc[i] = 0.0;

  double PDFAgeloc[pdfmap[5].size()];
  for (size_t i = 0; i < pdfmap[5].size(); ++i) PDFAgeloc[i] = 0.0;

  double PDFLdustloc[pdfmap[3].size()];
  for (size_t i = 0; i < pdfmap[3].size(); ++i) PDFLdustloc[i] = 0.0;

  double PDFcol1loc[pdfmap[6].size()];
  for (size_t i = 0; i < pdfmap[6].size(); ++i) PDFcol1loc[i] = 0.0;

  double PDFcol2loc[pdfmap[7].size()];
  for (size_t i = 0; i < pdfmap[7].size(); ++i) PDFcol2loc[i] = 0.0;

  double PDFmrefloc[pdfmap[8].size()];
  for (size_t i = 0; i < pdfmap[8].size(); ++i) PDFmrefloc[i] = 0.0;

  // Decide if the uncertainties on the rest-frame colors should be analysed
  bool colAnalysis;
  colAnalysis = ((fltColRF[0] >= 0) && (fltColRF[1] >= 0) &&
                 (fltColRF[2] >= 0) && (fltColRF[3] >= 0) && (fltREF >= 0));

  // parrallellize over each SED
#ifdef _OPENMP
  // double start = omp_get_wtime();
#pragma omp parallel private(pos, col1, col2, rfSED, prob) firstprivate( \
        thread_id, dimzg, number_threads) shared(locChi2, locInd, fulllib)
  {
    // Catch the name of the local thread in the parallelisation
    thread_id = omp_get_thread_num();
#pragma omp for schedule(static, 10000) reduction(                           \
        + : PDFzloc, PDFzqloc, PDFmassloc, PDFSFRloc, PDFsSFRloc, PDFAgeloc) \
    nowait
#endif
    // Loop over all SEDs, which is parallelized
    for (size_t i = 0; i < va.size(); i++) {
      size_t il = va[i];

      // Check that the model has a defined probability
      if (fulllib[il]->chi2 < 0.99e9) {
        // nlib, index z=0
        int nlibloc = fulllib[il]->nlib;

        // Marginalization for the galaxies
        if (nlibloc == 0) {
          // probability exp(-chi2/2), but multiplied by a common factor
          // exp(-chi2_min/2) for the same object Since the the PDF is
          // normalized later, this factor vanishes. It allows to compute
          // probabities with chi2>>1000 (not feasible with double type)
          prob = exp(-0.5 * (fulllib[il]->chi2 - chimin[0]));
          // photo-z PDF of galaxies
          // Sum the proba to marginalize
          pos = pdfmap[11].index(fulllib[il]->red);
          PDFzloc[pos] += prob;

          // If able to determine a normalisation and get a the mass (assume
          // that the other are feasible in this case)
          double dmloc = fulllib[il]->dm;
          double massloc = fulllib[il]->mass;
          if (dmloc > 0 && massloc > 0) {
            // 0:["MASS"] / 1:["SFR"] / 2:["SSFR"] / 3:["LDUST"] / 4:["LIR"] /
            // 5:["AGE"] / 6:["COL1"] / 7:["COL2"] / 8:["MREF"]/ 9:["MIN_ZG"] /
            // 10:["MIN_ZQ"] / 11:["BAY_ZG"] / 12:["BAY_ZQ"] stellar mass PDF of
            // galaxies
            pos = pdfmap[0].index(LOG10D(dmloc * massloc));
            PDFmassloc[pos] += prob;

            // SFR PDF of galaxies
            pos = pdfmap[1].index(LOG10D(fulllib[il]->sfr * dmloc));
            PDFSFRloc[pos] += prob;

            // sSFR PDF of galaxies
            pos = pdfmap[2].index(LOG10D(fulllib[il]->ssfr));
            PDFsSFRloc[pos] += prob;

            // Age PDF of galaxies
            pos = pdfmap[5].index(LOG10D(fulllib[il]->age));
            PDFAgeloc[pos] += prob;

            // Ldust PDF of galaxies, ltir already in log
            if (fulllib[il]->ltir >= 0) {
              pos = pdfmap[3].index(fulllib[il]->ltir + log10(dmloc));
              PDFLdustloc[pos] += prob;
            }
          }

          // Only if we need an analysis of the colors
          if (colAnalysis) {
            // Retrieved the corresponding SED at z=0
            rfSED = fulllib[fulllib[il]->index_z0];

            // First rest-frame color
            col1 = rfSED->mag[fltColRF[0] - 1] - rfSED->mag[fltColRF[1] - 1];
            pos = pdfmap[6].index(col1);
            PDFcol1loc[pos] += prob;

            // Second rest-frame color
            col2 = rfSED->mag[fltColRF[2] - 1] - rfSED->mag[fltColRF[3] - 1];
            pos = pdfmap[7].index(col2);
            PDFcol2loc[pos] += prob;

            // Reference absolute magnitude
            pos = pdfmap[8].index(rfSED->mag[fltREF - 1] - 2.5 * log10(dmloc));
            PDFmrefloc[pos] += prob;
          }

          // marginalization for the QSO
        } else if (nlibloc == 1) {
          prob = exp(-0.5 * (fulllib[il]->chi2 - chimin[1]));

          // photo-z PDF of QSO
          pos = pdfmap[12].index(fulllib[il]->red);
          PDFzqloc[pos] += prob;
        }
      }
    }

    // case with free redshift
    if (!zfix) {
      // Generate the PDF(z) from the chi2 stored in each SED
      // based on the chi2 stored in the SED class
      // parallelised loop over the library
#pragma omp for schedule(static, 10000)
      for (size_t i = 0; i < va.size(); i++) {
        size_t il = va[i];
        // Index of the considered redshift into the PDF
        double chi2loc = fulllib[il]->chi2;
        if (chi2loc < 0.99e9) {
          // 11: BAYZG
          int poszloc = pdfmap[11].index(fulllib[il]->red);
          int nlibloc = fulllib[il]->nlib;
          int indloc = fulllib[il]->index;
          // If local minimum inside the thread, store the new minimum for the
          // thread
          if (locChi2[nlibloc][poszloc][thread_id] > chi2loc) {
            locChi2[nlibloc][poszloc][thread_id] = chi2loc;
            locInd[nlibloc][poszloc][thread_id] = indloc;
          }
        }
      }

      // Find the minimum for each redshift step by collapsing the threads
#pragma omp for
      for (int i = 0; i < dimzg; i++) {
        for (int j = 0; j < number_threads; j++) {
          // 0:["MASS"] / 1:["SFR"] / 2:["SSFR"] / 3:["LDUST"] / 4:["LIR"] /
          // 5:["AGE"] / 6:["COL1"] / 7:["COL2"] / 8:["MREF"]/ 9:["MIN_ZG"] /
          // 10:["MIN_ZQ"] / 11:["BAY_ZG"] / 11:["BAY_ZQ"] look for the new
          // minimum among the threads / Galaxies
          if (pdfmap[9].chi2[i] > locChi2[0][i][j]) {
            pdfmap[9].chi2[i] = locChi2[0][i][j];
            pdfmap[9].ind[i] = locInd[0][i][j];
          }
          // look for the new minimum among the threads / AGN
          if (pdfmap[10].chi2[i] > locChi2[1][i][j]) {
            pdfmap[10].chi2[i] = locChi2[1][i][j];
            pdfmap[10].ind[i] = locInd[1][i][j];
          }
        }
      }
    }

#ifdef _OPENMP
  }
  // double end = omp_get_wtime();
  // cout << endl << "Times to run the OpenMP code in sec: " <<  end-start
  // <<endl;
#endif

  // 0:["MASS"] / 1:["SFR"] / 2:["SSFR"] / 3:["LDUST"] / 4:["LIR"] / 5:["AGE"] /
  // 6:["COL1"] / 7:["COL2"] / 8:["MREF"]/ 9:["MIN_ZG"] / 10:["MIN_ZQ"] /
  // 11:["BAY_ZG"] / 12:["BAY_ZQ"] Put back 1 dimension array into PDF objects
  for (size_t k = 0; k < pdfmap[11].size(); k++)
    pdfmap[11].vPDF[k] = PDFzloc[k];
  for (size_t k = 0; k < pdfmap[12].size(); k++)
    pdfmap[12].vPDF[k] = PDFzqloc[k];
  for (size_t k = 0; k < pdfmap[0].size(); k++)
    pdfmap[0].vPDF[k] = PDFmassloc[k];
  for (size_t k = 0; k < pdfmap[1].size(); k++)
    pdfmap[1].vPDF[k] = PDFSFRloc[k];
  for (size_t k = 0; k < pdfmap[2].size(); k++)
    pdfmap[2].vPDF[k] = PDFsSFRloc[k];
  for (size_t k = 0; k < pdfmap[5].size(); k++)
    pdfmap[5].vPDF[k] = PDFAgeloc[k];
  for (size_t k = 0; k < pdfmap[3].size(); k++)
    pdfmap[3].vPDF[k] = PDFLdustloc[k];
  for (size_t k = 0; k < pdfmap[6].size(); k++)
    pdfmap[6].vPDF[k] = PDFcol1loc[k];
  for (size_t k = 0; k < pdfmap[7].size(); k++)
    pdfmap[7].vPDF[k] = PDFcol2loc[k];
  for (size_t k = 0; k < pdfmap[8].size(); k++)
    pdfmap[8].vPDF[k] = PDFcol2loc[k];

  // Normalize the PDF

  for (const auto &key : {11, 12, 0, 1, 2, 5, 3, 6, 7, 8})
    pdfmap[key].normalization();

  // Convert minimum chi2 at a given redshift into proba
  //  9:["MIN_ZG"] / 10:["MIN_ZQ"]
  pdfmap[9].chi2toPDF();
  pdfmap[10].chi2toPDF();

  return;
}

/*
 Generate PDF marginalized over LIR
 based on the chi2 stored in the SED class
*/
void onesource::generatePDF_IR(vector<SED *> &fulllib) {
  // 4:["LIR"]
  double PDFlirloc[pdfmap[4].size()];
  for (size_t i = 0; i < pdfmap[4].size(); ++i) PDFlirloc[i] = 0.0;

// parrallellize over each SED
#ifdef _OPENMP
// double start = omp_get_wtime();
#pragma omp parallel
  {
#pragma omp for reduction(+ : PDFlirloc)
#endif

    // Loop over all SEDs
    for (vector<SED *>::iterator it = fulllib.begin(); it < fulllib.end();
         ++it) {
      double prob = exp(-0.5 * ((*it)->chi2 - chiminIR));

      // Check that the model has a defined probability
      if ((*it)->chi2 >= 0 && (*it)->chi2 < 0.99e9) {
        // If able to determine a normalisation
        if ((*it)->dm > 0) {
          // LIR PDF of galaxies, ltir already in log
          if ((*it)->ltir > 0) {
            // Index of the considered library age into the PDF
            int pos = pdfmap[4].index((*it)->ltir + log10((*it)->dm));
            // Sum the proba to marginalize
            PDFlirloc[pos] += prob;
          }
        }
      }
    }

#ifdef _OPENMP
  }
// double end = omp_get_wtime();
// cout << endl << "Times to run the OpenMP code in sec: " <<  end-start <<endl;
#endif

  for (size_t k = 0; k < pdfmap[4].size(); k++)
    pdfmap[4].vPDF[k] = PDFlirloc[k];

  // Normalize the PDF
  pdfmap[4].normalization();

  return;
}

/*
Compute the mode and the associated uncertainties based on the marginalized
error
 */
void onesource::mode() {
  // First element of zgmode is the mode of the marginalized PDF
  // Use the parabolic interpolation for galaxies and QSO
  zgmode.push_back(pdfmap[11].int_parabL());
  zqmode.push_back(pdfmap[12].int_parabL());

  // Define the confidence levels
  double confLev[3] = {68, 90, 99};

  // Add the uncertainties, lower and upper for 68%, 90% and 99%
  pair<double, double> interv;
  for (size_t k = 0; k < 3; k++) {
    // For galaxies
    interv = pdfmap[11].credible_interval(confLev[k], zgmode[0]);
    zgmode.push_back(interv.first);
    zgmode.push_back(interv.second);
    // For QSO
    interv = pdfmap[12].credible_interval(confLev[k], zqmode[0]);
    zqmode.push_back(interv.first);
    zqmode.push_back(interv.second);
  }
}

/*
 write the output file
*/
void onesource::write_out(vector<SED *> &fulllib, vector<SED *> &fulllibIR,
                          ofstream &stout, vector<string> outkeywords) {
  // Loop over each keyword
  for (const string &outkey : outkeywords) {
    // If a match a found, write the value into the stream
    if (outkey == "IDENT")
      stout << setw(15) << std::fixed << setprecision(0) << spec << " ";

    // GALAXY
    stout.unsetf(ios::fixed);  // There is a problem with Z_BEST: for rare
                               // pathological galaxies is =-1.e9
    if (outkey == "Z_BEST")
      stout << setw(16) << setprecision(4) << zgmin[0] << " ";
    stout.setf(ios::fixed, ios::floatfield);
    stout.precision(4);  // set float format and the n. of decimal digits
    if (outkey == "Z_BEST68_LOW") stout << setw(12) << zgmin[1] << " ";
    if (outkey == "Z_BEST68_HIGH") stout << setw(12) << zgmin[2] << " ";
    if (outkey == "Z_BEST90_LOW") stout << setw(12) << zgmin[3] << " ";
    if (outkey == "Z_BEST90_HIGH") stout << setw(12) << zgmin[4] << " ";
    if (outkey == "Z_BEST99_LOW") stout << setw(12) << zgmin[5] << " ";
    if (outkey == "Z_BEST99_HIGH") stout << setw(12) << zgmin[6] << " ";

    // median of the PDF. Keep old naming possible.
    if (outkey == "Z_MED" || outkey == "Z_ML")
      stout << setw(12) << zgmed[0] << " ";
    if (outkey == "Z_MED68_LOW" || outkey == "Z_ML68_LOW")
      stout << setw(12) << zgmed[1] << " ";
    if (outkey == "Z_MED68_HIGH" || outkey == "Z_ML68_HIGH")
      stout << setw(12) << zgmed[2] << " ";
    if (outkey == "Z_MED90_LOW" || outkey == "Z_ML90_LOW")
      stout << setw(12) << zgmed[3] << " ";
    if (outkey == "Z_MED90_HIGH" || outkey == "Z_ML90_HIGH")
      stout << setw(12) << zgmed[4] << " ";
    if (outkey == "Z_MED99_LOW" || outkey == "Z_ML99_LOW")
      stout << setw(12) << zgmed[5] << " ";
    if (outkey == "Z_MED99_HIGH" || outkey == "Z_ML99_HIGH")
      stout << setw(12) << zgmed[6] << " ";

    if (outkey == "Z_MODE") stout << setw(12) << zgmode[0] << " ";
    if (outkey == "Z_MODE68_LOW") stout << setw(12) << zgmode[1] << " ";
    if (outkey == "Z_MODE68_HIGH") stout << setw(12) << zgmode[2] << " ";
    if (outkey == "Z_MODE90_LOW") stout << setw(12) << zgmode[3] << " ";
    if (outkey == "Z_MODE90_HIGH") stout << setw(12) << zgmode[4] << " ";
    if (outkey == "Z_MODE99_LOW") stout << setw(12) << zgmode[5] << " ";
    if (outkey == "Z_MODE99_HIGH") stout << setw(12) << zgmode[6] << " ";

    // Integer
    stout.precision(
        0);  //  in fixed format, number of decimal digits in the output
    if (outkey == "MOD_BEST")
      stout << setw(6) << std::fixed << imasmin[0] << " ";
    if (outkey == "EXTLAW_BEST")
      stout << setw(6) << std::fixed << results["EXTLAW_BEST"] << " ";
    // Float
    stout.setf(ios::fixed, ios::floatfield);
    for (const auto &key :
         {"MASS_BEST", "LDUST_BEST", "LUM_TIR_BEST", "EBV_BEST", "LUM_NUV_BEST",
          "LUM_R_BEST", "LUM_K_BEST"}) {
      if (outkey == key)
        stout << setw(12) << setprecision(3) << results[key] << " ";
    }
    // Scientific
    stout.setf(ios::scientific, ios::floatfield);
    // note: PDZ_BEST has to be fixed
    if (outkey == "PDZ_BEST")
      stout << setw(14) << setprecision(5) << 99. << " ";
    if (outkey == "CHI_BEST")
      stout << setw(14) << setprecision(5) << chimin[0] << " ";
    if (outkey == "SCALE_BEST")
      stout << setw(14) << setprecision(5) << dmmin[0] << " ";
    for (const auto &key : {"SFR_BEST", "SSFR_BEST", "AGE_BEST"}) {
      if (outkey == key)
        stout << setw(12) << setprecision(3) << results[key] << " ";
    }

    // GENERAL
    if (outkey == "STRING_INPUT")
      stout << " " << str_inp << " ";  // string-like

    stout.setf(ios::fixed);
    stout.precision(0);  // Integer:
    if (outkey == "CONTEXT") stout << setw(16) << cont << " ";
    if (outkey == "NBAND_USED") stout << setw(5) << nbused << " ";
    if (outkey == "NBAND_ULIM") stout << setw(3) << nbul << " ";

    stout.setf(ios::fixed, ios::floatfield);  // Float:
    if (outkey == "ZSPEC") stout << setw(14) << setprecision(5) << zs << " ";

    if (outkey == "MAG_OBS()") {
      for (const auto &m : mab)
        stout << setw(10) << setprecision(3) << m << " ";
    }
    if (outkey == "ERR_MAG_OBS()") {
      for (const auto &m : msab)
        stout << setw(10) << setprecision(3) << m << " ";
    }
    if (outkey == "MAG_MOD()") {
      for (const auto &m : magm)
        stout << setw(10) << setprecision(3) << m << " ";
    }

    // QSO
    stout.setf(ios::fixed, ios::floatfield);
    stout.precision(4);  // Float:

    if (outkey == "ZQ_BEST") stout << setw(12) << zqmin[0] << " ";
    if (outkey == "ZQ_BEST68_LOW") stout << setw(12) << zqmin[1] << " ";
    if (outkey == "ZQ_BEST68_HIGH") stout << setw(12) << zqmin[2] << " ";
    if (outkey == "ZQ_BEST90_LOW") stout << setw(12) << zqmin[3] << " ";
    if (outkey == "ZQ_BEST90_HIGH") stout << setw(12) << zqmin[4] << " ";
    if (outkey == "ZQ_BEST99_LOW") stout << setw(12) << zqmin[5] << " ";
    if (outkey == "ZQ_BEST99_HIGH") stout << setw(12) << zqmin[6] << " ";

    if (outkey == "ZQ_MED") stout << setw(12) << zqmed[0] << " ";
    if (outkey == "ZQ_MED68_LOW") stout << setw(12) << zqmed[1] << " ";
    if (outkey == "ZQ_MED68_HIGH") stout << setw(12) << zqmed[2] << " ";
    if (outkey == "ZQ_MED90_LOW") stout << setw(12) << zqmed[3] << " ";
    if (outkey == "ZQ_MED90_HIGH") stout << setw(12) << zqmed[4] << " ";
    if (outkey == "ZQ_MED99_LOW") stout << setw(12) << zqmed[5] << " ";
    if (outkey == "ZQ_MED99_HIGH") stout << setw(12) << zqmed[6] << " ";

    if (outkey == "ZQ_MODE") stout << setw(12) << zqmode[0] << " ";
    if (outkey == "ZQ_MODE68_LOW") stout << setw(12) << zqmode[1] << " ";
    if (outkey == "ZQ_MODE68_HIGH") stout << setw(12) << zqmode[2] << " ";
    if (outkey == "ZQ_MODE90_LOW") stout << setw(12) << zqmode[3] << " ";
    if (outkey == "ZQ_MODE90_HIGH") stout << setw(12) << zqmode[4] << " ";
    if (outkey == "ZQ_MODE99_LOW") stout << setw(12) << zqmode[5] << " ";
    if (outkey == "ZQ_MODE99_HIGH") stout << setw(12) << zqmode[6] << " ";
    stout.setf(ios::scientific, ios::floatfield);
    stout.precision(5);  // Scientific:
    if (outkey == "CHI_QSO") stout << setw(14) << chimin[1] << " ";
    stout.setf(ios::fixed);
    stout.precision(0);  // Integer:
    if (outkey == "MOD_QSO") stout << setw(12) << imasmin[1] << " ";

    // STARS
    stout.setf(ios::fixed);
    stout.precision(0);  // Integer:
    if (outkey == "MOD_STAR") stout << setw(12) << imasmin[2] << " ";
    stout.setf(ios::scientific, ios::floatfield);
    stout.precision(5);  // Scientific:
    if (outkey == "CHI_STAR") stout << setw(14) << chimin[2] << " ";

    // SECONDARY GALAXY PEAK
    stout.setf(ios::fixed);
    stout.precision(0);  // Integer:
    if (outkey == "MOD_SEC") stout << setw(12) << zsecMod << " ";
    if (outkey == "EXTLAW_SEC") stout << setw(12) << zsecExtlaw << " ";
    stout.setf(ios::fixed, ios::floatfield);
    stout.precision(4);  // Float:
    if (outkey == "Z_SEC") stout << setw(12) << zsec << " ";
    if (outkey == "EBV_SEC") stout << setw(12) << zsecEbv << " ";
    stout.setf(ios::scientific, ios::floatfield);
    stout.precision(5);  // Scientific:
    if (outkey == "PDZ_SEC") stout << setw(14) << zsecProb << " ";
    if (outkey == "CHI_SEC") stout << setw(14) << zsecChi2 << " ";
    if (outkey == "SCALE_SEC") stout << setw(14) << zsecScale << " ";
    if (outkey == "AGE_SEC") stout << setw(14) << zsecAge << " ";

    // PHYSICAL PARAMETERS
    stout.setf(ios::fixed, ios::floatfield);
    stout.precision(5);  // Float:
    if (outkey == "AGE_MED") stout << setw(12) << agemed[0] << " ";
    if (outkey == "AGE_INF") stout << setw(12) << agemed[1] << " ";
    if (outkey == "AGE_SUP") stout << setw(12) << agemed[2] << " ";

    if (outkey == "LDUST_MED") stout << setw(12) << Ldustmed[0] << " ";
    if (outkey == "LDUST_INF") stout << setw(12) << Ldustmed[1] << " ";
    if (outkey == "LDUST_SUP") stout << setw(12) << Ldustmed[2] << " ";

    if (outkey == "LUM_TIR_MED") stout << setw(12) << LIRmed[0] << " ";
    if (outkey == "LUM_TIR_INF") stout << setw(12) << LIRmed[1] << " ";
    if (outkey == "LUM_TIR_SUP") stout << setw(12) << LIRmed[2] << " ";

    if (outkey == "MASS_MED") stout << setw(12) << massmed[0] << " ";
    if (outkey == "MASS_INF") stout << setw(12) << massmed[1] << " ";
    if (outkey == "MASS_SUP") stout << setw(12) << massmed[2] << " ";

    if (outkey == "SFR_MED") stout << setw(12) << SFRmed[0] << " ";
    if (outkey == "SFR_INF") stout << setw(12) << SFRmed[1] << " ";
    if (outkey == "SFR_SUP") stout << setw(12) << SFRmed[2] << " ";

    if (outkey == "SSFR_MED") stout << setw(12) << sSFRmed[0] << " ";
    if (outkey == "SSFR_INF") stout << setw(12) << sSFRmed[1] << " ";
    if (outkey == "SSFR_SUP") stout << setw(12) << sSFRmed[2] << " ";

    if (outkey == "COL1_MED") stout << setw(12) << col1med[0] << " ";
    if (outkey == "COL1_INF") stout << setw(12) << col1med[1] << " ";
    if (outkey == "COL1_SUP") stout << setw(12) << col1med[2] << " ";

    if (outkey == "COL2_MED") stout << setw(12) << col2med[0] << " ";
    if (outkey == "COL2_INF") stout << setw(12) << col2med[1] << " ";
    if (outkey == "COL2_SUP") stout << setw(12) << col2med[2] << " ";

    if (outkey == "MREF_MED") stout << setw(12) << Mrefmed[0] << " ";
    if (outkey == "MREF_INF") stout << setw(12) << Mrefmed[1] << " ";
    if (outkey == "MREF_SUP") stout << setw(12) << Mrefmed[2] << " ";

    // ABSOLUTE MAGNITUDES
    stout.setf(ios::fixed);
    stout.precision(0);  // Integer:
    if (outkey == "MABS_FILT()") {
      for (size_t l = 0; l < absfilt.size(); l++)
        stout << setw(4) << absfilt[l] << " ";
    }
    stout.setf(ios::fixed, ios::floatfield);
    stout.precision(3);  // Float:
    if (outkey == "K_COR()") {
      for (size_t l = 0; l < kap.size(); l++)
        stout << setw(13) << setprecision(5) << kap[l] << " ";
    }
    if (outkey == "MAG_ABS()") {
      for (size_t l = 0; l < mabs.size(); l++)
        stout << setw(12) << setprecision(5) << mabs[l] << " ";
    }
    if (outkey == "EMAG_ABS()") {
      for (size_t l = 0; l < emabs.size(); l++)
        stout << setw(13) << setprecision(5) << emabs[l] << " ";
    }

    // PREDICTED MAGNITUDES IN ADDITIONAL FILTERS
    stout.setf(ios::scientific, ios::floatfield);
    stout.precision(5);
    if (outkey == "MAG_PRED()") {
      for (size_t l = 0; l < magPred.size(); l++)
        stout << setw(15) << magPred[l] << " ";
    }
    if (outkey == "ABSMAG_PRED()") {
      for (size_t l = 0; l < absmagPred.size(); l++)
        stout << setw(15) << absmagPred[l] << " ";
    }

    // 8 MAIN EMISSION LINES LYA, OII, HB, OIIIA, OIIIB, HA, SIIIA, SIIIB
    stout.setf(ios::scientific, ios::floatfield);
    stout.precision(5);
    for (const auto &item : results_emission_lines) {
      if (outkey == item.first)
        stout << setw(15) << setprecision(5) << item.second << " ";
    }
    // Having all the flux in output
    if (outkey == "EM_FLUX()") {
      for (size_t l = 0; l < 65; l++) stout << setw(15) << fluxEL_SED[l] << " ";
    }

    // LIMITS FOR VMAX
    stout.setf(ios::fixed, ios::floatfield);
    stout.precision(4);  // Float:
    if (outkey == "LIMITS_ZMAX") stout << setw(12) << limits_zmax << " ";
    if (outkey == "LIMITS_MFAINT") stout << setw(12) << limits_Mfaint << " ";
  }
  stout << endl;

  return;
}

/*
 write the header of the PDF(z)
*/
void onesource::write_pdz_header(vector<string> pdztype,
                                 unordered_map<string, ofstream> &stpdz,
                                 const time_t &ti1) {
  // Loop over the PDF type wanted in output
  for (const auto &type : pdztype) {
    stpdz[type] << "# Creation date: " << asctime(localtime(&ti1));
    stpdz[type] << "# Probability associated to the following steps " << endl
                << "# Id ";
    for (const auto &xval : pdfmap[maptype[type]].xaxis)
      stpdz[type] << "P" << xval << " ";
    stpdz[type] << endl;
  }
  return;
}

/*
 write the PDF(z)
*/
void onesource::write_pdz(vector<string> pdztype,
                          unordered_map<string, ofstream> &stpdz) {
  // Loop over the PDF type wanted in output
  for (const auto &type : pdztype) {
    stpdz[type] << setw(15) << std::fixed << setprecision(4) << spec << " ";
    for (const auto &xval : pdfmap[maptype[type]].vPDF)
      stpdz[type] << setw(16) << std::scientific << xval << " ";
    stpdz[type] << endl;
  }
  return;
}

/*
 INTERPOLATE LINEARILY IN THE LIBRARY
 Do it only for GAL
*/
void onesource::interp_lib(vector<SED *> &fulllib, const int imagm,
                           cosmo lcdm) {
  magm.clear();

  // Take the value of zs to interpolate in the library
  if (indmin[0] >= 0) {
    // If zs is above or below the best fit bin
    int sens = int(copysign(1.0, consiz - (fulllib[indmin[0]]->red)));
    // Since galaxy templates are duplicated when add dispersion in emission
    // lines, need to go to the next template having a different redshift
    int deca = sens;
    if ((size_t)(indmin[0] + deca) <= fulllib.size() - 1) {
      while (fulllib[indmin[0]]->red == fulllib[indmin[0] + deca]->red) {
        deca += sens;
        if ((size_t)(indmin[0] + deca) >= fulllib.size() ||
            indmin[0] + deca < 0) {
          deca = 0;
          break;
        }
      }
    } else {
      deca = 0;
    }

    // Check if the next template is still the good one (could happen with BC03
    // if age close to the age of the Universe) If not, the redshift would be
    // zero, so no interpolation can be done if it happens
    if (fulllib[indmin[0] + deca]->red < 1.e-10) deca = 0;

    // Store in two SED around zs
    SED SEDa(*(fulllib[indmin[0]]));
    SED SEDb(*(fulllib[indmin[0] + deca]));

    // Check that the interpolation can be done
    // same model in the library around zs
    if (SEDa.is_same_model(SEDb) && deca != 0) {
      // Linear factor for interpolation
      double a = abs((consiz - SEDa.red) / (SEDa.red - SEDb.red));
      if (abs(SEDa.age - SEDb.age) > 1.)
        cout << " ALARM " << SEDa.age << " " << SEDb.age << endl;

      // Loop over the filters
      for (int k = 0; k < imagm; k++) {
        // Check if the interpolation if possible
        if (SEDb.mag[k] > 90) a = 0.0;
        // Interpolation of the predicted magnitudes.
        magm.push_back(a * SEDb.mag[k] + (1. - a) * SEDa.mag[k]);
        // Check if the interpolation if possible
        if (SEDb.kcorr[k] > 90) a = 0.0;
        // Interpolation of the k-corrections
        kap.push_back(a * SEDb.kcorr[k] + (1. - a) * SEDa.kcorr[k]);
      }

    } else {
      // No interpolation is possible
      // Loop over the filters
      for (int k = 0; k < imagm; k++) {
        magm.push_back(SEDa.mag[k]);
        kap.push_back(SEDa.kcorr[k]);
      }
    }

    // Retrieve the rest-frame magnitudes
    int index0 = indmin[0];
    // go down in the index as long as the redshift >0 and that the model is the
    // same
    while ((*(fulllib[index0])).red > 1.e-20 &&
           SEDa.is_same_model(*fulllib[index0])) {
      index0--;
    }
    if ((*(fulllib[index0])).red > 1.e-20)
      cout << " do not find redshift 0 for rest-frame colors. Problem." << endl;
    // Store the predicted magnitudes at z=0
    for (int k = 0; k < imagm; k++) {
      magm0.push_back((*(fulllib[index0])).mag[k]);
    }

    // Rescale the predicted magnitudes with the template scaling dm
    for (int k = 0; k < imagm; k++) {
      magm[k] = magm[k] - 2.5 * log10(dmmin[0]);
      magm0[k] = magm0[k] - 2.5 * log10(dmmin[0]);
    }
  } else {
    for (int k = 0; k < imagm; k++) {
      magm.push_back(-999);
      magm0.push_back(-999);
      kap.push_back(-999);
    }
  }
}

/*
 REDSHIFT INTERPOLATION  if  ZINTP=true
*/
void onesource::interp(const bool zfix, const bool zintp, cosmo lcdm) {
  // keep distance modulus before interpolation
  double distGbef = lcdm.distMod(zmin[0]);
  double distQbef = lcdm.distMod(zmin[1]);

  // If zfix=yes, use the spec-z as best value
  if (zfix) {
    // Modify the zmin for the galaxy
    zmin[0] = zs;
    zmin[1] = zs;
  } else {
    // If zfix=no but an interpolation is required
    if (zintp) {
      // The new zmin is coming from the PDF, the one with the minimum chi2
      zmin[0] = pdfmap[9].int_parab();   // Galaxies
      zmin[1] = pdfmap[10].int_parab();  // AGN
    }
  }

  // Need to rescale the normalisation accoring to the new redshift
  // distance modulus after interpolation
  double distGaft = lcdm.distMod(zmin[0]);
  double distQaft = lcdm.distMod(zmin[1]);
  dmmin[0] = dmmin[0] * pow(10., (0.4 * (distGaft - distGbef)));
  dmmin[1] = dmmin[1] * pow(10., (0.4 * (distQaft - distQbef)));
}

/*
 REDSHIFT UNCERTAINTIES (zmin,zmax) for dChi2=1,2.71,6.63
 derive the photo-z uncerstainties from the minimum chi2+level
*/
void onesource::uncertaintiesMin() {
  zgmin.push_back(zmin[0]);
  zqmin.push_back(zmin[1]);

  // Define the delta chi2
  double dchi[3] = {1., 2.71, 6.63};
  pair<double, double> interv;
  // Add the uncertainties, lower and upper for 68%, 90% and 99%
  for (size_t k = 0; k < 3; k++) {
    // Photo-z for galaxies
    interv = pdfmap[9].uncMin(dchi[k]);
    zgmin.push_back(interv.first);
    zgmin.push_back(interv.second);

    // Photo-z for QSO
    interv = pdfmap[10].uncMin(dchi[k]);
    zqmin.push_back(interv.first);
    zqmin.push_back(interv.second);
  }

  return;
}

/*
 REDSHIFT UNCERTAINTIES (zmin,zmax) for 68%, 90% and 99% of the total PDF area
*/
void onesource::uncertaintiesBay() {
  // 0:["MASS"] / 1:["SFR"] / 2:["SSFR"] / 3:["LDUST"] / 4:["LIR"] / 5:["AGE"] /
  // 6:["COL1"] / 7:["COL2"] / 8:["MREF"]/ 9:["MIN_ZG"] / 10:["MIN_ZQ"] /
  // 11:["BAY_ZG"] / 12:["BAY_ZQ"] First element is the median of the PDF
  zgmed.push_back(pdfmap[11].levelCumu2x(0.5));
  zqmed.push_back(pdfmap[12].levelCumu2x(0.5));
  massmed.push_back(pdfmap[0].levelCumu2x(0.5));
  SFRmed.push_back(pdfmap[1].levelCumu2x(0.5));
  sSFRmed.push_back(pdfmap[2].levelCumu2x(0.5));
  Ldustmed.push_back(pdfmap[3].levelCumu2x(0.5));
  agemed.push_back(pdfmap[5].levelCumu2x(0.5));
  col1med.push_back(pdfmap[6].levelCumu2x(0.5));
  col2med.push_back(pdfmap[7].levelCumu2x(0.5));
  Mrefmed.push_back(pdfmap[8].levelCumu2x(0.5));
  LIRmed.push_back(pdfmap[4].levelCumu2x(0.5));

  // Define the confidence levels
  double confLev[3] = {68, 90, 99};
  // Add the uncertainties, lower and upper for 68%, 90% and 99%
  pair<double, double> interv;
  for (size_t k = 0; k < 3; k++) {
    // Photo-z for galaxies
    interv = pdfmap[11].credible_interval(confLev[k], zgmed[0]);
    zgmed.push_back(interv.first);
    zgmed.push_back(interv.second);

    // Photo-z for QSO
    interv = pdfmap[11].credible_interval(confLev[k], zqmed[0]);
    zqmed.push_back(interv.first);
    zqmed.push_back(interv.second);

    // Mass for galaxies
    interv = pdfmap[0].credible_interval(confLev[k], massmed[0]);
    massmed.push_back(interv.first);
    massmed.push_back(interv.second);

    // SFR for galaxies
    interv = pdfmap[1].credible_interval(confLev[k], SFRmed[0]);
    SFRmed.push_back(interv.first);
    SFRmed.push_back(interv.second);

    // sSFR for galaxies
    interv = pdfmap[2].credible_interval(confLev[k], sSFRmed[0]);
    sSFRmed.push_back(interv.first);
    sSFRmed.push_back(interv.second);

    // age for galaxies
    interv = pdfmap[5].credible_interval(confLev[k], agemed[0]);
    agemed.push_back(interv.first);
    agemed.push_back(interv.second);

    // Ldust for galaxies
    interv = pdfmap[3].credible_interval(confLev[k], Ldustmed[0]);
    Ldustmed.push_back(interv.first);
    Ldustmed.push_back(interv.second);

    // First color rest-frame for galaxies
    interv = pdfmap[6].credible_interval(confLev[k], col1med[0]);
    col1med.push_back(interv.first);
    col1med.push_back(interv.second);

    // Second color rest-frame for galaxies
    interv = pdfmap[7].credible_interval(confLev[k], col2med[0]);
    col2med.push_back(interv.first);
    col2med.push_back(interv.second);

    // Reference absolute magnitude for galaxies
    interv = pdfmap[8].credible_interval(confLev[k], Mrefmed[0]);
    Mrefmed.push_back(interv.first);
    Mrefmed.push_back(interv.second);

    // Estimate the 68% uncertainties for the IR luminosity
    interv = pdfmap[5].credible_interval(confLev[k], LIRmed[0]);
    LIRmed.push_back(interv.first);
    LIRmed.push_back(interv.second);
  }

  return;
}

/*
     Search for the peaks in the profile likelihood function vs z
     And sort them from highest to smallest
     Keep the parameters corresponding to the second solution.
     As pdfmap[9]=="MIN_ZG" is harcoded here, this is only valid
     for the GAL solutions.
*/
void onesource::secondpeak(vector<SED *> &fulllib, const double dz_win,
                           const double min_thres) {
  // Detect the other maximum in the PDF
  pdfmap[9].secondMax(dz_win);
  // Default measurement
  zsec = -99.9;
  zsecChi2 = 1.e9;
  zsecEbv = -99.;
  zsecExtlaw = -99;
  zsecScale = -99.;
  zsecProb = -99.;
  zsecMod = -99;
  zsecAge = -99;
  indminSec = -99;
  // If the PDF is above the threshold, the second solution exist
  if (pdfmap[9].secondP.size() > 1) {
    if (pdfmap[9].secondP[1] > min_thres * pdfmap[9].secondP[0]) {
      // index of the corresponding SED
      int idx = pdfmap[9].secondInd[1];
      // Keep the information for all the secondary solution
      zsec = pdfmap[9].secondX[1];
      zsecChi2 = (fulllib[idx])->chi2;
      zsecEbv = (fulllib[idx])->ebv;
      zsecExtlaw = (fulllib[idx])->extlawId;
      zsecScale = (fulllib[idx])->dm;
      zsecProb = pdfmap[9].secondP[1];
      zsecMod = (fulllib[idx])->nummod;
      zsecAge = (fulllib[idx])->age;
      indminSec = idx;
    }
  }

  return;
}

/*
 ABSOLUTE MAGNITUDES
*/
void onesource::absmag(const vector<vector<int>> &bestFlt,
                       const vector<vector<double>> &maxkcolor, cosmo lcdm,
                       const vector<double> gridz) {
  int fobs;
  double errmagabs;

  // Find the index of the redshift into gridz
  // Increment as long as the redshift in the grid is lower than the considered
  // redshift
  int indexz = 0;
  while (consiz > gridz[indexz] && indexz < int(gridz.size() - 1)) indexz++;

  // loop over all filters in which the absolute magnitude needs to be computed
  for (size_t k = 0; k < mab.size(); ++k) {
    // Find the filter used for the corresponding apparent magnitude
    fobs = bestFlt[indexz][k];

    // If good filter defined and uncertainty in the filter not to high
    if (fobs >= 0 && msab[fobs] <= 0.3 && msab[fobs] >= 0) {
      // Absolute magnitude: mobs + color rest-frame (ref-fobs) -  k(fobs)
      // -DM(z)
      mabs.push_back(mab[fobs] + magm0[k] - magm0[fobs] - kap[fobs] -
                     lcdm.distMod(consiz));
      // uncertainty on the mag abs: sum the apparent magnitude uncertainties +
      // maximum uncertainties on the k-correction
      errmagabs = sqrt(pow(maxkcolor[indexz][k], 2.) + pow(sab[fobs], 2.));
      emabs.push_back(errmagabs);
      absfilt.push_back(fobs);  // filter used for apparent magnitude
    } else {
      // Directly from the template
      mabs.push_back(magm0[k]);
      emabs.push_back(-1);
      absfilt.push_back(-1);
    }
  }

  return;
}

/*
 Compute the emission lines to add them in the output later
*/
void onesource::computeEmFlux(vector<SED *> &fulllib, cosmo lcdm,
                              vector<opa> opaAll) {
  /// rest frame wavelengths [Lya, OII, Hb, OIIIa, OIIIb, Ha, SIIIa, SIIIb]
  array<double, 8> lambda_em = {1215.67, 3727.00, 4861.32, 4958.91,
                                5006.84, 6562.80, 9068.60, 9530.85};
  array<double, 8> flux_em = {0, 0, 0, 0, 0, 0, 0, 0};
  array<double, 8> EW_em = {-99, -99, -99, -99, -99, -99, -99, -99};

  int ind0;
  if (indmin[0] >= 0) {
    // Do it only if emission have been generated
    if (fulllib[indmin[0]]->has_emlines) {
      // Generate a SED based on the index at z=0: extract the rest-frame SED
      ind0 = (fulllib[indmin[0]])->index_z0;
      // Generate the continuum SED at z=0
      GalSED SEDz0_Gal(*(fulllib[ind0]));
      // Rescale SED continuum
      SEDz0_Gal.rescale(dmmin[0]);
      // Opacity applied in rest-frame, depending on the redshift of the source
      SEDz0_Gal.applyOpa(opaAll);

      // Create a SED of emission lines at z
      GalSED SEDz0_Em(*(fulllib[ind0]));
      // GalSED SEDz0_Em(*(fulllib[indmin[0]])); // To be re-activated if EL
      // ratio can change with z
      //// Generate the emission line spectra
      SEDz0_Em.generateEmSpectra(40);
      // Rescale the emission lines
      SEDz0_Em.rescale(dmmin[0]);
      // Opacity applied in rest-frame (lines are in rest-frame), depending on
      // the redshift of the source
      SEDz0_Em.applyOpa(opaAll);

      // Integrated flux divide by the filter area for emission lines. Reference
      // was at 10pc, which explain the 100 factor. The distance is in Mpc
      // explaining the 10^12
      double rescaleDist =
          100. / (pow(lcdm.distMet(consiz) * (1 + consiz), 2) * pow(10, 12));

      // Loop over each emission lines. Only the 8 main lines that we can have
      // in output (lyman alpha, OII, Halpha, etc)
      for (int k = 0; k < 8; k++) {
        // Create narrow filter, be careful to not encompass another line
        flt fltEm(lambda_em[k] - 10., lambda_em[k] + 10., 21);
        // Create narrow filter for the continuum
        flt fltEmC(lambda_em[k] - 30., lambda_em[k] + 30., 21);

        // Integrate the SED (emission lines & continuum) within the filter
        vector<double> emF, contiF;
        // If a sepctra exist
        if (SEDz0_Em.lamb_flux.size() > 0) {
          // int (F*T) dlambda
          // SEDz0_Em.warning_integrateSED(fltEm);
          emF = SEDz0_Em.integrateSED(fltEm);
          // int (F*T) dlambda / int (T) dlambda
          // SEDz0_Gal.warning_integrateSED(fltEmC);
          contiF = SEDz0_Gal.integrateSED(fltEmC);
          // Integrated flux divide by the filter area for emission lines.
          // Reference was at 10pc, which explain the 100 factor. The distance
          // is in Mpc explaining the 10^12
          flux_em[k] = emF[3] * rescaleDist;
          // Equivalent width
          // int (Fem*T) dlambda / [int (Fcont*T) dlambda / int (T) dlambda]
          EW_em[k] = emF[3] * contiF[0] / contiF[3];
          emF.clear();
          contiF.clear();
        } else {
          flux_em[k] = 0.;
          EW_em[k] = 0.;
        }
        // Keep all emission lines
        for (size_t l = 0; l < 65; l++)
          fluxEL_SED[l] = (SEDz0_Em.fac_line[l].val) * rescaleDist * dmmin[0];
        // Replace Lyman alpha by the integral of the line in the spectra, to
        // get the opacity of the intergalactic medium
        fluxEL_SED[0] = flux_em[0];
      }
    }
  }

  vector<string> line_names = {"LYA",   "OII", "HB",    "OIIIA",
                               "OIIIB", "HA",  "SIIIA", "SIIIB"};
  size_t count = 0;
  for (const string &name : line_names) {
    results_emission_lines["EM_FLUX_" + name] = flux_em[count];
    results_emission_lines["EM_EW_" + name] = EW_em[count];
    count++;
  }

  return;
}

/*
Compute predicted magnitude in new filters
*/
void onesource::computePredMag(vector<SED *> &fulllib, cosmo lcdm,
                               vector<opa> opaAll, vector<flt> allFltAdd) {
  double val;

  if (indmin[0] > 0) {
    // But we need first the fac_line at the right redshift
    // GalSED SED_emSave(*(fulllib[indmin[0]])); // To be re-activated if EL
    // ratio can change with z Extract the rest-frame spectra for the continuum
    // (spectra store only at z=0)
    int ind0 = (fulllib[indmin[0]])->index_z0;
    GalSED SED_gal(*(fulllib[ind0]));
    // Associate the right emission lines (fac_line depends on z)
    // SED_gal.fac_line=SED_emSave.fac_line; // To be re-activated if EL ratio
    // can change with z Generate the spectra at the right redshift, opacity,
    // emission lines
    if (consiz > 0) SED_gal.generate_spectra(consiz, dmmin[0], opaAll);
    // check that the filter cover a SED
    SED_gal.warning_integrateSED(allFltAdd, true);

    magPred.clear();
    // Loop over the filters
    for (size_t itf = 0; itf < allFltAdd.size(); ++itf) {
      // Derive the AB magnitudes in each filter
      vector<double> intFlux = SED_gal.integrateSED(allFltAdd[itf]);
      if (intFlux[3] > INVALID_FLUX) {
        if (intFlux[3] > NULL_FLUX) {
          val = flux2mag(intFlux[3] / intFlux[1], lcdm.distMod(consiz));
        } else
          val = HIGH_MAG;
      } else
        val = INVALID_MAG;
      // predicted mag stored into the vector of SED
      magPred.push_back(val);
    }

  } else {
    for (size_t itf = 0; itf < allFltAdd.size(); ++itf) magPred.push_back(-99.);
  }

  return;
}

/*
Compute absolute magnitudes in new filters
*/
void onesource::computePredAbsMag(vector<SED *> &fulllib, cosmo lcdm,
                                  vector<opa> opaAll, vector<flt> allFltAdd) {
  double val;

  if (indmin[0] > 0) {
    // But we need first the fac_line at the right redshift
    // GalSED SED_emSave(*(fulllib[indmin[0]])); // To be re-activated if EL
    // ratio can change with z Extract the rest-frame spectra for the continuum
    // (spectra store only at z=0)
    int ind0 = (fulllib[indmin[0]])->index_z0;
    GalSED SED_gal(*(fulllib[ind0]));
    // Associate the right emission lines (fac_line depends on z)
    // SED_gal.fac_line=SED_emSave.fac_line; // To be re-activated if EL ratio
    // can change with z Generate the spectra at the right redshift, opacity,
    // emission lines
    if (consiz > 0) SED_gal.generate_spectra(0., dmmin[0], opaAll);
    // check that the filter cover a SED
    SED_gal.warning_integrateSED(allFltAdd, true);

    absmagPred.clear();
    // Loop over the filters
    for (size_t itf = 0; itf < allFltAdd.size(); ++itf) {
      // Derive the AB magnitudes in each filter
      vector<double> intFlux = SED_gal.integrateSED(allFltAdd[itf]);
      if (intFlux[3] > 0) {
        val = flux2mag(intFlux[3] / intFlux[1]);
      } else
        val = HIGH_MAG;
      // predicted mag stored into the vector of SED
      absmagPred.push_back(val);
    }

  } else {
    for (size_t itf = 0; itf < allFltAdd.size(); ++itf)
      absmagPred.push_back(-99.);
  }

  return;
}

/*
Compute the maximum redshift at which a source can be observed given the
magnitude selection
*/
void onesource::limits(vector<SED *> &fulllib, vector<double> &limits_zbin,
                       int limits_ref, vector<int> &limits_sel,
                       vector<double> &limits_cut) {
  // Be careful with the convention of filter in the parameter file : start at 1
  if (limits_ref > 0) {
    limits_ref -= 1;
  }

  int deca = 1;
  int FltSel = 0;
  double MagSel = 90;
  double diff = 0;

  // Check if feasible
  if ((indmin[0] >= 0) && (size_t)(indmin[0]) < fulllib.size() - 1) {
    // Loop over the redshift bin for the selection filter
    FltSel = -1;
    // Loop over the redshift bins indicated in the selection vector
    for (int k = 0; k < int(limits_zbin.size()) - 1; k++) {
      // Check in which redshift bin the source falls
      if (consiz >= limits_zbin[k] && consiz < limits_zbin[k + 1]) {
        // Define the selection filter
        // Be careful with the convention of filter in the parameter file :
        // start at 1
        FltSel = limits_sel[k] - 1;
        MagSel = limits_cut[k];
      }
    }

    // check if the difference between the observed magnitude and magnitude
    // limit
    while (diff < (MagSel - mab[FltSel])) {
      // check if the difference between the observed magnitude at the minimum
      // redshift and maximum redshift
      diff = fulllib[indmin[0] + deca]->mag[FltSel] -
             fulllib[indmin[0]]->mag[FltSel];
      deca += 1;

      // Check that we are below the size of the library
      if ((size_t)(indmin[0] + deca) >= fulllib.size() - 1) {
        deca -= 1;
        break;
      }
      // Check that we didn't reach the maximum redshift of this model
      if ((fulllib[indmin[0] + deca]->red < fulllib[indmin[0]]->red)) {
        deca -= 1;
        break;
      }
    }
    // Assign the final z_max
    limits_zmax = (fulllib[indmin[0] + deca]->red);

    // Compute the absolute magnitude faint limit
    limits_Mfaint = -999.9;
    if (limits_ref >= 0 && (MagSel > mab[FltSel])) {
      limits_Mfaint = mabs[limits_ref] + MagSel - mab[FltSel];
    }

  } else {
    // cout << "The computation of z_max makes no sense " << endl;
  }

  return;
}

// Used for the python interface in the notebook
pair<vector<double>, vector<double>> onesource::best_spec_vec(
    short sol, vector<SED *> &fulllib, cosmo lcdm, vector<opa> opaAll,
    double minl, double maxl) {
  // sol=0 for GAL-1, sol=1 for GAL-2, sol=2 for FIR, sol=3 for QSO, sol=4 for
  // STAR
  int index = -1;     // index of the best fit
  int index_z0 = -1;  // index of the z=0 SED of the best fit
  double scale = 0.;  // normalisation of the best fit
  double red = 0.;    // redshift of the best fit
  if (sol == 0) {     // GAL-1
    index = indmin[0];
    if (index >= 0) {
      index_z0 = fulllib[index]->index_z0;
      scale = dmmin[0];
      red = consiz;
    }
  } else if (sol == 1) {  // GAL-2
    index = indminSec;
    if (index >= 0) {
      index_z0 = fulllib[index]->index_z0;
      scale = zsecScale;
      red = zsec;
    }
  } else if (sol == 2) {  // FIR
    index = indminIR;
    if (index >= 0) {
      index_z0 = fulllib[index]->index_z0;
      scale = dmminIR;
      red = zminIR;
    }
  } else if (sol == 3) {  // QSO
    index = indmin[1];
    if (index >= 0) {
      index_z0 = fulllib[index]->index_z0;
      scale = dmmin[1];
      red = zmin[1];
    }
  } else if (sol == 4) {  // STAR
    index = indmin[2];
    if (index >= 0) {
      index_z0 = index;
      scale = dmmin[2];
      red = zmin[2];  // this is 0 by construction for stars
    }
  }

  if (index >= 0) {
    SED tmpSED(*fulllib[index_z0]);
    // To be re-activated if EL ratio can change with z
    // if(sol==0 || sol == 1){
    //  SED emlinesSED(*fulllib[index]);
    //  tmpSED.fac_line = emlinesSED.fac_line;
    // }
    tmpSED.generate_spectra(red, scale, opaAll);
    // for STAR, red=zmin[2]=0 by construction
    return tmpSED.get_data_vector(minl, maxl, true, lcdm.distMod(red));
  }

  pair<vector<double>, vector<double>> p({}, {});
  return p;
}

/*
 WRITE .SPEC FILE
*/
void onesource::writeSpec(vector<SED *> &fulllib, vector<SED *> &fulllibIR,
                          cosmo lcdm, vector<opa> opaAll,
                          const vector<flt> &allflt, const string outspdir) {
  // open the output file
  ofstream stospec;
  string ospec = "Id" + string(spec) + ".spec";
  // if different from YES, include the file in a subdirectory indicated by the
  // users. In relative. If directory doesn't exist, create it.
  if (outspdir != "YES") {
    if (!(filesystem::exists(outspdir))) filesystem::create_directory(outspdir);
    ospec = outspdir + "/" + ospec;
  }
  stospec.open(ospec.c_str());

  /*
    Find the minimum and maximum lambda, to restrict the size of the SED in the
    .spec file
  */

  double minl = 1.e10;
  double maxl = 0;
  // Leff-2*width of the filter for the min, Leff+ 2*width for the max
  for (size_t k = 0; k < allflt.size(); k++) {
    if (minl > (allflt[k]).lmean - 2 * (allflt[k]).dwidth)
      minl = (allflt[k]).lmean - 2 * (allflt[k]).dwidth;
    if (maxl < (allflt[k]).lmean + 2 * (allflt[k]).dwidth)
      maxl = (allflt[k]).lmean + 2 * (allflt[k]).dwidth;
  }
  // In angstrom
  minl = minl * 10000.;
  maxl = maxl * 10000.;

  /*
    GALAXY CASE
  */
  vector<double> lG, mG;
  auto tmp = best_spec_vec(0, fulllib, lcdm, opaAll, minl, maxl);
  lG = tmp.first;
  mG = tmp.second;
  /*
    GALAXY CASE, SECOND SOLUTION
  */
  vector<double> lGsec, mGsec;
  tmp = best_spec_vec(1, fulllib, lcdm, opaAll, minl, maxl);
  lGsec = tmp.first;
  mGsec = tmp.second;
  /*
    GALAXY FIR CASE
  */
  vector<double> lIR, mIR;
  tmp = best_spec_vec(2, fulllibIR, lcdm, opaAll, minl, maxl);
  lIR = tmp.first;
  mIR = tmp.second;
  /*
    QSO CASE
  */
  vector<double> lQ, mQ;
  tmp = best_spec_vec(3, fulllib, lcdm, opaAll, minl, maxl);
  lQ = tmp.first;
  mQ = tmp.second;
  /*
    STAR CASE
  */
  vector<double> lS, mS;
  tmp = best_spec_vec(4, fulllib, lcdm, opaAll, minl, maxl);
  lS = tmp.first;
  mS = tmp.second;

  /*
  WRITE OUTPUT SPECTRA
  */

  // ... write Nfilt Zspec Zphot
  stospec << "# Ident Zspec Zphot " << endl;
  stospec << spec << " " << zs << " " << consiz << endl;

  stospec << "# Mag emag  Lbd_mean  Lbd_width Mag_gal  Mag_FIR  Mag_BCSTOCH "
          << endl;
  stospec << "FILTERS  " << allflt.size() << endl;

  stospec << "# Zstep  PDF " << endl;
  stospec << "PDF  " << pdfmap[11].size() << " " << pdfmap[9].size() << endl;

  stospec << "# Type Nline Model Library Nband  Zphot Zinf Zsup Chi2  PDF  "
             "Extlaw EB-V Lir Age  Mass SFR SSFR"
          << endl;

  // Galaxy
  if (indmin[0] > 0) {
    stospec << "GAL-1 " << lG.size() << " " << imasmin[0] << " 1 " << nbused
            << " " << consiz << " " << zgmin[0] << " " << zgmin[1] << " ";
    stospec << chimin[0] << " " << " -1" << " " << results["EXTLAW_BEST"] << " "
            << results["EBV_BEST"] << " " << results["LDUST_BEST"] << " "
            << results["AGE_BEST"] << " " << results["MASS_BEST"] << " "
            << results["SFR_BEST"] << " " << results["SSFR_BEST"] << endl;

  } else {
    stospec
        << "GAL-1 0 -1 -1 -1 -1. -1. -1. -1. -1. -1 -1. -1. -1. -1. -1. -1. "
        << endl;
  }

  // second solution
  stospec << "GAL-2 " << lGsec.size() << " " << zsecMod << " 1 " << nbused
          << " " << zsec << " -1 -1  " << zsecChi2 << " " << zsecProb
          << "  -1 -1. -1. -1. -1. -1. -1. " << endl;

  // Galaxy FIR
  if (indminIR > 0) {
    stospec << "GAL-FIR " << lIR.size() << " " << imasminIR << " 1 " << nbused
            << " " << zminIR << " -1 -1 ";
    stospec << chiminIR << " " << " 0 -1 -1 " << -999. << " -1 -1 -1 -1 "
            << endl;
  } else {
    stospec
        << "GAL-FIR 0 -1 -1 -1 -1. -1. -1. -1. -1. -1 -1. -1. -1. -1. -1. -1. "
        << endl;
  }

  // STOCH
  stospec
      << "GAL-STOCH 0 -1 -1 -1 -1. -1. -1. -1. -1. -1 -1. -1. -1. -1. -1. -1. "
      << endl;

  // QSO
  if (indmin[1] > 0) {
    stospec << "QSO " << lQ.size() << " " << imasmin[1] << " 2 " << nbused
            << " " << zmin[1] << " 0 0 ";
    stospec << chimin[1] << " 0. -1 -1. -1. -1. -1. -1. -1. " << endl;

  } else {
    stospec << "QSO 0 -1 -1 -1 -1. -1. -1. -1. -1. -1 -1. -1. -1. -1. -1. -1. "
            << endl;
  }

  // STAR
  if (indmin[2] > 0) {
    stospec << "STAR " << lS.size() << " " << imasmin[2] << " 3 " << nbused
            << " 0 0 0 ";
    stospec << chimin[2] << " 0. -1 -1. -1. -1. -1. -1. -1. " << endl;

  } else {
    stospec << "STAR 0 -1 -1 -1 -1. -1. -1. -1. -1. -1 -1. -1. -1. -1. -1. -1. "
            << endl;
  }

  // write mag obs + predicted + filters (rajouter flmoy(k),flwidth(k),
  // magb(k),magfirb(k),fluxphys(k)+bused
  for (size_t l = 0; l < mab.size(); l++)
    stospec << setw(13) << mab[l] << setw(13) << msab[l] << setw(10)
            << ((allflt[l]).lmean) * 10000 << setw(10)
            << ((allflt[l]).dwidth) * 10000 << " -1 -1 -1 " << busnorma[l]
            << " " << magm[l] << endl;

  // write PDF
  if (pdfmap[11].size() != pdfmap[9].vPDF.size())
    cout << " Problem PDF bay and min differ" << endl;
  for (size_t l = 0; l < pdfmap[11].size(); l++)
    stospec << setw(6) << pdfmap[11].xaxis[l] << setw(14) << pdfmap[11].vPDF[l]
            << setw(14) << pdfmap[9].vPDF[l] << endl;

  // GAL: Loop over the full vector with the two vectors concatenated
  for (size_t k = 0; k < lG.size(); ++k)
    stospec << fixed << setprecision(3) << lG[k] << " " << mG[k] << endl;

  // GAL second solution: Loop over the full vector with the two vectors
  // concatenated
  for (size_t k = 0; k < lGsec.size(); ++k)
    stospec << fixed << setprecision(3) << lGsec[k] << " " << mGsec[k] << endl;

  // GAL FIR : template from the IR library
  for (size_t k = 0; k < lIR.size(); ++k)
    stospec << fixed << setprecision(3) << lIR[k] << " " << mIR[k] << endl;

  // QSO: Loop over the full vector with the two vectors concatenated
  for (size_t k = 0; k < lQ.size(); ++k)
    stospec << fixed << setprecision(3) << lQ[k] << " " << mQ[k] << endl;

  // STAR: Loop over the full vector with the two vectors concatenated
  for (size_t k = 0; k < lS.size(); ++k)
    stospec << fixed << setprecision(3) << lS[k] << " " << mS[k] << endl;

  return;
}

/*
 WRITE FULL CHI2 FILE
*/
void onesource::writeFullChi(vector<SED *> &fulllib) {
  double sca;

  // open the output file
  ofstream stochi;
  string ochi = "Id" + string(spec) + ".chi";
  stochi.open(ochi.c_str());
  stochi << "# nlib z model Age Extlaw EB_V Ldust Luv Lr Lk Ldust2 Mo SFR Chi2"
         << endl;

  // Loop over all SEDs from the library
  for (size_t k = 0; k < fulllib.size(); k++) {
    // Normalisation
    sca = fulllib[k]->dm;
    // Write
    stochi << fulllib[k]->nlib << " ";
    stochi << fulllib[k]->red << " ";
    stochi << fulllib[k]->nummod << " ";
    stochi << fulllib[k]->age << " ";
    stochi << fulllib[k]->extlawId << " ";
    stochi << fulllib[k]->ebv << " ";
    // Check that the scaling and mass defined
    if (((fulllib[k])->mass > 0) && sca > 0) {
      stochi << (fulllib[k])->ltir + log10(sca) << " ";
      stochi << (fulllib[k])->luv + log10(3.e18 * 400 / pow(2300, 2)) -
                    log10(Lsol) + log10(sca)
             << " ";
      stochi << (fulllib[k])->lopt + log10(3.e18 * 1000 / pow(6000, 2)) -
                    log10(Lsol) + log10(sca)
             << " ";
      stochi << (fulllib[k])->lnir + log10(3.e18 * 2000 / pow(22000, 2)) -
                    log10(Lsol) + log10(sca)
             << " ";
      stochi << log10((fulllib[k])->mass * sca) << " ";
      stochi << log10((fulllib[k])->sfr * sca) << " ";
    } else {
      stochi << " -99 -99 -99 -99 -99 -99 ";
    }
    stochi << fulllib[k]->chi2 << endl;
  }
  return;
}
