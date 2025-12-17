/*

  17/11/14
  Implementation of the functions of the SED class

*/

#include "SED.h"

#include <algorithm>
#include <cmath>
#include <fstream>   // print output file
#include <iomanip>   // std::set precision
#include <iostream>  // print standard file
#include <sstream>
#include <string>
#include <vector>

#include "emission_lines.h"
#include "globals.h"
#include "oneElLambda.h"

using namespace std;

/*
  General construction of the basis class SED
*/
SED::SED(const string nameC, int nummodC, string type) {
  name = nameC;
  nummod = nummodC;  // number of the model in the list

  // Define nlib to have the type as an interger GAL=0, QSO=1, STAR=2
  nlib = string_to_object(type);

  has_emlines = false;
  idAge = 0;     // index of the age into the SED
  age = -999;    // Age  (yr)
  red = 0.;      // redshift considered for the SED
  index = 0;     // index of the model in the library
  index_z0 = 0;  // index of the model in the library at z=0
  mass = -999;   // mass (Mo)
  sfr = -999;    // SFR (Mo/yr)
  ltir = -999;   // Int_8um^1000um    L_lbda . dlbda    in Log unit Lo
  ebv = 0.;      // E(B-V) applied to the SED
  extlawId =
      0;  // index of the extinction law when dust attenuation has been applied
  chi2 = HIGH_CHI2;  // chi2 of the fit
  dm = -999.;        // Rescaling of the template
  distMod = 0;
  qi = {0., 0., 0., 0.};
}

/*
  extended constructor
*/
SED::SED(const string nameC, double tauC, double ageC, int nummodC,
         string typeC, int idAgeC)
    : SED(nameC, nummodC, typeC) {
  age = ageC;
  idAge = idAgeC;  // index of the age into the SED
}

/*
  destructor of the SEd class: remove the vectors
*/
SED::~SED() {
  mag.clear();
  lamb_flux.clear();
  kcorr.clear();
}

/*
  rescale the flux, valid for QSO/STARS/GAL
  input: rescaling factor
  output: void, modify the internal variable lamb_flux
*/
void SED::rescale(double scaleFac) {
  // Loop over the full vector with the two vectors concatenated
  for (vector<oneElLambda>::iterator it = lamb_flux.begin();
       it < lamb_flux.end(); ++it) {
    // multiply the flux by a given factor
    it->val = (it->val) * scaleFac;
  }
  return;
}

/*
  read the templates in ascii, ascii valid for QSO/STARS/GAL
  It will be overwritten in derived classes for more complex reading
*/
void SED::read(const string &sedFile) {
  ifstream ssed;
  string lit;

  // open the SED template file into a stream
  ssed.open(sedFile.c_str());
  if (!ssed) {
    throw invalid_argument(
        "Can't open the file with the list of SED to be used " + sedFile);
  }

  // Take the stream line by line
  while (getline(ssed, lit)) {
    // If the first character of the line is not #
    if (check_first_char(lit)) {
      // put the line into the stream ss again
      stringstream ss(lit);

      // fill the lambda/trans values of the SED
      double l, v;
      ss >> l;
      ss >> v;
      // if negative flux, put it at 0 because unphysical
      if (v < 0) v = 0.;
      lamb_flux.emplace_back(l, v, 1);
    }
  }

  if (lamb_flux.size() == 0) {
    throw runtime_error("SED::read(string sedFile) lam flux is zero");
  }

  // Close the stream
  ssed.close();

  // Sort the vector according to ascending lambda
  sort(lamb_flux.begin(), lamb_flux.end());

  return;
}

/*
  Write the sed in binary format with the basis format (used for STARS and QSO)
  That's the output of sedtolib
*/
void SED::writeSED(ofstream &ofsBin, ofstream &ofsPhys, ofstream &ofsDoc) {
  ofsBin.write((char *)&nummod, sizeof(int));

  long nbw = lamb_flux.size();               // Number of wavelength bins
  ofsBin.write((char *)&nbw, sizeof(long));  // Number of wavelength bins

  // Write the wavelength
  for (auto &oneEl : lamb_flux) {
    ofsBin.write((char *)&(oneEl.lamb), sizeof(double));
  }
  // Write the fluxes
  for (auto &oneEl : lamb_flux) {
    ofsBin.write((char *)&(oneEl.val), sizeof(double));
  }

  // Write the documentation
  ofsDoc << "MOD_" << nummod << " " << setw(6) << nummod << " " << setw(6)
         << nbw << " " << setw(10) << " " << name << endl;

  return;
}

/*
  read the sed library in binary format in the basis case (Stars and QSO)
*/
void SED::readSEDBin(ifstream &ins) {
  ins.read((char *)&nummod, sizeof(int));
  name = "MOD_" + to_string(nummod);

  long nbw;
  // second line with the flux
  ins.read((char *)&nbw, sizeof(long));

  if (nbw <= 0) {
    throw runtime_error("SED::readSEDBin(ifstream& ins) lamb flux is zero");
  }

  lamb_flux.resize(nbw, oneElLambda(-99, -99, 1));
  // Read the wavelength
  for (auto &oneEl : lamb_flux) {
    ins.read((char *)&oneEl.lamb, sizeof(double));
  }
  // Read the flux
  for (auto &oneEl : lamb_flux) {
    ins.read((char *)&oneEl.val, sizeof(double));
  }

  return;
}

/*
  Check that we can integrate within a filter. If a risk, issue a warning, then
  extrapolate
*/
void SED::warning_integrateSED(const vector<flt> &filters, bool verbose) {
  // Loop over the filters
  for (const auto &filter : filters) {
    if (((lamb_flux.begin())->lamb) * (1. + red) > filter.lmin()) {
      // if(verbose){
      // cout << "A problem could occur since minimum of SED " <<
      // (lamb_flux.begin())->lamb << " above min of the filter " <<
      // filter.lmin() ; cout << " with filters bluer than " << filter.name << "
      // and SED " << name << " and z " << red << "." ; cout << " Add lambda=0 ;
      // flux=0 to extralolate in blue.  " << endl;
      //}
      // Put the extreme value at 0, in order to define the SED in the blue part
      lamb_flux.emplace(lamb_flux.begin(), 0, 0, 1);
    }

    if (((lamb_flux.end() - 1)->lamb) * (1. + red) < filter.lmax()) {
      if (verbose && (red == 0)) {
        cout << "A problem could occur since maximum of SED "
             << lamb_flux.back().lamb << " below max of the filter "
             << filter.lmax();
        cout << " with filters redder than " << filter.name << " and SED "
             << name << " and z " << red << ".";
        cout << " Add lambda=1.e8 ; flux=0 to extralolate in red. Really "
                "risky: check templates. linear extrapolation not physical. "
             << endl;
      }
      // Put the extreme value at lambda=1e8 and flux=0, in order to define the
      // SED in the red part This is a linear extrapolation from the last point
      // defined in the SED. The extrapolation should be done in the template
      // itself, with a physical meaning. Need to avoid such situation.
      lamb_flux.emplace_back(1.e8 * (1. + red), 0, 1);
    }
  }
}

double SED::integrate(const double lmin, const double lmax) {
  // restrict to cases where the SED is defined over the
  // whole range
  if (lamb_flux.front().lamb > lmin || lamb_flux.back().lamb < lmax) {
    return INVALID_VAL;
  }

  auto up =
      lower_bound(lamb_flux.begin(), lamb_flux.end(), oneElLambda(lmin, 1., 1));
  size_t j = std::distance(lamb_flux.begin(), up) - 1;

  // integrate from lmin to lamb_flux[j+1]
  double x1 = lamb_flux[j].lamb;
  double x2 = lamb_flux[j + 1].lamb;
  double y1 = lamb_flux[j].val;
  double y2 = lamb_flux[j + 1].val;

  double slope = (y2 - y1) / (x2 - x1);
  double interp = y1 + slope * (lmin - x1);
  double res = (y2 + interp) * 0.5 * (x2 - lmin);
  size_t lastidx = j + 1;

  // #pragma omp parallel for reduction(+:res)
  for (size_t i = j + 1; i < lamb_flux.size() - 1; i++) {
    // if(lamb_flux[i].lamb<lmin) continue;
    if (lamb_flux[i + 1].lamb >= lmax) {
      lastidx = i;
      break;
    }
    double fmean = (lamb_flux[i].val + lamb_flux[i + 1].val) / 2;
    double dlbd = (lamb_flux[i + 1].lamb - lamb_flux[i].lamb);
    res += dlbd * fmean;
  }
  // integrate from lamb_flux[lastidx] to lmax
  x1 = lamb_flux[lastidx].lamb;
  x2 = lamb_flux[lastidx + 1].lamb;
  y1 = lamb_flux[lastidx].val;
  y2 = lamb_flux[lastidx + 1].val;

  slope = (y2 - y1) / (x2 - x1);
  interp = y1 + slope * (lmax - x1);
  res += (y1 + interp) * 0.5 * (lmax - x1);

  return res;
}

// Integral of lamb_flux with the trapezoidal method
double SED::trapzd() {
  double s = 0.0;
  for (size_t k = 0; k < lamb_flux.size() - 1; k++) {
    s += (lamb_flux[k + 1].val + lamb_flux[k].val) * 0.5 *
         (lamb_flux[k + 1].lamb - lamb_flux[k].lamb);
  }
  return s;
}

vector<oneElLambda> SED::resample(vector<oneElLambda> &lamb_all,
                                  const int origine, const double lmin,
                                  const double lmax) {
  vector<oneElLambda> lamb_interp;
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
      if (prevEl.lamb >= 0) {
        // Increase the iterator to reach the next item with the right origin
        for (vector<oneElLambda>::iterator kt = it; kt < lamb_interp.end();
             ++kt) {
          // Store the upper value for the interpolation
          if (kt->ori == origine) {
            nextEl = *kt;
            break;
          }

          // Case when the final position in the vector is examined without
          // finding any item with the right origin Not possible to find an
          // upper value for the interpolation
          if (kt == lamb_interp.end() - 1) nextEl.lamb = -999;
        }

        // Change the vector as the result of the interpolation
        // Mofify the origin variable
        it->interp(prevEl, nextEl);
        // If the interpolation was possible
        if (it->ori >= 0) it->ori = origine;
      } else {
        // If the interpolation wasn't possible
        it->val = -99;
        it->ori = -99;
      }
    }
  }
  return lamb_interp;
}

/*
  Create a SED corresponding to the calibration studied, used for filter
*/
void SED::generateCalib(double lmin, double lmax, int Nsteps, int calib) {
  // Define a lambda range of Nsteps steps between lambda min and max
  double dl = (lmax - lmin) / double(Nsteps);
  double lamb, flux;

  lamb_flux.clear();

  // Define the SED depending on the keyword "CALIB"
  SED calibSED("calib");
  // Loop over the 1000 steps
  for (int k = 0; k <= Nsteps; ++k) {
    // Step in lambda
    lamb = lmin + dl * double(k);

    switch (calib) {
      case 1:
        // reference function: 1/lambda (nuBnu=cst)
        flux = 1. / lamb;
        break;
      case 2:
        // reference function: 1/lambda^3 (Bnu=nu)
        flux = 1. / pow(lamb, 3.);
        break;
      case 3:
        // reference function: black body at 10000K
        flux = blackbody(10000, lamb);
        break;
      case 4:
        // reference function: black body at 10000K
        flux = blackbody(10000, lamb);
        break;
      case 5:
        // reference function: 1/lambda^3 (Bnu=nu)
        flux = 1. / pow(lamb, 3.);
        break;
      default:
        // reference function: 1/lambda^2 (Bnu=cst)
        flux = 1. / pow(lamb, 2.);  // Standard case 0 and others >5
        break;
    }

    // fill the lambda/trans values of the object
    lamb_flux.emplace_back(lamb, flux, 1);
  }

  return;
}

void SED::sumSpectra(SED addSED, const double rescal) {
  // no need to work if condition applies
  if (addSED.lamb_flux.empty() || rescal == 0.) return;

  // we need to interpolate both SED, as sed1+rescal*sed2
  // is supposed to be defined on the union of the sed1 and
  // sed2 domains, with sed1 or sed2=0 where they are not defined
  auto [x1, y1] = to_tuple(lamb_flux);
  auto [x2, y2] = to_tuple(addSED.lamb_flux);
  // addSED sampled at SED position
  auto newy2 = fast_interpolate(x2, y2, x1, 0);
  // SED sampled at addSED position
  auto newy1 = fast_interpolate(x1, y1, x2, 0);
  oneElVector allpoints;
  allpoints.reserve(x1.size() + x2.size());
  for (size_t k = 0; k < x1.size(); k++) {
    allpoints.emplace_back(x1[k], y1[k] + rescal * newy2[k], 1);
  }
  for (size_t k = 0; k < x2.size(); k++) {
    allpoints.emplace_back(x2[k], newy1[k] + rescal * y2[k], 1);
  }
  std::sort(allpoints.begin(), allpoints.end());
  allpoints.erase(std::unique(allpoints.begin(), allpoints.end()),
                  allpoints.end());

  lamb_flux.clear();
  lamb_flux = allpoints;
  return;
}

/*
  Limit the size of the data in memory
*/
void SED::reduce_memory(vector<flt> allFlt) {
  // Only in the case z>0 and lamb_flux is defined
  if ((red >= 1.e-20) && (lamb_flux.size() > 0)) {
    vector<bool> flags;
    flags.resize(lamb_flux.size(), false);

    // Loop over the filters, and set the origin at -1 when encompassed in the
    // filter curve
    for (const auto &flt : allFlt) {
      // Loop over all the lambda of the SED
      for (size_t k = 1; k < lamb_flux.size() - 1; ++k) {
        // Set the flag if the current SED point is within the filter support,
        // Keep the value just before the filter
        if ((lamb_flux[k].lamb <= flt.lmin()) &&
            (lamb_flux[k + 1].lamb >= flt.lmin()))
          flags[k] = true;
        // Keep the values inside the filter
        if ((lamb_flux[k].lamb >= flt.lmin()) &&
            (lamb_flux[k].lamb <= flt.lmax()))
          flags[k] = true;
        // Keep the value just after the filter
        if ((lamb_flux[k - 1].lamb <= flt.lmax()) &&
            (lamb_flux[k].lamb >= flt.lmax()))
          flags[k] = true;
      }
    }

    // Keep the two extreme for interpolation
    flags.front() = true;
    flags.back() = true;

    // Loop over the SED and push in reducedSED only the part below the filter
    // curves
    vector<oneElLambda> new_vec;
    for (size_t k = 0; k < lamb_flux.size(); ++k) {
      if (flags[k]) {
        // Remove non flagged entries
        new_vec.push_back(lamb_flux[k]);
      }
    }
    lamb_flux.assign(new_vec.begin(), new_vec.end());
  }
  return;
}

void SED::generate_spectra(double zin, double dmin, vector<opa> opaAll) {
  if (nlib == 0) {
    // GAL
    red = zin;
    // Create a SED of emission lines
    GalSED SEDz0_Em(*this);
    // Don't keep the original lamb_flux, replace it with emission lines.
    SEDz0_Em.lamb_flux.clear();
    SEDz0_Em.generateEmSpectra(40);
    // Add them to the continuum
    sumSpectra(SEDz0_Em, 1.0);
    // Opacity applied in rest-frame, depending on the redshift of the source
    applyOpa(opaAll);
    // Redshift the SED at the redshift extracted from the minimum chi2
    redshift();
    // Rescale la SED
    rescale(dmin);
  } else if (nlib == 1) {
    // QSO
    red = zin;
    // Rescale la SED
    rescale(dmin);
    // Opacity applied in rest-frame, depending on the redshift of the source
    applyOpa(opaAll);
    //// Redshift the SED at the redshift extracted from the minimum chi2
    redshift();
  } else if (nlib == 2) {
    // STAR
    //  Rescale la SED
    rescale(dmin);
  }
  return;
}

/*
  redshift the SED
  input: rescaling factor
  output: void, modify the internal variable lamb_flux
*/
void SED::redshift() {
  // Loop over the full vector with the two vectors concatenated
  for (vector<oneElLambda>::iterator it = lamb_flux.begin();
       it < lamb_flux.end(); ++it) {
    // lambda(1+z) and flux/(1+z)
    it->val = (it->val) / (1 + red);
    it->lamb = (it->lamb) * (1 + red);
  }

  return;
}

void SED::apply_extinction(const double ebv, const ext &oneext) {
  (*this).ebv = ebv;
  // No need to loose time if E(b-V)~0
  if (ebv <= 1.e-20) return;
  // if empty spectrum return
  if (lamb_flux.empty()) return;

  // the extinction is the (x,y) curve
  auto [x, y] = to_tuple(oneext.lamb_ext);
  auto [z, t] = to_tuple(lamb_flux);
  // interpolate the extinction curve at the position of the SED lambdas
  // extrapolations set to 0 so that the exponential in the for loop is 1
  auto newext_val = fast_interpolate(x, y, z, 0);

  for (size_t k = 0; k < z.size(); ++k) {
    double val = t[k];
    double lamb = z[k];
    double ext_val = newext_val[k];
    // Change the value of the flux according to the dust extinction
    lamb_flux[k].val *= pow(10., -0.4 * ebv * ext_val);
  }

  // Indicate what is the value of the ebv and the index of the extinction law
  extlawId = oneext.numext;
  return;
}

/*
  Apply dust extinction on the emission lines (fac_line)
  Only for galaxies and QSO
*/
void SED::applyExtLines2(double ebv, const ext &oneext) {
  // If you want to use the redshift dependency of the attenuation line versus
  // continuum F=F(z=0)+a*z If such option is re-activated, absolutely need to
  // change the code in read_lib and onesource.cpp
  // double a=0.2;
  // No redshift dependency
  double a = 0.;

  // needed to establish the escape fraction.
  // Relation by Hayes et al 2011 :
  // https://iopscience.iop.org/article/10.1088/0004-637X/730/1/8
  double CLya = 0.445;
  double kLya = 13.8;

  if (ebv <= 1.e-20) return;

  // Interpolate the extinction curve to the line positions
  auto [extx, exty] = to_tuple(oneext.lamb_ext);
  for (size_t i = 0; i < NUM_EMISSION_LINES; i++) {
    double line_lamb = fac_line[i].lamb;
    double line_val = fac_line[i].val;
    double attenuation = interp_linear_point(extx, exty, line_lamb);

    double f = (ebvFac[i] + a * red);
    if (f > 1) f = 1.;

    line_val *= pow(10., -0.4 * ebv / f * attenuation);
    // For the Lya line one must also account for escape.
    //  Relation by Hayes et al 2011 to derive the escape fraction
    if (line_lamb > 1215. && line_lamb < 1216.) {
      // we could condition on i=0 as Lyalpha is the first line
      double f_esc_Lya = CLya * pow(10., (-0.4 * ebv / f * kLya));
      line_val *= f_esc_Lya;
    }
    fac_line[i].val = line_val;
  }
  return;
}

/*
  Apply dust extinction on the emission lines (fac_line)
  Only for galaxies and QSO
*/
void SED::applyExtLines(double ebv, const ext &oneext) {
  // Local value from Calzetti
  // double ebvFac[65] =
  // {0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44,0.44
  // }; use the same attenuation for stellar continuum and emission lines
  double ebvFac[65] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                       1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  // Check size
  int longEL = sizeof(ebvFac) / sizeof(ebvFac[0]);
  if (longEL != 65) cout << " Error with size of ebv lines " << longEL << endl;
  // If you want to use the redshift dependency of the attenuation line versus
  // continuum F=F(z=0)+a*z If such option is re-activated, absolutely need to
  // change the code in read_lib and onesource.cpp
  // double a=0.2;
  // No redshift dependency
  double a = 0.;

  // needed to establish the escape fraction.
  // Relation by Hayes et al 2011 :
  // https://iopscience.iop.org/article/10.1088/0004-637X/730/1/8
  double CLya = 0.445;
  double kLya = 13.8;

  if (ebv <= 1.e-20) return;

  // work with fac_line as spectra
  vector<oneElLambda> line_all = fac_line;

  // Concatenate two vectors composed of "oneElLambda" including the
  // extinction law and the SED
  line_all.insert(line_all.end(), oneext.lamb_ext.begin(),
                  oneext.lamb_ext.end());

  // Sort the vector in increasing lambda (vector including the extinction law
  // and the SED)
  sort(line_all.begin(), line_all.end());

  // Resample the extinction in a lambda range combining the SED and the
  // extinction law
  vector<oneElLambda> new_ext = resample(line_all, 2, 0., 1e20);

  // Loop over the SED and extinction curves using the concatenate ext+SED
  // vector
  vector<oneElLambda> lamb_new;
  int l = 0;
  for (size_t k = 0; k < line_all.size(); ++k) {
    // If we are well looking at the original lines
    if ((line_all[k]).ori == 5) {
      double lamb = line_all[k].lamb;
      double val = line_all[k].val;
      if (oneext.lmin < lamb && oneext.lmax > lamb) {
        // Change the value of the flux according to the duts extinction
        double f = (ebvFac[l] + a * red);
        if (f > 1) f = 1.;
        // Check that the extinction is well defined, other put 0
        if (new_ext[k].ori < 0) new_ext[k].val = 0.;
        val = (line_all[k]).val * pow(10., (-0.4 * ebv / f * (new_ext[k]).val));
        // Relation by Hayes et al 2011 to derive the escape fraction
        if (line_all[k].lamb > 1215. && line_all[k].lamb < 1216.) {
          double f_esc_Lya = CLya * pow(10., (-0.4 * ebv / f * kLya));
          val = val * f_esc_Lya;
        }
      }
      // Count what is the line according to the original fac_line
      l++;
      // Store the result into a new lamb_flux vector
      lamb_new.emplace_back(lamb, val, 1);
    }
  }

  // Replace the old lambda flux vector with the new one including extinction
  fac_line.clear();
  fac_line = lamb_new;
  lamb_new.clear();
  new_ext.clear();

  return;
}

/*
  product of the sed with the opacity of the IGM
  only for galaxies and QSO
*/
void SED::applyOpa(const vector<opa> &opaAll) {
  // Select the right opacity file according to the redshift of the source
  // Do not use simply 0.1 -> need to have the same as the fortran version and
  // some rounding problem (e.g. z=0.15 with step=0.01)
  int ind = lround(red / 0.100000000000001);
  // In lepharedir/opa/ there are opacity files up to z=8 (file #80),
  // use the last one (transmission=0 at all lambda<1216) also when z>8
  if (ind > 80) {
    ind = 80;
  }
  auto opa = opaAll[ind];

  // No need to do anything above the high limit of the intergalactic
  // attenuation curves
  auto limit = lower_bound(lamb_flux.begin(), lamb_flux.end(), 1215.67);

  oneElVector temp(lamb_flux.begin(), limit);
  // the opacity is the (x,y) curve
  auto [x, y] = to_tuple(opa.lamb_opa);
  auto [z, t] = to_tuple(temp);
  // interpolate the opacity curve at the position of the SED lambdas
  // extrapolations set to 1 so that the input spectrum is left untouched
  // by multiplication
  auto newopa = fast_interpolate(x, y, z, 1);

  // directly modify the lamb_fluxes values for lambda<1215.67,
  // which is z.size()
  for (size_t k = 0; k < z.size(); ++k) {
    lamb_flux[k].val *= newopa[k];
  }
  return;
}

/*
  GALAXIES
*/

/*
  Constructor
*/
GalSED::GalSED(const string nameC, int nummodC) : SED(nameC, nummodC, "GAL") {
  format = 'A';  // ascii by default

  d4000 = -999;  // D4000 = Int4050-4250 Fl / Int 3750_3950 dL   dimensionless
  tau = -999;    // Tau (yr)   only for SFH as exp-t/Tau

  zmet = -999;  // Metallicity
  lnir = -999;  // Int_2.1um^2.3um L_lbda.dlbda        in Log unit erg/s/Hz
  luv = -999;   // Int_0.21um^0.25um L_lbda . dlbda    in Log unit erg/s/Hz
  lopt = -999;  // Int_0.55um^0.65um L_lbda . dlbda    in Log unit erg/s/Hz

  fracEm = 1.;  // fraction of the emmission line considered
}

/*
  extended constructor
*/
GalSED::GalSED(const string nameC, double tauC, double ageC, string formatC,
               int nummodC, int idAgeC)
    : SED(nameC, tauC, ageC, nummodC, "GAL", idAgeC) {
  format = formatC;
  tau = tauC;
  d4000 = -999;
  zmet = -999;
  lnir = -999;
  luv = -999;
  lopt = -999;
  fracEm = 1.;  // fraction of the emmission line considered
}

GalSED GalSED::generateEmSED(const string &emtype) {
  // Only if emission lines
  GalSED oneEm("");
  if (emtype[0] == 'P') {
    // new method to include emission lines, with physical recipes
    auto nebular_contribution = add_neb_cont(qi[2]);  // Compute the continuum
    oneEm.generateEmPhys(zmet, qi[2]);  // Generate the emission lines
  } else if (emtype.compare("EMP_UV") == 0) {
    // Empirical method for emission lines
    // Need the SED M(NUV) and NUV_R before applying extinction to get the
    // emission lines
    double MNUV_int = -2.5 * luv + 51.59;
    double nuvr = -2.5 * (luv - lopt);
    oneEm.generateEmEmpUV(MNUV_int, nuvr);  // Generate the emission lines
  } else if (emtype.compare("EMP_SFR") == 0) {
    // Empirical method for emission lines
    // Need the SED SFR
    double nuvr = -2.5 * (luv - lopt);
    oneEm.generateEmEmpSFR(sfr, nuvr);  // Generate the emission lines
  }
  return oneEm;
}

/*
  Calculate number of continuum photons shortward of ionising edges.
*/
void GalSED::calc_ph() {
  // ################ ionisation edges  ################################

  // Definition of ionisations edges

  // Sorted by increasing order for the integration
  // hc = const [eV.A] (globals.cpp);
  double wedge[4];        // tab de 4 edges
  wedge[0] = hc / 54.42;  // HeII edge in Angstroem
  wedge[1] = hc / 24.59;  // HeI edge in Angstroem
  wedge[2] = hc / 13.60;  // H edge in Angstroem
  wedge[3] = 1108.7;  // H_2 excitation from ground state (B-X) ############ A
                      // vérifier (11.18 eV)############

  // Definition :
  // qi[4] define in SED.h = ionising fluxes  for 4 elements [#photons.cm-2/s]

  double qsum = 0;  // sum var
  int iedge = 0;    // incrementation var of wedge[]
  int i = 1;        // we need (i-1) to integrate

  // Start the integration (trapezium method)
  do {
    // Integrate and multiply by lambda/hc energy (1.6022e-12 to convert from
    // e.V to erg)
    qsum = qsum + (lamb_flux[i].val + lamb_flux[i - 1].val) * 0.5 *
                      (lamb_flux[i].lamb - lamb_flux[i - 1].lamb) *
                      lamb_flux[i].lamb / (hc * 1.6022e-12);
    // if wave is <= to the edge wave
    if (lamb_flux[i].lamb <= wedge[iedge]) {
      qi[iedge] = qsum;  // fill qi
    } else {
      iedge++;  // next qi
    }
    i++;  // next wavelength
  } while (lamb_flux[i].lamb <= wedge[3]);  // limite

  // cout << "nombre de photons ionisants : \n pour HeII : " << qi[0]
  //      << "\n pour HeI : " << qi[1]
  //      << "\n pour H (Lyman-continuum flux) : " << qi[2]
  //      << "\n pour H_2 : " << qi[3] << endl;

  return;
}

/*
  Add the continuum from the nebular regions
  Work done by Cedric Dubois
*/
vector<double> GalSED::add_neb_cont(double qi) {
  /* we assume that the emitting gas has an electron temperature of Te = 10000
     K, an electron density N = 100 cm-3 (low density limite), and a helium
     abundance of 10% by number relative to hydrogen.
  */

  oneElVector ga;
  for (size_t i = 0; i < 71; i++) {
    double val = ga_total[i];
    val = val <= 0 ? -100 : log10(val);
    ga.emplace_back(ga_lamb[i], val, 4);
  }

  ga.insert(ga.end(), lamb_flux.begin(), lamb_flux.end());
  sort(ga.begin(), ga.end());
  auto ga_interp = resample(ga, 4, 0, 1.6e+6);
  for (size_t i = 0; i < ga_interp.size(); i++)
    ga_interp[i].val = pow(10, ga_interp[i].val);

  // Take the vector lamb_flux and add the nebular continu in .val
  double flux_neb;
  vector<double> neb_contrib;
  double cst = c * 1e-40 * qi * f_ga / alpha_B;
  for (size_t i = 0; i < lamb_flux.size(); i++) {
    double x = lamb_flux[i].lamb;
    double v = lamb_flux[i].val;
    // c/lambda^2*gamma/alpha_B * number of ionizing photons * fraction of
    // absorbed photons 1e-40 since c in A/s and gamma in 10^-40 erg
    flux_neb = cst / (x * x) * ga_interp[i].val;
    neb_contrib.push_back(flux_neb);
    // Sum the nebular flux to the original flux from stellar population
    lamb_flux[i].val = v + flux_neb;
  }

  return neb_contrib;
}

/*
  Add the continuum from the nebular regions
  Work done by Cedric Dubois
*/
vector<double> GalSED::add_neb_cont2(double qi) {
  /* we assume that the emitting gas has an electron temperature of Te = 10000
     K, an electron density N = 100 cm-3 (low density limite), and a helium
     abundance of 10% by number relative to hydrogen.
  */

  oneElVector ga;
  // note: this is done repeatedly for each template
  // of the initial set (sedtolib output), while it does not depend on it
  // So this could/should be improved
  for (size_t i = 0; i < NUM_NEBULAR_DATA; i++) {
    double val = ga_total[i];
    //    val = val <= 0 ? -100 : log10(val);
    ga.emplace_back(ga_lamb[i], val, 4);
  }
  auto [x1, y1] = to_tuple(lamb_flux);
  auto [x2, y2] = to_tuple(ga);

  auto x3 = x1;
  x3.insert(x3.end(), x2.begin(), x2.end());
  std::sort(x3.begin(),
            x3.end());  // no unique as the nebular continuum has vertical drops

  lamb_flux.clear();
  lamb_flux.reserve(x3.size());

  double scale = c * 1.e-40 / alpha_B * qi * f_ga;
  double flux_neb;
  vector<double> neb_contrib;
  for (double x : x3) {
    double v1 = interp_linear_point(x1, y1, x);
    double v2 = interp_linear_point(x2, y2, x);
    //    flux_neb = scale / (x * x) * pow(10, v2);
    flux_neb = scale / (x * x) * v2;
    neb_contrib.push_back(flux_neb);
    lamb_flux.emplace_back(x, v1 + flux_neb, 1);
  }
  // return for control
  return neb_contrib;
}

/*
  Generate the emission lines in rest-frame without renormaliation or extinction
  based on the physical prescriptions
  Work done by Cedric Dubois
*/
void GalSED::generateEmPhys(double zmet, double qi) {
  // reference luminosity (H_beta in case B) :
  // L_ref = f_ga*qi*h*c/lamb_Hbeta*alpha_Hbeta/alpha_B
  // Osterbrock data (H. Krueger et al, 1995)
  double L_ref =
      f_ga * qi * 4.780e-13;  // [erg.s-1.cm-2], different value than Schearer
                              // because of the different alpha_B value

  // Select the right table depending on the metallicity Zmet
  // Multiply by the reference luminosity
  fac_line.clear();
  for (size_t i = 0; i < 65; i++) {
    // Find the line ratio
    if (zmet / 0.02 <= 0.03) {
      fac_line.emplace_back(emission_lines[i], L_ref * Z1_line[i], 5);
    } else if (zmet / 0.02 > 0.03 && zmet / 0.02 <= 0.3) {
      fac_line.emplace_back(emission_lines[i], L_ref * Z2_line[i], 5);

    } else if (0.3 < zmet / 0.02) {
      fac_line.emplace_back(emission_lines[i], L_ref * Z3_line[i], 5);
    } else {
      cout << "enter a valid value of Z metalicity " << endl;
    }
  }

  // sort the vector by wavelength
  sort(fac_line.begin(), fac_line.end());

  return;
}

//---------------------------------------------------------------------------------------------------------------------

/*
  Generate the emission lines in rest-frame without renormaliation or extinction
  based on the Kenicutt relation and the UV light
*/
void GalSED::generateEmEmpUV(double MNUV_int, double NUVR) {
  // Rescale the lines according to the Kennicutt 2012 relation, only for blue
  // galaxies Using Halpha as a reference.
  // 2.85 is the intrinsic Halpha to Hbeta flux line ratio
  double scaleFac = 0;
  // correct only for blue enough galaxies
  if (NUVR < 4) scaleFac = pow(10., -0.4 * MNUV_int - 6.224) / 2.85;

  // Create the emission line and rescale them
  fac_line.clear();
  for (size_t i = 0; i < NUM_EMISSION_LINES; i++) {
    fac_line.emplace_back(emission_lines[i], scaleFac * empirical_ratio[i], 5);
  }
  // sort the vector by wavelength
  sort(fac_line.begin(), fac_line.end());
}

/*
  Generate the emission lines in rest-frame without renormaliation or extinction
  based on the Kenicutt relation and the SFR
*/
void GalSED::generateEmEmpSFR(double sfr, double NUVR) {
  // Rescale the lines according to the Kennicutt 2012 relation, only for blue
  // galaxies Using Halpha as a reference
  // 2.85 is the intrinsic Halpha to Hbeta flux line ratio
  double scaleFac = 0;
  // correct only for blue enough galaxies
  if (NUVR < 4)
    scaleFac = pow(10., log10(sfr) + 41.27 -
                            log10(4 * pi * 100 * pow(3.08568, 2)) - 36) /
               2.85;

  // Create the emission line and rescale them
  fac_line.clear();
  for (size_t i = 0; i < 65; i++) {
    fac_line.emplace_back(emission_lines[i], scaleFac * empirical_ratio2[i], 5);
  }
  // sort the vector by wavelength
  sort(fac_line.begin(), fac_line.end());
}

/*
  Generate the spectra based on emission lines
  Adapt the step with the sigma, but emission lines shouldn't be blended
*/
void GalSED::generateEmSpectra(int nstep) {
  // Number of sigma to generate the gaussian.
  double nbsigma = 6;

  lamb_flux.clear();

  // Define sigma of the lines
  double sigma[65];
  for (size_t j = 0; j < fac_line.size(); j++) {
    // Distribute the flux through an effective sigma width
    // no real idea of width to take (depends on inclination, how the HII
    // regions are distributed, etc) Take 100 km/s
    sigma[j] = 100. * fac_line[j].lamb / c / 2.355 * 1.e13;
  }

  // Loop over each emission line considered here in order to define the lambda
  // array in which to distribute the EL line
  for (size_t j = 0; j < fac_line.size(); j++) {
    // Generate the line only of necessary
    if (fac_line[j].val > 0) {
      // step in lambda
      double step = nbsigma * 2. * sigma[j] / double(nstep);
      // Loop to cover the full line
      for (int i = 0; i < nstep; i++) {
        // Lambda corresponding to this index
        double lb = fac_line[j].lamb - nbsigma * sigma[j] + double(i) * step;
        // Add it to the emission line spectra
        // Start with a -nbsigma*sigma put at 0
        lamb_flux.emplace_back(lb, 0., 1);
      }
    }
  }
  // Sort the lambda, necessary in case of overlapp between lines
  sort(lamb_flux.begin(), lamb_flux.end());

  // Loop over each emission line considered here in order to estimate the flux
  // and sum it to the other lines
  for (size_t j = 0; j < fac_line.size(); j++) {
    // Generate the line only of necessary
    if (fac_line[j].val > 0) {
      // Loop over the original lambda-flux
      for (size_t k = 0; k < lamb_flux.size(); k++) {
        double lineFl = (fac_line[j].val) / (sigma[j] * sqrt(2. * pi)) *
                        exp(-pow((lamb_flux[k].lamb - fac_line[j].lamb), 2.) /
                            (2 * pow(sigma[j], 2.)));
        lamb_flux[k].val += lineFl;
      }
    }
  }

  // Add 0 before/after each line
  if (lamb_flux.size() > 0) {
    vector<oneElLambda> insert_lamb;
    // Loop over all the line elements
    for (size_t k = 1; k < lamb_flux.size() - 1; k++) {
      // if nothing at more than 1A after the considered lambda, it is the end
      // of one line
      if ((lamb_flux[k].lamb + 1.) < lamb_flux[k + 1].lamb) {
        insert_lamb.emplace_back(lamb_flux[k].lamb + 0.1, 0., 1);
      }
      // if nothing at more than 1A below the considered lambda, it is the
      // beginning of a line
      if ((lamb_flux[k - 1].lamb + 1.) < lamb_flux[k].lamb) {
        insert_lamb.emplace_back(lamb_flux[k].lamb - 0.1, 0., 1);
      }
    }
    // Add val=0 at extreme wavelengths
    insert_lamb.emplace_back(0, 0., 1);
    insert_lamb.emplace_back(10000000., 0., 1);

    // Add this vector to the first one
    lamb_flux.insert(lamb_flux.end(), insert_lamb.begin(), insert_lamb.end());
    // Sort the lambda, necessary in case of overlapp between lines
    sort(lamb_flux.begin(), lamb_flux.end());
  }

  return;
}

pair<vector<double>, vector<double>> SED::get_data_vector(double minl,
                                                          double maxl, bool mag,
                                                          double offset) {
  vector<double> lambs, vals;
  double lamb, val;
  for (const auto &onel : lamb_flux) {
    lamb = onel.lamb;
    val = onel.val;
    if (lamb <= minl || lamb >= maxl) continue;

    lambs.push_back(lamb);

    if (mag) {
      val = val <= 0.0 ? HIGH_MAG : flux2mag(val * lamb * lamb / c, offset);
    }

    vals.push_back(val);
  }
  return make_pair(lambs, vals);
}

void GalSED::compute_luminosities() {
  double fluxconv = (4 * pi * 100 * pow(pc, 2));

  // construct a heaviside filter between 0.21 and 0.25 micron
  // Integrate the SED in NUV : between 0.21 and 0.25 micron
  luv = this->integrate(2100., 2500.);
  if (luv > 0) luv = log10(luv * pow(2300, 2) / 400 / c * fluxconv);

  // Integrate the SED in optical : between 0.55 and 0.65 micron
  lopt = this->integrate(5500., 6500.);
  if (lopt > 0) lopt = log10(lopt * pow(6000, 2.) / 1000. / c * fluxconv);

  // Integrate the SED in NUV : between 2.1 and 2.3 micron
  lnir = this->integrate(21000., 23000.);
  if (lnir > 0) lnir = log10(lnir * pow(22000, 2) / 2000 / c * fluxconv);

  // Integrate the SED before and after the Balmer break
  double lBalm1 = this->integrate(3750., 3950.);
  double lBalm2 = this->integrate(4050., 4250.);
  if (lBalm1 > 0 && lBalm2 > 0) d4000 = lBalm2 / lBalm1;

  // Integrate the SED in IR : between 8 and 1000 micron
  ltir = this->integrate(80000., 10000000.);
  if (ltir > 0) ltir = log10(ltir / Lsol * fluxconv);

  // compute the number of available ionizing photons at edge of
  // HeII HeI H and H_2
  this->calc_ph();

  return;
}

/*
  Write the galaxy sed in binary format, the physical parameters linked to the
  sed, the doc
*/
void GalSED::writeSED(ofstream &ofsBin, ofstream &ofsPhys, ofstream &ofsDoc) {
  SED::writeSED(ofsBin, ofsPhys, ofsDoc);

  // Write in the binary file the physical information
  ofsBin.write((char *)&luv, sizeof(double));
  ofsBin.write((char *)&lopt, sizeof(double));
  ofsBin.write((char *)&lnir, sizeof(double));
  ofsBin.write((char *)&ltir, sizeof(double));
  ofsBin.write((char *)&mass, sizeof(double));
  ofsBin.write((char *)&sfr, sizeof(double));
  ofsBin.write((char *)&zmet, sizeof(double));
  ofsBin.write((char *)&tau, sizeof(double));
  ofsBin.write((char *)&d4000, sizeof(double));
  ofsBin.write((char *)&qi[2], sizeof(double));
  ofsBin.write((char *)&age, sizeof(double));

  // Physical parameters in the ascii file
  ofsPhys << setw(12) << age << " ";
  ofsPhys << setw(12) << luv << " ";
  ofsPhys << setw(12) << lopt << " ";
  ofsPhys << setw(12) << lnir << " ";
  ofsPhys << setw(12) << ltir << " ";
  ofsPhys << setw(12) << mass << " ";
  ofsPhys << setw(12) << sfr << " ";
  ofsPhys << setw(12) << zmet << " ";
  ofsPhys << setw(12) << tau << " ";
  ofsPhys << setw(12) << d4000 << " ";
  ofsPhys << setw(12) << qi[2] << endl;

  return;
}

/*
  Read the sed in binary format in the galaxy case
*/
void GalSED::readSEDBin(ifstream &ins) {
  SED::readSEDBin(ins);

  ins.read((char *)&luv, sizeof(double));
  ins.read((char *)&lopt, sizeof(double));
  ins.read((char *)&lnir, sizeof(double));
  ins.read((char *)&ltir, sizeof(double));
  ins.read((char *)&mass, sizeof(double));
  ins.read((char *)&sfr, sizeof(double));
  ins.read((char *)&zmet, sizeof(double));
  ins.read((char *)&tau, sizeof(double));
  ins.read((char *)&d4000, sizeof(double));
  ins.read((char *)&qi[2], sizeof(double));
  ins.read((char *)&age, sizeof(double));

  return;
}

/*
  Write the predicted magnitudes in binary format for the galaxy library
  Result of mag_gal
*/
void GalSED::writeMag(bool outasc, ofstream &ofsBin, ofstream &ofsDat,
                      vector<flt> allFilters, string magtyp) const {
  // number of filters
  int nbFlt = mag.size();

  ofsBin.write((char *)&luv, sizeof(double));
  ofsBin.write((char *)&lopt, sizeof(double));
  ofsBin.write((char *)&lnir, sizeof(double));
  ofsBin.write((char *)&ltir, sizeof(double));
  ofsBin.write((char *)&mass, sizeof(double));
  ofsBin.write((char *)&sfr, sizeof(double));
  ofsBin.write((char *)&zmet, sizeof(double));
  ofsBin.write((char *)&tau, sizeof(double));
  ofsBin.write((char *)&d4000, sizeof(double));

  // Write the library in a binary file
  ofsBin.write((char *)&nummod, sizeof(int));      // Index du template
  ofsBin.write((char *)&extlawId, sizeof(int));    // type du model
  ofsBin.write((char *)&ebv, sizeof(double));      // E(B-V)
  ofsBin.write((char *)&fracEm, sizeof(double));   // fracEm index
  ofsBin.write((char *)&red, sizeof(double));      // redshift
  ofsBin.write((char *)&distMod, sizeof(double));  // distance modulus
  ofsBin.write((char *)&age, sizeof(double));      // age
  ofsBin.write((char *)&nbFlt, sizeof(int));       // Number of filters

  // Write the predicted magnitudes
  for (int k = 0; k < nbFlt; k++)
    ofsBin.write((char *)&(mag[k]), sizeof(double));
  // Write the k-correction
  for (int k = 0; k < nbFlt; k++)
    ofsBin.write((char *)&(kcorr[k]), sizeof(double));
  // Write the emission lines
  if (has_emlines) {
    // Write the flux predicted into the filters
    for (int k = 0; k < nbFlt; k++)
      ofsBin.write((char *)&(flEm[k]), sizeof(double));
    // Don't write the flux at z>0, to save RAM. With such cleaning, no z
    // dependency on line ratio can be implemented.
    if (red < 1.e-20) {
      // Write the line fluxes
      int nbEm = fac_line.size();
      ofsBin.write((char *)&(nbEm), sizeof(int));
      for (int k = 0; k < nbEm; k++)
        ofsBin.write((char *)&(fac_line[k].lamb), sizeof(double));
      for (int k = 0; k < nbEm; k++)
        ofsBin.write((char *)&(fac_line[k].val), sizeof(double));
    }
  }

  // Write the spectra only if the redshift is 0, for the output .spec
  if (red < 1.e-20) {
    // Write the continuum spectra
    int nbLamb = lamb_flux.size();
    ofsBin.write((char *)&(nbLamb), sizeof(int));
    for (int k = 0; k < nbLamb; k++)
      ofsBin.write((char *)&(lamb_flux[k].lamb), sizeof(double));
    for (int k = 0; k < nbLamb; k++)
      ofsBin.write((char *)&(lamb_flux[k].val), sizeof(double));
  }

  // Case with the ascii file
  if (outasc) {
    // Write output
    ofsDat << setw(6) << nummod << " ";
    // start the numbering of attenuation curves at 1 in output
    ofsDat << setw(3) << extlawId + 1 << " ";
    ofsDat << setw(3) << ebv << " ";
    ofsDat << setw(12) << fracEm << " ";
    ofsDat << setw(5) << red << " ";
    ofsDat << setw(12) << distMod << " ";
    ofsDat << setw(12) << age << " ";
    ofsDat << setw(4) << nbFlt << " ";
    if (magtyp[0] != 'V') {  // AB
      for (int k = 0; k < nbFlt; k++) {
        ofsDat << setw(6) << mag[k] << " ";
      }
    } else {  // VEGA
      for (int k = 0; k < nbFlt; k++) {
        ofsDat << setw(6) << mag[k] + allFilters[k].ab << " ";
      }
    }
    for (int k = 0; k < nbFlt; k++) {
      ofsDat << setw(6) << kcorr[k] << " ";
    }
    if (has_emlines) {
      for (int k = 0; k < nbFlt; k++) {
        ofsDat << setw(6) << flEm[k] << " ";
      }
    }
    ofsDat << endl;
  }

  return;
}

/*
  Read the predicted magnitudes in binary format
  Case for the galaxies
  Used in zphota
*/
void GalSED::readMagBin(ifstream &ins) {
  int nbFlt;

  ins.read((char *)&luv, sizeof(double));
  ins.read((char *)&lopt, sizeof(double));
  ins.read((char *)&lnir, sizeof(double));
  ins.read((char *)&ltir, sizeof(double));
  ins.read((char *)&mass, sizeof(double));
  ins.read((char *)&sfr, sizeof(double));
  ins.read((char *)&zmet, sizeof(double));
  ins.read((char *)&tau, sizeof(double));
  ins.read((char *)&d4000, sizeof(double));

  // Read information conserning the SED which is read
  ins.read((char *)&nummod, sizeof(int));
  ins.read((char *)&extlawId, sizeof(int));
  ins.read((char *)&ebv, sizeof(double));
  ins.read((char *)&fracEm, sizeof(double));
  ins.read((char *)&red, sizeof(double));
  ins.read((char *)&distMod, sizeof(double));
  ins.read((char *)&age, sizeof(double));
  ins.read((char *)&nbFlt, sizeof(int));

  // define the ssfr
  if (mass > 0) {
    ssfr = sfr / mass;
  } else {
    ssfr = -999.;
  }

  // Read the magnitudes
  mag.resize(nbFlt, 99);
  for (auto &m : mag) {
    ins.read((char *)&m, sizeof(double));
  }
  kcorr.resize(nbFlt, 0);
  for (auto &k : kcorr) {
    ins.read((char *)&k, sizeof(double));
  }

  // Read the emission lines fluxes integrated into filters
  if (has_emlines) {
    flEm.resize(nbFlt, 0.);
    for (auto &flem : flEm) {
      ins.read((char *)&flem, sizeof(double));
    }
    // Don't read the lines at z>0, to save RAM. With such cleaning, no z
    // dependency on line ratio can be implemented.
    if (red < 1.e-20) {
      // Read the emission lines fluxes
      int nbEm;
      ins.read((char *)&nbEm, sizeof(int));
      fac_line.resize(nbEm, oneElLambda(-999, -999, 1));
      for (auto &eml : fac_line) {
        ins.read((char *)&eml.lamb, sizeof(double));
      }
      for (auto &eml : fac_line) {
        ins.read((char *)&eml.val, sizeof(double));
      }
    }
  }

  // read the spectra only if the redshift is 0
  if (red < 1.e-20) {
    int nblamb;
    ins.read((char *)&nblamb, sizeof(int));
    // put ori=1 because it's a SED
    lamb_flux.resize(nblamb, oneElLambda(-999, -999, 1));
    for (auto &oneEl : lamb_flux) {
      ins.read((char *)&oneEl.lamb, sizeof(double));
    }
    for (auto &oneEl : lamb_flux) {
      ins.read((char *)&oneEl.val, sizeof(double));
    }
  }

  return;
}

/*
  Sum the emission line flux and the continuum flux
*/
void GalSED::sumEmLines() {
  // Loop over each filter.
  for (size_t k = 0; k < mag.size(); k++) {
    // sum it to the original flux, before having applied any emission line
    double fluxinter = pow(10., -0.4 * (mag[k] + 48.6)) + flEm[k];
    if (fluxinter > 0) {
      mag[k] = -2.5 * log10(fluxinter) - 48.6;
    } else {
      mag[k] = 999.9;
    }
  }
  // Clean the emission line from the SED after its used. Save memory.
  flEm.clear();

  return;
}

/*
  Recompute k-corr. Necessary because of the emission lines contribution not
  included in mag_gal
*/
void GalSED::kcorrec(const vector<double> &magz0) {
  // Loop over each filter.
  for (size_t k = 0; k < mag.size(); k++) {
    // Substract the magnitude at z and the one at z=0
    // kcorr[k]=mag[k]-mag_z0[k]-distMod;
    kcorr[k] = mag[k] - magz0[k] - distMod;
  }

  return;
}

/*
  Rescale the emission lines
*/
void GalSED::rescaleEmLines() {
  // All lines
  double adjust_line[65] = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
                            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
  // Adjust only the OIII doublet
  // double adjust_line[65] =
  // {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  // double adjust_line[65] =
  // {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  // Check
  int longEL = sizeof(adjust_line) / sizeof(adjust_line[0]);
  if (longEL != 65)
    cout << " Error with size of adjust lines " << longEL << endl;

  // In case the redshift is 0, so fac_line is defined and we apply the same
  // fraction
  if (fac_line.size() > 0) {
    for (size_t k = 0; k < fac_line.size(); k++) {
      if (adjust_line[k] == 1) {
        fac_line[k].val = fac_line[k].val * fracEm;
      }
    }
  }

  return;
}

/*
  Apply a z dependence of the emission line ratio
  Choose if the reference is z=0 or z=2.12 with flag
*/
void GalSED::zdepEmLines(int flag) {
  double zdep;

  // multiply the OIII flux by a factor depending on redshift.
  if (flag == 0) {
    // reference at z=0
    // Take a value of 1.5*z+0.9 for the OIIIb/Hb ratio
    if (red < 4) {
      zdep = 1.5 * red + 0.9;
    } else {
      zdep = 6.9;
    }
  } else {
    // We need to fix the reference at z=2.12 to be able to use the physical
    // recipes. We fix OIII/Hbeta at 4.081 at z=2.12 (4.081−0.9)/1.5
    if (red < 4) {
      zdep = (0.3675 * (red - 2.12) + 1);
    } else {
      zdep = 1.69;
    }
  }

  // Change the OIIIa & OIIIb emission lines
  fac_line[21].val = fac_line[21].val * zdep;
  fac_line[22].val = fac_line[22].val * zdep;

  return;
}

/*

  QSO

*/

/*
  Write the magnitudes in binary format for the QSO library
*/
void QSOSED::writeMag(bool outasc, ofstream &ofsBin, ofstream &ofsDat,
                      vector<flt> allFilters, string magtyp) const {
  // Info screen
  // cout << "Compute magnitudes for QSO SEDs " << name << " z " <<  setw(6) <<
  // red << " Ext law " <<  extlawId <<  "  E(B-V) "  << ebv; cout <<  "  \r "
  // << flush;

  // number of filters
  int nbFlt = mag.size();
  // Write the library in a binary file
  ofsBin.write((char *)&nummod, sizeof(int));      // Index du template
  ofsBin.write((char *)&extlawId, sizeof(int));    // type du model
  ofsBin.write((char *)&ebv, sizeof(double));      // E(B-V)
  ofsBin.write((char *)&red, sizeof(double));      // redshift
  ofsBin.write((char *)&distMod, sizeof(double));  // distance modulus
  ofsBin.write((char *)&nbFlt, sizeof(int));       // Number of filters

  // Write the predicted magnitudes
  for (int k = 0; k < nbFlt; k++)
    ofsBin.write((char *)&(mag[k]), sizeof(double));
  // Write the k-corrections
  for (int k = 0; k < nbFlt; k++)
    ofsBin.write((char *)&(kcorr[k]), sizeof(double));

  // Write the spectra only if the redshift is 0
  if (red < 1.e-20) {
    int nbLamb = lamb_flux.size();
    ofsBin.write((char *)&(nbLamb), sizeof(int));
    for (size_t k = 0; k < lamb_flux.size(); k++)
      ofsBin.write((char *)&(lamb_flux[k].lamb), sizeof(double));
    for (size_t k = 0; k < lamb_flux.size(); k++)
      ofsBin.write((char *)&(lamb_flux[k].val), sizeof(double));
  }

  // Case with the ascii file
  if (outasc) {
    // Write output
    ofsDat << setw(6) << nummod << " ";
    // start the numbering of attenuation curves at 1 in output
    ofsDat << setw(3) << extlawId + 1 << " ";
    ofsDat << setw(3) << ebv << " ";
    ofsDat << setw(5) << red << " ";
    ofsDat << setw(12) << distMod << " ";
    ofsDat << setw(4) << nbFlt << " ";
    if (magtyp[0] != 'V') {  // AB
      for (int k = 0; k < nbFlt; k++) {
        ofsDat << setw(6) << mag[k] << " ";
      }
    } else {  // VEGA
      for (int k = 0; k < nbFlt; k++) {
        ofsDat << setw(6) << mag[k] + allFilters[k].ab << " ";
      }
    }
    for (int k = 0; k < nbFlt; k++) {
      ofsDat << setw(6) << kcorr[k] << " ";
    }
    ofsDat << endl;
  }

  return;
}

/*
  Read the predicted magnitudes in binary format
  Case for the galaxies
  Used in zphota
*/
void QSOSED::readMagBin(ifstream &ins) {
  int nbFlt;

  // Read information conserning the SED which is read
  ins.read((char *)&nummod, sizeof(int));
  ins.read((char *)&extlawId, sizeof(int));
  ins.read((char *)&ebv, sizeof(double));
  ins.read((char *)&red, sizeof(double));
  ins.read((char *)&distMod, sizeof(double));
  ins.read((char *)&nbFlt, sizeof(int));

  // Read the magnitudes
  mag.assign(nbFlt, HIGH_MAG);
  for (auto &m : mag) {
    ins.read((char *)&m, sizeof(double));
  }
  kcorr.assign(nbFlt, 0);
  for (auto &k : kcorr) {
    ins.read((char *)&k, sizeof(double));
  }

  // read the spectra only if the redshift is 0
  if (red < 1.e-20) {
    int nblamb;
    ins.read((char *)&nblamb, sizeof(int));
    // put ori=1 because it's a SED
    lamb_flux.resize(nblamb, oneElLambda(-999, -999, 1));
    for (auto &oneEl : lamb_flux) {
      ins.read((char *)&oneEl.lamb, sizeof(double));
    }
    for (auto &oneEl : lamb_flux) {
      ins.read((char *)&oneEl.val, sizeof(double));
    }
  }

  return;
}

/*

  STARS

*/

/*
  Read the predicted magnitudes in binary format
  Case for the stars
  Used in zphota
*/
void StarSED::readMagBin(ifstream &ins) {
  int nbFlt;
  // Read information conserning the SED which is read
  ins.read((char *)&nummod, sizeof(int));
  ins.read((char *)&nbFlt, sizeof(int));

  // Read the magnitudes
  mag.resize(nbFlt, 99);
  for (auto &m : mag) {
    ins.read((char *)&m, sizeof(double));
  }

  // read the spectra only if the redshift is 0
  int nblamb;
  ins.read((char *)&nblamb, sizeof(int));
  // put ori=1 because it's a SED
  lamb_flux.resize(nblamb, oneElLambda(-999, -999, 1));
  for (auto &oneEl : lamb_flux) {
    ins.read((char *)&oneEl.lamb, sizeof(double));
  }
  for (auto &oneEl : lamb_flux) {
    ins.read((char *)&oneEl.val, sizeof(double));
  }

  return;
}

/*
  Write the magnitudes in binary format for the stars
*/
void StarSED::writeMag(bool outasc, ofstream &ofsBin, ofstream &ofsDat,
                       vector<flt> allFilters, string magtyp) const {
  // Info screen
  // cout << "Compute magnitudes for SED " << name  << "  \r " << flush;

  // number of filters
  int nbFlt = mag.size();
  // Write the library in a binary file
  ofsBin.write((char *)&nummod, sizeof(int));  // Index du template
  ofsBin.write((char *)&nbFlt, sizeof(int));   // Number of filters

  // Write the predicted magnitudes
  for (int k = 0; k < nbFlt; k++)
    ofsBin.write((char *)&(mag[k]), sizeof(double));

  // Write the spectra
  int nbLamb = lamb_flux.size();
  ofsBin.write((char *)&(nbLamb), sizeof(int));
  for (size_t k = 0; k < lamb_flux.size(); k++)
    ofsBin.write((char *)&(lamb_flux[k].lamb), sizeof(double));
  for (size_t k = 0; k < lamb_flux.size(); k++)
    ofsBin.write((char *)&(lamb_flux[k].val), sizeof(double));

  // Case with the ascii file
  if (outasc) {
    // Write output
    ofsDat << setw(6) << nummod << " ";
    ofsDat << setw(4) << nbFlt << " ";
    if (magtyp[0] != 'V') {  // AB
      for (int k = 0; k < nbFlt; k++) {
        ofsDat << setw(6) << mag[k] << " ";
      }
    } else {  // VEGA
      for (int k = 0; k < nbFlt; k++) {
        ofsDat << setw(6) << mag[k] + allFilters[k].ab << " ";
      }
    }
    ofsDat << endl;
  }

  return;
}

/*
   compute the magnitudes in each of the filters
*/
void SED::compute_magnitudes(const vector<flt> &filters) {
  double val;
  for (const auto &filter : filters) {
    // Derive the AB magnitudes in each filter
    vector<double> intFlux;
    intFlux = integrateSED(filter);
    if (intFlux[3] != INVALID_VAL) {
      if (intFlux[3] > 0.0) {
        val = -2.5 * LOG10D(intFlux[3] / intFlux[1] * filter.fcorr) - 48.6 +
              distMod;
      } else
        val = HIGH_MAG;
    } else
      val = INVALID_MAG;
    // mag stored into the vector of SED
    mag.push_back(val);
  }
}

vector<double> SED::compute_fluxes(const vector<flt> &filters) {
  size_t imagm = filters.size();
  vector<double> result(imagm, NULL_FLUX);
  // check that the SED is defined
  for (size_t k = 0; k < imagm; k++) {
    vector<double> intFlux;
    intFlux = integrateSED(filters[k]);
    result[k] =
        (intFlux[3] == INVALID_VAL) ? INVALID_FLUX : intFlux[3] / intFlux[1];
  }
  return result;
}

vector<double> SED::integrateSED(const flt &filter) {
  vector<double> results(6, 0.);

  if (lamb_flux.front().lamb > filter.lmin() ||
      lamb_flux.back().lamb < filter.lmax()) {
    return vector<double>(6, INVALID_VAL);
  }

  auto [x, newsed, newflt] =
      restricted_resampling(lamb_flux, filter.lamb_trans, -1);

  double r0, r1, r2, r3, r4, r5;
#pragma omp parallel for reduction(+ : r0, r1, r2, r3, r4, r5)
  for (size_t i = 0; i < x.size() - 1; i++) {
    double lmean = (x[i] + x[i + 1]) / 2;
    double tmean = (newflt[i] + newflt[i + 1]) / 2;
    double fmean = (newsed[i] + newsed[i + 1]) / 2;
    double dlbd = (x[i + 1] - x[i]);
    double tmp = tmean * dlbd;
    r0 += tmp;
    r2 += tmp * lmean;
    r3 += tmp * fmean;
    r4 += tmp * fmean * lmean;
    r5 += tmp / pow(lmean, 2);
  }
  r1 = r5 * c;
  return vector<double>{r0, r1, r2, r3, r4, r5};
}
