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
  char t = toupper(type[0]);
  if (t == 'S') {
    nlib = STAR;
  } else if (t == 'Q') {
    nlib = QSO;
  } else if (t == 'G') {
    nlib = GAL;
  } else {
    throw invalid_argument("Object type not recognized: " + type);
  }

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
  chi2 = 1.e9;  // chi2 of the fit
  dm = -999.;   // Rescaling of the template
  distMod = 0;
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
  Add a shift to the predicted magnitudes, valid for QSO/STARS/GAL
*/
void SED::applyShift(const vector<double> &shifts, const int imagm) {
  // Check that the expected number of filters in the library matches the size
  // of the vector shifts
  if (imagm == int(shifts.size())) {
    // Loop over the filters
    for (int k = 0; k < imagm; k++) {
      mag[k] = mag[k] + shifts[k];
    }
  } else {
    if (shifts.size() > 1)
      cout << "Not the same number of shifts (apply_shift) than the number of "
              "filters in the library "
           << imagm << " " << shifts.size() << endl;
  }

  return;
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
    throw invalid_argument("Can't open SED file " + sedFile);
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

    if ((((lamb_flux.end() - 1)->lamb) * (1. + red) < filter.lmax()) &&
        (red == 0)) {
      if (verbose) {
        cout << "A problem could occur since maximum of SED "
             << (lamb_flux.end() - 1)->lamb << " below max of the filter "
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
      lamb_flux.emplace_back(1.e8, 0, 1);
    }
  }
}

/*
  Integrate within a filter, valid in the three cases GAL/QSO/STARS
*/
vector<double> SED::integrateSED(const flt &filter) {
  vector<oneElLambda> lamb_all;
  vector<double> mag;
  double area = 0.;
  double arean = 0.;
  double areaCorr = 0.;
  double fmel = 0.;
  double lmoy = 0.;
  double lnum = 0.;

  // Check that the full SED is defined over the filter
  if (lamb_flux.front().lamb < filter.lmin() &&
      lamb_flux.back().lamb > filter.lmax()) {
    // Check that there is a minimum coverage between the SED and the filter
    // if( (lamb_flux.begin())->lamb < filter.lmax() &&
    // (lamb_flux.end()-1)->lamb > filter.lmin()){

    // Initialise with a flux=0 at lambda=0
    lamb_all.emplace_back(0, 0, 1);
    // Loop over the SED
    for (vector<oneElLambda>::iterator it = lamb_flux.begin();
         it < lamb_flux.end(); ++it) {
      // Keep the first element before lmin (replace as long as you don't meet
      // the filter).
      if (it->lamb < filter.lmin()) lamb_all[0] = *it;
      // keep the SED when well included in the filter lambda range
      if ((it->lamb) >= filter.lmin() && (it->lamb) <= filter.lmax())
        lamb_all.push_back(*it);
      // keep only the last element after lmax
      if (it->lamb > filter.lmax()) {
        lamb_all.push_back(*it);
        break;
      }
    }

    // Concatenate two vectors composed of the filter and the SED
    lamb_all.insert(lamb_all.end(), (filter.lamb_trans).begin(),
                    (filter.lamb_trans).end());

    // Sort the vector in increasing lambda (vector including the filter and the
    // SED)
    sort(lamb_all.begin(), lamb_all.end());

    // be sure that the first and the last element of the concatenate vector is
    // the SED, not the filter If not, add an element with a flux=0 before or
    // after
    if (lamb_all.begin()->ori == 0) {
      lamb_all.emplace(lamb_all.begin(), lamb_all.begin()->lamb - 1., 0, 1);
    }
    if (lamb_all.back().ori == 0) {
      lamb_all.emplace_back(lamb_all.back().lamb + 1., 0, 1);
    }

    // Resample in two vectors with a common lambda range (the combination of
    // the filter and SED lambda) filter
    vector<oneElLambda> new_trans;
    resample(lamb_all, new_trans, 0, 0, 1.e50);
    // sed
    vector<oneElLambda> new_sed;
    resample(lamb_all, new_sed, 1, 0, 1.e50);

    // Loop to find the mean lambda, mean flux, mean trans
    // sum to get the area, arean, the integrated flux, the integrated flux
    // accros the filter
    for (int k = 0; k < int(new_sed.size()) - 1; ++k) {
      // In case the resampling was impossible, put the value at 0 (means no
      // transmission or no SED emission)
      if (new_sed[k].ori < 0) new_sed[k].val = 0;
      if (new_sed[k + 1].ori < 0) new_sed[k + 1].val = 0;
      if (new_trans[k].ori < 0) new_trans[k].val = 0;
      if (new_trans[k + 1].ori < 0) new_trans[k + 1].val = 0;
      // lambda mean
      double lmean = ((new_sed[k]).lamb + (new_sed[k + 1]).lamb) / 2.;
      // trans mean
      double tmean = ((new_trans[k]).val + (new_trans[k + 1]).val) / 2.;
      // flux mean
      double fmean = ((new_sed[k]).val + (new_sed[k + 1]).val) / 2.;
      // delta lambda
      double dlbd = ((new_sed[k + 1]).lamb - (new_sed[k]).lamb);
      // conv
      double conv = c / pow(lmean, 2.);
      // integral T dlambda
      area = area + tmean * dlbd;
      // integral (T*c/lambda**2) dlambda
      arean = arean + tmean * conv * dlbd;
      // integral (T*lambda) dlambda
      lmoy = lmoy + lmean * tmean * dlbd;
      // integral (F*T) dlambda
      fmel = fmel + fmean * tmean * dlbd;
      // integral (F*T*lambda) dlambda
      lnum = lnum + fmean * tmean * lmean * dlbd;
      // integral (T/lambda^2) dlambda
      areaCorr = areaCorr + tmean / pow(lmean, 2.) * dlbd;
    }

  } else {
    // SED not defined over the full filter
    area = INVALID_VAL;
    arean = INVALID_VAL;
    lmoy = INVALID_VAL;
    fmel = INVALID_VAL;
    lnum = INVALID_VAL;
    areaCorr = INVALID_VAL;
  }

  // 0el: int (T) dlambda
  // 1el: int (T*c/lambda^2) dlambda
  // 2el: int (T*lambda) dlambda
  // 3el: int (F*T) dlambda
  // 4el: int (F*T*lambda) dlambda
  // 5el: int (T/lambda^2) dlambda
  mag.push_back(area);
  mag.push_back(arean);
  mag.push_back(lmoy);
  mag.push_back(fmel);
  mag.push_back(lnum);
  mag.push_back(areaCorr);

  lamb_all.clear();

  return mag;
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

/*
  resample the vector

  INPUT
  lamb_all= all elements concatenated (filter+SED)
  origine=  indicate how the vector needs to be filled (0 for filter, 1 for SED,
  2 for attenuation laws, 3 opacity) lmin= lambda min to be considered lmax=
  lambda max to be considered

  OUTPUT
  lamb_interp= vector with a value in each lambda element, obtained with linear
  interpolation
*/
void SED::resample(vector<oneElLambda> &lamb_all,
                   vector<oneElLambda> &lamb_interp, const int origine,
                   const double lmin, const double lmax) {
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
      if (prevEl.lamb > 0) {
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

/*
  Sum an additional contribution to the existing spectra
*/
void SED::sumSpectra(SED addSED, const double rescal) {
  // Change the origin of lamb_all
  for (size_t i = 0; i < addSED.lamb_flux.size(); i++) {
    addSED.lamb_flux[i].ori = 1;
  }
  for (size_t i = 0; i < lamb_flux.size(); i++) {
    lamb_flux[i].ori = 0;
  }

  // Concatenate two vectors composed of "oneElLambda" including this spectra
  // and the one to be added
  lamb_flux.insert(lamb_flux.end(), addSED.lamb_flux.begin(),
                   addSED.lamb_flux.end());

  // Sort the vector in increasing lambda
  sort(lamb_flux.begin(), lamb_flux.end());

  // Resample the first spectra into a common lambda range
  vector<oneElLambda> spectra_ori;
  resample(lamb_flux, spectra_ori, 0, 0., 1e20);
  // Resample the second spectra into a common lambda range
  vector<oneElLambda> spectra_add;
  resample(lamb_flux, spectra_add, 1, 0., 1e20);

  // Sum the two spectra and generate a new lamb_flux
  lamb_flux.clear();
  for (size_t i = 0; i < spectra_ori.size(); i++) {
    // Check that the two have been well resample before
    if (spectra_add[i].ori < 0) spectra_add[i].val = 0.;
    if (spectra_ori[i].ori < 0) spectra_ori[i].val = 0.;
    lamb_flux.emplace_back(spectra_add[i].lamb,
                           (spectra_add[i].val * rescal) + spectra_ori[i].val,
                           1);
  }

  // clean
  spectra_add.clear();
  spectra_ori.clear();
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

/*
  Apply dust extinction on the original flux of the sed
  Only for galaxies and QSO
*/
void SED::applyExt(const double ebv, const ext &oneext) {
  (*this).ebv = ebv;
  // if needed because E(b-V)>0
  if (ebv > 1.e-20) {
    // work with the original lamb_flux
    vector<oneElLambda> lamb_all = lamb_flux;
    vector<oneElLambda> new_ext;

    // Concatenate two vectors composed of "oneElLambda" including the
    // extinction law and the SED
    lamb_all.insert(lamb_all.end(), oneext.lamb_ext.begin(),
                    oneext.lamb_ext.end());

    // Sort the vector in increasing lambda (vector including the extinction law
    // and the SED)
    sort(lamb_all.begin(), lamb_all.end());
    // Resample the extinction in a lambda range combining the SED and the
    // extinction law
    resample(lamb_all, new_ext, 2, 0., 1e20);

    // Loop over the SED and extinction curves using the concatenate ext+SED
    // vector
    vector<oneElLambda> lamb_new;
    for (size_t k = 0; k < lamb_all.size(); ++k) {
      // If we are well looking at the original SED
      if (lamb_all[k].ori == 1) {
        double val = lamb_all[k].val;
        double lamb = lamb_all[k].lamb;
        if (oneext.lmin < lamb && oneext.lmax > lamb) {
          // If extinction not defined (resampling failed), put it at 0
          if (new_ext[k].ori < 0) new_ext[k].val = 0;
          // Change the value of the flux according to the dust extinction
          val = (lamb_all[k]).val * pow(10., (-0.4 * ebv * (new_ext[k]).val));
        }
        // Store the result into a new lamb_flux vector
        lamb_new.emplace_back(lamb, val, 1);
      }
    }

    // Replace the old lambda flux vector with the new one including extinction
    lamb_flux = lamb_new;

    // Indicate what is the value of the ebv and the index of the extinction law
    extlawId = oneext.numext;

    // Clean
    lamb_all.clear();
    lamb_new.clear();
    new_ext.clear();
  }

  return;
}

/*
  Apply dust extinction on the emission lines (fac_line)
  Only for galaxies and QSO
*/
void SED::applyExtLines(const ext &oneext) {
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
  // change the code in read_lib and onesource.cpp double a=0.2; No redshift
  // dependency
  double a = 0.;

  // if needed because E(b-V)>0
  if (ebv > 1.e-20) {
    // work with fac_line as spectra
    vector<oneElLambda> line_all = fac_line;
    vector<oneElLambda> new_ext;

    // Concatenate two vectors composed of "oneElLambda" including the
    // extinction law and the SED
    line_all.insert(line_all.end(), oneext.lamb_ext.begin(),
                    oneext.lamb_ext.end());

    // Sort the vector in increasing lambda (vector including the extinction law
    // and the SED)
    sort(line_all.begin(), line_all.end());

    // Resample the extinction in a lambda range combining the SED and the
    // extinction law
    resample(line_all, new_ext, 2, 0., 1e20);

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
          val =
              (line_all[k]).val * pow(10., (-0.4 * ebv / f * (new_ext[k]).val));
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
  }

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

  vector<oneElLambda> lamb_all, lamb_end;
  // Split in two the SED lambda (below and above 1216A)
  for (size_t k = 0; k < lamb_flux.size(); ++k) {
    // If below 1215.67, fill lamb_all, otherwise keep that in lamb_end
    if (lamb_flux[k].lamb < 1215.67) {
      lamb_all.push_back(lamb_flux[k]);
    } else {
      lamb_end.push_back(lamb_flux[k]);
    }
  }
  // Concatenate two vectors composed of "oneElLambda" including the IGM opacity
  // and the SED below 1216
  lamb_all.insert(lamb_all.end(), (opaAll[ind].lamb_opa).begin(),
                  (opaAll[ind].lamb_opa).end());

  // Sort the vector in increasing lambda (vector including the extinction law
  // and the SED)
  sort(lamb_all.begin(), lamb_all.end());

  vector<oneElLambda> new_opa;
  // Resample the IGM opacity in a lambda range combining the SED and the
  // opacity
  resample(lamb_all, new_opa, 3, 0., 1.e20);

  // Loop over the SED and opacity using the concatenate opa+SED vector
  vector<oneElLambda> lamb_new;
  for (size_t k = 0; k < lamb_all.size(); ++k) {
    // If we are well looking at the original SED
    if ((lamb_all[k]).ori == 1) {
      double lamb = lamb_all[k].lamb;
      double val = lamb_all[k].val;
      if (opaAll[ind].lmin < lamb && opaAll[ind].lmax > lamb) {
        // Put the opacity at 1 if the resample is not done (do nothing on the
        // flux)
        if ((new_opa[k]).ori < 0) (new_opa[k]).val = 1;
        // Multiply the flux with the opacity
        val = (lamb_all[k]).val * (new_opa[k]).val;
      }
      lamb_new.emplace_back(lamb, val, 1);
    }
  }

  // Replace the old lambda flux vector with the new one including IGM
  lamb_flux = lamb_new;
  // Add the SED at lambda > 1216
  lamb_flux.insert(lamb_flux.end(), lamb_end.begin(), lamb_end.end());
  lamb_new.clear();
  lamb_all.clear();

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
               int nummodC, string typeC, int idAgeC)
    : SED(nameC, tauC, ageC, nummodC, typeC, idAgeC) {
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
    calc_ph();       // compute the number of ionizing photons
    add_neb_cont();  // Compute the continuum
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
                      // vérifier ############

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

  // cout << "nombre de photons ionisants : \n pour HeII : "  << qi[0] << "\n
  // pour HeI : "  << qi[1] << "\n pour H (Lyman-continuum flux) : " << qi[2] <<
  // "\n pour H_2 : " << qi[3] << endl;

  return;
}

/*
  Add the continuum from the nebular regions
  Work done by Cedric Dubois
*/
void GalSED::add_neb_cont() {
  /* we assume that the emitting gas has an electron temperature of Te = 10000
     K, an electron density N = 100 cm-3 (low density limite), and a helium
     abundance of 10% by number relative to hydrogen.
  */

  // Atomic data :
  // double alpha_B = 2.59e-13; // [cm^3 s^-1] : total recombination coeff for
  // hydrogen in case B (except to groundstate), for Te = 10kK
  //  Different from Schearer, use Osterbrock
  double n_heII =
      0.1;  // proportion of Helium compared to hydrogen = n(HeII)/n(HI)

  // Wavelength grid for wavelength & emission coefficient :
  // Points from Aller (1984)   for lambda=1000 to 10000.1 (54 pts)
  // Points from Ferland (1980) for lambda > 10000.1     (17 pts)

  // Wavelength
  double ga_lamb[71] = {
      1000.,  1200.,  1300.,  1400.,  1500.,  1600.,   1800.,   2051.,  2053.,
      2200.,  2400.1, 2599.4, 2600.8, 2725.4, 2855.2,  2997.9,  3121.4, 3122.,
      3331.,  3421.4, 3422.,  3527.,  3642.,  3648.,   3679.2,  3679.9, 4000.,
      4282.8, 4499.9, 4996.5, 5096.,  5450.8, 5695.8,  5700.,   5995.9, 6633.8,
      6635.8, 6999.9, 7438.7, 7441.1, 7848.,  7850.4,  8193.3,  8196.2, 8196.7,
      8198.5, 8207.,  8209.,  8265.4, 8268.1, 8499.9,  9000.,   9500.1, 10000.1,
      14591., 14592., 22799., 22800., 32805., 32806.,  44705.,  44706., 58460.,
      58461., 74145., 74146., 91199., 91200., 109878., 109879., 132173.};

  /*  Emission coefficients from Aller (1984, Table 4-9, p.102)
      for lambda <=10000.1 Ang,
      and Ferland (1980, PASP 92, 506, Table 1) for lambda >=14591 Ang
      all for Te = 10000 K:

      gamma_H in [1.d-40 ergs cm^3 s^-1 Hz^-1]   */
  // Hydrogen recombination
  double ga_H_val[71] = {
      0.001,  0.009,  0.023,  0.051,  0.100,  0.181,  0.486,  1.272,  1.280,
      2.026,  3.453,  5.404,  5.419,  6.928,  8.740,  11.020, 13.231, 13.242,
      17.477, 19.490, 19.505, 21.976, 24.841, 1.387,  1.434,  1.435,  1.950,
      2.461,  2.883,  3.929,  4.150,  4.959,  5.534,  5.544,  6.252,  7.801,
      7.806,  8.693,  9.754,  9.760,  10.729, 10.735, 11.537, 11.544, 11.545,
      11.549, 11.569, 4.317,  4.368,  4.371,  4.578,  5.020,  5.449,  5.867,
      9.04,   5.93,   8.51,   6.90,   8.50,   7.56,   8.66,   8.06,   8.87,
      8.47,   9.11,   8.82,   9.34,   9.14,   9.58,   9.42,   9.80};

  /*   gamma_2q in [1.d-40 ergs cm^3 s^-1 Hz^-1]
       assume ga_2q=0 for lambda >=10000 A    */
  // To be check
  double ga_2q_val[71] = {
      0.000, 0.000, 5.611, 8.236, 9.210, 9.478, 9.196, 8.433, 8.426,
      7.947, 7.314, 6.747, 6.743, 6.423, 6.096, 5.770, 5.508, 5.507,
      5.206, 4.947, 4.946, 4.771, 4.588, 4.579, 4.532, 4.531, 4.088,
      3.746, 3.512, 3.055, 2.975, 2.715, 2.554, 2.552, 2.375, 2.064,
      2.063, 1.913, 1.751, 1.750, 1.616, 1.615, 1.515, 1.514, 1.514,
      1.514, 1.511, 1.511, 1.496, 1.496, 1.438, 1.324, 1.221, 1.129,
      0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.,
      0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.};

  /*  gamma_HeI in [1.d-40 ergs cm^3 s^-1 Hz^-1]
      assume ga_hei=ga_h for lambda >=10000 A  */
  // Helium recombination
  double ga_HeI_val[71] = {
      0.046,  0.079,  0.126,  0.210,  0.400,  0.869,  1.650,  3.895,  3.926,
      6.087,  9.014,  11.861, 7.156,  9.976,  12.651, 15.326, 17.443, 15.873,
      20.979, 22.989, 4.321,  4.902,  5.501,  5.532,  5.687,  1.450,  2.070,
      2.807,  3.310,  4.943,  5.231,  6.176,  6.759,  6.771,  7.583,  9.086,
      8.661,  9.391,  10.176, 10.000, 10.607, 8.081,  8.497,  5.880,  5.880,
      5.041,  5.048,  5.050,  5.098,  4.360,  4.596,  5.065,  5.483,  5.860,
      9.04,   5.93,   8.51,   6.90,   8.50,   7.56,   8.66,   8.06,   8.87,
      8.47,   9.11,   8.82,   9.34,   9.14,   9.58,   9.42,   9.80};

  // store wavelength, emission coefficients &  type in vector <oneElLambda> for
  // H, 2q(2 ph decay) & HeI
  vector<oneElLambda> ga_H;
  vector<oneElLambda> ga_2q;
  vector<oneElLambda> ga_HeI;
  int i;

  // fill the vectors gamma
  // notation 4 for the gama values
  for (i = 0; i < 71; i++) {
    // H :
    ga_H.emplace_back(ga_lamb[i], ga_H_val[i], 4);
    // 2q :
    ga_2q.emplace_back(ga_lamb[i], ga_2q_val[i], 4);
    // HeI :
    ga_HeI.emplace_back(ga_lamb[i], ga_HeI_val[i], 4);
  }

  // on utilise une fonction pour concatener lamb_flux et ga_H/2q/HeI (avec les
  // .val en log pour ensuite interpoler en log NB : pas besoin de prendre le
  // log des valeurs de lamb_flux.val car on va les écraser)

  // Concatenate lamb_flux & ga_H/2q/HeI in order have gama in the same lambda
  // as the SED. Interpolate in log scale.
  vector<oneElLambda> ga_H_all, ga_2q_all, ga_HeI_all;

  // For H :
  ga_H_all = lamb_flux;
  ga_H_all.insert(ga_H_all.end(), ga_H.begin(), ga_H.end());
  // Convert the value in log
  for (size_t i = 0; i < ga_H_all.size(); i++) {
    ga_H_all[i].val = log10(ga_H_all[i].val);
  }

  // For 2q :
  ga_2q_all = lamb_flux;
  ga_2q_all.insert(ga_2q_all.end(), ga_2q.begin(), ga_2q.end());
  for (size_t i = 0; i < ga_2q_all.size(); i++) {
    if (ga_2q_all[i].val > 0) {  // Some value of gamma are at 0. Need to deal
                                 // with them to use the log
      ga_2q_all[i].val = log10(ga_2q_all[i].val);  // convert in log
    } else {
      ga_2q_all[i].val =
          -100;  // Use an extreme value in log, since the original value is 0
    }
  }

  // For HeI :
  ga_HeI_all = lamb_flux;
  ga_HeI_all.insert(ga_HeI_all.end(), ga_HeI.begin(), ga_HeI.end());
  for (size_t i = 0; i < ga_HeI_all.size(); i++) {
    ga_HeI_all[i].val = log10(ga_HeI_all[i].val);
  }

  // Log interpolation of the 3 gamma functions using "resample"
  vector<oneElLambda> ga_H_interp, ga_2q_interp, ga_HeI_interp;

  // Sort by increasing order :
  sort(ga_H_all.begin(), ga_H_all.end());
  sort(ga_2q_all.begin(), ga_2q_all.end());
  sort(ga_HeI_all.begin(), ga_HeI_all.end());

  // Log interpolation
  resample(ga_H_all, ga_H_interp, 4, 0, 1.6e+6);
  resample(ga_2q_all, ga_2q_interp, 4, 0, 1.6e+6);
  resample(ga_HeI_all, ga_HeI_interp, 4, 0, 1.6e+6);

  // Remove the log used only for interpolation
  // H :
  for (size_t i = 0; i < ga_H_interp.size(); i++)
    ga_H_interp[i].val = pow(10, ga_H_interp[i].val);
  // 2q :
  for (size_t i = 0; i < ga_2q_interp.size(); i++)
    ga_2q_interp[i].val = pow(10, ga_2q_interp[i].val);
  // HeI :
  for (size_t i = 0; i < ga_HeI_interp.size(); i++)
    ga_HeI_interp[i].val = pow(10, ga_HeI_interp[i].val);

  // Sum the three gamma values
  vector<oneElLambda> ga_tot;
  for (size_t i = 0; i < ga_H_interp.size(); i++) {
    double val = 0;
    if (ga_H_interp[i].lamb > 1000)
      val = (ga_H_interp[i].val + ga_2q_interp[i].val +
             ga_HeI_interp[i].val * n_heII);  // Sum the gamma values
    ga_tot.emplace_back(ga_H_interp[i].lamb, val, 4);
  }

  return;
}

/*
  Generate the emission lines in rest-frame without renormaliation or extinction
  based on the physical prescriptions
  Work done by Cedric Dubois
*/
void GalSED::generateEmPhys(double zmet, double qi) {
  // reference luminosity (H_beta in case B) :
  // L_ref = f_ga*qi*h*c/lamb_Hbeta*alpha_Hbeta/alpha_B   // Osterbrock data (H.
  // Krueger et al, 1995)
  double L_ref =
      f_ga * qi * 4.780e-13;  // [erg.s-1.cm-2], different value than Schearer
                              // because of the different alpha_B value

  /*
    Line ratios for 3 different metallicities Z1 (Z < 0.02), Z2 (0.02 < Z <=
    0.2), Z3-Z5 (0.2 < Z) BC03 convention: 62 solar 0.02, 52 subsolar 0.008, 42
    subsolar 0.004 'H_Lya', '[CII]', '[OIII]', '[CIII]', '[NII]', '[CII]',
    '[MgII]', '[OII]a', '[OII]b', '[NeIII]', 'H_zeta+[HeI]', 'H_eps+[NeIII]',
    '[HeI]', '[SII]', '[SII]', 'H_delta', '[OIII]', 'H_gamma', '[HeI]',
    '[ArIV]+[HeI]', 'H_beta', '[OIII]', '[OIII]', '[NI]', '[NII]', '[HeI]',
    '[OI]', '[SIII]', '[NII]', 'H_alpha', '[NII]', '[HeI]', '[SII]', '[SII]',
    '[HeI]', '[ArIII]', '[OII]', '[OII]', '[ArIII]', '[SIII]', '[SIII]',
    '[SII]', '[SII]', '[SII]', '[HeI]', 'Pa_Pa', 'Pb_Pa', 'Pg_Pa', 'Pd_Pa',
    'Pe_Pa', 'P6_HI', 'P7_HI', 'P8_HI', 'P9_HI', 'P10_HI', 'P11_HI', 'P12_HI',
    'Bra', 'Brb', 'Brg', 'Brd', 'Bre', 'Br6', 'Br7', 'Br8'
  */

  // Wavelength air
  double lb_line[65] = {
      1215.67,   1335.708,  1663.48,   1908.73,   2140.896,  2325.726,
      2799.116,  3726.04,   3728.8,    3869.063,  3888.856,  3968.771,
      4026.191,  4068.601,  4076.349,  4101.734,  4363.21,   4340.472,
      4471.48,   4712.249,  4861.352,  4958.911,  5006.843,  5199.079,
      5754.594,  5875.62,   6300.304,  6312.064,  6548.041,  6562.787,
      6583.461,  6678.151,  6716.44,   6730.812,  7065.19,   7135.803,
      7319.458,  7330.196,  7751.057,  9068.611,  9530.586,  10286.731,
      10320.491, 10336.407, 10827.334, 18750.976, 12818.072, 10938.174,
      10049.368, 9545.969,  9229.667,  9014.909,  8862.782,  8750.473,
      8665.018,  8598.392,  8545.383,  40511.742, 26251.517, 21655.267,
      19445.581, 18174.237, 17362.143, 16806.509, 16407.21};

  // Check
  int longEL = sizeof(lb_line) / sizeof(lb_line[0]);
  if (longEL != 65) cout << " Error with size of lines " << longEL << endl;

  // For Z < 0.03 in solar unit, Z1=0.0004
  double Z1_line[65] = {
      22.2,    0.0,    0.0,    0.0,    0.0,    0.0,    0.0,     0.245,
      0.245,   0.295,  0.203,  0.27,   0.015,  0.005,  0.002,   0.256,
      0.109,   0.466,  0.036,  0.01,   1.0,    1.097,  3.159,   0.003,
      0.0,     0.096,  0.008,  0.009,  0.005,  2.87,   0.015,   0.026,
      0.037,   0.029,  0.028,  0.027,  0.012,  0.007,  0.067,   0.0,
      0.0,     0.0,    0.0,    0.0,    0.232,  0.352,  0.165,   0.0906,
      0.0368,  0.0366, 0.0254, 0.0184, 0.0138, 0.0106, 0.0185,  0.00666,
      0.00541, 0.0884, 0.0471, 0.0281, 0.0186, 0.0127, 0.00914, 0.00683,
      0.00524};
  // check
  longEL = sizeof(Z1_line) / sizeof(Z1_line[0]);
  if (longEL != 65) cout << " Error with size of lines Z1 " << longEL << endl;

  // 0.03 < Z <= 0.3 in solar unit, Z2=0.004
  double Z2_line[65] = {
      22.2,    0.0,    0.058,  0.0,    0.0,    0.0,    0.31,    0.896,
      0.896,   0.416,  0.192,  0.283,  0.015,  0.017,  0.007,   0.256,
      0.066,   0.466,  0.036,  0.014,  1.0,    1.617,  4.752,   0.01,
      0.0,     0.108,  0.041,  0.017,  0.059,  2.87,   0.175,   0.03,
      0.188,   0.138,  0.023,  0.071,  0.027,  0.014,  0.176,   0.2,
      0.51,    0.0,    0.0,    0.0,    0.225,  0.352,  0.165,   0.0906,
      0.0368,  0.0366, 0.0254, 0.0184, 0.0138, 0.0106, 0.0185,  0.00666,
      0.00541, 0.0884, 0.0471, 0.0281, 0.0186, 0.0127, 0.00914, 0.00683,
      0.00524};
  // Check
  longEL = sizeof(Z2_line) / sizeof(Z2_line[0]);
  if (longEL != 65) cout << " Error with size of lines Z2 " << longEL << endl;

  // 0.3 < Z in solar unit, Z3-Z5=0.008-0.05
  double Z3_line[65] = {
      22.2,    0.11,   0.01,   0.18,   0.01,   0.29,   0.07,    1.505,
      1.505,   0.3,    0.107,  0.159,  0.015,  0.029,  0.011,   0.256,
      0.01,    0.466,  0.05,   0.0,    1.0,    1.399,  4.081,   0.03,
      0.01,    0.14,   0.13,   0.03,   0.136,  2.87,   0.404,   0.03,
      0.3,     0.21,   0.04,   0.035,  0.026,  0.014,  0.086,   0.365,
      0.945,   0.048,  0.058,  0.054,  0.229,  0.352,  0.165,   0.0906,
      0.0368,  0.0366, 0.0254, 0.0184, 0.0138, 0.0106, 0.0185,  0.00666,
      0.00541, 0.0884, 0.0471, 0.0281, 0.0186, 0.0127, 0.00914, 0.00683,
      0.00524};
  // Check
  longEL = sizeof(Z3_line) / sizeof(Z3_line[0]);
  if (longEL != 65) cout << " Error with size of lines Z3 " << longEL << endl;

  // Select the right table depending on the metallicity Zmet
  // Multiply by the reference luminosity
  fac_line.clear();
  for (int i = 0; i < 65; i++) {
    // Find the line ratio
    if (zmet / 0.02 <= 0.03) {
      fac_line.emplace_back(lb_line[i], L_ref * Z1_line[i], 5);
    } else if (zmet / 0.02 > 0.03 && zmet / 0.02 <= 0.3) {
      fac_line.emplace_back(lb_line[i], L_ref * Z2_line[i], 5);

    } else if (0.3 < zmet / 0.02) {
      fac_line.emplace_back(lb_line[i], L_ref * Z3_line[i], 5);
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
  // Wavelength air
  double lb_line[65] = {
      1215.67,   1335.708,  1663.48,   1908.73,   2140.896,  2325.726,
      2799.116,  3726.04,   3728.8,    3869.063,  3888.856,  3968.771,
      4026.191,  4068.601,  4076.349,  4101.734,  4363.21,   4340.472,
      4471.48,   4712.249,  4861.352,  4958.911,  5006.843,  5199.079,
      5754.594,  5875.62,   6300.304,  6312.064,  6548.041,  6562.787,
      6583.461,  6678.151,  6716.44,   6730.812,  7065.19,   7135.803,
      7319.458,  7330.196,  7751.057,  9068.611,  9530.586,  10286.731,
      10320.491, 10336.407, 10827.334, 18750.976, 12818.072, 10938.174,
      10049.368, 9545.969,  9229.667,  9014.909,  8862.782,  8750.473,
      8665.018,  8598.392,  8545.383,  40511.742, 26251.517, 21655.267,
      19445.581, 18174.237, 17362.143, 16806.509, 16407.21};
  // Check
  int longEL = sizeof(lb_line) / sizeof(lb_line[0]);
  if (longEL != 65) cout << " Error with size of lines " << longEL << endl;

  // conversion factor between Halpa and other lines
  // Use a ratio 1 for Halpha/OII, according to Iary. Drop the 1.77 since
  // attenuation seems to be alreay included I took OIII5007/OIII4959=3 I
  // took 4.081 betwen 5007 and Hbeta, as the physical recipes. This ratio will
  // be modified later as a function of redshift. I took 0.3 for NII/Halpha
  // according to the BPT diagram 6 for Lyman_alpha, but no real idea
  // double Z0_line[65] =
  // {6,0,0,0,0,0,0,1.425,1.425,0,0,0,0,0,0,0,0,0,0,0,1,1.36,4.081,0,0,0,0,0,0,2.85,0.86,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  // Original ratio used in 2009 paper
  double Z0_line[65] = {1.62, 0, 0, 0,    0,  0, 0, 0.81, 0.81, 0,    0, 0, 0,
                        0,    0, 0, 0,    0,  0, 0, 1,    0.21, 0.59, 0, 0, 0,
                        0,    0, 0, 2.90, 0., 0, 0, 0,    0,    0,    0, 0, 0,
                        0,    0, 0, 0,    0,  0, 0, 0,    0,    0,    0, 0, 0,
                        0,    0, 0, 0,    0,  0, 0, 0,    0,    0,    0, 0, 0};
  // check
  longEL = sizeof(Z0_line) / sizeof(Z0_line[0]);
  if (longEL != 65) cout << " Error with size of lines Z0 " << longEL << endl;

  // Rescale the lines according to the Kennicutt 2012 relation, only for blue
  // galaxies Using Halpha as a reference
  double scaleFac;
  if (NUVR < 4)
    scaleFac = pow(10., -0.4 * MNUV_int - 6.224);
  else
    scaleFac = 0;
  // Using Hbeta as a reference
  scaleFac = scaleFac / 2.85;

  // Create the emission line and rescale them
  fac_line.clear();
  for (int i = 0; i < 65; i++) {
    fac_line.emplace_back(lb_line[i], scaleFac * Z0_line[i], 5);
  }
  // sort the vector by wavelength
  sort(fac_line.begin(), fac_line.end());
}

/*
  Generate the emission lines in rest-frame without renormaliation or extinction
  based on the Kenicutt relation and the SFR
*/
void GalSED::generateEmEmpSFR(double sfr, double NUVR) {
  // Wavelength air
  double lb_line[65] = {
      1215.67,   1335.708,  1663.48,   1908.73,   2140.896,  2325.726,
      2799.116,  3726.04,   3728.8,    3869.063,  3888.856,  3968.771,
      4026.191,  4068.601,  4076.349,  4101.734,  4363.21,   4340.472,
      4471.48,   4712.249,  4861.352,  4958.911,  5006.843,  5199.079,
      5754.594,  5875.62,   6300.304,  6312.064,  6548.041,  6562.787,
      6583.461,  6678.151,  6716.44,   6730.812,  7065.19,   7135.803,
      7319.458,  7330.196,  7751.057,  9068.611,  9530.586,  10286.731,
      10320.491, 10336.407, 10827.334, 18750.976, 12818.072, 10938.174,
      10049.368, 9545.969,  9229.667,  9014.909,  8862.782,  8750.473,
      8665.018,  8598.392,  8545.383,  40511.742, 26251.517, 21655.267,
      19445.581, 18174.237, 17362.143, 16806.509, 16407.21};
  int longEL = sizeof(lb_line) / sizeof(lb_line[0]);
  if (longEL != 65) cout << " Error with size of lines " << longEL << endl;

  // conversion factor between Halpa and other lines
  // Use a ratio 1 for Halpha/OII, according to Iary. Drop the 1.77 since
  // attenuation seems to be alreay included I took OIII5007/OIII4959=3 I
  // took 4.081 betwen 5007 and Hbeta, as the physical recipes. This ratio will
  // be modified later as a function of redshift. I took 0.3 for NII/Halpha
  // according to the BPT diagram 6 for Lyman_alpha, but no real idea
  double Z0_line[65] = {6, 0,    0,    0, 0, 0, 0, 1.425, 1.425, 0, 0, 0, 0, 0,
                        0, 0,    0,    0, 0, 0, 1, 1.36,  4.081, 0, 0, 0, 0, 0,
                        0, 2.85, 0.86, 0, 0, 0, 0, 0,     0,     0, 0, 0, 0, 0,
                        0, 0,    0,    0, 0, 0, 0, 0,     0,     0, 0, 0, 0, 0,
                        0, 0,    0,    0, 0, 0, 0, 0,     0};
  // Checks
  longEL = sizeof(Z0_line) / sizeof(Z0_line[0]);
  if (longEL != 65) cout << " Error with size of lines Z0 " << longEL << endl;

  // Rescale the lines according to the Kennicutt 2012 relation, only for blue
  // galaxies Using Halpha as a reference
  double scaleFac;
  if (NUVR < 4)
    scaleFac = pow(
        10., log10(sfr) + 41.27 - log10(4 * pi * 100 * pow(3.08568, 2)) - 36);
  else
    scaleFac = 0;
  // Using Hbeta as a reference
  scaleFac = scaleFac / 2.85;

  // Create the emission line and rescale them
  fac_line.clear();
  for (int i = 0; i < 65; i++) {
    fac_line.emplace_back(lb_line[i], scaleFac * Z0_line[i], 5);
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

/*
  Fill the variable with the SED properties like D4000, or luv
*/
void GalSED::SEDproperties() {
  double fluxconv = 1. / (4 * pi * 100 * pow(pc, 2));

  // construct a heaviside filter between 0.21 and 0.25 micron
  flt fltUV(2100., 2500., 100);
  // Integrate the SED within this filter in NUV
  luv = (this->integrateSED(fltUV))[3];
  if (luv > 0) luv = log10(luv * pow(2300, 2) / 400 / c / fluxconv);
  // construct a heaviside filter between 0.55 and 0.65 micron
  flt fltOpt(5500., 6500., 100);
  // Integrate the SED within this filter in optical
  lopt = (this->integrateSED(fltOpt))[3];
  if (lopt > 0) lopt = log10(lopt * pow(6000, 2.) / 1000. / c / fluxconv);

  // construct a heaviside filter between 2.1 and 2.3 micron
  flt fltNIR(21000., 23000., 100.);
  // Integrate the SED within this filter in NUV
  lnir = (this->integrateSED(fltNIR))[3];
  if (lnir > 0) lnir = log10(lnir * pow(22000, 2) / 2000 / c / fluxconv);

  // construct two heaviside filters between 0.375 and 0.425 micron to measure
  // the balmer break
  flt fltBalm1(3750., 3950., 100);
  flt fltBalm2(4050., 4250., 100);
  // Integrate the SED within this filter in optical
  double lBalm1 = (this->integrateSED(fltBalm1))[3];
  double lBalm2 = (this->integrateSED(fltBalm2))[3];
  if (lBalm1 > 0 && lBalm2 > 0) d4000 = lBalm2 / lBalm1;

  // construct a heaviside filter between 8 and 1000 micron
  flt fltIR(80000., 10000000., 100);
  // Integrate the SED within this filter in IR
  ltir = (this->integrateSED(fltIR))[3];
  if (ltir > 0) ltir = log10(ltir / Lsol / fluxconv);

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
  ofsBin.write((char *)&ltir, sizeof(double));     // LTIR
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
    ofsDat << setw(3) << extlawId << " ";
    ofsDat << setw(3) << ebv << " ";
    ofsDat << setw(12) << ltir << " ";
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
  ins.read((char *)&ltir, sizeof(double));
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
    ofsDat << setw(3) << extlawId << " ";
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
    vector<double> intFlux = integrateSED(filter);
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
    vector<double> intFlux = integrateSED(filters[k]);
    result[k] =
        (intFlux[3] == INVALID_VAL) ? INVALID_FLUX : intFlux[3] / intFlux[1];
  }
  return result;
}
