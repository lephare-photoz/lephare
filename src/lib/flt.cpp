/*

  10/11/14
  Implementation of function for the "flt" class

*/

#include <unistd.h>  //posix interface, prooviding path access method

#include <algorithm>  // sort
#include <cmath>      // for the log
#include <fstream>    // print output file
#include <iomanip>    // std::setprecision
#include <iostream>   // print standard file
#include <sstream>
#include <string>
#include <vector>

// La Phare class
#include "SED.h"
#include "flt.h"
#include "globals.h"
#include "oneElLambda.h"

using namespace std;

/*
 READ THE FILTER GIVEN THE NAME OF THE FILE (for filter.cpp)
*/
void flt::read(const string &fltFile) {
  ifstream sflt;
  string lit;

  // open the filter file into a stream
  sflt.open(fltFile.c_str());
  if (!(sflt)) {
    throw invalid_argument("Can't open file " + fltFile);
  }

  // Take the stream line by line
  while (getline(sflt, lit)) {
    // If the first character of the line is not #
    if (check_first_char(lit)) {
      // put the line into the stream ss again
      stringstream ss(lit);

      // fill the lambda/trans values into one element
      double l, v;
      ss >> l;  // 1 col: lambda
      ss >> v;  // 2 column: transmission
      // origin at 0 to indicate that this element refers to a filter
      lamb_trans.emplace_back(l, v, 0);
    }
  }

  // Close the stream
  sflt.close();

  // Sort the filter, remove the low values on the edge, be sure that the
  // extreme points are ending with 0
  clean();

  // transform it to have all the filters with a T as a function of energy
  // trans needs to be after clean because lambdaMean computation supposes that
  // lambda_trans be ordered.
  trans();

  // need to call this to compute dwidth
  width();

  return;
}

/*
 Read the filter given a stream (for mag_gal.cpp)
*/
void flt::read(ifstream &sfiltIn) {
  int bid, nbLines;
  char firstChar;

  // read the line with the comments
  // The informations on the name, the number of wavelength are stored there
  sfiltIn >> firstChar >> nbLines >> name >> calibtyp >> id;

  lamb_trans.clear();
  if (firstChar != '#') {
    cout << "Don t read the expected character in the filter file " << endl;
  } else {
    // Loop over each lambda, according to the number of lines indicated in
    // comment
    for (int i = 0; i < nbLines; i++) {
      // fill the lambda/trans values into one element
      double l, v;
      sfiltIn >> l;    // 1 col: lambda
      sfiltIn >> v;    // 2 column: transmission
      sfiltIn >> bid;  // consume comment string in file
      lamb_trans.emplace_back(l, v, 0);
    }
  }

  return;
}

/*
  CLEAN THE FILTER
*/
// Sort the filter in lambda, remove the values with a low transmission on the
// edge, be sure that the extreme points are ending with 0
void flt::clean() {
  // Derive the highest value of the transmission
  tpeak = peak();

  // Sort the vector according to ascending lambda
  sort(lamb_trans.begin(), lamb_trans.end());

  // Remove the extreme elements of the filter as long as the transmission is
  // below 1% of the maximum transmission From the end of the filter to the
  // first significant transmission
  vector<oneElLambda>::iterator it = (lamb_trans.end() - 1);
  while ((*it).val < tpeak * 0.01) it--;
  lamb_trans.erase(it + 1, lamb_trans.end());
  // From the beginning to the first significant transmission
  it = lamb_trans.begin();
  while ((*it).val < tpeak * 0.01) it++;
  lamb_trans.erase(lamb_trans.begin(), it);

  // Add 0 value at the two extreme ends of the filter, 10 Angtrom below and
  // above
  oneElLambda addFirstEl((*lamb_trans.begin()).lamb - 10, 0, 0);
  lamb_trans.emplace(lamb_trans.begin(), lamb_trans.front().lamb - 10, 0, 0);
  oneElLambda addLastEl((*(lamb_trans.end() - 1)).lamb + 10, 0, 0);
  lamb_trans.emplace(lamb_trans.end(), lamb_trans.back().lamb + 10, 0, 0);

  return;
}

/*
  TRANSFORM THE FILTER IN TRANS/ENERGY
*/
void flt::trans() {
  lmean = lambdaMean();
  // if transtype=0, no need to do anything
  if (transtyp == 1 && lmean > 0) {
    // Loop over all the lambda of this filter
    for (vector<oneElLambda>::iterator it = lamb_trans.begin();
         it < lamb_trans.end(); ++it) {
      // Transform the filter according to trans*(lbd/lbmean)**transtyp (if
      // transtype=0 -> nothing)
      (it->val) = (it->val) * (it->lamb) / lmean;
    }
  }

  return;
}

/*
 DERIVED THE MAXIMUM TRANSMISSION OF THE FILTER
*/
double flt::peak() {
  // Duplicate the filter
  vector<oneElLambda> sortVec = lamb_trans;
  // Sort the vector according to descending transmission
  sort(sortVec.begin(), sortVec.end(),
       [](oneElLambda i, oneElLambda j) { return (i.val > j.val); });

  // First element of the sorted vector -> Maximum transmission
  tpeak = (sortVec[0]).val;

  return tpeak;
}

/*
 DERIVED LAMBDA MEAN OF THE FILTER
*/
double flt::lambdaMean() {
  double area = 0.;
  double lbmean = 0.;

  // Loop over all the lambda of this filter
  int nb = lamb_trans.size();
  for (int k = 0; k < nb - 1; ++k) {
    double trans = ((lamb_trans[k]).val + (lamb_trans[k + 1]).val) /
                   2.;  // averaged transmission
    double diffLamb =
        (lamb_trans[k + 1]).lamb - (lamb_trans[k]).lamb;  // step in lambda
    double avLamb = ((lamb_trans[k + 1]).lamb + (lamb_trans[k]).lamb) /
                    2.;                           // average lambda
    area = area + trans * diffLamb;               // int (T) dlambda
    lbmean = lbmean + trans * diffLamb * avLamb;  // int(T*lambda) dlambda
  }
  // int(T*lambda) dlambda  / int (T) dlambda
  lmean = lbmean / area;

  return lmean;
}

/*
 DERIVED THE PEAK OF THE FILTER AND ITS WIDTH
*/
double flt::width() {
  if (name == "Heavy") {
    return lmax() - lmin() - 2;
  }

  // Width defined when the transmission is above 0.5*maximum transmission
  double treshold = 0.5 * peak();

  double lmin = -999.;
  double lmax = -999.;
  double slope;

  // Iterate over each element in increasing lambda -> find the lower limit in
  // lambda at T=0.5
  int nb = lamb_trans.size();
  for (int k = 0; k < nb; ++k) {
    oneElLambda elP0(lamb_trans[k]);
    oneElLambda elP1(lamb_trans[k + 1]);
    // If the transmission of the P0 element below the threshold while the P1
    // element is above
    if (elP0.val < treshold && elP1.val > treshold) {
      // interpolation
      slope = (elP1.lamb - elP0.lamb) / (elP1.val - elP0.val);
      // Minimum lambda refined with interpolation
      lmin = elP0.lamb + (treshold - elP0.val) * slope;
      break;
    }
  }

  // Iterate over each element in decreasing lambda -> find the higher limit in
  // lambda at T=0.5
  for (int k = nb - 1; k > 0; --k) {
    oneElLambda elP0(lamb_trans[k - 1]);
    oneElLambda elP1(lamb_trans[k]);
    // If the transmission of the P0 element above the threshold while the P1
    // element is below
    if (elP0.val > treshold && elP1.val < treshold) {
      // interpolation
      slope = (elP1.lamb - elP0.lamb) / (elP1.val - elP0.val);
      // Maximum lambda refined with interpolation
      lmax = elP0.lamb + (treshold - elP0.val) * slope;
      break;
    }
  }

  // width of the filter
  dwidth = lmax - lmin;

  return dwidth;
}

/*
  EFFECTIVE WAVELENGTH
  Derived the effective wevelength of the filters
  Effective wavelenth based on the Vega reference spectra
*/
double flt::lambdaEff() {
  // create an SED object based on the Vega spectra
  SED vegSED("Vega");
  string vegaFile = lepharedir + "/vega/VegaLCB.sed";
  vegSED.read(vegaFile);
  // integrate the vega spectrum over the filter considered in this instance of
  // the object Various information on the integral is stored within magVeg
  vector<double> magVeg = vegSED.integrateSED(*this);

  // int (Fvega T lambda) dLambda /  int (Fvega T) dLambda
  if (magVeg[4] > 0) leff = magVeg[4] / magVeg[3];

  return leff;
}

/*
 DERIVED THE AB-VEGA CORRECTION
*/
double flt::abcorr() {
  // create an SED object based on the Vega spectra
  SED vegSED("Vega");
  string vegaFile = lepharedir + "/vega/VegaLCB.sed";
  vegSED.read(vegaFile);
  // integrate the vega spectrum over the filter considered in this instance of
  // the object Various information on the integral is stored within magVeg
  vector<double> magVeg = vegSED.integrateSED(*this);

  // int (Fvega T) dLambda /  int ( T*c/lambda^2) dLambda
  if (magVeg[3] > 0) ab = -2.5 * log10(magVeg[3] / magVeg[1]) - 48.6;

  return ab;
}

/*
 DERIVED THE THUAN GUNN CORRECTION
*/
double flt::tgcorr() {
  // create an SED object based on the Vega spectra
  SED vegSED("Vega");
  string vegaFile = lepharedir + "/vega/VegaLCB.sed";
  vegSED.read(vegaFile);
  // integrate the vega spectrum over the filter considered in this instance of
  // the object Various information on the integral is stored within magVeg
  vector<double> magVeg = vegSED.integrateSED(*this);

  // create an SED object based on the BD17 spectra
  SED tgSED("TG");
  string tgFile = lepharedir + "/vega/BD+17o4708.sed";
  tgSED.read(tgFile);
  // integrate the BD17 spectrum over the filter considered in this instance of
  // the object Various information on the integral is stored within magTG
  vector<double> magTG = tgSED.integrateSED(*this);

  // mtg(*) = mvega(*) + TGcor : TGcor=2.5log(Int[F(BD)TdL]/Int[F(Vega)TdL])
  // +9.50-0.03
  if (magTG[3] > 0)
    tg = 2.5 * log10(magTG[3] / magVeg[3]) + 9.5 - 0.03;
  else
    tg = -99;

  return tg;
}

/*
 DERIVED THE VEGA MAGNITUDE
*/
double flt::vega() {
  // create an SED object based on the Vega spectra
  SED vegSED("Vega");
  string vegaFile = lepharedir + "/vega/VegaLCB.sed";
  vegSED.read(vegaFile);
  vector<double> magVeg = vegSED.integrateSED(*this);

  // int (Fvega T) dLambda /  int ( T) dLambda
  if (magVeg[3] > 0) veg = 2.5 * log10(magVeg[3] / magVeg[0]);

  return veg;
}

/*
 DERIVED THE ABSOLUTE MAGNITUDES OF THE SUN
*/
double flt::magsun() {
  // create an SED object based on the Sun spectra
  SED sunSED("Sun");
  string sunFile = lepharedir + "/vega/SunLCB.sed";
  sunSED.read(sunFile);
  vector<double> magSun = sunSED.integrateSED(*this);

  /*
  mo-Mo= 5log D -5 , with D=1U.A. expressed in pc :  mo-Mo= -31.572
  fo = Lo / 4Pi A^2   with A=1U.A. express in cm
  mo = -2.5 lg(Lo,nu / 4Pi A^2) -48.59
  int (Fsun T) dLambda /  int ( T*c/lambda^2) dLambda
  */
  if (magSun[3] > 0)
    msun = -2.5 * log10(magSun[3] / magSun[1]) + 68.6227 - 48.6 + 31.572;

  return msun;
}

/*
 DERIVED THE FLUX CORRECTION DEPENDING ON THE CALIB
 fcorr= Fcalib(leff) * leff^2 *  (int (T / lambda^2) dLambda /  int ( Fcalib T )
 dLambda )
*/
double flt::fcorrec() {
  // Generate a SED which is used as a reference (Fcalib) depending on the
  // keyword "CALIB" Define a lambda range of 1000 steps between lambda min and
  // max of the filter
  SED calibSED("calib");
  calibSED.generateCalib(lmin() - 10, lmax() + 10, 1000, calibtyp);

  // Define the second part of the correction (Fcalib is the SED used as a
  // reference) int (T / lambda^2) dLambda /  int ( Fcalib T ) dLambda
  vector<double> magCalib = calibSED.integrateSED(*this);
  fcorr = magCalib[5] / magCalib[3];

  // If calib=4 or 5, take nuFnu=cst (case calib 1) to compute lambda_eff ->
  // change the calibration SED
  if (calibtyp == 4 || calibtyp == 5) {
    calibSED.generateCalib(lmin() - 10, lmax() + 10, 1000, 1);
    magCalib = calibSED.integrateSED(*this);
  }

  // int (Fcalib T lambda) dLambda /  int ( Fcalib T ) dLambda
  double leff2 = magCalib[4] / magCalib[3];

  // Compute Fcalib at leff2
  double B0;
  switch (calibtyp) {
    case 1:
      // reference function: 1/lambda (nuBnu=cst)
      B0 = 1. / leff2;
      break;
    case 2:
      // reference function: 1/lambda^3 (Bnu=nu)
      B0 = 1. / pow(leff2, 3.);
      break;
    case 3:
      // reference function: black body at 10000K
      B0 = blackbody(10000, leff2);
      break;
    case 4:
      // reference function: black body at 10000K
      B0 = blackbody(10000, leff2);
      break;
    case 5:
      // reference function: 1/lambda^3 (Bnu=nu)
      B0 = 1. / pow(leff2, 3.);
      break;
    default:
      // reference function: 1/lambda^2 (Bnu=cst)
      B0 = 1. / pow(leff2, 2.);  // Standard case 0 and others >5
      break;
  }
  fcorr = B0 * fcorr * pow(leff2, 2.);

  return fcorr;
}

/*
 EFFECTIVE WEVELENGTH COMPUTED WITH A DIFFERENT CALIBRATION SED
*/
double flt::lambdaEff2() {
  // Generate a SED which is used as a reference (Fcalib) depending on the
  // keyword "CALIB" Define a lambda range of 1000 steps between lambda min and
  // max of the filter
  SED calibSED("calib");
  calibSED.generateCalib(lmin() - 10, lmax() + 10, 1000, calibtyp);

  // If calib=4 or 5, take nuFnu=cst (case calib 1) to compute lambda_eff ->
  // change the calibration SED
  if (calibtyp == 4 || calibtyp == 5)
    calibSED.generateCalib(lmin() - 10, lmax() + 10, 1000, 1);

  vector<double> magCalib = calibSED.integrateSED(*this);
  // int (Fcalib T lambda) dLambda /  int ( Fcalib T ) dLambda
  double leff2 = magCalib[4] / magCalib[3];

  return leff2;
}

/*
write_output_filter
original: 11/11/2014

Write the output of the filter file

Input: vector of the filters
Output: the files
*/
void write_output_filter(string &filtfile, string &filtdoc,
                         vector<flt> vecFlt) {
  // create a stream with the two output files
  ofstream stfiltfile, stfiltdoc;
  stfiltfile.open(filtfile.c_str());
  // Check if file is opened
  if (!stfiltfile) {
    throw invalid_argument("Can't open file " + filtfile);
  }
  stfiltdoc.open(filtdoc.c_str());
  // Check if file is opened
  if (!stfiltdoc) {
    throw invalid_argument("Can't open file " + filtdoc);
  }

  // documentation on the screen
  cout << std::fixed;
  cout << setw(30) << left << "# NAME";
  cout << setw(8) << left << "IDENT";
  cout << setw(12) << right << "Lbda_mean";
  cout << setw(12) << right << "Lbeff(Vega)";
  cout << setw(12) << right << "FWHM ";
  cout << setw(10) << right << "AB-cor";
  cout << setw(10) << right << "TG-cor";
  cout << setw(10) << right << "VEGA";
  cout << setw(10) << right << "M_sun(AB)";
  cout << setw(8) << right << "CALIB";
  cout << setw(12) << right << "Lb_eff";
  cout << setw(12) << right << "Fac_corr" << endl;

  // Write the number of filters
  stfiltfile << "# " << vecFlt.size() << endl;

  // header for the .doc
  stfiltdoc << "# Name id Lambda_mean Lambda_eff width AB Vega fcorr transtype "
               "calib Nb_lines "
            << endl;

  // Loop over all the filters
  for (vector<flt>::iterator it = vecFlt.begin(); it < vecFlt.end(); ++it) {
    // derive filter properties
    double lambm = (it->lambdaMean()) / 10000.;
    double lambeff = (it->lambdaEff()) / 10000.;
    double wid = ((it->width()) / 10000.);
    double ab = (it->abcorr());
    double veg = (it->vega());
    double fcorr = (it->fcorrec());

    // OUTPUT SCREEN
    string shortName = (it->name).substr(((it->name).rfind("/") + 1));
    int sizeName = shortName.size() + 2;
    cout << setw(max(30, sizeName)) << left << shortName;
    cout << setw(8) << left << (it->id + 1);
    cout << setw(12) << right << setprecision(4) << lambm;
    cout << setw(12) << right << setprecision(4) << lambeff;
    cout << setw(12) << right << setprecision(4) << wid;
    cout << setw(10) << right << setprecision(4) << ab;
    cout << setw(10) << right << setprecision(4) << (it->tgcorr());
    cout << setw(10) << right << setprecision(4) << veg;
    cout << setw(10) << right << setprecision(4) << (it->magsun());
    cout << setw(8) << right << setprecision(2) << (it->calibtyp);
    cout << setw(12) << right << setprecision(4) << (it->lambdaEff2()) / 10000.;
    cout << setw(12) << right << setprecision(4) << fcorr << endl;

    // OUTPUT DOC
    stfiltdoc << setw(max(30, sizeName)) << left << shortName;
    stfiltdoc << setw(8) << left << (it->id + 1);
    stfiltdoc << setw(12) << right << setprecision(4) << lambm;
    stfiltdoc << setw(12) << right << setprecision(4) << lambeff;
    stfiltdoc << setw(10) << right << setprecision(4) << wid;
    stfiltdoc << setw(10) << right << setprecision(4) << ab;
    stfiltdoc << setw(10) << right << setprecision(4) << veg;
    stfiltdoc << setw(10) << right << setprecision(4) << fcorr;
    stfiltdoc << setw(8) << right << (it->transtyp);
    stfiltdoc << setw(8) << right << (it->calibtyp);
    stfiltdoc << setw(8) << right << (it->lamb_trans.size()) << endl;

    // OUTPUT FILTER
    stfiltfile << "# " << (it->lamb_trans).size() << " " << it->name << " "
               << it->calibtyp << " " << ((it->id) + 1) << endl;

    // Loop over all the lambda of this filter
    vector<oneElLambda> oneFlt = (it->lamb_trans);
    for (vector<oneElLambda>::iterator it2 = oneFlt.begin(); it2 < oneFlt.end();
         ++it2) {
      stfiltfile << it2->lamb << "  " << it2->val << "  " << ((it->id) + 1)
                 << endl;
    }
  }
  return;
}

/*
read_doc_filters
original: 04/05/2015

create a vector of filters from the documentation file including all
informations as the AB correction, etc.

Input: name of the input file
Output: vector including all filters, and a boolan to indicate if the
documentatio file is present
*/
vector<flt> read_doc_filters(const string filtFile) {
  // the output
  vector<flt> allFilt;
  string line;

  // Open the documentation of the filter file
  string fltdoc = lepharework + "/filt/" + filtFile + ".doc";

  // Create a stream with the doc file
  ifstream stfltdoc;
  stfltdoc.open(fltdoc.c_str());
  // Check that the file exist.
  bool Fexist = stfltdoc.good();
  // Check that the file exist. Stop if not the case
  if (!(Fexist)) {
    throw invalid_argument(
        "Filter documentation does not exist in LEPHAREWORK/filt/" + filtFile);
  }
  // Open the filter file
  string fltfile = lepharework + "/filt/" + filtFile + ".dat";
  ifstream sfltIn;
  sfltIn.open(fltfile.c_str());

  // Read each line of the doc file
  while (getline(stfltdoc, line)) {
    // If the first character of the line is not #
    if (check_first_char(line)) {
      // Transfer the line content to line_stream.
      stringstream line_stream(line);

      // Generate one object "flt"
      flt oneFilt;

      // read informations on this filter and fill the object flt
      int nbLines;
      line_stream >> (oneFilt.name) >> (oneFilt.id) >> (oneFilt.lmean) >>
          (oneFilt.leff) >> (oneFilt.dwidth) >> (oneFilt.ab) >> (oneFilt.veg) >>
          (oneFilt.fcorr) >> (oneFilt.transtyp) >> (oneFilt.calibtyp) >>
          nbLines;

      // store all filters in one vector
      oneFilt.lamb_trans.resize(nbLines, oneElLambda(-999, -999, -999));
      allFilt.push_back(oneFilt);
    }
  }

  // Read the filter
  getline(sfltIn, line);
  // Loop over each filter
  for (size_t k = 0; k < allFilt.size(); k++) {
    // read the lambda/transmission
    allFilt[k].read(sfltIn);
  }

  /* Outout on the screen */

  cout << std::fixed;
  cout << setw(30) << left << "# NAME";
  cout << setw(8) << left << "IDENT";
  cout << setw(12) << right << "Lbda_mean";
  cout << setw(12) << right << "Lbeff(Vega)";
  cout << setw(12) << right << "FWHM ";
  cout << setw(10) << right << "AB-cor";
  cout << setw(10) << right << "VEGA";
  cout << setw(8) << right << "CALIB";
  cout << setw(12) << right << "Fac_corr" << endl;

  // Loop over all the filters
  for (vector<flt>::iterator it = allFilt.begin(); it < allFilt.end(); ++it) {
    cout << setw(30) << left << (it->name).substr(((it->name).rfind("/") + 1));
    cout << setw(8) << left << it->id;
    cout << setw(12) << right << setprecision(4) << it->lmean;
    cout << setw(12) << right << setprecision(4) << it->leff;
    cout << setw(12) << right << setprecision(4) << it->dwidth;
    cout << setw(10) << right << setprecision(4) << it->ab;
    cout << setw(10) << right << setprecision(4) << it->veg;
    cout << setw(8) << right << setprecision(2) << it->calibtyp;
    cout << setw(12) << right << setprecision(4) << it->fcorr << endl;
  }

  // Return the vector containing all filters
  return allFilt;
}

// Compute the various quantities related to the filter
void flt::compute_all() {
  // Compute fcorr, useful for FIR filters
  fcorrec();
  // Compute the AB correction
  abcorr();
}
