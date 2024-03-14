/*

  27/11/14
  Implementation of the functions to read the various SED files depending on
  their format The output should be a vector of "SED" (one SED for a simple
  ascci file and a vector of SED when several templates are included (e.g. BC03)

*/

#include <cmath>
#include <fstream>   // print output file
#include <iostream>  // print standard file
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "SED.h"
#include "globals.h"
#include "oneElLambda.h"

using namespace std;

// prototypes
vector<GalSED> readBC03(string sedFile, int nummod, string type,
                        vector<double> &ageSel);
vector<GalSED> readPEGASE(string sedFile, int nummod, string type,
                          vector<double> &ageSel);
vector<bool> closeAge(vector<double> ageSel, vector<double> age);

/*
READ THE SED GENERATED WITH BRUZUAL & CHARLOT
*/
vector<GalSED> readBC03(string sedFile, int nummod, string type,
                        vector<double> &ageSel) {
  ifstream ssed;
  string lit, id, id2, id3;
  vector<GalSED> outSED;
  vector<double> age, w, h, f;
  double dage, dw, dh, df;
  int nages;

  // conversion  from Lum {Lo/A} to flux {erg/cm2/s/A}: fluxconv = 3.197e-07
  double convol = Lsol / (4 * pi * 100 * pow(pc, 2));

  // open the filter file into a stream
  ssed.open(sedFile.c_str());
  // Check if file is opened
  if (!ssed) {
    throw invalid_argument("Can't open file " + sedFile);
  }

  // Number of ages and read the ages
  ssed >> nages;
  for (int k = 0; k < nages; k++) {
    ssed >> dage;
    age.push_back(dage);
  }
  // Define the vector with the age to be selected
  vector<bool> useAge = closeAge(ageSel, age);
  // Few infos used later
  int iseg, jo, iop, inw, ik, ix;
  double ml, mu, xx, lm, um, baux, cn, cc, totm, totn, avs, tau, tau1, tau2,
      tau3, tau4;
  char stelib;
  ssed >> ml >> mu >> iseg;
  for (int k = 0; k < iseg; k++) {
    ssed >> xx >> lm >> um >> baux >> cn >> cc;
  }
  ssed >> totm >> totn >> avs >> jo >> tau >> tau1 >> tau2 >> tau3 >> tau4 >>
      iop >> stelib;
  // Skip the 3 lines
  getline(ssed, id);  // add this one but it should not be ... Probably the end
                      // of the previous line
  getline(ssed, id);
  double pos, zmett;
  string zme;
  pos = id.find("Z=");  // find the position of "Z="
  zme = id.substr(
      pos + 2, 10);  // read in id the string starting in pos+2 during 10 char.
  zmett = stod(zme, 0);  // to convert the string in double
  // Still two lines to be ignored
  getline(ssed, id2);
  getline(ssed, id3);

  // Number of wavelengths and read lambda axis
  ssed >> inw;
  for (int k = 0; k < inw; k++) {
    ssed >> dw;
    w.push_back(dw);
  }
  // Loop over the ages
  // Create an object SED and store it in the output vector
  for (int i = 0; i < nages; i++) {
    // Construct the object SED and store the various relevant informations
    GalSED oneSED(sedFile, tau, age[i], "BC03", nummod, type, i);

    // Attribute the metallicity to the SED
    oneSED.zmet = zmett;
    // Number of flux values
    ssed >> ik;
    if (ik != inw) {
      throw runtime_error(
          "Number of wavelength different from the number of flux in BC03 ");
    }
    // read the flux and generate lambda/flux of the SED
    for (int k = 0; k < ik; k++) {
      // read the flux and construct one lambda/flux element
      ssed >> dh;
      // convolve the flux by Lsun/(4pi*10pc^2)
      dh = dh * convol;
      oneElLambda litOne(w[k], dh, 1);
      // fill the lambda/trans values of the object SED
      (oneSED.lamb_flux).push_back(litOne);
    }
    oneElLambda litOne(10000000., 0., 1);
    (oneSED.lamb_flux).push_back(litOne);
    // Number of values of ?? and read them
    ssed >> ix;
    for (int k = 0; k < ix; k++) {
      ssed >> df;
      f.push_back(df);
    }
    // The core of the SED is created (its name, age, tau and lambda/flux),
    // physical parameters added later
    if (useAge[i]) {
      outSED.push_back(oneSED);
    }
  }

  /*
   Final physical information stored at the end of the BC03 file, include the
   relevant one in each SED ix represent the number of ages
  */
  double dbol, dmstr, dsf, devf, dsnr, dpnr, dbh, dsn, dwd, drm;
  vector<double> vmass, vsfr;
  ssed >> ix;
  for (int k = 0; k < ix; k++) {
    ssed >> dbol;
  }
  ssed >> ix;
  for (int k = 0; k < ix; k++) {
    ssed >> dmstr;
    if (useAge[k]) vmass.push_back(dmstr);
  }  // Mass-to-light ratio
  ssed >> ix;
  for (int k = 0; k < ix; k++) {
    ssed >> dsf;
    if (useAge[k]) vsfr.push_back(dsf);
  }  // star formation history
  ssed >> ix;
  for (int k = 0; k < ix; k++) {
    ssed >> devf;
  }
  ssed >> ix;
  for (int k = 0; k < ix; k++) {
    ssed >> dsnr;
  }
  ssed >> ix;
  for (int k = 0; k < ix; k++) {
    ssed >> dpnr;
  }
  ssed >> ix;
  for (int k = 0; k < ix; k++) {
    ssed >> dbh;
  }
  ssed >> ix;
  for (int k = 0; k < ix; k++) {
    ssed >> dsn;
  }
  ssed >> ix;
  for (int k = 0; k < ix; k++) {
    ssed >> dwd;
  }
  ssed >> ix;
  for (int k = 0; k < ix; k++) {
    ssed >> drm;
  }

  // Fill the mass and sfr values corresponding to each age
  for (unsigned i = 0; i < vmass.size(); i++) {
    outSED[i].mass = vmass[i];
    outSED[i].sfr = vsfr[i];
  }

  return outSED;
}

/*
READ THE SED GENERATED WITH PEGASE
*/
vector<GalSED> readPEGASE(string sedFile, int nummod, string type,
                          vector<double> &ageSel) {
  vector<GalSED> outSED;

  cout << " Need to implement the PEGASE format in SED reading " << endl;

  return outSED;
}

/*
SELECT THE AGE WHICH CAN BE KEPT
*/
vector<bool> closeAge(vector<double> ageSel, vector<double> age) {
  // Initialise the output vector
  vector<bool> outAge;
  // if the age for selection is not defined, include all ages
  if (ageSel.size() == 2) {
    for (vector<double>::iterator it = age.begin(); it < age.end(); ++it)
      outAge.push_back(true);
  } else {
    // Initialise the output vector at 0 (age not used)
    for (vector<double>::iterator it = age.begin(); it < age.end(); ++it)
      outAge.push_back(false);
  }

  // range in ages defined by the parameter AGE_RANGE
  double agemin = ageSel[0];
  double agemax = ageSel[1];

  // Loop over the ages to be selected
  for (vector<double>::iterator iti = ageSel.begin() + 2; iti < ageSel.end();
       ++iti) {
    double dage = 10.e10, diff;
    int k = 0, kmin = -99;
    // Loop over the ages to be examined
    for (vector<double>::iterator itj = age.begin(); itj < age.end(); ++itj) {
      // difference between selection and SED ages
      diff = abs(*itj - *iti);
      // keep the one which minimize the difference
      if (diff < dage && diff < 2.e9 && *itj >= agemin &&
          (*itj <= agemax || agemax <= 0)) {
        dage = diff;
        kmin = k;
      }
      k++;
    }
    // The age which minimizes the difference can be used
    if (kmin > 0) outAge[kmin] = true;
  }

  return outAge;
}
