/*

  16/12/14
  Implementation of the functions of the cosmo class

*/

#include "cosmology.h"

#include <algorithm>
#include <cmath>
#include <fstream>   // print output file
#include <iostream>  // print standard file
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "globals.h"

using namespace std;

/*
  COMPUTE THE DISTANCE MODULUS
*/
double cosmo::distMod(double z) const {
  double funz;

  if (z < 1.e-10) {
    funz = 0.;
  } else {
    funz = 5. * log10((1 + z) * distMet(z)) + 25;
  }

  return funz;
}

/*
  Compute the metric distance dmet in Mpc
  dlum = dmet*(1+z)
  dang = dmet/(1+z) = dlum/(1+z)^2
*/
double cosmo::distMet(double z) const {
  double dmet, ao;

  // case without the cosmological constant
  if (l0 == 0) {
    //  ao = c/(h0*sqrt(ABS(1-omt)))
    //  in fact we use x = ao * x(z) with x(z) from eq 8 of
    //  Moscardini et al.  So we don't need to compute ao
    if (om0 > 0) {
      ao = 1.;
      dmet = om0 * z - (om0 - 2) * (1 - sqrt(1 + om0 * z));
      dmet = 2 * ckms / (ao * h0 * om0 * om0 * (1 + z)) * dmet;
    } else {
      ao = 1.;
      dmet = ckms * z * (1. + z / 2.) / (h0 * (1 + z));
    }

  } else if (om0 < 1 && l0 != 0) {
    ao = 1.;
    double sum = 0.;
    double dz = z / 50.;
    for (int i = 0; i < 50; i++) {
      double zi = (double(i) + 0.5) * dz;
      double Ez = sqrt(om0 * pow((1. + zi), 3.) +
                       (1 - om0 - l0) * pow((1. + zi), 2.) + l0);
      sum = sum + dz / Ez;
    }
    dmet = ckms / (h0 * ao) * sum;

  } else {
    throw runtime_error("Cosmology not included : h0=" + to_string(h0) +
                        " Om0=" + to_string(om0) + " l0=" + to_string(l0));
  }

  return dmet;
}

/*
  Compute cosmological time from z=infinty  to z
  as a function of cosmology.  Age given in year !!
  Note : use lambda0 non zero if Omega_o+lambda_o=1
*/
double cosmo::time(double z) const {
  double timy = 0., val;
  double hy = h0 * 1.0224e-12;

  if (abs(om0 - 1) < 1.e-6 && l0 == 0) {
    timy = 2. * pow((1 + z), -1.5) / (3 * hy);

  } else if (om0 == 0 && l0 == 0) {
    timy = 1. / (hy * (1 + z));

  } else if (om0 < 1 && om0 > 0 && l0 == 0) {
    val = (om0 * z - om0 + 2.) / (om0 * (1 + z));
    timy = 2. * sqrt((1 - om0) * (om0 * z + 1)) / (om0 * (1 + z));
    timy = timy - log10(val + sqrt(val * val - 1));
    timy = timy * om0 / (2. * hy * pow((1 - om0), 1.5));

  } else if (om0 > 1 && l0 == 0) {
    timy = acos((om0 * z - om0 + 2.) / (om0 * (1 + z)));
    timy = timy - 2 * sqrt((om0 - 1) * (om0 * z + 1)) / (om0 * (1 + z));
    timy = timy * om0 / (2 * hy * pow((om0 - 1), 1.5));

  } else if (om0 < 1 && abs(om0 + l0 - 1) < 1.e-5) {
    val = sqrt(1 - om0) / (sqrt(om0) * pow((1 + z), 1.5));
    timy = log(val + sqrt(val * val + 1));
    timy = timy * 2. / (3. * hy * sqrt(1 - om0));

  } else {
    throw runtime_error(" Not the right cosmology to derive the time ");
  }

  return timy;
}

// Two possible grid in redshift : linear or in (1+z)
// Possible now to define a minimum redshift to be considered
vector<double> zgrid(int gridType, double dz, double zmin, double &zmax) {
  if (zmax < zmin) {
    throw invalid_argument(
        "You are probably using the old parametrisation of "
        "Z_STEP since Z MIN > Z MAX in Z_STEP.");
  }

  vector<double> z;

  // first redshift at 0
  z.push_back(0.);

  // Start at zmin
  if (zmin > 0) z.push_back(zmin);

  unsigned int count = 1;
  double zinter = zmin;
  // Define a vector with the redshift grid according to the given type
  switch (gridType) {
      // grid in dz*(1+z)
    case 1:
      while (zinter < zmax) {
        // Step in dz*(1+z)
        zinter = zinter + (1. + zinter) * dz;
        // keep only in the redshift range zmin-zmax defined in Z_STEP
        if (zinter > zmin && zinter < zmax) z.push_back(zinter);
      }
      break;
      // Linear grid in redshift
    default:
      while (zinter < zmax) {
        // Step in dz
        zinter = zmin + count * dz;
        // keep only in the redshift range zmin-zmax defined in Z_STEP
        if (zinter > zmin && zinter < zmax) {
          z.push_back(zinter);
        }
        count++;
      }
      break;
  }
  z.push_back(zmax);

  return z;
}

int indexz(const double red, const vector<double> &gridz) {
  // gridz is assumed to be sorted
  // case 1 : red <= gridz[0] : return index 0
  if (red <= gridz.front()) return 0;
  // case 2 : red >= gridz[size-1] : return index size-1
  if (red >= gridz.back()) return gridz.size() - 1;
  // case 3 : red = gridz[k] : return k
  auto val = find(gridz.begin(), gridz.end(), red);
  //  auto val = binary_search(gridz.begin(), gridz.end(), red);//if gridz
  //  guaranteed sorted
  if (val != gridz.end()) return val - gridz.begin();
  // case 4 : gridz[0] < red < gridz[size-1] : find closest match in gridz and
  // return its index
  auto up = lower_bound(gridz.begin(), gridz.end(), red);
  auto low = up - 1;
  if (abs(*up - red) <= abs(red - *low))
    return up - gridz.begin();
  else
    return low - gridz.begin();
}
