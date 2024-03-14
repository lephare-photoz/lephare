/*

  15/05/2015
  Implementation of the functions of the PDF class

*/

#include "PDF.h"

#include <algorithm>  // sort
#include <cmath>
#include <fstream>   // print output file
#include <iostream>  // print standard file
#include <sstream>
#include <string>
#include <vector>

#include "globals.h"

using namespace std;

/*
 CONSTRUCTOR OF THE PDF
*/
PDF::PDF(const double min, const double max, const size_t size) {
  // by construction
  scaleLinear = 1;
  scaleMin = min;
  scaleMax = max;
  vsize = size;
  xaxis.assign(size, 0.0);
  scaleStep = (max - min) / (vsize - 1);
  for (size_t k = 0; k < vsize; k++) {
    xaxis[k] = min + double(k) * scaleStep;
  }
  vPDF.assign(vsize, 0.0);
  chi2.assign(vsize, HIGH_CHI2);
  // Index of the SED corresponding to the minima or mode
  ind.assign(vsize, 0.0);
  // Store the inverse of scale step for computational time
  invScaleStep = 1.0 / scaleStep;

  return;
}

/*
 NORMALIZATION
*/
void PDF::normalization() {
  double somme = 0.;

  // Integral of the PDF
  for (size_t k = 0; k < vPDF.size() - 1; k++) {
    somme += (vPDF[k] + vPDF[k + 1]) / 2. * (xaxis[k + 1] - xaxis[k]);
  }

  // Renormalize the PDF
  if (somme > 0) {
    for (size_t k = 0; k < vsize; k++) {
      vPDF[k] = vPDF[k] / somme;
    }
  }

  return;
}

/*
 CONVERT CHI2 INTO PDF
*/
void PDF::chi2toPDF() {
  // Normalisation
  for (size_t k = 0; k < chi2.size(); k++) {
    // transform the chi2 into PDF
    vPDF[k] = exp(-0.5 * chi2[k]);
  }

  return;
}

/*
 INDEX OF THE MINIMUM CHI2 VALUE
*/
int PDF::chi2mini() {
  double chi2min = 1.e50;
  double indexMin = -99;

  // Loop over each value of the xaxis
  for (size_t k = 0; k < xaxis.size(); k++) {
    // if the chi2 is lower, store the new chi2 and the new index
    if (chi2[k] < chi2min) {
      chi2min = chi2[k];
      indexMin = k;
    }
  }

  return indexMin;
}

/*
 Interpolation into the chi2 distribution
 Description : Parabolic interpolation for minmum Chi2
 Origine     : Bevington book
*/
double PDF::int_parab() {
  double xmin = -99;
  double dxb, dxa, num, den;

  // Find the index of the minimum chi2
  int ib = this->chi2mini();
  // if the index is at the limit of the xaxis, can't do any interpolation
  if (ib == 0 || (size_t)ib == (xaxis.size() - 1)) {
    // keep the original value
    xmin = xaxis[ib];
  } else {
    // if one of the chi2 necessary for the interpolation k-1 -> k+1 too high,
    // can't interpolate
    if (chi2[ib - 1] > 1.e9 || chi2[ib] > 1.e9 || chi2[ib + 1] > 1.e9) {
      // keep the original value
      xmin = xaxis[ib];
    } else {
      // compute the best value from the parabolic interpolation of chi2
      // define the step for the interpolation (before and after the considered
      // point)
      dxb = xaxis[ib] - xaxis[ib - 1];
      dxa = xaxis[ib + 1] - xaxis[ib];
      num = (chi2[ib - 1] - chi2[ib]) * pow(dxa, 2.) -
            (chi2[ib + 1] - chi2[ib]) * pow(dxb, 2.);
      den = (chi2[ib - 1] - chi2[ib]) * dxa + (chi2[ib + 1] - chi2[ib]) * dxb;
      if (den == 0.e0 || abs(dxa) > 0.5 || abs(dxb) > 0.5) {
        xmin = xaxis[ib];
      } else {
        // parabolic interpolation possible
        // old version. True only if the x step is constant
        // xmin=xaxis[ib+1]-dxb*( (chi2[ib+1]-chi2[ib]) /
        // (chi2[ib+1]-2*chi2[ib]+chi2[ib-1]) +0.5); New generic version with
        // flexible step
        xmin = xaxis[ib] + 0.5 * num / den;
      }
    }
  }

  return xmin;
}

/*
 Find the values of xaxis crossing chiMin+dchi
*/
pair<double, double> PDF::uncMin(double dchi) {
  double chiLim;  // chi2 limit to be considered

  // Initialize the uncertainties with the range considered
  pair<double, double> result;
  result.first = xaxis[0];
  result.second = xaxis[xaxis.size() - 1];

  // Find the index of the minimum chi2
  int ib = this->chi2mini();

  // Define chi2 lim = best chi2 + delta chi2
  chiLim = chi2[ib] + dchi;
  // cout << " Start " << ib << "  " << chi2[ib] << "   " << xaxis[ib] << "   "
  // << chiLim << "   " << endl;

  // Loop over each point of the PDF to find the minimum error
  for (int k = 0; k < ib; k++) {
    // When the considered chi2 is below the chi2 limit and the following one
    // above
    if (chi2[k] >= chiLim && chi2[k + 1] < chiLim) {
      // Linear interpolation
      double slope = (chiLim - chi2[k]) / (chi2[k + 1] - chi2[k]);
      // store the minimum error
      result.first = xaxis[k] + slope * (xaxis[k + 1] - xaxis[k]);
      break;
    }
  }

  // Loop over each point of the PDF to find the maximum error
  for (size_t k = ib; k < xaxis.size(); k++) {
    // When the considered chi2 is above the chi2 limit and the following one
    // below
    if (chi2[k] <= chiLim && chi2[k + 1] > chiLim) {
      // Linear interpolation
      double slope = (chiLim - chi2[k]) / (chi2[k + 1] - chi2[k]);
      // store the minimum error
      result.second = xaxis[k] + slope * (xaxis[k + 1] - xaxis[k]);
      break;
    }
  }

  return result;
}

/*
 Derive the mode of the PDF, as well as the uncertainties
 Deal with asymetric PDF for the uncertainties
 First argument is the confidence level to consider, and the second argument is
 the considered value
 */
pair<double, double> PDF::credible_interval(float level, double val) {
  pair<double, double> result;
  vector<double>::iterator bound_val, bound_left_id, bound_right_id, boundpos,
      boundneg;

  // if levels given in percentage
  if (level > 1.) level /= 100.;

  // Find the index corrresponding to the value given in input
  bound_val = upper_bound(xaxis.begin(), xaxis.end(), val);
  size_t maxid = bound_val - xaxis.begin();
  // Take the index of the closest value
  if (maxid > 0) {
    if ((xaxis[maxid] - val) > (val - xaxis[maxid - 1])) maxid = maxid - 1;
  }

  // Compute the full cumulant
  vector<double> cumulant;
  double tmp = 0;
  cumulant.push_back(0.0);
  for (size_t k = 0; k < xaxis.size() - 1; k++) {
    tmp += (xaxis[k + 1] - xaxis[k]) * (vPDF[k + 1] + vPDF[k]) / 2.;
    cumulant.push_back(tmp);
  }

  // If cumulative defined
  if (cumulant.back() > 0) {
    // Normalization
    for (auto &c : cumulant) c /= cumulant.back();

    // - Use cumulant to compute integral from minimum x value to the mode
    double cumul_left = cumulant[maxid];
    // - Use cumulant to compute integral from the mode to the end
    double cumul_right = 1. - cumul_left;

    // Define the limits, taken into account the asymetry of the PDF
    double lowerLevel, upperLevel;
    // If integral of the PDF above the mode is higher than xval/2 and lower
    // than 0.5, enough area on each side
    if (cumul_right >= level / 2. && cumul_left >= level / 2.) {
      upperLevel = cumul_left + level / 2.;
      lowerLevel = cumul_left - level / 2.;
      // Case with the intergral between the min and the mode which is lower
      // than level/2
    } else if (cumul_right <= level / 2.) {
      // Integrate on the right in order to encompass 99.9% of the PDF
      upperLevel = 0.999;
      lowerLevel = 0.999 - level;
      // Case with the intergral between the min and the mode which is lower
      // than level/2
    } else {
      // Integrate on the left, leaving only 0.1% below the lower limit
      upperLevel = level - 0.001;
      lowerLevel = 0.001;
    }

    // Iterator pointing on the item with a value for the upper-level integral
    bound_right_id = upper_bound(cumulant.begin(), cumulant.end(), upperLevel);
    size_t indR = bound_right_id - cumulant.begin();
    // Linear interpolation - right case
    if (xaxis[indR - 1] > val) {
      // args in reverse order y, y1, x1, y2, x2, in order to get x instead of y
      result.second =
          linear_interp(upperLevel, cumulant[indR - 1], xaxis[indR - 1],
                        cumulant[indR], xaxis[indR]);
    } else {
      result.second = xaxis[indR];
    }
    // Linear interpolation - left case
    bound_left_id = upper_bound(cumulant.begin(), cumulant.end(), lowerLevel);
    size_t indL = bound_left_id - cumulant.begin() - 1;
    if (xaxis[indL + 1] < val) {
      result.first = linear_interp(lowerLevel, cumulant[indL], xaxis[indL],
                                   cumulant[indL + 1], xaxis[indL + 1]);
    } else {
      result.first = xaxis[indL];
    }

  } else {
    result.first = xaxis[0];
    result.second = xaxis.back();
  }
  return result;
}

/*
find the xaxis value corresponding to a level in the cumulative function
 */
double PDF::levelCumu2x(float level) {
  double result = -99.9;
  vector<double>::iterator bound_right_id;

  // if levels given in percentage
  if (level > 1.) level /= 100.;

  // Compute the full cumulant
  vector<double> cumulant;
  double tmp = 0;
  cumulant.push_back(0.0);
  for (size_t k = 0; k < xaxis.size() - 1; k++) {
    tmp += (xaxis[k + 1] - xaxis[k]) * (vPDF[k + 1] + vPDF[k]) / 2.;
    cumulant.push_back(tmp);
  }
  // If cumulative defined
  if (cumulant.back() > 0) {
    // Normalization
    for (auto &c : cumulant) c /= cumulant.back();

    // Iterator pointing on the item with a value for the upper-level integral
    bound_right_id = upper_bound(cumulant.begin(), cumulant.end(), level);
    size_t indR = bound_right_id - cumulant.begin();
    // Linear interpolation - right case
    if (indR > 0) {
      // args in reverse order y, y1, x1, y2, x2, in order to get x instead of y
      result = linear_interp(level, cumulant[indR - 1], xaxis[indR - 1],
                             cumulant[indR], xaxis[indR]);
    } else {
      result = xaxis[0];
    }
  }

  return result;
}

//
// Find the index in the PDF of the considered value
// If out of range, output a negative value
// Only grid with constant step accepted
size_t PDF::index(const double inVal) const {
  size_t res = 0;
  // If the value is well with the allowed scale
  if (inVal >= scaleMin && inVal <= scaleMax) {
    res = int(round((inVal - scaleMin) * invScaleStep));
  } else {
    // cout << " Value is out of range for the PDF:" << inVal << " " << scaleMin
    // << " " <<scaleMax<< endl;
    if (inVal < scaleMin) {
      res = 0;
    } else {
      res = vsize - 1;
    }
  }

  return res;
}

//
// Search for high peaks in The ML function vs xaxis
// And sort them from highest to smallest peaks in ML
void PDF::secondMax(const double win) {
  // generate an array of index
  // search for possible maximums
  // and define the windows
  int j = 0;
  vector<int> idxLim, idx, poss, possInd;
  vector<double> possX, possP;

  // Check for the possible maxima
  // At redshift 0
  if (vPDF[0] > vPDF[1]) {
    idx.push_back(0);
    poss.push_back(1);
    possX.push_back(xaxis[0]);
    possP.push_back(vPDF[0]);
    possInd.push_back(ind[0]);
    j++;
  }

  // Check for the possible maxima
  for (size_t i = 1; i < vsize - 1; ++i) {
    // Check if it could be a local maximum, the two adjacent bin should be
    // lower
    if (vPDF[i] > vPDF[i - 1] && vPDF[i] > vPDF[i + 1]) {
      idx.push_back(j);
      poss.push_back(1);
      possX.push_back(xaxis[i]);
      possP.push_back(vPDF[i]);
      possInd.push_back(ind[i]);
      j++;
    }
  }

  // Check for the possible maxima
  // At the maximum redshift
  if (vPDF[vsize - 1] > vPDF[vsize - 2]) {
    idx.push_back(j);
    poss.push_back(1);
    possX.push_back(xaxis[vsize - 1]);
    possP.push_back(vPDF[vsize - 1]);
    possInd.push_back(ind[vsize - 1]);
    j++;
  }

  // Sort according to the probability within the window
  sort(idx.begin(), idx.end(),
       [&](int i, int j) { return possP[i] > possP[j]; });

  // Loop over the potential maximums
  secondX.clear();
  secondP.clear();
  secondInd.clear();
  for (size_t k = 0; k < idx.size(); k++) {
    // If not in a window of previous solution, keep
    // Check in the previous best solution
    int already = 0;
    for (int j = k - 1; j >= 0; j--) {
      if (abs(possX[idx[j]] - possX[idx[k]]) < win && poss[idx[j]] == 1) {
        already = 1;
        poss[idx[k]] = 0;
      }
    }
    if (already == 0) {
      // Save the local maximum
      secondX.push_back(possX[idx[k]]);
      secondP.push_back(possP[idx[k]]);
      secondInd.push_back(possInd[idx[k]]);
    }
  }

  possX.clear();
  possP.clear();
  possInd.clear();

  return;
}
