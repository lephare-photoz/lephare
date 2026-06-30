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

pair<double, double> quadratic_extremum(double x1, double x2, double x3,
                                        double y1, double y2, double y3) {
  double x23 = x2 - x3;
  double x31 = x3 - x1;
  double x12 = x1 - x2;
  double delta = x23 * x1 * x1 + x31 * x2 * x2 + x12 * x3 * x3;
  double a = y1 * x23 + y2 * x31 + y3 * x12;
  if (a == 0. || delta == 0.) {  // possibly fragile
    throw invalid_argument("The three points are aligned");
  }
  double b = y1 * (x2 * x2 - x3 * x3) - y2 * (x1 * x1 - x3 * x3) +
             y3 * (x1 * x1 - x2 * x2);
  double c = y1 * x2 * x3 * x23 + y2 * x1 * x3 * x31 + y3 * x1 * x2 * x12;
  double ba = b / a;
  double xmin = 0.5 * ba;
  double ymin = -0.25 * ba * b / delta + c / delta;
  return make_pair(xmin, ymin);
}

/*
 CONSTRUCTOR OF THE PDF
*/
PDF::PDF(const double min, const double max, const size_t size) {
  if (size <= 0 || min >= max)
    throw std::invalid_argument(
        "Invalid argument : size must be greater than 0, and max must be "
        "greater than min.");

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
double PDF::normalization() {
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

  return somme;
}

vector<double> PDF::cumulant() {
  vector<double> cumulant;
  double tmp = 0.;
  cumulant.push_back(0.0);
  for (size_t k = 0; k < vsize - 1; k++) {
    tmp += (xaxis[k + 1] - xaxis[k]) * (vPDF[k + 1] + vPDF[k]) / 2.;
    cumulant.push_back(tmp);
  }
  return cumulant;
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

pair<double, double> PDF::improve_extremum(bool is_chi2) const {
  size_t id;
  double x[3];
  double y[3];

  // the extremum the precision of which is to be improved
  // is either the chi2 min or the PDF max
  id = is_chi2 ? min_element(chi2.begin(), chi2.end()) - chi2.begin()
               : max_element(vPDF.begin(), vPDF.end()) - vPDF.begin();

  // no quadratic search possible at boundaries
  if (id == 0 || id == vsize - 1)
    return is_chi2 ? make_pair(xaxis[id], chi2[id])
                   : make_pair(xaxis[id], vPDF[id]);

  for (size_t k = 0; k < 3; k++) {
    x[k] = xaxis[id - 1 + k];
    y[k] = is_chi2 ? chi2[id - 1 + k] : vPDF[id - 1 + k];
    // give up quadratic search if one of the chi2 is invalid
    if (is_chi2 && y[k] == HIGH_CHI2) return make_pair(xaxis[id], chi2[id]);
    // in the case of vPDF, it is a bit more problematic, but let's assume that
    // 0. is indeed a non acceptable value for the parabolic algorithm
    if (!is_chi2 && y[k] == 0.) return make_pair(xaxis[id], vPDF[id]);
  }

  // limit the separation between the three consecutive points
  // this can happen if the zgrid does not start at 0, and 0 is
  // added by the code; then the first bin in the zgrid can be too large.
  if (abs(x[1] - x[0]) > 0.5 || abs(x[2] - x[1]) > 0.5) {
    return make_pair(x[1], y[1]);
  }
  return quadratic_extremum(x[0], x[1], x[2], y[0], y[1], y[2]);
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
  for (size_t k = ib; k < xaxis.size() - 1; k++) {
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

pair<double, double> PDF::confidence_interval(float dchi) {
  // Initialize the uncertainties with the range considered
  pair<double, double> result = make_pair(xaxis[0], xaxis[vsize - 1]);

  // get the chi2 minimum with quadratic improvement
  auto extremum = improve_extremum(true);
  double zmin = extremum.first;
  double chi2min = extremum.second;

  size_t min_ib = index(zmin);
  double chiLim = chi2min + dchi;

  // Now we search for the position where the chi2 curve goes above chiLim,
  // right and left from index min_ib, which remains reasonable compared to zmin
  // position Loop over each point of the PDF to find the minimum error
  for (int k = 0; k < min_ib; k++) {
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
  for (size_t k = min_ib; k < xaxis.size() - 1; k++) {
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
 Derive the credible interval of the PDF, possibly asymetric
 First argument is the confidence level to consider, and the second argument is
 the considered value centering the interval (typically the mode in the redshift
 distribution)
 */
pair<double, double> PDF::credible_interval(float level, double val) {
  pair<double, double> result;
  vector<double>::iterator bound_val, bound_left_id, bound_right_id, boundpos,
      boundneg;

  // val is the centered value of the credible interval; it normally has to be
  // located inside the xaxis. But we have pdfmaps that will never be correct,
  // e.g. when BC03 is not used (no physical parameter estimation), and thus we
  // return an empty interval.
  // Likewise when level <= 0, which is not a meaningful value.
  // When the value is not meaningful (e.g. -99), put the same value for the
  // interval
  if (val < scaleMin || val >= scaleMax || level <= 0.) {
    return make_pair(val, val);
  }

  // if levels given in percentage
  if (level > 1.) level /= 100.;
  // we assume that level>1 means level provided as a percentage
  // but this means that now level should be <1

  // If level still not consistent with the expected value, put the full
  // range as credible interval
  if (level > 1.) return make_pair(scaleMin, scaleMax);

  // Compute the full cumulant
  vector<double> cumulant;
  double tmp = 0;
  // put the first item at 0
  cumulant.push_back(0.0);
  for (size_t k = 0; k < vsize - 1; k++) {
    tmp += (xaxis[k + 1] - xaxis[k]) * (vPDF[k + 1] + vPDF[k]) / 2.;
    cumulant.push_back(tmp);
  }

  // normalize if cumulative defined
  if (cumulant.back() <= 0) {
    return make_pair(xaxis.front(), xaxis.back());
  } else {
    for (auto &c : cumulant) c /= cumulant.back();
  }

  // PDF x axis is by construction uniformly sampled, thus so is the cumulant
  // vector. Then the closest to but lower than val index is immediately
  // computable. val<scaleMax => idx_lo < vsize-1 so that idx_hi is at most
  // vsize-1
  int idx_lo = (val - scaleMin) / scaleStep;
  double cum_lo = cumulant[idx_lo];
  int idx_hi = idx_lo + 1;
  double cum_hi = cumulant[idx_hi];
  // interpolated value of cumulant at xaxis value val
  double cum_val = cum_lo + (cum_hi - cum_lo) /
                                (xaxis[idx_hi] - xaxis[idx_lo]) *
                                (val - xaxis[idx_lo]);

  // - Use cumulant to compute integral from minimum x value to the mode
  double cumul_left = cum_val;
  // - Use cumulant to compute integral from the mode to the end
  double cumul_right = 1. - cumul_left;

  // Define the limits, taking into account the asymetry of the PDF
  double lowerLevel, upperLevel;

  // If integral of the PDF above the mode is higher than level/2 and lower
  // than 0.5, enough area on each side
  if (cumul_right >= level / 2. && cumul_left >= level / 2.) {
    upperLevel = cumul_left + level / 2.;
    lowerLevel = cumul_left - level / 2.;
    // Case with the integral between the min and the mode which is lower
    // than level/2
  } else if (cumul_right <= level / 2.) {
    // Integrate on the right in order to encompass the full PDF
    upperLevel = 1.;
    lowerLevel = 1. - level;
    // Case with the integral between the min and the mode which is lower
    // than level/2
  } else {
    // Integrate on the left, leaving only 0% below the lower limit
    upperLevel = level;
    lowerLevel = 0.0;
  }

  // Iterator pointing on the first item that has upperLevel < item
  bound_right_id = upper_bound(cumulant.begin(), cumulant.end(), upperLevel);
  size_t indR = bound_right_id - cumulant.begin();
  // Linear interpolation - right case
  if (cumulant[indR - 1] < upperLevel) {
    // args in reverse order y, y1, x1, y2, x2, in order to get x instead of y
    result.second = linear_interp(upperLevel, cumulant[indR - 1],
                                  xaxis[indR - 1], cumulant[indR], xaxis[indR]);
  } else {
    result.second = xaxis[indR - 1];
  }
  // Linear interpolation - left case
  bound_left_id = upper_bound(cumulant.begin(), cumulant.end(), lowerLevel);
  size_t indL = bound_left_id - cumulant.begin();
  if (cumulant[indL - 1] < lowerLevel) {
    result.first = linear_interp(lowerLevel, cumulant[indL - 1],
                                 xaxis[indL - 1], cumulant[indL], xaxis[indL]);
  } else {
    result.first = xaxis[indL - 1];
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
// Only grid with constant step accepted
size_t PDF::index(const double inVal) const {
  size_t res = 0;
  // If the value is well with the allowed scale
  if (inVal >= scaleMin && inVal <= scaleMax) {
    res = int(round((inVal - scaleMin) * invScaleStep));
  } else {
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
