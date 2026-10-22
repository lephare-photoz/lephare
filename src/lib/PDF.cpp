/*

  15/05/2015
  Implementation of the functions of the PDF class

*/

#include "PDF.h"

#include <algorithm>  // sort
#include <cctype>
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
    for (auto& c : cumulant) c /= cumulant.back();
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
    for (auto& c : cumulant) c /= cumulant.back();

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

/* ------------------------------------------------------------------------ */
/*                       Redshift quality flag                               */
/* ------------------------------------------------------------------------ */

namespace {

string trim_lower(const string& s, bool lower) {
  size_t b = 0, e = s.size();
  while (b < e && isspace(static_cast<unsigned char>(s[b]))) b++;
  while (e > b && isspace(static_cast<unsigned char>(s[e - 1]))) e--;
  string out = s.substr(b, e - b);
  if (lower)
    for (auto& c : out)
      c = static_cast<char>(tolower(static_cast<unsigned char>(c)));
  return out;
}

double parse_number(const string& name, const string& val) {
  size_t pos = 0;
  double d;
  try {
    d = stod(val, &pos);
  } catch (const exception&) {
    throw invalid_argument("Z_FLAG_PAR: cannot read a number for '" + name +
                           "' from '" + val + "'");
  }
  if (pos != val.size() || !isfinite(d))
    throw invalid_argument("Z_FLAG_PAR: cannot read a number for '" + name +
                           "' from '" + val + "'");
  return d;
}

void check_range(const string& name, bool ok, const string& expect) {
  if (!ok)
    throw invalid_argument("Z_FLAG_PAR: '" + name + "' must be " + expect);
}

}  // namespace

FlagParams FlagParams::from_string(const string& str) {
  FlagParams p;
  const string all = trim_lower(str, true);
  if (all.empty() || all == "none") return p;

  vector<string> seen;
  size_t start = 0;
  while (start <= all.size()) {
    size_t comma = all.find(',', start);
    if (comma == string::npos) comma = all.size();
    const string item = trim_lower(all.substr(start, comma - start), true);
    start = comma + 1;
    if (item.empty()) continue;
    const size_t eq = item.find('=');
    if (eq == string::npos)
      throw invalid_argument("Z_FLAG_PAR: '" + item +
                             "' is not of the form name=value");
    const string name = trim_lower(item.substr(0, eq), true);
    const string val = trim_lower(item.substr(eq + 1), false);
    if (find(seen.begin(), seen.end(), name) != seen.end())
      throw invalid_argument("Z_FLAG_PAR: '" + name + "' is given twice");
    seen.push_back(name);
    const double d = parse_number(name, val);
    if (name == "nb_peak") {
      check_range(name, d >= 1 && d == floor(d), "an integer >= 1");
      p.nb_peak = static_cast<int>(d);
    } else if (name == "height") {
      check_range(name, d > 0 && d <= 1, "in ]0,1]");
      p.height = d;
    } else if (name == "peak_dz") {
      check_range(name, d >= 0, ">= 0");
      p.peak_dz = d;
    } else if (name == "tail") {
      check_range(name, d >= 0 && d <= 1, "in [0,1]");
      p.tail = d;
    } else if (name == "n_window") {
      check_range(name, d > 0, "> 0");
      p.n_window = d;
    } else if (name == "peak_ratio") {
      check_range(name, d >= 0, ">= 0");
      p.peak_ratio = d;
    } else if (name == "error") {
      check_range(name, d >= 0, ">= 0");
      p.error = d;
    } else if (name == "fwhm_level") {
      check_range(name, d > 0 && d < 1, "in ]0,1[");
      p.fwhm_level = d;
    } else {
      throw invalid_argument(
          "Z_FLAG_PAR: unknown parameter '" + name +
          "'. Valid names: nb_peak, height, peak_dz, tail, n_window, "
          "peak_ratio, error, fwhm_level");
    }
  }
  return p;
}

string FlagParams::to_string() const {
  ostringstream os;
  os.precision(10);
  os << "nb_peak=" << nb_peak << ",height=" << height << ",peak_dz=" << peak_dz
     << ",tail=" << tail << ",n_window=" << n_window
     << ",peak_ratio=" << peak_ratio << ",error=" << error
     << ",fwhm_level=" << fwhm_level;
  return os.str();
}

namespace {

double trapezoid(const double* y, const double* x, size_t n) {
  double s = 0.;
  for (size_t i = 1; i < n; i++)
    s += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
  return s;
}

/* Number of local maxima of v higher than hmin and at least `distance` samples
   apart. Same rules as scipy.signal.find_peaks(v, height=hmin,
   distance=distance): plateaus count once (middle point), the first and last
   samples are never peaks, and among too close peaks the highest is kept. */
int count_peaks(const double* v, size_t n, double hmin, size_t distance) {
  vector<size_t> peaks;
  size_t i = 1;
  const size_t i_max = n - 1;
  while (i < i_max) {
    if (v[i - 1] < v[i]) {
      size_t ahead = i + 1;
      while (ahead < i_max && v[ahead] == v[i]) ahead++;
      if (v[ahead] < v[i]) {
        const size_t mid = (i + ahead - 1) / 2;
        if (v[mid] >= hmin) peaks.push_back(mid);
        i = ahead;
      }
    }
    i++;
  }
  if (distance <= 1 || peaks.size() < 2) return static_cast<int>(peaks.size());

  const size_t np = peaks.size();
  vector<size_t> order(np);
  for (size_t k = 0; k < np; k++) order[k] = k;
  stable_sort(order.begin(), order.end(),
              [&](size_t a, size_t b) { return v[peaks[a]] < v[peaks[b]]; });
  vector<char> keep(np, 1);
  for (size_t k = np; k-- > 0;) {
    const size_t j = order[k];
    if (!keep[j]) continue;
    for (size_t m = j; m-- > 0 && peaks[j] - peaks[m] < distance;) keep[m] = 0;
    for (size_t m = j + 1; m < np && peaks[m] - peaks[j] < distance; m++)
      keep[m] = 0;
  }
  int count = 0;
  for (char c : keep) count += c;
  return count;
}

}  // namespace

ZFlagResult zflag(const double* x, const double* v, size_t n,
                  const FlagParams& p, double estimate) {
  ZFlagResult res;
  if (n < 3) return res;

  // The PDF must be a non negative, finite, non null function
  double total = 0.;
  for (size_t i = 0; i < n; i++) {
    if (!(v[i] >= 0.) || !isfinite(v[i])) return res;
    total += v[i];
  }
  const double area = trapezoid(v, x, n);
  if (!(total > 0.) || !(area > 0.) || !isfinite(area)) return res;

  vector<double> q(n), work(n);
  for (size_t i = 0; i < n; i++) q[i] = v[i] / area;
  const size_t imax = max_element(q.begin(), q.end()) - q.begin();
  const double qmax = q[imax];
  const double est = isfinite(estimate) ? estimate : x[imax];

  // Pseudo standard deviation around the estimate
  for (size_t i = 0; i < n; i++) work[i] = (x[i] - est) * (x[i] - est) * q[i];
  const double error = sqrt(trapezoid(work.data(), x, n));

  // Gaussian equivalent width from the width of the main peak at the level
  // fwhm_level*max: for a Gaussian width = 2 sigma sqrt(2 ln(1/level))
  const double level = p.fwhm_level * qmax;
  double xl = x[0], xr = x[n - 1];
  {
    size_t i = imax;
    while (i > 0 && q[i] >= level) i--;
    if (q[i] < level)
      xl = x[i] + (level - q[i]) / (q[i + 1] - q[i]) * (x[i + 1] - x[i]);
    i = imax;
    while (i < n - 1 && q[i] >= level) i++;
    if (q[i] < level)
      xr =
          x[i - 1] + (q[i - 1] - level) / (q[i - 1] - q[i]) * (x[i] - x[i - 1]);
  }
  const double dx = x[1] - x[0];
  const double sigma =
      max(xr - xl, dx) / (2. * sqrt(2. * log(1. / p.fwhm_level)));

  // Probability mass outside +- n_window sigma around the estimate
  double tail = 0.;
  if (sigma > dx) {
    const double bound = p.n_window * sigma;
    for (size_t i = 0; i < n; i++)
      work[i] = (x[i] < est - bound || x[i] > est + bound) ? q[i] : 0.;
    tail = trapezoid(work.data(), x, n);
  }

  // Mean to max ratio, and number of peaks
  double mean = 0.;
  for (size_t i = 0; i < n; i++) mean += q[i];
  mean /= static_cast<double>(n);
  const double ratio = mean / qmax;
  const size_t dist =
      max<size_t>(1, static_cast<size_t>(ceil(p.peak_dz / dx - 1.e-6)));
  const int npeaks = count_peaks(q.data(), n, p.height * qmax, dist);

  // Assemble the flag
  const double range = x[n - 1] - x[0];
  int flag = 0;
  if (error > p.error && error < range / 2) flag += 1;
  if (ratio > p.peak_ratio) flag += 2;
  if (tail > p.tail) flag += 4;
  if (npeaks >= p.nb_peak) flag += 8;
  if (error >= range / 2) flag += 16;

  res.flag = flag;
  res.estimate = est;
  res.error = error;
  res.peak_ratio = ratio;
  res.tail_mass = tail;
  res.n_peaks = npeaks;
  res.sigma = sigma;
  return res;
}
