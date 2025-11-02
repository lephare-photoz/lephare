/*

  18/11/14
  Implementation of functions for the onElLambda class

*/

#include "oneElLambda.h"

// #include <fstream>   // print output file
#include <omp.h>

#include <algorithm>
#include <iostream>  // print standard file
#include <set>
#include <stdexcept>
#include <vector>

#include "globals.h"

using namespace std;

/*
  Linear interpolation in lambda
*/
void oneElLambda::interp(const oneElLambda& previousEl,
                         const oneElLambda& nextEl) {
  // normal linear interpolation possible (the two lambda are positive and the
  // next value is above the previous one
  if (previousEl.lamb >= 0 && nextEl.lamb > previousEl.lamb) {
    double slope =
        (nextEl.val - previousEl.val) / (nextEl.lamb - previousEl.lamb);
    val = previousEl.val + (lamb - previousEl.lamb) * slope;
  } else {
    // Interpolation is not possible -> put the value at the latest one and
    // origin at -99
    val = -99.;
    ori = -99;
  }

  return;
}

oneElVector interpLinearOpenMP_restrict(const oneElVector& in1,
                                        const oneElVector& in2)
// const std::vector<double>& x,
//      const std::vector<double>& y,
//      const std::vector<double>& z)
{
  // We want to interpolate in1 at the lambda values of in2

  // First we need to define the start
  //  if (x.size() != y.size() || x.size() < 2)
  //      throw std::invalid_argument("x et y doivent avoir la même taille (>=
  //      2).");
  //  if (!std::is_sorted(x.begin(), x.end()))
  //      throw std::invalid_argument("x doit être trié par ordre croissant.");

  const size_t m = in2.size();
  oneElVector result;

#pragma omp parallel for schedule(static)
  for (size_t j = 0; j < in2.size(); ++j) {
    double zi = in2[j].lamb;

    // restrict
    if (zi <= in1.front().lamb) {
      continue;
    }
    if (zi >= in1.back().lamb) {
      continue;
    }
    //        if (zi <= x.front()) { result[j] = y.front(); continue; }
    //        if (zi >= x.back())  { result[j] = y.back();  continue; }

    auto it = std::lower_bound(in1.begin(), in1.end(), in2[j]);
    size_t i = std::distance(in1.begin(), it) - 1;

    double t = (zi - in1[i].lamb) / (in1[i + 1].lamb - in1[i].lamb);
    double interpval = in1[i].val + t * (in1[i + 1].val - in1[i].val);
    result.emplace_back(zi, interpval, 0);
  }

  return result;
}

// --- interpolation linéaire (x croissant, sans extrapolation) ---
static inline double interp_linear_point(const std::vector<double>& x,
                                         const std::vector<double>& y,
                                         double xi) {
  if (xi <= x.front()) return y.front();
  if (xi >= x.back()) return y.back();

  auto it = std::lower_bound(x.begin(), x.end(), xi);
  size_t idx = std::distance(x.begin(), it);
  if (idx == 0) return y[0];
  if (idx >= x.size()) return y.back();

  double x0 = x[idx - 1], x1 = x[idx];
  double y0 = y[idx - 1], y1v = y[idx];
  double t = (xi - x0) / (x1 - x0);
  return y0 + t * (y1v - y0);
}

// --- interpolation vectorisée (parallèle OpenMP) ---
std::vector<double> interp_linear_vec(const std::vector<double>& x,
                                      const std::vector<double>& y,
                                      const std::vector<double>& q) {
  std::vector<double> out(q.size());
#pragma omp parallel for
  for (long long i = 0; i < (long long)q.size(); ++i)
    out[i] = interp_linear_point(x, y, q[i]);
  return out;
}

// --- création d’une grille régulière ---
std::vector<double> make_regular_grid(double lo, double hi, double dx) {
  if (dx <= 0) throw std::runtime_error("dx must be positive");
  size_t n = static_cast<size_t>((hi - lo) / dx) + 1;
  std::vector<double> grid(n);
#pragma omp parallel for
  for (long long i = 0; i < (long long)n; ++i) grid[i] = lo + i * dx;
  grid.back() = hi;  // pour assurer une fin exacte
  return grid;
}

// --- création de la grille issue de l’union de x1 et x2 ---
std::vector<double> make_union_grid(const std::vector<double>& x1,
                                    const std::vector<double>& x2, double lo,
                                    double hi) {
  std::set<double> s;
  for (double v : x1)
    if (v >= lo && v <= hi) s.insert(v);
  for (double v : x2)
    if (v >= lo && v <= hi) s.insert(v);
  return {s.begin(), s.end()};
}

// --- fonction principale unifiée ---
// dx > 0  → grille régulière
// dx < 0  → union triée de x1 et x2 dans l’intervalle commun
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
common_interpolate_combined(const std::vector<double>& x1,
                            const std::vector<double>& y1,
                            const std::vector<double>& x2,
                            const std::vector<double>& y2, double dx) {
  if (x1.size() != y1.size() || x2.size() != y2.size())
    throw std::runtime_error("x/y size mismatch");

  // Déterminer intervalle commun
  double lo = std::max(x1.front(), x2.front());
  double hi = std::min(x1.back(), x2.back());
  if (lo >= hi) return {{}, {}, {}};  // pas de recouvrement

  // Construire x' selon dx
  std::vector<double> x_common;
  if (dx > 0)
    x_common = make_regular_grid(lo, hi, dx);
  else
    x_common = make_union_grid(x1, x2, lo, hi);

  // Interpolation parallèle
  std::vector<double> y1_interp = interp_linear_vec(x1, y1, x_common);
  std::vector<double> y2_interp = interp_linear_vec(x2, y2, x_common);

  return {x_common, y1_interp, y2_interp};
}
