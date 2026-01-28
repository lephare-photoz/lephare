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

double interp_linear_point(const std::vector<double>& x,
                           const std::vector<double>& y, double xi) {
  if (xi < x.front() || xi > x.back()) return 0;
  if (xi == x.front()) return y.front();
  if (xi == x.back()) return y.back();

  auto it = std::lower_bound(x.begin(), x.end(), xi);
  size_t idx = std::distance(x.begin(), it);
  if (idx == 0) return y[0];
  // à tester le cas > ... Ca ne devrait jamais arriver....
  if (idx >= x.size()) return y.back();

  double x0 = x[idx - 1], x1 = x[idx];
  double y0 = y[idx - 1], y1v = y[idx];
  double t = (xi - x0) / (x1 - x0);
  return y0 + t * (y1v - y0);
}

// --- création d’une grille régulière ---
// J'ai défini un code qui généralise l'ancien resample et peut
// creer une grille uniforme si besoin, plutôt que de rassembler
// les points des deux vecteur d'input.
// Cela peut disparaitre à terme
std::vector<double> make_regular_grid(double lo, double hi, double dx) {
  if (dx <= 0) throw std::runtime_error("dx must be positive");
  size_t n = static_cast<size_t>((hi - lo) / dx) + 1;
  std::vector<double> grid(n);
  // #pragma omp parallel for
  for (long long i = 0; i < (long long)n; ++i) grid[i] = lo + i * dx;
  grid.back() = hi;  // ensure strict ending at hi
  return grid;
}

// --- création de la grille issue de l’union de x1 et x2 ---
// Ca c'est le cas resample classique de lephare. C'est chatgpt
// qui m'a sorti l'usage de std::set. Je pense que son utilisation
// impose de ne pas paralléliser avec OMP, mais je n'ai pas creusé
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

std::tuple<std::vector<double>, std::vector<double>, std::vector<double>>
common_interpolate_combined(const std::vector<double>& x1,
                            const std::vector<double>& y1,
                            const std::vector<double>& x2,
                            const std::vector<double>& y2, double dx) {
  if (x1.size() != y1.size() || x2.size() != y2.size())
    throw std::runtime_error("x/y size mismatch");

  // determine common interval
  double lo = std::max(x1.front(), x2.front());
  double hi = std::min(x1.back(), x2.back());
  if (lo >= hi) return {{}, {}, {}};  // pas de recouvrement

  // build the grid
  std::vector<double> x_common;
  if (dx > 0)
    x_common = make_regular_grid(lo, hi, dx);
  else
    x_common = make_union_grid(x1, x2, lo, hi);

  // interpolate
  std::vector<double> y1_interp = interp_linear_vec(x1, y1, x_common);
  std::vector<double> y2_interp = interp_linear_vec(x2, y2, x_common);

  return {x_common, y1_interp, y2_interp};
}
