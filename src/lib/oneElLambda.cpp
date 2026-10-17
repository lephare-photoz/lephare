/*

  18/11/14
  Implementation of functions for the onElLambda class

*/

#include "oneElLambda.h"

#include <algorithm>
#include <vector>

using namespace std;

std::pair<double, size_t> interp_linear_point(const std::vector<double>& x,
                                              const std::vector<double>& y,
                                              double xi, size_t x_hint) {
  if (xi < x.front() || xi > x.back()) return {0, x_hint};
  if (xi == x.front()) return {y.front(), x_hint};
  if (xi == x.back()) return {y.back(), x_hint};

  auto it = std::lower_bound(x.begin() + x_hint, x.end(), xi);
  size_t idx = std::distance(x.begin(), it);
  if (idx == 0) return {y[0], x_hint};
  // à tester le cas > ... Ca ne devrait jamais arriver....
  if (idx >= x.size()) return {y.back(), idx};

  double x0 = x[idx - 1], x1 = x[idx];
  double y0 = y[idx - 1], y1v = y[idx];
  double t = (xi - x0) / (x1 - x0);
  return {y0 + t * (y1v - y0), idx};
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
