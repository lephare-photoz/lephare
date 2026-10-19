/*

  18/11/14
  Implementation of functions for the onElLambda class

*/

#include "oneElLambda.h"

#include <vector>

using namespace std;

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
