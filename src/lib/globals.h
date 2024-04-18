#ifndef GLOBAL_H
#define GLOBAL_H

#include <bitset>
#include <cmath>
#include <string>

#define MAX_CONTEXT 1024
#define INVALID_FLUX -9999.0
#define INVALID_VAL -9999.0
#define INVALID_MAG -9999.0
#define INVALID_PHYS -999.0
#define NULL_FLUX 0.0
#define HIGH_MAG 1000.0
#define HIGH_CHI2 1.0e9
#define EPS_Z 1.0e-10

using namespace std;

// extern tells the compiler this variable is declared elsewhere
extern string lepharedir, lepharework;
extern string fakeString;
extern const string nonestring;
extern int fakeInt;
extern const double c, ckms;
extern const double hplanck;
extern const double kboltzmann;
extern const double Lsol;
extern const double pc;
extern const double pi;
extern const double hc;
extern const double f_ga;

void get_lephare_env();

bool check_first_char(string maligne);

double blackbody(double T, double lambda);

int bdincl(int n, long cont, int max);

inline string bool2string(const bool &b) {
  string sb;
  b ? sb = "YES" : sb = "NO";
  return sb;
}

inline bool CHECK_CONTEXT_BIT(unsigned long context, unsigned int n) {
  return std::bitset<MAX_CONTEXT>(context).test(n);
}
inline double POW10D(double x) { return exp(2.302585092994046 * x); }
inline double POW10DSLOW(double x) { return pow(10.0, x); }
// This is a fast approximation to log2()
// Y = C[0]*F*F*F + C[1]*F*F + C[2]*F + C[3] + E;
inline double log2f_approx(double X) {
  double Y, F;
  int E;
  F = frexp(fabs(X), &E);
  Y = 1.23149591368684;
  Y *= F;
  Y += -4.11852516267426;
  Y *= F;
  Y += 6.02197014179219;
  Y *= F;
  Y += -3.13396450166353;
  Y += E;
  return (Y);
}
#define LOG10D LOG10D_SLOW
inline double LOG10D_SLOW(double x) { return log10(x); }
inline double LOG10D_FAST(double x) {
  return log2f_approx(x) * 0.3010299956639812f;
}
inline double mag2flux(double x, double zp = 0.0) {
  return POW10D(-0.4 * (x + 48.6 - zp));
}
inline double flux2mag(double x, double offset = 0.0) {
  return -2.5 * LOG10D(x) - 48.6 + offset;
}
#endif
