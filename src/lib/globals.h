#ifndef GLOBAL_H
#define GLOBAL_H

#include <string>
#include <cmath>

#include <bitset>

#define MAX_CONTEXT 1024
#define INVALID_FLUX -9999.d
#define INVALID_VAL -9999.d
#define INVALID_MAG -9999.d
#define INVALID_PHYS -999.d
#define NULL_FLUX 0.d
#define HIGH_MAG 1000.d
#define HIGH_CHI2 1.e9d
#define EPS_Z 1.e-10d

using namespace std;

// extern tells the compiler this variable is declared elsewhere
extern string lepharedir,lepharework;
extern string fakeString;
extern const string nonestring;
extern int fakeInt;
extern const double c,ckms;
extern const double hplanck;
extern const double kboltzmann;
extern const double Lsol;
extern const double pc;
extern const double pi;
extern const double hc;
extern const double f_ga;


void get_lephare_env();

bool test_first_char(string maligne);

double blackbody(double T, double lambda);

int  bdincl(int n,long cont,int max);

inline string bool2string(const bool & b) {string sb; b?sb="YES": sb="NO"; return sb;}

inline bool CHECK_CONTEXT_BIT(unsigned long context, unsigned int n) {return std::bitset<MAX_CONTEXT>(context).test(n);}
inline double POW10D(double x) {return exp(2.302585092994046d*x);}
inline double POW10DSLOW(double x) {return pow(10.d,x);}
// This is a fast approximation to log2()
// Y = C[0]*F*F*F + C[1]*F*F + C[2]*F + C[3] + E;
inline double log2f_approx(double X) {
  double Y, F;
  int E;
  F = frexpf(fabsf(X), &E);
  Y = 1.23149591368684f;
  Y *= F;
  Y += -4.11852516267426f;
  Y *= F;
  Y += 6.02197014179219f;
  Y *= F;
  Y += -3.13396450166353f;
  Y += E;
  return(Y);
}
#define LOG10D LOG10D_SLOW
inline double LOG10D_SLOW(double x) {return log10(x);}
inline double LOG10D_FAST(double x) {return log2f_approx(x)*0.3010299956639812f;}
inline double mag2flux(double x, double zp=0.d) {return POW10D(-0.4d*(x + 48.6d - zp));}
inline double flux2mag(double x, double offset=0.d) {return -2.5d*LOG10D(x) - 48.6d + offset;}
#endif
