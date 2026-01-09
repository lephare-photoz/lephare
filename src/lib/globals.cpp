/*
declaration of global variables
*/
#include "globals.h"

#include <dirent.h>  // deal with directory

#include <cmath>
#include <cstring>
#include <fstream>   // print output file
#include <iostream>  // print standard file
#include <string>

using namespace std;

string fakeString;
int fakeInt;
// lepharedir and lepharework as global variables
string lepharedir, lepharework;
// NONE
const string nonestring = "NONE";
// pi
const double pi = 3.14159265359;
// c
const double c = 2.99792458e18;  // en A/s
const double ckms = 2.99792458e5;
// h
const double hplanck = 6.62606957e-34;
// k Boltzmann
const double kboltzmann = 1.3806488e-23;
// L solar in erg/s
const double Lsol = 3.826e33;
// pc en cm
const double pc = 3.086e18;
// hc from Cedric
const double hc = 12398.42;  // [eV.A]
// f_ga from Cedric
const double f_ga = 1;

static const string LEPHAREWORK = "LEPHAREWORK";
static const string LEPHAREDIR = "LEPHAREDIR";

/*
10/11/2014
Get the environment variable LEPHAREWORK and LEPHAREDIR
Stop the code if LEPHAREDIR not defined
*/
// LCOV_EXCL_START
void get_lephare_env() {
  char const *temp;
  temp = getenv(LEPHAREDIR.c_str());
  if (temp != NULL) {
    lepharedir = string(temp);
  } else {
    cout << "Environment variable LEPHAREDIR not defined, need to stop.";
  }
  temp = getenv(LEPHAREWORK.c_str());
  if (temp != NULL) {
    lepharework = string(temp);
  } else {
    cout
        << "Environment variable LEPHAREWORK not defined, use $LEPHAREDIR/work."
        << endl;
    lepharework = lepharedir + "/work";
  }
}
// LCOV_EXCL_STOP

/*
13/12/2013
Function to test if the first character of a string corresponds to a comment
String starting with # or ! or - are considered as comment
If the string has a length of 0, return also false
input: string to analyse
output: boolean
*/
bool check_first_char(const string &maligne) {
  for (string::const_iterator it = maligne.begin(); it != maligne.end(); it++) {
    if (*it == ' ' || *it == '\t') {
      continue;
    } else if (*it != '#' && *it != '!') {
      return true;
    } else {
      return false;
    }
  }
  return false;
}

double blackbody(double T, double lambda) {
  double hckt = hplanck * c / kboltzmann / T;
  double val;
  if ((exp(hckt / lambda) - 1) < 1.e-5) {
    val = 1. / (hckt / lambda) / pow(lambda, 5.);
  } else {
    val = 1. / (exp(hckt / lambda) - 1.) / pow(lambda, 5.);
  }

  return val;
}

vector<size_t> indexes_in_vec(const double &value, const vector<double> &vec,
                              const float &precision) {
  vector<size_t> result;
  for (size_t i = 0; i < vec.size(); i++) {
    if (abs(vec[i] - value) < precision) result.push_back(i);
  }

  return result;
}  // LCOV_EXCL_LINE

vector<double> fast_interpolate(const std::vector<double> &x,
                                const std::vector<double> &y,
                                const std::vector<double> &z, double d) {
  std::vector<double> out;
  out.reserve(z.size());

  std::size_t i = 0;  // pointer in x,y

  // Loop through all target points z[k]
  for (double zk : z) {
    // Out of bounds → return d
    if (zk < x.front() || zk > x.back()) {
      out.push_back(d);
      continue;
    }
    // Advance i until x[i] <= zk <= x[i+1]
    while (i + 1 < x.size() && x[i + 1] < zk) {
      ++i;
    }
    // Now interpolate between (x[i], y[i]) and (x[i+1], y[i+1])
    double x0 = x[i], x1 = x[i + 1];
    double y0 = y[i], y1 = y[i + 1];
    double t = (zk - x0) / (x1 - x0);
    out.push_back(y0 + t * (y1 - y0));
  }
  return out;
}  // LCOV_EXCL_LINE
