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

/*
13/12/2013
Function to test if the first character of a string corresponds to a comment
String starting with # or ! or - are considered as comment
If the string has a length of 0, return also false
input: string to analyse
output: boolean
*/
bool check_first_char(string maligne) {
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

/*

  19/11/14
  Implementation of the black body function

*/
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

// DEPRECATED FOR CHECK_CONTEXT_BIT
int bdincl(int n, long cont, int max) {
  int res = 1;
  long sum = cont;

  // If the context is positive, otherwise accept the filter
  if (sum > 0) {
    // Start from the last filter and go down
    for (int i = max - 1; i >= n; i--) {
      // substract the context corresponding to the filter
      long summ = sum - pow(2, i);
      if (summ >= 0 &&
          i > n) {  // The filter is included into the context (summ>0) and the
                    // considered filter is not reached
        sum = summ;
      } else if (i == n) {  // Reached the considered filter

        if (summ >= 0) {  // The sum is still positive -> include the filter
          res = 1;
        } else {  // The sum is negative -> do not include the filter
          res = 0;
        }
        break;
      }
    }
  }
  return res;
}
