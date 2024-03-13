/*
original: 03/01/2014
updated:

Read the information file, i.e. the output of the non-para LF estimator

input LFst: stream from the input file
input/output LFall: vector of objects from the class "LFdata"
input zmoy: average redshift of the bin
*/
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "LFdata.h"

using namespace std;

// Declare propototypes
bool test_first_char(string maligne);

void read_LFdata(istream &LFst, vector<LFdata> &LFall, const double zmoy,
                 const int k, const double norma, const double normaerr) {
  string lit;
  LFdata LFone(
      zmoy, k, norma,
      normaerr);  // use the constructor which instances the redshift bin

  // Take the stream line by line
  while (getline(LFst, lit)) {
    // If the first character of the line is not #
    if (test_first_char(lit)) {
      // put the line into the stream ss again
      stringstream ss(lit);

      // fill the "LFdata" object
      ss >> LFone.M >> LFone.value >> LFone.errp >> LFone.errm;

      // Multiply everything by 1000. to get a better precision on phistar
      LFone.value = LFone.value + 3;

      // store all the strings into a vector
      LFall.push_back(LFone);
    }
  }
}
