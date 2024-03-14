/*

  18/11/14
  Implementation of functions for the onElLambda class

*/

#include "oneElLambda.h"

#include <fstream>   // print output file
#include <iostream>  // print standard file

#include "globals.h"

using namespace std;

/*
  Linear interpolation in lambda
*/
void oneElLambda::interp(const oneElLambda &previousEl,
                         const oneElLambda &nextEl) {
  // normal linear interpolation possible (the two lambda are positive and the
  // next value is above the previous one
  if (previousEl.lamb > 0 && nextEl.lamb > previousEl.lamb) {
    double slope =
        (nextEl.val - previousEl.val) / (nextEl.lamb - previousEl.lamb);
    val = previousEl.val + (lamb - previousEl.lamb) * slope;
  } else {
    // Interpolation is not possible -> put the value at the latest one and
    // origin at -99
    val = previousEl.val;
    ori = -99;
  }

  return;
}
