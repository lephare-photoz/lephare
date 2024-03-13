/*
19/12/2013
convert string to double
*/
#include <sstream>
#include <string>
using namespace std;

double strtodouble(const string &inputstring) {
  istringstream instr(inputstring);
  double val;
  instr >> val;
  return val;
}
