/*
original: 18/12/2013
updated:

Create a vector including the double elements of a c_string splitted acording to
delim

input s: string to be parsed
input delim: char used to parse the string
output: vector of double
*/
#include <sstream>
#include <string>
#include <vector>

using namespace std;

double strtodouble(const string &inputstring);

vector<double> split_double(const string s, char delim) {
  vector<double> elems;
  string item;
  stringstream ss(s);
  while (getline(ss, item, delim)) {
    elems.push_back(strtodouble(item));
  }
  return elems;
}
