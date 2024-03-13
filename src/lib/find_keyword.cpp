/*
original: 18/12/2013
updated:

Match the object of the class keyword found in the command line or in the config
file with the reference list of keywords

input: object of type keyword
output: position of the keyword in the reference list
*/
#include <cstring>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "keyword.h"

using namespace std;

int find_keyword(const keyword onekey, const string *list_keywords,
                 const int nb_ref_key) {
  // Find if the keyword correspond to anything in the keyword reference list
  int iarg = -1;

  // Need to use a c_string with the keyword name
  char *strname = new char[((onekey).name).length() + 1];
  std::strcpy(strname, ((onekey).name).c_str());

  for (int k = 0; k < nb_ref_key; k++) {
    // Need to use a c_string with the keyword from the reference list
    char *strkey = new char[(list_keywords[k]).length() + 1];
    strcpy(strkey, (list_keywords[k]).c_str());

    // store in iarg the position of the keyword in the reference list
    if (strcmp(strname, strkey) == 0) iarg = k;
  }
  return iarg;
}
