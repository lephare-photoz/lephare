/*
  19/12/2013
  Structure and class related to the keyword
*/

// avoid multiple def of the same class
#ifndef KEYWORD_H  // check that this keyword has been set already
#define KEYWORD_H  // define the keyword to be checked

#include <algorithm>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

// read a keyword and convert it into a format given by T
template <typename T>
T stoT(const string val) {
  T value;
  stringstream ss(val);
  ss >> value;
  return value;
}

/// class keyword to store name and value of kweywords from the config file or
/// the command line
class keyword {
 public:
  string name,  ///< keyword name
      value;    ///< keyword value

  int def;  ///< indicate if the keyword is provided in the config or user input

  /* MP: added constructor and destructor*/
  keyword() {
    name = "No name";
    value = "";
    def = 0;
  };
  keyword(string n, string v) {
    name = n;
    value = v;
    def = 0;
    expand_path();
  }
  ~keyword(){};

  /// expand $HOME and $LEPHAREDIR if they are found in string value
  void expand_path();

  /// return default_val or value depending on def equal to or different from 0;
  template <typename T>
  T one(T default_val) {
    T result;

    // Case in which the keyword is not defined
    if (def == 0) {
      // replace its value by the default one
      result = default_val;
    } else {
      result = value;
    }

    return result;
  }

  template <typename T>
  vector<T> split(string default_val, int nbItems) {
    // split first in string
    vector<string> skey = split_string(default_val, nbItems);
    // Convert all the item string into double
    vector<T> dkey;
    dkey.reserve(skey.size());

    transform(skey.begin(), skey.end(), back_inserter(dkey),
              [](const string &val) { return stoT<T>(val); });

    return dkey;
  }

  /* Template functions are added */
  vector<string> split_string(string default_val, int nbItems);
  vector<int> split_int(string default_val, int nbItems);
  vector<long> split_long(string default_val, int nbItems);
  vector<double> split_double(string default_val, int nbItems);
  vector<bool> split_bool(string default_val, int nbItems);
  // template<class Te> vector<Te> split(string default_val, int nbItems);
};

typedef map<string, keyword> keymap;
keymap read_command(int argc, char *argv[]);
keymap read_config(const string &configfile);
keymap analyse_keywords(int argc, char *argv[], const string *list_keywords,
                        int nb_ref_key);

#endif
