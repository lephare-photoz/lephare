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

/// class keyword to store name and value of kweywords from the config file or
/// the command line
class keyword {
 public:
  string name,  ///< keyword name
      value;    ///< keyword value

  /* MP: added constructor and destructor*/
  keyword() {
    name = "No name";
    value = "";
  };
  keyword(string n, string v) {
    name = n;
    value = v;
    expand_path();
  }
  ~keyword(){};

  /// expand $HOME and $LEPHAREDIR if they are found in string value
  void expand_path();

  /*! split keyword value into an array of strings
   *
   * \param default_val: default value to use, if provided
   * \param nbItems: size of the expected array
   *
   * The delimiter is the comma.
   * If nbItems is negative, no check on the size of the
   * output array is made. If nbItems is positive and does not match the array
   obtained
   * from the input, a warning is printed out and an array of size
   * nbItems filled with default_val is returned.
   * If the input is a singleton, the returned array will be of size
   * nbItems filled with the singleton value
   */
  vector<string> split_string(string default_val, int nbItems);
  /// split into an array of strings, then convert to integer
  vector<int> split_int(string default_val, int nbItems);
  /// split into an array of strings, then convert to long integer
  vector<long> split_long(string default_val, int nbItems);
  /// split into an array of strings, then convert to double
  vector<double> split_double(string default_val, int nbItems);
  /// split into an array of strings, then convert to boolean
  vector<bool> split_bool(string default_val, int nbItems);
};

typedef map<string, keyword> keymap;
keymap read_command(int argc, char *argv[]);
keymap read_config(const string &configfile);
keymap analyse_keywords(int argc, char *argv[], const string *list_keywords,
                        int nb_ref_key);

#endif
