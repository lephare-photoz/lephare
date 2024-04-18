/*

  10/11/14
  Implementation of function for the keyword class (split in a vector of double,
  etc)

*/

#include "keyword.h"

#include <algorithm>
#include <fstream>   // print output file
#include <iostream>  // print standard file
#include <sstream>
#include <string>
#include <vector>

#include "globals.h"

using namespace std;

// Create a vector of double from the keyword value
// If the keyword is not defined, use the default value
vector<double> keyword::split_double(string default_val, int nbItems) {
  // split first in string
  vector<string> skey = split_string(default_val, nbItems);
  // Convert all the item string into double
  vector<double> dkey;
  dkey.reserve(skey.size());
  transform(
      skey.begin(), skey.end(), back_inserter(dkey), [](const string &val) {
        try {
          return stod(val);
        } catch (exception const &e) {
          throw invalid_argument(
              "Bad input parameter in a keyword. Waiting double rather than: " +
              val);
        }
      });

  return dkey;
}

// Create a vector of long from the keyword value
// If the keyword is not defined, use the default value
vector<long> keyword::split_long(string default_val, int nbItems) {
  // split first in string
  vector<string> skey = split_string(default_val, nbItems);
  // Convert all the item string into integer
  vector<long> ikey;
  ikey.reserve(skey.size());
  static string nm = name;
  transform(skey.begin(), skey.end(), back_inserter(ikey),
            [](const string &val) {
              try {
                return stol(val);
              } catch (exception const &e) {
                throw invalid_argument("Bad input parameter in keyword " + nm +
                                       ". Waiting long rather than: " + val);
              }
            });

  return ikey;
}

// Create a vector of integer from the keyword value
// If the keyword is not defined, use the default value
vector<int> keyword::split_int(string default_val, int nbItems) {
  // split first in string
  vector<string> skey = split_string(default_val, nbItems);
  // Convert all the item string into integer
  vector<int> ikey;
  ikey.reserve(skey.size());
  // Transform each string into integer, checking that it is feasible
  transform(skey.begin(), skey.end(), back_inserter(ikey),
            [](const string &val) {
              try {
                return stoi(val);
              } catch (exception const &e) {
                throw invalid_argument(
                    "Bad input parameter in a keyword. "
                    "Waiting integer rather than: " +
                    val);
              }
            });

  return ikey;
}

// Create a vector of string from the keyword value
// If the keyword is not defined, use the default value
vector<string> keyword::split_string(string default_val, int nbItems) {
  vector<string> elems;
  string item;

  // separator in a list of values
  char delim = ',';

  // Case in which the keyword is not defined
  if (value.empty()) {
    value = default_val;
  }
  stringstream ss(value);

  // Split it according to the delimitor
  while (getline(ss, item, delim)) {
    item.erase(std::remove(item.begin(), item.end(), '\t'), item.end());
    item.erase(std::remove(item.begin(), item.end(), ' '), item.end());
    elems.push_back(item);
  }

  // If only one element read in the keyword while we expect a larger number
  int nbElems = (int)elems.size();
  if (nbItems > 1 && nbElems == 1) {
    // Fill the vector with the missing ones by duplicating the first value
    for (int k = 1; k < nbItems; k++) {
      elems.push_back(elems[0]);
    };
  }

  // Check if the number of expected values correspond to the number of read
  // values
  nbElems = (int)elems.size();
  if (nbItems > 0 && nbElems != nbItems) {
    cout << "WRONG NUMBER OF ARGUMENTS FOR OPTION " << name << endl;
    cout << "We have " << nbElems << " instead of " << nbItems << endl;
    cout << "Use default value " << default_val << " for all filters " << endl;
    elems.clear();
    for (int k = 0; k < nbItems; k++) {
      elems.push_back(default_val);
    };
  }

  return elems;
}

vector<bool> keyword::split_bool(string default_val, int nbItems) {
  vector<bool> result;
  vector<string> strings = split_string(default_val, nbItems);
  for (int k = 0; k < nbItems; k++) {
    string s = strings[k];
    char t = toupper(s[0]);
    if (t == 'Y' || t == '1' || t == 'T') {
      result.push_back(true);
    } else {
      result.push_back(false);
    }
  }
  return result;
}

// Expand the environment variable if present in the keyword value
// For instance, if a keyword is defined as $HOME/bid.dat, this function replace
// HOME by its actual value
void keyword::expand_path() {
  static const string LEPHAREDIR = "LEPHAREDIR";
  static const string HOME = "HOME";

  // Check if the LEPHAREDIR environment variable is present at the beginning of
  // the string
  string str = string("LEPHAREDIR");
  // Search the string
  size_t found = value.find(str);
  // if found : replace by using the environment variable
  if (found != string::npos)
    value.replace(value.find(str) - 1, str.length() + 1,
                  getenv(LEPHAREDIR.c_str()));

  // Check if the HOME environment variable is present at the beginning of the
  // string
  str = string("HOME");
  // Search the string
  found = value.find(str);
  // if found : replace by using the environment variable
  if (found != string::npos)
    value.replace(value.find(str) - 1, str.length() + 1, getenv(HOME.c_str()));

  return;
}

keymap read_command(int argc, char *argv[]) {
  keymap result;
  result.clear();
  string lit, configfile;
  int nb_keywords(0);
  // Loop of each argument of the command line
  //  Expect one keyword, followed by a space and a string
  for (int k = 1; k < argc; k = k + 2) {
    // look for the config file
    if (strcmp(argv[k], "-c") == 0) {
      if (k + 1 <= argc) {
        configfile = argv[k + 1];
      } else {
        cout << "Miss the name of the para file after -c " << endl;
      }
    } else {
      // Normal procedure if it is not the config file
      // Look for the option line character "-"
      lit = argv[k];
      if (lit.at(0) == '-') {
        if (k + 1 < argc) {
          // Start at the second characters
          string name(lit.substr(1));
          if (name.at(0) == '-') {
            // make room for the case of args prepended with -- instead of -
            name = name.substr(1);
          }
          string value(argv[k + 1]);
          // store all the strings into a vector
          result[name] = keyword(name, value);
          nb_keywords++;
        } else {
          throw runtime_error(
              "Problem in format of the command line. Use one "
              "keyword followed by one value ");
        }
      }
    }
  }

  cout << "Number of keywords read at the command line (excluding -c config): "
       << nb_keywords << endl;
  // Stop if no config file at the command line
  if (configfile.length() == 0) {
    throw runtime_error("No configuration file indicated with the -c option.");
  } else {
    // store in the first keyword the name of the config file
    result["c"] = keyword("c", configfile);
  }

  return result;
}

// Read the configuration file (.para)
keymap read_config(const string &configfile) {
  // define a vector of keyword
  keymap result;
  string lit;

  // Stop if no config file at the command line
  ifstream configst;
  configst.open(configfile.c_str());
  if (!configst) {
    throw invalid_argument("ERROR: can't open config file " + configfile);
  }

  // Do that as long as the code can find a new line in the stream
  unsigned int nb_keywords(0);
  while (getline(configst, lit)) {
    // If the first character of the line is not #
    bool test = check_first_char(lit);
    if (test) {
      // put the line into the stream ss again
      stringstream ss(lit);

      // fill the keyword object
      string name(" "), value(" ");
      ss >> name >> value;
      if (value.back() == ',') {
        string value2;
        while (!ss.eof()) {
          ss >> value2;
          if (check_first_char(value2)) {
            value += value2;
          } else {
            break;
          }
        }
      }
      // store all the strings into a vector
      result[name] = keyword(name, value);
      nb_keywords++;
    }
  }
  cout << "Number of keywords read in the config file: " << nb_keywords << endl;
  configst.close();
  return result;
}

keymap analyse_keywords(int argc, char *argv[], const string *list_keywords,
                        int nb_ref_key) {
  keymap result, key_config;

  // Read the command line
  result = read_command(argc, argv);

  // Read the configuration file
  string configfile = result["c"].value;
  cout << "Reading keywords from " << configfile << endl;
  key_config = read_config(configfile);

  // add config into result. If key is present in command
  //  it is skipped, so this guarantees that command takes precedence
  result.insert(key_config.begin(), key_config.end());

  // finally, check the list_keywords and default keyword if they have not
  // been provided by the config or command
  for (int k = 0; k < nb_ref_key; k++) {
    if (result.find(list_keywords[k]) == result.end()) {
      cout << "Keyword " << list_keywords[k] << " not provided " << endl;
      result[list_keywords[k]] = keyword(list_keywords[k], "");
    }
  }
  return result;
}
