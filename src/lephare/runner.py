import argparse
import os
import sys
import time
from contextlib import suppress

from ._lephare import get_lephare_env, keyword

__all__ = [
    "Runner",
]


class Runner:
    """

    Runner is the base class of the classes Filter, Sedtolib, MagGal, and Zphota,
    responsible to drive the execution of the filter, sedtolib, mag_gal, and
    `zphota` scripts, respectively.

    Configuration of the classes uses either a configuration file, a dictionary of key names
    and corresponding keyword objects, or set of key/values passed as arguments to the constructor

    The order of precedence is : key/values pairs overwrite config_keymap entries, which
    themselves overwrite config_file entries.

    Arguments
    ----------
         config_keys : `list` or `None`, optional
             List of all admissible configuration keys, provided by the inheriting classes
         config_file : `string` or `None`, optional
             Path to a configuration file in LePHARE .para format
         config_keymap : `dict` or `None`, optional
             Dictionary of configuration values provided as 'key/keyword object' pairs.
    """

    def __init__(self, config_keys=None, config_file=None, config_keymap=None, **kwargs):
        # set the LEPHAREDIR and LEPHAREWORK env variable
        get_lephare_env()
        self.keymap = {}
        self.config = ""
        self.verbose = False
        self.typ = None
        self.timer = False

        # Check that the relevant keyword names are defined
        if config_keys is None:
            raise RuntimeError("Runner is a base class and cannot be initialized")
        self.config_keys = config_keys

        # a config_file is passed as argument to the Python constructor
        # note : VERBOSE can be in the config_file, but note type
        if config_file is not None:
            self.keymap = self.parse_config_file(config_file)

        if config_keymap is not None:
            # check validity of config_keymap entries
            self.validate_config_dict(config_keymap)
            # if config_file is not provided but config_keympa is, then use
            # config_keymap. If both are provided, merge them, with
            # config_keymap entries overriding in case of duplicates
            self.keymap = self.keymap | config_keymap
            if "VERBOSE" in self.keymap:
                self.verbose = bool(self.keymap["VERBOSE"].value == "YES")

        # Finally, take the directly provided arguments, which supersede
        # both config_file and config_keymap inputs
        self.validate_config_dict(kwargs, no_raise=False)
        for k, v in kwargs.items():
            uk = k.upper()
            if uk == "TYP":
                self.typ = v.upper()
            # allow verbose to be passed as a bool arg, more pythonesque
            if uk == "VERBOSE":
                if v.__class__ is bool:
                    self.verbose = v
                    kwargs[k] = "YES" if v else "NO"
                else:
                    self.verbose = bool(v.upper() == "YES")

            self.keymap[uk] = keyword(uk, str(kwargs[k]))

        # If self.keymap is still empty, it either means that
        # 1/ an instance of Runner or its inherited classes was called from an executable,
        # 2/ or that it was instantiated in a python session without argument
        if self.keymap == {} and not hasattr(sys, "ps1"):
            self.args = self.config_parser()

        self.update_help()
        return

    # This function takes the config file, read it line per linem and output a keyword map
    def parse_config_file(self, filename):
        """Load config file and set config values.

        Parameters
        ----------
        filename : `string`
            Path to config file
        """
        if not os.path.exists(filename):
            raise RuntimeError("File %s not found" % filename)
        self.config = filename
        keymap = {}

        def config_reader(filename):
            for row in open(filename, "r"):  # noqa: SIM115
                yield row

        config = config_reader(filename)
        for line in config:
            if line[0] != "#" and line != "\n" and not line.isspace():
                splits = line.split()
                splits[0].lstrip()  # remove leading spaces
                if splits[0] != "#" and splits[0] in self.config_keys:
                    try:
                        keymap[splits[0]] = keyword(splits[0], splits[1])
                    except:  # noqa: E722
                        keymap[splits[0]] = keyword(splits[0], "")

                    if splits[0] in self.config_keys and hasattr(self, "parser"):  # pragma no cover
                        self.parser.set_defaults(
                            **{
                                splits[0]: splits[1],
                            }
                        )
        return keymap

    def config_parser(self):
        """Create command line config parser from list of keys"""

        self.parser = argparse.ArgumentParser()
        # these two arguments are not part of the config keys: config is provided directly in python,
        # and timer is dealt with differently in a python ecosystem.
        self.parser.add_argument("--timer", help="switch timer on to time execution", action="store_true")
        self.parser.add_argument("-c", "--config", type=str, default="", help="Path to config file.")
        # add config keys as authorized command line args:
        self.add_authorized_keys()
        args, unknown = self.parser.parse_known_args()

        keymap = {}
        if args.config != "":
            keymap = self.parse_config_file(args.config)

        # capture specific args if they is passed as script argument
        with suppress(Exception):
            self.typ = args.typ

        with suppress(Exception):
            self.timer = args.timer

        # copy the args keywords back into the keymap,
        # setting value to "" for verbose, timer and the like
        for key in self.config_keys:
            if key == "verbose":
                self.verbose = args.verbose
                keymap[key] = keyword(key, "YES") if self.verbose else keyword(key, "NO")
                continue
            try:
                keymap[key] = keyword(key, getattr(args, key))
            except:  # noqa: E722
                if key not in keymap:
                    keymap[key] = keyword(key, "")

        self.keymap = keymap
        return args

    def run(self, **kwargs):
        if self.timer:
            self.start = time.time()

        # update keymap and verbosity based on call arguments
        # this is only when the code is called from python session
        v = kwargs.pop("verbose", self.verbose)
        if v.__class__ is bool:
            self.verbose = v
            self.keymap["VERBOSE"] = keyword("VERBOSE", "YES") if v else keyword("VERBOSE", "NO")
        if v.__class__ is str:
            if v.upper() == "YES":
                self.keymap["VERBOSE"] = keyword("VERBOSE", "YES")
                self.verbose = True
            elif v.upper() == "NO":
                self.keymap["VERBOSE"] = keyword("VERBOSE", "NO")
                self.verbose = False
            else:
                raise RuntimeError("verbose is to be set to False/True or YES/yes//NO/no")

        if "typ" in kwargs:
            self.typ = kwargs.pop("typ", self.typ)
            self.typ = self.typ.upper()
        for k, v in kwargs.items():
            if k.upper() in self.keymap:
                self.keymap[k.upper()] = keyword(k.upper(), str(v))

    def end(self):
        self.stop = time.time()
        if self.timer:
            print("execution time: %.4g" % (self.stop - self.start))

    def add_authorized_keys(self):
        """Add authorized keys in config keys to the parser"""
        for k, v in self.config_keys.items():
            if k == "typ":
                self.parser.add_argument("-t", "--typ", type=str, default="", help=v, required=True)
            elif k == "verbose":
                self.parser.add_argument("--verbose", help="increase onscreen verbosity", action="store_true")
            else:
                self.parser.add_argument("--%s" % k, type=str, metavar="", help=v)

    def validate_config_dict(self, input_dict, no_raise=True):
        """Check that the input dictionary match the licit arguments"""
        for key in input_dict.copy():
            if key.upper() not in self.config_keys:
                if no_raise:
                    pass
                else:
                    raise RuntimeError(f"{key} is not a recognized argument of {self.__class__.__name__}.")

    def update_help(self):
        """Method to be overloaded by ineriting classes"""
        doc = "Runner is a base class that needs to be inherited from"
        with suppress(Exception):
            self.parser.usage = doc
        self.__doc__ = doc


if __name__ == "__main__":
    raise TypeError("Runner is a base class, not an executable")  # pragma: no cover
