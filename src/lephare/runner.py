import argparse
import os
import time

from lephare._lephare import get_lephare_env, keyword

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

        # Check that the relevant keyword names are defined
        if config_keys is None:
            raise RuntimeError("Runner is a base class and cannot be initialized")
        self.config_keys = config_keys

        if config_file is not None:
            # this only happens if the code is called from python
            # Read the config file and creates a keyword map
            self.parse_config_file(config_file)
        if config_keymap is not None:
            # merge the config_file and config_keymap, keeping the config_keymap in case of duplicate
            self.keymap = self.keymap | config_keymap

        for key in kwargs:
            ukey = key.upper()
            if ukey in self.config_keys:
                if key == "type":
                    self.typ = kwargs[key]
                else:
                    self.keymap[ukey] = keyword(ukey, str(kwargs[key]))
            else:
                raise RuntimeError(f"{key} is not a recognized argument of {self.__class__.__name__}.")

        if config_keymap is None and config_file is None and kwargs == {}:
            # this only happens if the code is called as an executed script
            # Consider the keywords given in the line command
            self.args = self.config_parser()
            if self.timer:
                self.start = time.time()
        # set verbosity. check keymap is not set on the commandline.
        if not self.verbose and "VERBOSE" in self.keymap:
            self.verbose = self.keymap["VERBOSE"].split_bool("NO", 1)[0]

    # This function take the config file, read it line per linem and output a keyword map
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
                if splits[0] != "#":
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

        self.keymap = keymap

    def config_parser(self):
        """Create command line config parser from list of keys"""

        self.parser = argparse.ArgumentParser()
        self.parser.add_argument("--verbose", help="increase onscreen verbosity", action="store_true")
        self.parser.add_argument("--timer", help="switch timer on to time execution", action="store_true")
        self.parser.add_argument("-c", "--config", type=str, default="", help="Path to config file.")
        if self.__class__.__name__ in ["Sedtolib", "MagGal"]:
            self.parser.add_argument(
                "-t",
                "--typ",
                type=str,
                default="",
                help="Type of SED object : GAL QSO or STAR",
                required=True,
            )

        # add authorized command line args:
        self.add_authorized_keys()
        args, unknown = self.parser.parse_known_args()

        if args.config != "":
            self.parse_config_file(args.config)
        else:
            self.keymap = {}
        try:  # noqa: SIM105
            # capture the type if it is passed as script argument
            self.typ = args.typ
        except:  # noqa: E722
            pass

        self.verbose = args.verbose
        self.timer = args.timer
        # copy the args keywords back into the keymap,
        # setting value to "" for verbose, timer and the like
        for key in self.config_keys:
            try:
                self.keymap[key] = keyword(key, getattr(args, key))
            except:  # noqa: E722
                if key not in self.keymap:
                    self.keymap[key] = keyword(key, "")

        return args

    def run(self):
        raise Exception("runner.py is an abstract class")

    def end(self):
        if self.args.timer:
            print("execution time: %.4g" % (time.time() - self.start))

    def add_authorized_keys(self):
        """Add authorized keys in config keys to the parser"""
        for key in self.config_keys:
            self.parser.add_argument("--%s" % key, type=str, metavar="", help=self.config_keys[key])


if __name__ == "__main__":
    raise TypeError("Runner is a base class, not an executable")  # pragma: no cover
