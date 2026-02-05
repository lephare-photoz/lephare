# runner.py
import time
from pathlib import Path

from lephare import get_lephare_env, keyword

__all__ = ["Runner"]

DEBUG = False


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

    def __init__(self, config_keys=None, config_file="", config_keymap=None, **kwargs):
        get_lephare_env()

        if config_keys is None:
            raise RuntimeError("Runner is a base class and cannot be initialized directly")

        self.config_keys = config_keys
        self.verbose = False
        self.typ = None
        self.timer = False

        # 1) config file
        self.config_file = config_file
        self.keymap = self._load_config_file()

        # 2) explicit keymap
        self._load_keymap(config_keymap)

        # 3) kwargs
        self._load_kwargs(**kwargs)

        self.update_help()

    # ---------------- config loading ----------------

    def _load_config_file(self):
        """
        Load defaults from ASCII file with 3 columns: NAME VALUE COMMENT
        Lines starting with # are ignored.
        Only keep authorized keys.
        """
        defaults = {}
        if self.config_file:
            path = Path(self.config_file)

            if not path.exists():
                raise RuntimeError(f"File {self.config_file} not found")

            with path.open() as f:
                with open(self.config_file) as f:
                    for line in f:
                        line = line.strip()
                        if not line or line.startswith("#"):
                            continue
                        parts = line.split(None, 2)
                        if len(parts) < 2:
                            continue
                        name, value = parts[0], parts[1]
                        if name in self.config_keys:
                            defaults[name] = keyword(name, value)

        for k in self.config_keys:
            if k not in defaults:
                defaults[k] = keyword(k, "")

        return defaults

    def _load_keymap(self, config_keymap):
        if not config_keymap:
            return
        self.validate_config_dict(config_keymap)
        self.keymap |= config_keymap

        if "VERBOSE" in self.keymap:
            self.verbose = self.keymap["VERBOSE"].value.upper() == "YES"

    def _load_kwargs(self, **kwargs):
        self.validate_config_dict(kwargs, no_raise=False)
        for k, v in kwargs.items():
            uk = k.upper()
            if uk == "VERBOSE":
                self.verbose = bool(v)
                v = "YES" if v else "NO"
            self.keymap[uk] = keyword(uk, str(v))

    # ---------------- runtime ----------------

    def run(self, **kwargs):
        if DEBUG:  # pragma no cover
            print("#######################################")
            print(f"# Running {self.name} with the following arguments")
            print(f"# Config file : {self.config_file}")
        if self.timer:
            self.start = time.time()
        if DEBUG:  # pragma no cover
            print(f"# timer {self.timer}")

        if "VERBOSE" in kwargs:
            val = kwargs.pop("VERBOSE")
            self.verbose = bool(val.upper() == "YES")
        if DEBUG:  # pragma no cover
            print(f"# verbose {self.verbose}")

        if "typ" in self.config_keys:
            if "typ" in kwargs:
                self.typ = kwargs.pop("typ", self.typ)
                self.typ = self.typ.upper()
            if DEBUG:  # pragma no cover
                print(f"# typ {self.typ}")

        for k, v in kwargs.items():
            if k.upper() in self.config_keys:
                self.keymap[k.upper()] = keyword(k.upper(), str(v))
        # for debugging only
        if DEBUG:  # pragma no cover
            for k in self.config_keys:
                if k.name not in ["VERBOSE", "typ"]:
                    key = self.keymap[k.upper()]
                    print(f"# {key.name} : {key.value}")

    def end(self):
        if self.timer:
            print(f"execution time: {time.time() - self.start:.4g}")

    # ---------------- utils ----------------

    def validate_config_dict(self, input_dict, no_raise=True):
        for key in input_dict:
            if key == "TYP" and "typ" in self.config_keys:
                continue
            if key.upper() not in self.config_keys and not no_raise:
                raise RuntimeError(f"{key} is not a recognized argument")

    def update_help(self):
        self.__doc__ = "Runner base class"
