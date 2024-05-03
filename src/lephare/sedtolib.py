import time

from lephare._lephare import GalSEDLib, QSOSEDLib, StarSEDLib, keyword
from lephare.runner import Runner

__all__ = [
    "Sedtolib",
]

global sedtolib_config_keys
# List of keywords associated to setolib
sedtolib_config_keys = [
    "GAL_SED",
    "GAL_FSCALE",
    "GAL_LIB",
    "SEL_AGE",
    "AGE_RANGE",
    "QSO_SED",
    "QSO_FSCALE",
    "QSO_LIB",
    "STAR_SED",
    "STAR_LIB",
    "STAR_FSCALE",
    "LIB_ASCII",
]


class Sedtolib(Runner):
    """Read a configurable set of SED, compute extinction corrections, and store the
    results into a binary library for later use.

    The run method is equivalent to the terminal sedtolib command.

    Parameters
    ----------
    config_file : `string` or `None`, optional
        Path to config file in LePHARE .para format
    config_keymap : `dict` or `None`, optional
        Dictionary of all config values as alternative to config file.
    """

    # Initalisation
    def __init__(self, config_file=None, config_keymap=None):
        # Class heritates from runner. So, __init__ in runner.py
        super().__init__(sedtolib_config_keys, config_file, config_keymap)

    def run(self, **kwargs):
        """Take keymap and set SED library as class variable.

        Must be run independently for stars, galaxies, and QSO.
        """
        # update keymap and verbosity based on call arguments
        # this is only when the code is called from python session
        self.verbose = kwargs.pop("verbose", self.verbose)

        # loop over a dictionnary passed in argument (key and value)
        for k, v in kwargs.items():
            if k == "typ":
                self.typ = v.upper()
                continue
            # if k == "config":
            #     self.config = config_file
            #     continue
            # if the key in argument is in the keymap, update the keymap with the argument
            if k.upper() in self.keymap:
                self.keymap[k.upper()] = keyword(k.upper(), str(v))

        if self.typ[0] == "G":
            sed_library = GalSEDLib(self.keymap, self.config, self.typ)
        elif self.typ[0] == "Q":
            sed_library = QSOSEDLib(self.keymap, self.config, self.typ)
        elif self.typ[0] == "S":
            sed_library = StarSEDLib(self.keymap, self.config, self.typ)
        else:
            raise KeyError("-t arg must start with G/g Q/q or S/s for Galaxy QSO and Star respectively.")

        sed_library.print_info()
        sed_library.read_model_list()
        sed_library.write_SED_lib()
        sed_library.print_time_tofile(int(time.time()))
        # we need to call the close method here because run can
        # be called within a python session that stays alive afterwards
        sed_library.close_output_files()

        self.SEDLib = sed_library
        return


def main():
    runner = Sedtolib()
    runner.run()
    runner.end()


if __name__ == "__main__":
    main()
