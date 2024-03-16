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
]


class Sedtolib(Runner):
    # Initalisation
    def __init__(self, config_file=None, config_keymap=None):
        # Class heritates from runner. So, __init__ in runner.py
        super().__init__(sedtolib_config_keys, config_file, config_keymap)

    def run(self, **kwargs):
        # update keymap and verbosity based on call arguments
        # this is only when the code is called from python session
        self.verbose = kwargs.pop("verbose", self.verbose)

        # loop over a dictionnary passed in argument (key and value)
        for k, v in kwargs.items():
            if k == "typ":
                self.typ = v.upper()
                continue
            if k == "config":
                self.config = config_file
                continue
            # if the key in argument is in the keymap, update the keymap with the argument
            if k.upper() in self.keymap.keys():
                self.keymap[k.upper()] = keyword(k.upper(), str(v))

        if self.typ[0] == "G":
            SEDLibrary = GalSEDLib(self.keymap, self.config, self.typ)
        elif self.typ[0] == "Q":
            SEDLibrary = QSOSEDLib(self.keymap, self.config, self.typ)
        elif self.typ[0] == "S":
            SEDLibrary = StarSEDLib(self.keymap, self.config, self.typ)
        else:
            raise KeyError("-t arg must start with G/g Q/q or S/s for Galaxy QSO and Star respectively.")

        SEDLibrary.print_info()
        SEDLibrary.read_model_list()
        SEDLibrary.write_SED_lib()
        SEDLibrary.print_time_tofile(int(time.time()))
        # we need to call the close method here because run can
        # be called within a python session that stays alive afterwards
        SEDLibrary.close_output_files()

        self.SEDLib = SEDLibrary
        return


if __name__ == "__main__":
    runner = Sedtolib()
    runner.run()
    runner.end()
