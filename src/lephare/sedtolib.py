import time

from lephare._lephare import GalSEDLib, QSOSEDLib, StarSEDLib, keyword
from lephare.runner import Runner

__all__ = [
    "Sedtolib",
]

# List of keywords associated to setolib
config_keys = {
    "type": "define what kind of objects these SED belong to : GAL, QSO, or STAR",
    "GAL_SED": "file listing the galaxy SEDs to be used",
    "GAL_FSCALE": "arbitrary Flux scale for galaxy templates",
    "GAL_LIB": "name of the output binary SED file for the galaxies (relative to $LEPHAREWORK/lib_bin/)",
    "SEL_AGE": "file listing the different galaxy ages to consider",
    "AGE_RANGE": "minimal and maximal age in year to consider",
    "QSO_SED": "same for QSO/AGN templates",
    "QSO_FSCALE": "same for QSO/AGN templates",
    "QSO_LIB": "same for QSO/AGN templates",
    "STAR_SED": "same for STAR templates",
    "STAR_LIB": "same for STAR templates",
    "STAR_FSCALE": "same for STAR templates",
    "LIB_ASCII": "if set to YES, also provide the output in ascii",
}


class Sedtolib(Runner):
    """
    The specific arguments to the Sedtolib class are

    type:
           define what kind of objects these SED belong to : GAL, QSO, or STAR
    GAL_SED:
           file listing the galaxy SEDs to be used
    GAL_FSCALE":
           arbitrary Flux scale for galaxy templates
    GAL_LIB:
           name of the output binary SED file for the galaxies (relative to $ZPHOTWORK/lib_bin/)
    SEL_AGE:
           file listing the different galaxy ages to consider
    AGE_RANGE:
           minimal and maximal age in year to consider
    QSO_SED, QSO_FSCALE, QSO_LIB, STAR_SED, STAR_LIB, STAR_FSCALE :
           same for QSO/AGN and STAR SED types
    LIB_ASCII:
           if set to YES, also provide the output in ascii
    """

    def add_authorized_keys(self):
        """Add the specific Sedtolib arguments to the argument parser"""
        for key in config_keys:
            self.parser.add_argument("--%s" % key, type=str, metavar="", help=config_keys[key])
        self.parser.usage = "Build the LePHARE internal representation of the set of SED templates to be used"

    def __init__(self, config_file=None, config_keymap=None, **kwargs):
        super().__init__(config_keys, config_file, config_keymap, **kwargs)

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
