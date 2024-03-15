from lephare._lephare import GalMag, QSOMag, StarMag, keyword
from lephare.runner import Runner

__all__ = [
    "MagGal",
]

global mag_gal_config_keys
mag_gal_config_keys = [
    "COSMOLOGY",
    "FILTER_FILE",
    "MAGTYPE",
    "EXTINC_LAW",
    "EB_V",
    "MOD_EXTINC",
    "ZGRID_TYPE",
    "Z_STEP",
    "GAL_LIB_IN",
    "QSO_LIB_IN",
    "STAR_LIB_IN",
    "GAL_LIB_OUT",
    "QSO_LIB_OUT",
    "STAR_LIB_OUT",
    "LIB_ASCII",
    "EM_LINES",
    "EM_DISPERSION",
    "ADD_DUSTEM",
]


class MagGal(Runner):
    def __init__(self, config_file=None, config_keymap=None):
        super().__init__(mag_gal_config_keys, config_file, config_keymap)

    def run(self, **kwargs):
        # update keymap and verbosity based on call arguments
        # this is only when the code is called from python session
        self.verbose = kwargs.pop("verbose", self.verbose)
        for k, v in kwargs.items():
            if k == "typ":
                self.typ = v.upper()
                continue
            if k == "config":
                self.config = config_file
                continue
            if k.upper() in self.keymap.keys():
                self.keymap[k.upper()] = keyword(k.upper(), str(v))

        # needed to comply with the C++ interface
        # Define the type (Galaxy, QSO, Stars)
        self.keymap["t"] = keyword("t", self.typ)
        # Parameter file
        self.keymap["c"] = keyword("c", self.config)

        if self.typ[0] == "G":
            Mag = GalMag(self.keymap)
        elif self.typ[0] == "Q":
            Mag = QSOMag(self.keymap)
        elif self.typ[0] == "S":
            Mag = StarMag(self.keymap)
        else:
            raise KeyError("-t arg must start with G/g Q/q or S/s for Galaxy QSO and Star respectively.")
        Mag.open_files()
        Mag.print_info()
        # Read dust extinction laws
        Mag.read_ext()
        # Define the redshift grid
        Mag.def_zgrid()
        # Read B12 templates to add dust emission to BC03
        Mag.read_B12()
        # Read sed, apply extinction and IGM opacity
        Mag.read_SED()
        Mag.write_doc()
        # we need to call the close method here because run can
        # be called within a python session that stays alive afterwards
        Mag.close_files()

        self.Mag = Mag
        return


if __name__ == "__main__":
    runner = MagGal()
    runner.run()
    runner.end()
