from contextlib import suppress

from ._lephare import GalMag, QSOMag, StarMag, keyword
from .runner import Runner

__all__ = [
    "MagGal",
]

config_keys = {
    "typ": "define what kind of objects these SED belong to : GAL, QSO, or STAR",
    "verbose": "increase onscreen verbosity",
    "COSMOLOGY": "fiducial cosmology used for absolute magnitudes evaluations",
    "FILTER_FILE": "filter file provided by filter script or the Filter class",
    "MAGTYPE": "AB or VEGA system",
    "EXTINC_LAW": "list of extinction files used for extinction correction\
    (if relative, in $LEPHAREDIR/ext/)",
    "MOD_EXTINC": "list of SED template ranges for which each extinction law is going to be applied",
    "EB_V": "E(B-V) grid of values to be applied to the SEDs (comma separated)",
    "Z_STEP": "dz, zmin, zmax defining the redshift grid",
    "GAL_LIB_IN": "input SED library file, as output by sedtolib script or SEDtolib class\
    (in $LEPHAREWORK/lib_bin/)",
    "GAL_LIB_OUT": "output file of the synthetic magnitudes (in $LEPHAREWORK/lib_mag/)",
    "QSO_LIB_IN": "same for QSO/AGN templates",
    "QSO_LIB_OUT": "same for QSO/AGN templates",
    "STAR_LIB_IN": "same for STAR templates",
    "STAR_LIB_OUT": "same for STAR templates",
    "LIB_ASCII": "if set to YES, also provide the output in ascii",
    "EM_LINES": "[NO/EMP_UV/EMP_SFR/PHYS] choice of prescription for emission line computation",
    "EM_DISPERSION": "rescaling values for the emission lines",
    "ADD_DUSTEM": "add the dust emission in templates when missing",
}


class MagGal(Runner):
    """The specific arguments to the MagGal class are

    typ:
           define what kind of objects these SED belong to : GAL, QSO, or STAR
    verbose:
           increase onscreen verbosity
    COSMOLOGY:
           fiducial h0, Omega_m0, and LambdaO used to define a flat LCDM cosmology
           used for absolute magnitudes evaluations
    FILTER_FILE:
           filter file provided by filter script or the Filter class
    MAGTYPE:
           AB or VEGA system
    EXTINC_LAW:
           list of extinction files used for extinction correction
           (if relative, in $LEPHAREDIR/ext/)
    MOD_EXTINC:
           list of SED template ranges for which each extinction law is going to be applied
    EB_V:
           E(B-V) grid of values to be applied to the SEDs (comma separated)
    Z_STEP:
           dz, zmin, zmax defining the redshift grid
    GAL_LIB_IN:
           input SED library file, as output by sedtolib script or SEDtolib class
           (in $LEPHAREWORK/lib_bin/)
    GAL_LIB_OUT:
           output file of the synthetic magnitudes (in $LEPHAREWORK/lib_mag/)
    QSO_LIB_IN, QSO_LIB_OUT, STAR_LIB_IN, STAR_LIB_OUT
           same for STAR templates
    LIB_ASCII:
           if set to YES, also provide the output in ascii
    EM_LINES:
           [NO/EMP_UV/EMP_SFR/PHYS] choice of prescription for emission line computation
    EM_DISPERSION:
           possible rescaling values for the emission lines
    ADD_DUSTEM:
           add the dust emission in templates when missing
    VERBOSE:
           add verbosity
    """

    def update_help(self):
        """Add the specific MagGal help"""
        doc = "Build the LePHARE internal representation of the set of SED templates to be used"
        with suppress(Exception):
            self.parser.usage = "Build the LePHARE synthetic magnitudes"
        self.__doc__ = doc + "\n"  # + inspect.getdoc(MagGal)

    def __init__(self, config_file=None, config_keymap=None, **kwargs):
        super().__init__(config_keys, config_file, config_keymap, **kwargs)

    def run(self, **kwargs):
        """Compute the model magnitudes across the redshift grid.

        Returns
        -------
        None
        """
        super().run(**kwargs)

        # needed to comply with the C++ interface
        # Define the type (Galaxy, QSO, Stars)
        self.keymap["t"] = keyword("t", self.typ)
        # Parameter file
        self.keymap["c"] = keyword("c", self.config)

        if self.typ[0] == "G":
            mag = GalMag(self.keymap)
        elif self.typ[0] == "Q":
            mag = QSOMag(self.keymap)
        elif self.typ[0] == "S":
            mag = StarMag(self.keymap)
        else:
            raise KeyError("-t arg must start with G/g Q/q or S/s for Galaxy QSO and Star respectively.")
        mag.open_files()
        mag.print_info()
        # Read dust extinction laws
        mag.read_ext()
        # Define the redshift grid
        mag.def_zgrid()
        # Read B12 templates to add dust emission to BC03
        add_dust = self.keymap["ADD_DUSTEM"].split_bool("NO", 1)
        if add_dust:
            mag.read_B12()
        # Read sed, apply extinction and IGM opacity
        mag.read_SED()
        mag.write_doc()
        # we need to call the close method here because run can
        # be called within a python session that stays alive afterwards
        mag.close_files()

        self.Mag = mag
        return


def main():  # pragma no cover
    runner = MagGal()
    runner.run()
    runner.end()


if __name__ == "__main__":  # pragma no cover
    main()
