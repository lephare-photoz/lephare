# ruff: noqa: N999

from lephare import GalMag, QSOMag, StarMag, keyword, read_config

__all__ = [
    "MagSvc",
]


class MagSvc:
    """Magnitude service to return magnitudes of a single SED."""

    @classmethod
    def from_config(cls, objtype, config_file):
        """Load config and return magnitudes

        Parameters
        ----------
        objtype : `str`
            The type of object galaxy, star, or qso based on case insensitive first
            letter of string.
        config_file : `str`
            Path to config file.
        """
        keywords = [
            "COSMOLOGY",
            "FILTER_FILE",
            "MAGTYPE",
            "EXTINC_LAW",
            "EB_V",
            "MOD_EXTINC",
            "ZGRID_TYPE",
            "Z_STEP",
            "LIB_ASCII",
            "VERBOSE",
            "ADD_DUSTEM",
        ]
        if objtype[0].upper() == "G":
            keywords += ["GAL_LIB_IN", "GAL_LIB_OUT", "EM_LINES", "EM_DISPERSION", "ADD_DUSTEM"]
            instance = GalMag
        elif objtype[0].upper() == "Q":
            keywords += ["QSO_LIB_IN", "QSO_LIB_OUT"]
            instance = QSOMag
        elif objtype[0].upper() == "S":
            keywords += ["STAR_LIB_IN", "STAR_LIB_OUT"]
            instance = StarMag
        else:
            raise ValueError(f"Unknown objtype: {objtype}")

        keymap = read_config(config_file)
        keymap["t"] = keyword("t", objtype)
        keymap["c"] = keyword("c", config_file)
        return instance(keymap)
