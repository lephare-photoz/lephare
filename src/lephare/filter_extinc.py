import os
from contextlib import suppress

import numpy as np

import lephare as lp

from ._lephare import GalMag, compute_filter_extinction, ext
from .runner import Runner

__all__ = [
    "FiltExt",
    "calculate_extinction_values",
]

config_keys = {
    "verbose": "increase onscreen verbosity",
    "FILTER_FILE": "path to filter file on which to compute extinction,\
    output of the `filter` execution, to be found in $LEPHAREWORK/filt",
    "EXT_CURVE": "extinction law to use, to be searched in $LEPHAREDIR/ext if relative",
    "GAL_CURVE": "extinction curve in the galaxy, to be searched in $LEPHAREDIR/ext if relative",
    "OUTPUT": "output file name",
}


class FiltExt(Runner):
    """
    The specific arguments to the Filter class are

    verbose
        increase onscreen verbosity
    FILTER_FILE
        Path to filter file on which to compute extinction, output of the
        `filter` execution, to be found in $LEPHAREWORK/filt
    EXT_CURVE
        Extinction law to use, to be searched in $LEPHAREDIR/ext if relative
    GAL_CURVE
        Extinction curve in the galaxy, to be searched in $LEPHAREDIR/ext if relative
    OUTPUT
        Output file name
    """

    def update_help(self):
        """Add specific help information"""
        doc = "Compute atmospheric extinction for a given filter\n"
        with suppress(Exception):
            self.parser.usage = doc
        self.__doc__ = doc + "\n"

    def __init__(self, config_file=None, config_keymap=None, **kwargs):
        super().__init__(config_keys, config_file, config_keymap, **kwargs)

    def run(self, **kwargs):
        """Update keymap and verbosity based on call arguments.
        This is only when the code is called from python session.
        """
        super().run(**kwargs)
        keymap = self.keymap
        # Get the parameters
        filters = keymap["FILTER_FILE"].split_string("unknown", 1)[0]
        if not os.path.isabs(filters):
            filters = os.path.join(os.environ["LEPHAREWORK"], "filt", filters)
        atmec = keymap["EXT_CURVE"].split_string("NONE", 1)[0]
        output = keymap["OUTPUT"].split_string("filter_extinc.dat", 1)[0]
        galec = keymap["GAL_CURVE"].split_string("CARDELLI", 1)[0]

        all_filters, aint, albdav, albd = calculate_extinction_values(
            filters, atmec, galec, verbose=self.verbose
        )

        if self.verbose:
            print("#######################################")
            print("# Computing ATMOSPHERIC AND GALACTIC EXTINCTION ")
            print("# with the following options:")
            print(f"# Config file: {self.config}")
            print(f"# FILTER_FILE: {filters}")
            print(f"# EXT_CURVE (Atmospheric extinction curve): {atmec}")
            print(f"# GAL_CURVE (Galactic extinction curve: {galec}")
            print(f"# Output file: {output}")
            print("#######################################")
            print(" Filters Ext(mag/airmass) Albda/Av Albda/E(B-V) ")
            for k, f in enumerate(all_filters):
                name = os.path.basename(f.name)
                print(f"{name:20} {aint[k]:<20} {albdav[k]:<20} {albd[k]:<20}")

        with open(output, "w") as out:
            out.write("#######################################\n")
            out.write("# Computing ATMOSPHERIC AND GALACTIC EXTINCTION \n")
            out.write("# with the following options:\n")
            out.write(f"# Config file: {self.config}\n")
            out.write(f"# FILTER_FILE: {filters}\n")
            out.write(f"# EXT_CURVE (Atmospheric extinction curve): {atmec}\n")
            out.write(f"# GAL_CURVE (Galactic extinction curve: {galec}\n")
            out.write(f"# Output file: {output}\n")
            out.write("#######################################\n")
            out.write(" Filters Ext(mag/airmass) Albda/Av Albda/E(B-V) \n")
            for k, f in enumerate(all_filters):
                name = os.path.basename(f.name)
                out.write(f"{name:20} {aint[k]:<20} {albdav[k]:<20} {albd[k]:<20}\n")

        return all_filters, aint, albdav, albd


def calculate_extinction_values(filters, atmec, galec, verbose=False):
    """Calculate the extinction values for a set of filters

    Parameters
    ==========

    filters :  str
        The file containing the lephare compiled filter list
    atmec : str
        Extinction law to use, to be searched in $LEPHAREDIR/ext if relative
    galec : str
        Extinction curve in the galaxy, to be searched in $LEPHAREDIR/ext if relative
    verbose : bool
        Increase onscreen verbosity

    Returns
    =======

    all_filters : list of lephare.GalMag
        The list of filter objects
    aint : np.array
        Atmospheric extinction in each filter (mag/airmass)
    albdav : np.array
        Galactic extinction in each filter A(lbd)/Av
    albd : np.array
        Galactic extinction in each filter A(lbd)/E(B-V)
    """
    all_filters = GalMag.read_flt(filters)
    if atmec == "NONE":
        aint = np.full(len(all_filters), 99.0).tolist()
    else:
        if not os.path.isabs(atmec):
            atmec = os.path.join(os.environ["LEPHAREDIR"], "ext", atmec)
        atmospheric_ext = ext(atmec, 0)
        atmospheric_ext.read(atmec)
        aint = [compute_filter_extinction(f, atmospheric_ext) for f in all_filters]

    if galec == "CARDELLI":
        # If cardelli (hardcoded) with Rv=3.1 by  default
        # output A(lbd)/Av->A(lbd)/(Rv*E(B-V))->A(lbd)/E(B-V)=Rv*A(lbd)/Av
        albdav = np.array([lp.cardelli_ext(f) for f in all_filters])
        albd = albdav * 3.1
    else:
        if not os.path.isabs(galec):
            galec = os.path.join(os.environ["LEPHAREDIR"], "ext", galec)
        galactic_ext = ext(galec, 1)
        galactic_ext.read(galec)

        #  Rv=3.1 except for Calzetti law (4.05) and SMC Prevot (2.72)
        rv = 3.1
        if "SMC_prevot" in galec:
            rv = 2.72
            if "calzetti" in galec:
                rv = 4.05
        if verbose:
            print(f"assuming Rv={rv} for this Extinction law {galec}")
        # galactic curves given in k(lbda) (=A(lbda)/E(B-V))
        #  -> A(lbda)/Av = A(lbda)/E(B-V) / Rv)
        albd = np.array([compute_filter_extinction(f, galactic_ext) for f in all_filters])
        albdav = albd / rv
    return all_filters, aint, albdav, albd


def main():  # pragma no cover
    runner = FiltExt()
    runner.run()
    runner.end()


if __name__ == "__main__":  # pragma no cover
    main()
