import os
from contextlib import suppress
import numpy as np

from ._lephare import GalMag, compute_filter_extinction, ext, flt
from .runner import Runner

__all__ = [
    "FiltExt",
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
                path to filter file on which to compute extinction,
                output of the `filter` execution, to be found in $LEPHAREWORK/filt
    EXT_CURVE
                extinction law to use, to be searched in $LEPHAREDIR/ext if relative
    GAL_CURVE
                extinction curve in the galaxy, to be searched in $LEPHAREDIR/ext if relative
    OUTPUT
                output file name
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

        filters = keymap["FILTER_FILE"].split_string("unknown", 1)[0]
        if not os.path.isabs(filters):
            filters = os.path.join(os.environ["LEPHAREWORK"], "filt", filters)
        all_filters = GalMag.read_flt(filters)
        
        atmec = keymap["EXT_CURVE"].split_string("NONE", 1)[0]
        if atmec is "NONE":
            aint = 99.*np.ones(len(all_filters)).tolist()
        else:
            if not os.path.isabs(atmec):
                atmec = os.path.join(os.environ["LEPHAREDIR"], "ext", atmec)
            atmospheric_ext = ext(atmec, 0)
            atmospheric_ext.read(atmec)
            aint = [compute_filter_extinction(f, atmospheric_ext) for f in all_filters]

        
        output = keymap["OUTPUT"].split_string("filter_extinc.dat", 1)[0]

        galec = keymap["GAL_CURVE"].split_string("CARDELLI", 1)[0]

        if galec == "CARDELLI":
            # If cardelli (hardcoded) with Rv=3.1 by  default
            # output A(lbd)/Av->A(lbd)/(Rv*E(B-V))->A(lbd)/E(B-V)=Rv*A(lbd)/Av
            albdav = np.array([cardelli_ext(f) for f in all_filters])
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
            if self.verbose:
                print(f"assuming Rv={rv} for this Extinction law {galec}")
            # galactic curves given in k(lbda) (=A(lbda)/E(B-V))
            #  -> A(lbda)/Av = A(lbda)/E(B-V) / Rv)
            albd = np.array([compute_filter_extinction(f, galactic_ext) for f in all_filters])
            albdav = albd / rv

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
                
        return


# compute galactic extinction in the filter based on Cardelli et al., 1989, ApJ 345
def cardelli_ext(filt: flt):
    # Define the limits of this filter
    lmin = filt.lmin()
    lmax = filt.lmax()
    one_ext = ext("CARDELLI", 2)

    # computes the galactic extinction
    dlbd = (lmax - lmin) / 400.0
    for i in range(402):
        lextg = lmin + (i - 1) * dlbd
        extg = cardelli_law(lextg)
        one_ext.add_element(lextg, extg, 2)

    return compute_filter_extinction(filt, one_ext)


#  compute albd/av at a given lambda (A) for the Cardelli law
def cardelli_law(lb):
    rv = 3.1
    x = 10000.0 / lb
    y = x - 1.82

    if x <= 1.1:
        f1 = 0.574 * pow(x, 1.61)
        f2 = -0.527 * pow(x, 1.61)
    elif x > 1.1 and x < 3.3:
        f1 = 1 + 0.17699 * y - 0.50447 * y * y - 0.02427 * y * y * y + 0.72085 * y * y * y * y
        f1 = f1 + 0.01979 * pow(y, 5) - 0.77530 * pow(y, 6) + 0.32999 * pow(y, 7)
        f2 = 1.41338 * y + 2.28305 * y * y + 1.07233 * y * y * y
        f2 = f2 - 5.38434 * pow(y, 4) - 0.62251 * pow(y, 5) + 5.30260 * pow(y, 6) - 2.09002 * pow(y, 7)
    elif x >= 3.3 and x < 5.9:
        f1 = 1.752 - 0.316 * x - 0.104 / ((x - 0.467) * (x - 0.467) + 0.341)
        f2 = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) * (x - 4.62) + 0.262)
    elif x >= 5.9 and x < 8:
        fa = -0.04473 * (x - 5.9) * (x - 5.9) - 0.009779 * (x - 5.9) * (x - 5.9) * (x - 5.9)
        fb = 0.2130 * (x - 5.9) * (x - 5.9) + 0.1207 * (x - 5.9) * (x - 5.9) * (x - 5.9)
        f1 = 1.752 - 0.316 * x - 0.104 / ((x - 0.467) * (x - 0.467) + 0.341) + fa
        f2 = -3.090 + 1.825 * x + 1.206 / ((x - 4.62) * (x - 4.62) + 0.262) + fb
    else:
        f1 = -1.073 - 0.628 * (x - 8) + 0.137 * (x - 8) * (x - 8) - 0.070 * (x - 8) * (x - 8) * (x - 8)
        f2 = 13.670 + 4.257 * (x - 8) - 0.420 * (x - 8) * (x - 8) + 0.374 * (x - 8) * (x - 8) * (x - 8)

    return f1 + f2 / rv


def main():  # pragma no cover
    runner = FiltExt()
    runner.run()
    runner.end()


if __name__ == "__main__":  # pragma no cover
    main()
