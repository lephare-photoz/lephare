import os
from contextlib import suppress

from ._lephare import flt, write_output_filter
from .runner import Runner

__all__ = [
    "Filter",
]

config_keys = {
    "verbose": "increase onscreen verbosity",
    "FILTER_REP": "path to repository where filter files are searched",
    "FILTER_LIST": "list of filter files, to be searched in FILTER_REP if relative",
    "TRANS_TYPE": "Transmission curve type: 0[def] for energy, 1 for photon nb",
    "FILTER_CALIB": "calibration system for the filter:\
    0[def]: fnu=cst 1: nu.fnu=cst 2: fnu=nu 3: fnu=Black Body @ T=10000K\
    4: for MIPS (leff with nu fnu=ctt and flux with BB @ 10000K",
    "FILTER_FILE": "output filter filename, will be saved in $LEPHAREWORK/filt/",
}


class Filter(Runner):
    """
    The specific arguments to the Filter class are

    verbose
                increase onscreen verbosity
    FILTER_REP
                path to repository where filter files are searched
    FILTER_LIST
                list of filter files, to be searched in FILTER_REP if relative
    TRANS_TYPE
                Transmission curve type: 0[def] for energy, 1 for photon nb
    FILTER_CALIB
                calibration system for the filter:
                  0[def]: fnu=cst
                  1: nu.fnu=cst
                  2: fnu=nu
                  3: fnu=Black Body @ T=10000K
                  4: for MIPS (leff with nu fnu=ctt and flux with BB @ 10000K
    FILTER_FILE
                output filter filename, will be saved in $LEPHAREWORK/filt/
    """

    def update_help(self):
        """Add specific help information"""
        doc = "Build the LePHARE internal representation of the set of filters to be used\n"
        with suppress(Exception):
            self.parser.usage = doc
        self.__doc__ = doc + "\n"  # + inspect.getdoc(Filter)

    def __init__(self, config_file=None, config_keymap=None, **kwargs):
        super().__init__(config_keys, config_file, config_keymap, **kwargs)

    def run(self, **kwargs):
        """Update keymap and verbosity based on call arguments.

        This is only when the code is called from python session.
        """

        super().run(**kwargs)

        keymap = self.keymap
        flt_rep = keymap["FILTER_REP"].split_string(os.path.join(os.environ["LEPHAREDIR"], "filt/"), 1)[0]
        flt_files = keymap["FILTER_LIST"].split_string("flt.pb", -99)
        ref_size = len(flt_files)
        # Transmission in energy or photons
        transtyp = (keymap["TRANS_TYPE"]).split_int("0", ref_size)
        # calibration depending on the instrument
        calibtyp = (keymap["FILTER_CALIB"]).split_int("0", ref_size)
        # Output file
        output_name = keymap["FILTER_FILE"].split_string("filters", 1)[0]
        filtfile = os.path.join(os.environ["LEPHAREWORK"], "filt", output_name + ".dat")
        filtdoc = os.path.join(os.environ["LEPHAREWORK"], "filt", output_name + ".doc")

        if self.verbose:
            print("#######################################")
            print("# Build the filter file with the following options: ")
            print(f"# Config file: {self.config}")
            print(f"# FILTER_REP: {flt_rep}")
            print(f"# FILTER_LIST: {flt_files}")
            print(f"# TRANS_TYPE: {transtyp}")
            print(f"# FILTER_CALIB: {calibtyp}")
            print(f"# FILTER_FILE: {filtfile}")
            print(f"# FILTER_FILE.doc: {filtdoc}")
            print("#######################################")

        vec_flt = []
        for k, (f, t, c) in enumerate(zip(flt_files, transtyp, calibtyp)):
            flt_file = f if os.path.isabs(f) else os.path.join(flt_rep, f)
            one_filt = flt(k, flt_file, t, c)
            vec_flt.append(one_filt)

        write_output_filter(filtfile, filtdoc, vec_flt)
        return vec_flt


def main():  # pragma no cover
    runner = Filter()
    _ = runner.run()
    runner.end()


if __name__ == "__main__":  # pragma no cover
    main()
