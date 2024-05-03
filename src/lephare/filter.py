import os

from lephare._lephare import flt, keyword, write_output_filter
from lephare.runner import Runner

__all__ = [
    "Filter",
]

global filter_config_keys
filter_config_keys = ["FILTER_REP", "FILTER_LIST", "TRANS_TYPE", "FILTER_CALIB", "FILTER_FILE"]


class Filter(Runner):
    """Build the filter files based on config

    Parameters
    ----------
    config_file : `string`
    """

    def __init__(self, config_file=None, config_keymap=None):
        super().__init__(filter_config_keys, config_file, config_keymap)

    def run(self, **kwargs):
        """Update keymap and verbosity based on call arguments.

        This is only when the code is called from python session.
        """
        self.verbose = kwargs.pop("verbose", self.verbose)
        for k, v in kwargs.items():
            if k.upper() in self.keymap:
                self.keymap[k.upper()] = keyword(k.upper(), str(v))

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
            flt_file = os.path.join(flt_rep, f)
            one_filt = flt(k, f, t, c)
            one_filt.read(flt_file)
            vec_flt.append(one_filt)

        write_output_filter(filtfile, filtdoc, vec_flt)


def main():
    runner = Filter()
    runner.run()
    runner.end()


if __name__ == "__main__":
    main()
