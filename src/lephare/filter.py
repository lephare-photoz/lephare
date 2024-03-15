import os

from lephare._lephare import flt, keyword, write_output_filter
from lephare.runner import Runner

__all__ = [
    "Filter",
]

global filter_config_keys
filter_config_keys = ["FILTER_REP", "FILTER_LIST", "TRANS_TYPE", "FILTER_CALIB", "FILTER_FILE"]


class Filter(Runner):
    def __init__(self, config_file=None, config_keymap=None):
        super().__init__(filter_config_keys, config_file, config_keymap)

    def run(self, **kwargs):
        # update keymap and verbosity based on call arguments
        # this is only when the code is called from python session
        self.verbose = kwargs.pop("verbose", self.verbose)
        for k, v in kwargs.items():
            if k.upper() in self.keymap.keys():
                self.keymap[k.upper()] = keyword(k.upper(), str(v))

        keymap = self.keymap
        fltRep = keymap["FILTER_REP"].split_string(os.path.join(os.environ["LEPHAREDIR"], "filt/"), 1)[0]
        fltFiles = keymap["FILTER_LIST"].split_string("flt.pb", -99)
        ref_size = len(fltFiles)
        # Transmission in energy or photons
        transtyp = (keymap["TRANS_TYPE"]).split_int("0", ref_size)
        # calibration depending on the instrument
        calibtyp = (keymap["FILTER_CALIB"]).split_int("0", ref_size)
        # Output file
        outputName = keymap["FILTER_FILE"].split_string("filters", 1)[0]
        filtfile = os.path.join(os.environ["LEPHAREWORK"], "filt", outputName + ".dat")
        filtdoc = os.path.join(os.environ["LEPHAREWORK"], "filt", outputName + ".doc")

        if self.verbose:
            print("#######################################")
            print("# Build the filter file with the following options: ")
            print(f"# Config file: {self.config}")
            print(f"# FILTER_REP: {fltRep}")
            print(f"# FILTER_LIST: {fltFiles}")
            print(f"# TRANS_TYPE: {transtyp}")
            print(f"# FILTER_CALIB: {calibtyp}")
            print(f"# FILTER_FILE: {filtfile}")
            print(f"# FILTER_FILE.doc: {filtdoc}")
            print("#######################################")

        vecFlt = []
        for k, (f, t, c) in enumerate(zip(fltFiles, transtyp, calibtyp)):
            fltFile = os.path.join(fltRep, f)
            oneFilt = flt(k, f, t, c)
            oneFilt.read(fltFile)
            vecFlt.append(oneFilt)

        write_output_filter(filtfile, filtdoc, vecFlt)


if __name__ == "__main__":
    runner = Filter()
    runner.run()
    runner.end()
