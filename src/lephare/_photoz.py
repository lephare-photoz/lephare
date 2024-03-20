import os

import numpy as np
from astropy.table import Column, Table

from ._lephare import PhotoZ, readOutKeywords
from ._utils import continueClass

__all__ = [
    "PhotoZ",
]


@continueClass
class PhotoZ:  # noqa: F811
    pdfdict = {
        "MASS": 0,
        "SFR": 1,
        "SSFR": 2,
        "LDUST": 3,
        "LIR": 4,
        "AGE": 5,
        "COL1": 6,
        "COL2": 7,
        "MREF": 8,
        "MIN_ZG": 9,
        "MIN_ZQ": 10,
        "BAY_ZG": 11,
        "BAY_ZQ": 12,
    }

    def build_output_tables(self, srclist, para_out=None, filename=None):
        # BUILD THE TABLE OF THE OUTPUT PARAMETERS
        d = np.loadtxt(os.path.join(os.environ["LEPHAREDIR"], "alloutputkeys.txt"), dtype="str")
        allkeys = {}
        for count, label in enumerate(d[:, 0]):
            allkeys[label] = (d[count, 1], d[count, 2])
        outkeys = readOutKeywords(self.outpara) if para_out is None else readOutKeywords(para_out)
        outputs = {}
        t = Table()
        for key in outkeys:
            typ, attr = allkeys[key]
            is_array = False
            if "[" in attr:
                tmp = attr.split("[")
                attr = tmp[0]
                try:
                    index = int(tmp[1][0])
                except ValueError:
                    index = tmp[1][1:-2]
                is_array = True
            outputs[key] = []
            for src in srclist:
                if is_array is False:
                    outputs[key].append(getattr(src, attr))
                else:
                    outputs[key].append(getattr(src, attr)[index])
            t.add_column(Column(name=key, dtype=typ, data=outputs[key]))
        # NOW THE PDFS, IN THE SAME TABLE FOR NOW
        for typ in self.pdftype:
            arr = []
            key = self.pdfdict[typ]
            for src in srclist:
                pdf = src.pdfmap[key]
                arr.append(pdf.vPDF)
            t.add_column(Column(name=typ, dtype=float, data=arr))
            # use the last src object to get the x-axis values
            t.meta[typ] = " ".join(str(e) for e in pdf.xaxis)

        if filename is not None:
            self.save_table(t, filename)
        return t

    def save_table(self, table, filename, fmt="fits", overwrite=True):
        if fmt == "fits":
            table.write(filename, "fits", overwrite)
