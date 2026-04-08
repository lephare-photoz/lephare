import os

import lephare as lp
import numpy as np

# explicit import of a private function
from lephare._lephare import _get_opa_vector

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATADIR = os.path.join(TESTDIR, "../data")


def test_opa():
    opalist = _get_opa_vector()
    assert len(opalist) == 81
    sed = lp.GalSED("", 0)
    sed.read(os.path.join(TESTDATADIR, "sed/GAL/COSMOS_SED/Ell1_A_0.sed"))
    sed.red = 0
    x, y = sed.data()
    sed.applyOpa(opalist)
    x1, y1 = sed.data()
    assert np.array_equal(x, x1)
    print(y, y1)
    mask = np.where(y1 != y)[0]

    opax = [ll.lamb for ll in opalist[0].lamb_opa]
    opay = [ll.val for ll in opalist[0].lamb_opa]

    newopa = np.interp(x[mask], opax, opay)
    assert np.allclose(y1[mask] / y[mask], newopa)
