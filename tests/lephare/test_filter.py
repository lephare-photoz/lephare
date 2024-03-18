import os
import unittest
import numpy as np
from lephare import LEPHAREDIR, flt
from lephare.filterSvc import FilterSvc

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATADIR = os.path.join(TESTDIR, "../data")


class filter_tests(unittest.TestCase):

    def test_flt_class(self):
        tophat = flt(100.0, 200.0, 50)
        self.assertEqual(tophat.name, "Heavy")
        self.assertEqual(tophat.width(), 100.0)
        #
        filename = os.path.join(TESTDATADIR, "filt/IB527.pb")
        f = flt(-1, filename, trans=1, calib=1)
        self.assertAlmostEqual(f.width(), 241.9479, 4)
        self.assertAlmostEqual(f.lambdaMean(), 5262.2831, 4)

    def test_FilterSvc(self):
        filename = os.path.join(TESTDATADIR, "filt/IB527.pb")
        f = FilterSvc.from_file(filename, trans=1, calib=0)
        self.assertAlmostEqual(f.width(), 241.9479, 4)
        self.assertAlmostEqual(f.lambdaMean(), 5262.2831, 4)

        g = FilterSvc.from_svo(-1, "Subaru/Suprime.IB527", "AB")
        if g is not None:
            self.assertTrue(np.allclose(f.data()[1], g.data()[1]))
        #
        fltVec = FilterSvc.from_config("../examples/COSMOS.para")
        self.assertEqual(len(fltVec), 30)


if __name__ == "__main__":
    unittest.main()
