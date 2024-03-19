import os
import tempfile
import unittest

import numpy as np
from lephare import QSOSED, SED, GalSED, StarSED, flt

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATADIR = os.path.join(TESTDIR, "../data")


class sed_tests(unittest.TestCase):
    def test_sed_constructors(self):
        sed = SED("toto", 10, "GAL")
        self.assertEqual(sed.name, "toto")
        self.assertEqual(sed.nummod, 10)
        self.assertEqual(sed.is_gal(), True)
        sed2 = StarSED(sed)
        self.assertEqual(sed2.is_star(), True)
        self.assertEqual(sed2.name, "toto")
        self.assertEqual(sed2.nummod, 10)
        sed2 = QSOSED(sed)
        self.assertEqual(sed2.is_qso(), True)
        self.assertEqual(sed2.name, "toto")
        self.assertEqual(sed2.nummod, 10)
        sed2 = GalSED(sed)
        self.assertEqual(sed2.is_gal(), True)
        self.assertEqual(sed2.name, "toto")
        self.assertEqual(sed2.nummod, 10)

    def test_sed_fromfile(self):
        sed = SED("toto", 10, "GAL")
        sed.read(os.path.join(TESTDATADIR, "sed/o5v.sed.ext"))
        self.assertEqual(sed.size(), 4860)
        d = sed.data()
        x = d[0]
        y = d[1]
        sed.rescale(2.0)
        y2 = sed.data()[1]
        self.assertTrue(np.allclose(y2, 2 * y))

    def test_sed_readwrite(self):
        star = StarSED("star_test")
        star.read(os.path.join(TESTDATADIR, "sed/WDgd71.sed.ext"))

        with tempfile.TemporaryDirectory() as dir_name:
            file_path = os.path.join(dir_name, "test_star.fits")
            star.writeSED(file_path + ".bin", file_path + ".phys", file_path + ".doc")

            newsed = StarSED("newsed")
            newsed.readSEDBin(file_path + ".bin")
            np.testing.assert_array_equal(newsed.data()[0], star.data()[0])
            np.testing.assert_array_equal(newsed.data()[1], star.data()[1])

    def test_integration(self):
        sed = SED("toto", 10, "GAL")
        # Tophat filter
        lmin = 200
        lmax = 300
        hat = flt(lmin, lmax, 100)
        # flat SED ov val=1
        x = np.linspace(100, 500, 1000)
        sed.set_vector(x, np.ones_like(x))
        result = sed.integrateSED(hat)
        c = 2.99792458e18
        true_res = [
            lmax - lmin,
            c * (1 / lmin - 1 / lmax),
            0.5 * (lmax**2 - lmin**2),
            lmax - lmin,
            0.5 * (lmax**2 - lmin**2),
            1 / lmin - 1 / lmax,
        ]
        print(result)
        print(true_res)
        relative_error = (np.array(result) - np.array(true_res)) / np.array(true_res)
        self.assertTrue(np.allclose(np.array(result), np.array(true_res), 1.1e-2))


if __name__ == "__main__":
    unittest.main()
