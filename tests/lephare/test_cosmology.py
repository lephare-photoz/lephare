import os
import unittest

import numpy as np
from lephare import *

TESTDIR = os.path.abspath(os.path.dirname(__file__))
LEPHAREDIR = os.path.join(TESTDIR, "..")


class globals(unittest.TestCase):
    def test_cosmology(self):
        # test default arguments
        cosmo_test = cosmo()
        # test setting arguments
        cosmo_test = cosmo(h0=67, om0=0.27, l0=0.73)

    def test_zgrid(self):
        # linear grid
        zmin = 0.000
        zmax = 1
        dz = 0.1
        # test linear grid
        grid = zgrid(0, dz, zmin, zmax)
        dummy = np.arange(zmin, zmax, dz)
        if zmin > 0:
            dummy = np.insert(dummy, 0, 0)
        if dummy[-1] != zmax:
            pygrid = np.insert(dummy, len(dummy), zmax)
        self.assertIsNone(np.testing.assert_almost_equal(grid, pygrid))
        # test (1+z)dz grid
        grid = zgrid(1, dz, zmin, zmax)
        pygrid = np.array([])
        z = zmin
        while z < zmax:
            pygrid = np.append(pygrid, z)
            z = z + (1 + z) * dz
        pygrid = np.append(pygrid, zmax)
        if zmin > 0:
            pygrid = np.insert(pygrid, 0, 0)
        self.assertIsNone(np.testing.assert_almost_equal(grid, pygrid))

    def test_indexz(self):
        grid = np.arange(1, 5, 1)
        # test value outside of grid on the lower side
        self.assertEqual(indexz(0, grid), 0)
        # test value outside of grid on the upper side
        self.assertEqual(indexz(6, grid), grid.size - 1)
        # test exact valuesin the grid
        self.assertEqual(indexz(1.0, grid), 0)
        self.assertEqual(indexz(2.0, grid), 1)
        self.assertEqual(indexz(3.0, grid), 2)
        self.assertEqual(indexz(4.0, grid), 3)
        # test value greater and close to a grid value
        self.assertEqual(indexz(2.1, grid), 1)
        # test value half to two grid values :
        # convention is to return the value above.
        self.assertEqual(indexz(2.5, grid), 2)
        # test value smaller and close to a grid value
        self.assertEqual(indexz(2.9, grid), 2)


if __name__ == "__main__":
    unittest.main()
