import os
import unittest

from lephare import _lephare

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATADIR = os.path.join(TESTDIR, "../data")


class oneElLambda(unittest.TestCase):
    def test_constructor(self):
        el1 = _lephare.oneElLambda(10, 10, 0)
        # test that 1 + 1 = 2
        self.assertEqual(el1.lamb, 10)
        self.assertEqual(el1.val, 10)
        self.assertEqual(el1.ori, 0)
        el2 = _lephare.oneElLambda(el1)
        self.assertEqual(el1.lamb, el2.lamb)
        self.assertEqual(el1.val, el2.val)
        self.assertEqual(el1.ori, el2.ori)

    def test_interp(self):
        el1 = _lephare.oneElLambda(10, 10, 0)
        el2 = _lephare.oneElLambda(20, 20, 0)
        el3 = _lephare.oneElLambda(15, 0, 0)
        el3.interp(el1, el2)
        self.assertEqual(el3.val, 15)
        # next fails as interp expect ordering
        # which could be made to work
        # el3.interp(el2,el1)
        # self.assertEqual(el3.val, 15)
        pass


class ext(unittest.TestCase):
    def test_constructor(self):
        ext = _lephare.ext(name="test", numext=10)
        self.assertEqual(ext.name, "test")
        self.assertEqual(ext.numext, 10)
        self.assertEqual(len(ext.lamb_ext), 0)

    def test_read(self):
        filename = os.path.join(TESTDATADIR, "calzetti.dat")
        ext = _lephare.ext(name="test", numext=10)
        ext.read(filename)
        self.assertEqual(ext.lmin(), 400.0)
        self.assertEqual(ext.lmax(), 40000.0)


class ext(unittest.TestCase):
    def test_constructor(self):
        ext = _lephare.ext(name="test", numext=10)
        self.assertEqual(ext.name, "test")
        self.assertEqual(ext.numext, 10)
        self.assertEqual(len(ext.lamb_ext), 0)

    def test_read(self):
        filename = os.path.join(TESTDATADIR, "tau10.out")
        ext = _lephare.ext(name="test", numext=10)
        ext.read(filename)


if __name__ == "__main__":
    unittest.main()
