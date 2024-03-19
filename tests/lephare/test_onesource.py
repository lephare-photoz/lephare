import unittest

from lephare import onesource

# TESTDIR = os.path.abspath(os.path.dirname(__file__))
# LEPHAREDIR = os.path.join(TESTDIR,"..")


class onesource_tests(unittest.TestCase):
    def test_creation(self):
        src = onesource()
        self.assertEqual(src.spec, "1")
        self.assertEqual(src.zs, -99.9)
        self.assertEqual(src.cont, 0)
        self.assertEqual(src.closest_red, 0)
        self.assertEqual(src.zmin, [-99.9, -99.9, -99.9])
        self.assertEqual(src.chimin, [1.0e9, 1.0e9, 1.0e9])
        self.assertEqual(src.indmin, [-99, -99, -99])
        self.assertEqual(src.imasmin, [-99, -99, -99])

        src = onesource(31, [0, 0.5, 1])
        self.assertEqual(src.pos, 31)
        self.assertEqual(src.pdfmap[9].size(), 3)
        self.assertEqual(src.spec, "1")
        self.assertEqual(src.zs, -99.9)
        self.assertEqual(src.cont, 0)
        self.assertEqual(src.closest_red, 0)
        self.assertEqual(src.zmin, [-99.9, -99.9, -99.9])
        self.assertEqual(src.chimin, [1.0e9, 1.0e9, 1.0e9])
        self.assertEqual(src.indmin, [-99, -99, -99])
        self.assertEqual(src.imasmin, [-99, -99, -99])

    def test_setPriors(self):
        src = onesource()
        src.setPriors([0.0, 1000.0], [0, 1000])
        self.assertEqual(src.priorLib, [0.0, 0.0, 1000.0, 1000.0])


if __name__ == "__main__":
    unittest.main()
