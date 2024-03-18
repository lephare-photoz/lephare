import os
import unittest
from lephare._lephare import get_lephare_env, test_first_char, blackbody

TESTDIR = os.path.abspath(os.path.dirname(__file__))
LEPHAREDIR = os.path.join(TESTDIR, "..")


class globals(unittest.TestCase):
    def test_first_char(self):
        self.assertEqual(test_first_char("test"), True)
        self.assertEqual(test_first_char("   test"), True)
        self.assertEqual(test_first_char(""), False)
        self.assertEqual(test_first_char(" "), False)
        self.assertEqual(test_first_char("!test"), False)
        self.assertEqual(test_first_char("#test"), False)
        self.assertEqual(test_first_char("\t#"), False)

    def test_blackbody(self):
        self.assertAlmostEqual(blackbody(10000, 500), 1.018807e-26)


if __name__ == "__main__":
    unittest.main()
