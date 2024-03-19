import os
import unittest

from lephare import keyword, read_command, read_config

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATADIR = os.path.join(TESTDIR, "../data")
LEPHAREDIR = os.path.join(TESTDIR, "..")


class filter_tests(unittest.TestCase):
    def test_basic(self):
        # test that default is taken when value is empty
        k = keyword("k1", "")
        self.assertEqual(k.split_string("DEFAULT", 1)[0], "DEFAULT")
        # check basic assignment
        k = keyword("k2", "1")
        self.assertEqual(k.name, "k2")
        self.assertEqual(k.value, "1")
        # test mismatch in size
        self.assertEqual(k.split_string("0", 2), ["1", "1"])
        k = keyword("k3", "1,1")
        self.assertEqual(k.split_string("0", 3), ["0", "0", "0"])

    def test_array(self):
        k = keyword("k", "1,1")
        self.assertEqual(k.split_string("0", 2), ["1", "1"])
        self.assertEqual(k.split_int("0", 2), [1, 1])
        self.assertEqual(k.split_long("0", 2), [1, 1])
        self.assertEqual(k.split_bool("0", 2), [True, True])
        self.assertEqual(k.split_double("0", 2), [1.0, 1.0])

    def test_bool(self):
        k = keyword("k", "1")
        self.assertEqual(k.split_bool("0", 1)[0], True)
        k = keyword("k", "YES")
        self.assertEqual(k.split_bool("0", 1)[0], True)
        k = keyword("k", "yes")
        self.assertEqual(k.split_bool("0", 1)[0], True)
        k = keyword("k", "True")
        self.assertEqual(k.split_bool("0", 1)[0], True)
        k = keyword("k", "true")
        self.assertEqual(k.split_bool("0", 1)[0], True)
        k = keyword("k", "ANYTHING ELSE")
        self.assertEqual(k.split_bool("0", 1)[0], False)

    def test_expand_path(self):
        lepharedir = os.environ["LEPHAREDIR"]
        home = os.environ["HOME"]
        k = keyword("k", "$LEPHAREDIR/unittest")
        self.assertEqual(k.value, os.path.join(lepharedir, "unittest"))
        k = keyword("k", "$HOME/unittest")
        self.assertEqual(k.value, os.path.join(home, "unittest"))

    def test_read_config(self):
        config_file = os.path.join(TESTDATADIR, "examples/COSMOS.para")
        km = read_config(config_file)
        for k in km.keys():
            self.assertNotEqual(k[0], "#")
            self.assertNotEqual(k[0], " ")
            self.assertNotEqual(k[0], "\t")

    def test_read_command(self):
        km = read_command(["exec", "-c", "config", "--arg1", "val1", "-arg2", "val2"])
        for k in km.keys():
            if k == "c":
                self.assertEqual(km[k].value, "config")
            if k == "arg1":
                self.assertEqual(km[k].value, "val1")
            if k == "arg2":
                self.assertEqual(km[k].value, "val2")

    def test_misformed_kwd(self):
        # normal
        k = keyword("k", "1,2,3")
        self.assertEqual(k.split_string("0", -1), ["1", "2", "3"])
        # space
        k = keyword("k", "1,2, 3")
        self.assertEqual(k.split_string("0", -1), ["1", "2", "3"])
        # tab
        k = keyword("k", "1,2,\t3")
        self.assertEqual(k.split_string("0", -1), ["1", "2", "3"])
        with open("tmp", "w") as f:
            f.writelines(
                [
                    "FILTER_CALIB    0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,\t0",
                    "\nFILTER_FILE     filter_cosmos",
                ]
            )
        km = read_config("tmp")
        self.assertEqual(km["FILTER_CALIB"].value, "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0")
        self.assertEqual(km["FILTER_FILE"].value, "filter_cosmos")
        os.remove("tmp")


if __name__ == "__main__":
    unittest.main()
