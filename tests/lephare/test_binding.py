import os
import unittest

os.environ["LEPHAREDIR"] = "."
os.environ["LEPHAREWORK"] = "."

from lephare import zphota


class binding_tests(unittest.TestCase):
    def test_zphota(self):
        _ = zphota.time.CLOCK_MONOTONIC
