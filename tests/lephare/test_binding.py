import os

os.environ["LEPHAREDIR"] = "."
os.environ["LEPHAREWORK"] = "."

from lephare import zphota


def test_zphota():
    _ = zphota.time.CLOCK_MONOTONIC
