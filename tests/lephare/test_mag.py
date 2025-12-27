import os

from lephare import GalMag, QSOMag, StarMag

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATADIR = os.path.join(TESTDIR, "../data")


def test_basic():
    # empty constructors
    mag = StarMag()
    mag = QSOMag()


def test_readb12():
    mag = GalMag()
    mag.set_zgrid(0.5, 0, 1)
    mag.read_B12()
