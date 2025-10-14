import os
import tempfile

import lephare as lp
import numpy as np
import pytest
from lephare import StarMag, GalMag, QSOMag

TESTDIR = os.path.abspath(os.path.dirname(__file__))
TESTDATADIR = os.path.join(TESTDIR, "../data")


def test_readB12():
    mag = GalMag()
    mag.set_zgrid(0.5,0,1)
    mag.read_B12()
    
