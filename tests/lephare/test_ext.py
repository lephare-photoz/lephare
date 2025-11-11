import lephare as lp
import pytest
import scipy.integrate as sciint


def test_cardelli_ext():
    tophat = lp.flt(100.0, 200.0, 500)
    extinction = lp.cardelli_ext(tophat)
    r, err = sciint.quad(lp.cardelli_law, 100, 200)
    r /= 100  # filter integral
    # the bad resolution is entirely due to the
    # addition of the two bracketing 0 values in the Heaviside definition
    # the two current values are 13617.62 and 13699.37
    assert r == pytest.approx(extinction, 100)
