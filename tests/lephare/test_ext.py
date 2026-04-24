import lephare as lp
import numpy as np
import pytest
import scipy.integrate as sciint


def test_cardelli_ext():
    tophat = lp.flt(100.0, 1000.0, 500)
    extinction = lp.cardelli_ext(tophat)

    # almost complete extinction below 200
    def cardelli_factor(lamb):
        return pow(10.0, -0.4 * lp.cardelli_law(lamb))

    r, err = sciint.quad(cardelli_factor, 100, 1000)
    r /= 900  # filter integral
    # the bad resolution is entirely due to the
    # addition of the two bracketing 0 values in the Heaviside definition
    # the two current values are 13617.62 and 13699.37
    assert -2.5 * np.log10(r) == pytest.approx(extinction, 1.0e-4)


def test_basics():
    e = lp.ext("", 0)
    # set_vector with unmatched size
    with pytest.raises(RuntimeError):
        e.set_vector([0, 1, 2], [0, 2])
    # basic set_vector unit testing
    e.set_vector([0, 1], [0, 1])
    assert len(e.lamb_ext) == 2
    assert e.lamb_ext[0].lamb == 0
    assert e.lamb_ext[0].val == 0
    assert e.lamb_ext[-1].lamb == 1
    assert e.lamb_ext[-1].val == 1
    assert e.lmin == 0
    assert e.lmax == 1
    # add_element check
    e.add_element(2, 3)
    assert e.lmax == 2
    assert e.lamb_ext[-1].val == 3
