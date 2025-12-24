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


def test_basics():
    e = lp.ext("", 0)
    # set_vector with inmatched size
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
