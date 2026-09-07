import lephare


def test_version():
    """Check to see that we can get the package version"""
    assert lephare.__version__ is not None
