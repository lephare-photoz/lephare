def test_version():
    """Case 1: _version.py exists — normal import works."""
    # Create a fake _version module
    import lephare

    assert isinstance(lephare.__version__, str)
