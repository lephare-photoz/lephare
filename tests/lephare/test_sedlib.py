import pytest
import os
import lephare as lp 

def test_read_ages_from_file_valid_data(tmp_path):
    """
    Test that valid ages are read, converted to 1e9, 
    and prefixed with agemin/agemax.
    """
    # 1. Create a temporary file on the fly
    d = tmp_path / "sub"
    d.mkdir()
    file_path = d / "test_ages.txt"
    
    # Write mock data: a comment, a valid age, and another valid age
    content = "# This is a comment\n13.8\n10.5\n"
    file_path.write_text(content)

    # 2. Call the C++ function via Python
    agemin = 0.1
    agemax = 15.0
    # Convert path object to string for the C++ interface
    result = lp._read_ages_from_file(str(file_path), agemin, agemax)

    # 3. Assertions
    assert len(result) == 4
    assert result[0] == agemin
    assert result[1] == agemax
    assert result[2] == 13.8e9  # Check conversion to years
    assert result[3] == 10.5e9

def test_read_ages_from_file_none():
    """Test behavior when filename is 'none'."""
    result = lp._read_ages_from_file("none", 1.0, 2.0)
    assert result == [1.0, 2.0]

def test_read_ages_missing_file(capfd):
    """Test behavior when file does not exist (should print to cerr)."""
    result = lp._read_ages_from_file("non_existent.txt", 0.0, 5.0)
    
    # Capture the output to stderr (cerr in C++)
    captured = capfd.readouterr()
    assert "Can't open file" in captured.err
    assert result == [0.0, 5.0]

def test_closeAge_select_all():
    """If ageSel only contains 2 elements (min/max), all ages should be True."""
    age_library = [1.0e6, 2.0e6, 3.0e6]
    # Only bounds provided, no target ages
    age_sel = [0.0, 1.0e12] 
    
    result = lp._closeAge(age_sel, age_library)
    assert result == [True, True, True]

def test_closeAge_exact_match():
    """Test that target ages match the closest library entry within range."""
    age_library = [1.0e9, 2.0e9, 3.0e9]
    # min=0, max=5e9, targets=[2.0e9]
    age_sel = [0.0, 5.0e9, 2.0e9]
    
    result = lp._closeAge(age_sel, age_library)
    # Only the second element (2.0e9) should be True
    assert result == [False, True, False]

def test_closeAge_out_of_range_or_tolerance():
    """Test that ages outside the [min, max] or tolerance are rejected."""
    age_library = [1.0e9, 4.0e9]
    
    # Target 4.0e9 is valid, but agemax is set to 3.0e9
    age_sel_out_of_range = [0.0, 3.0e9, 4.0e9]
    res1 = lp._closeAge(age_sel_out_of_range, age_library)
    assert res1 == [False, False]

    # Target 1.0001e9 is too far from 1.0e9 (diff > 1e-5 tolerance)
    age_sel_tolerance = [0.0, 5.0e9, 1.0001e9]
    res2 = lp._closeAge(age_sel_tolerance, age_library)
    assert res2 == [False, False]

def test_closeAge_multiple_targets():
    """Test matching multiple targets to the library."""
    age_library = [1.0, 2.0, 3.0, 4.0, 5.0]
    # Targets 1.0 and 4.0
    age_sel = [0.0, 10.0, 1.0, 4.0]
    
    result = lp._closeAge(age_sel, age_library)
    assert result == [True, False, False, True, False]
