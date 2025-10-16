import os
import shutil

import lephare as lp


def test_filter_extinc():
    """Test the filter_extinc runner makes a file and has good values"""
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    os.makedirs(os.path.join(test_dir, "../tmp", "filt"), exist_ok=True)
    # Move filter file into the tmp WORK directory
    src = os.path.join(
        test_dir, "../data", "filt", "LSST_FILTERS.dat"
    )  # replace 'myfile.txt' with your filename
    dst = os.path.join(test_dir, "../tmp", "filt", "LSST_FILTERS.dat")
    shutil.copy(src, dst)
    out_file = os.path.join(test_dir, "tmp", "filter_extinc.dat")
    options = {
        "FILTER_FILE": "LSST_FILTERS.dat",
        "EXT_CURVE": "SB_calzetti.dat",
        "GAL_CURVE": "CARDELLI",
        "OUTPUT": out_file,
    }
    runner = lp.FiltExt(config_keymap=lp.all_types_to_keymap(lp.default_cosmos_config), **options)
    runner.run()
    # Check it made the file
    assert os.path.exists(out_file)
    with open(out_file, "r") as f:
        contents = f.read()
    assert contents.split()[-1] == "1.3039147701114893"
