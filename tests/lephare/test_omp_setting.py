import importlib
import multiprocessing
import os


def test_omp_setting():
    max_cpus = multiprocessing.cpu_count()
    os.environ["OMP_NUM_THREADS"] = "1"
    import lephare as lp

    assert os.environ["OMP_NUM_THREADS"] == "1"

    os.environ["OMP_NUM_THREADS"] = "-5"
    importlib.reload(lp)
    assert os.environ["OMP_NUM_THREADS"] == "1"

    os.environ["OMP_NUM_THREADS"] = "50"
    importlib.reload(lp)
    assert os.environ["OMP_NUM_THREADS"] == str(max(1, max_cpus - 1))
