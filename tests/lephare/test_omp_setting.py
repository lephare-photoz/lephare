import os
import importlib
import multiprocessing

def test_omp_setting():
    max_cpus = multiprocessing.cpu_count()
    os.environ["OMP_NUM_THREADS"] = "5"
    import lephare as lp
    assert os.environ["OMP_NUM_THREADS"] == "5"

    os.environ["OMP_NUM_THREADS"] = "-5"
    importlib.reload(lp)
    assert os.environ["OMP_NUM_THREADS"] == "1" 

    os.environ["OMP_NUM_THREADS"] = "50"
    importlib.reload(lp)
    assert os.environ["OMP_NUM_THREADS"] == str(max_cpus - 1)
