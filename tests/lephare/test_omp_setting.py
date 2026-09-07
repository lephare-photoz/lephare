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


def test_omp_setting_unset(monkeypatch, capsys):
    """When OMP_NUM_THREADS is unset it defaults to one fewer than the CPU count."""
    from lephare._set_omp_num_threads import _set_omp_num_threads

    max_cpus = multiprocessing.cpu_count()
    monkeypatch.delenv("OMP_NUM_THREADS", raising=False)
    _set_omp_num_threads()

    assert os.environ["OMP_NUM_THREADS"] == str(max(1, max_cpus - 1))
    assert f"max available CPUs: {max_cpus}" in capsys.readouterr().out


def test_omp_setting_not_an_integer(monkeypatch):
    """A non-numeric OMP_NUM_THREADS is replaced rather than allowed to raise."""
    from lephare._set_omp_num_threads import _set_omp_num_threads

    max_cpus = multiprocessing.cpu_count()
    monkeypatch.setenv("OMP_NUM_THREADS", "not-a-number")
    _set_omp_num_threads()

    assert os.environ["OMP_NUM_THREADS"] == str(max(1, max_cpus - 1))


def test_omp_setting_in_range_is_respected(monkeypatch):
    """A value the machine can honour is left untouched."""
    from lephare._set_omp_num_threads import _set_omp_num_threads

    monkeypatch.setenv("OMP_NUM_THREADS", "1")
    _set_omp_num_threads()
    assert os.environ["OMP_NUM_THREADS"] == "1"


def test_omp_setting_zero_is_clamped(monkeypatch):
    """Zero is below the minimum and is clamped up to one."""
    from lephare._set_omp_num_threads import _set_omp_num_threads

    monkeypatch.setenv("OMP_NUM_THREADS", "0")
    _set_omp_num_threads()
    assert os.environ["OMP_NUM_THREADS"] == "1"
