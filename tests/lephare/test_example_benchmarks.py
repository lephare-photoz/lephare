"""Tests for the sample functions exercised by the asv benchmark suite."""

from lephare import example_benchmarks


def test_runtime_computation_sleeps(monkeypatch):
    """runtime_computation sleeps for a bounded, non-negative interval."""
    slept = []
    monkeypatch.setattr(example_benchmarks.time, "sleep", slept.append)

    assert example_benchmarks.runtime_computation() is None

    (duration,) = slept
    assert 0 <= duration <= 5


def test_memory_computation_length_is_bounded():
    """memory_computation returns a zero-filled list of at most 512 entries."""
    result = example_benchmarks.memory_computation()

    assert isinstance(result, list)
    assert 0 <= len(result) <= 512
    assert set(result) <= {0}
