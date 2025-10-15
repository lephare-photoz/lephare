import os
from pathlib import Path
from unittest.mock import mock_open, patch

import lephare as lp
import pooch
import pytest
from lephare.data_retrieval import (
    DEFAULT_BASE_DATA_URL,
    DEFAULT_LOCAL_DATA_PATH,
    MAX_RETRY_ATTEMPTS,
    _check_downloaded_files,
    _create_directories_from_files,
    download_all_files,
    download_file,
    filter_files_by_prefix,
    make_default_retriever,
    make_retriever,
    read_list_file,
)


def test_filter_file_by_prefix(test_data_dir: str):
    file_path = os.path.join(test_data_dir, "test_file_names.list")
    target_prefixes = ["prefix1", "prefix2"]
    expected_lines = ["prefix1_file1", "prefix2_file2"]
    assert filter_files_by_prefix(file_path, target_prefixes) == expected_lines


def test_read_list_file(test_data_dir: str):
    file_path = os.path.join(test_data_dir, "test_file_names.list")
    expected_files = ["prefix1_file1", "prefix2_file2"]
    assert read_list_file(file_path) == expected_files


@patch("requests.get")
def test_read_list_file_remote(mock_get):
    mock_get.return_value.status_code = 200
    mock_get.return_value.text = "file1\tA\nfile2\tA\nfile3\tA"

    expected_files = ["file1", "file2", "file3"]
    assert read_list_file("http://example.com") == expected_files


def test_make_default_retriever(data_registry_file: str):
    mock_registry_content = "file1.txt hash1\nfile2.txt hash2"

    # Create a mock_open instance
    m_open = mock_open(read_data=mock_registry_content)

    # Mock __enter__ and __exit__ for context management
    m_open.return_value.__enter__ = lambda self: self
    m_open.return_value.__exit__ = lambda self, exc_type, exc_val, exc_tb: None

    # Patch the built-in open method with the mocked version
    with patch("builtins.open", m_open):
        retriever = make_default_retriever()
        assert retriever.base_url == DEFAULT_BASE_DATA_URL
        assert retriever.path == Path(DEFAULT_LOCAL_DATA_PATH)
        assert retriever.registry["file1.txt"] == "hash1"
        assert retriever.registry["file2.txt"] == "hash2"


def test_make_retriever(data_registry_file: str):
    base_url = "http://example.com/"
    data_path = "/tmp"
    registry_file = data_registry_file
    retriever = make_retriever(base_url=base_url, data_path=data_path, registry_file=registry_file)
    assert retriever.base_url == base_url
    assert retriever.path == Path(data_path)
    assert retriever.registry["file1.txt"] == "hash1"
    assert retriever.registry["file2.txt"] == "hash2"
    with pytest.raises(KeyError):
        _ = retriever.registry["file3.txt"]


@patch("os.makedirs")
def test_create_directories_from_files(mock_makedirs):
    file_names = ["dir1/file1.txt", "dir2/file2.txt", "dir1/file3.txt"]
    _create_directories_from_files(file_names)
    mock_makedirs.assert_any_call("dir1")
    mock_makedirs.assert_any_call("dir2")
    assert mock_makedirs.call_count == 2


@patch("os.path.getsize")
def test_check_downloaded_files_success(mock_getsize):
    file_names = ["/tmp/file1.txt", "/tmp/file2.txt"]
    downloaded_files = ["/tmp/file1.txt", "/tmp/file2.txt"]
    mock_getsize.return_value = 10
    assert _check_downloaded_files(file_names, downloaded_files)


@patch("os.path.getsize")
def test_check_downloaded_files_missing(mock_getsize):
    file_names = ["file1.txt", "file2.txt"]
    downloaded_files = ["/tmp/file1.txt"]
    mock_getsize.return_value = 10
    assert not _check_downloaded_files(file_names, downloaded_files)


@patch("os.path.getsize")
def test_check_downloaded_files_empty(mock_getsize):
    file_names = ["file1.txt", "file2.txt"]
    downloaded_files = ["/tmp/file1.txt", "/tmp/file2.txt"]
    mock_getsize.side_effect = [10, 0]
    assert not _check_downloaded_files(file_names, downloaded_files)


def test_download_single_file(data_registry_file: str):
    retriever = make_retriever(registry_file=data_registry_file)
    downloader = pooch.HTTPDownloader(headers={"User-Agent": "LePHARE"})
    file_name = "file1.txt"
    with patch.object(retriever, "fetch", return_value="/path/to/downloaded/file1.txt") as mock_fetch:
        download_file(retriever, file_name, downloader=downloader)
        mock_fetch.assert_called_once_with(file_name, downloader=downloader)


def test_download_single_file_without_downloader(data_registry_file: str):
    """Simple test just to check it works without a downloader specified."""
    retriever = make_retriever(registry_file=data_registry_file)
    file_name = "file1.txt"
    with patch.object(retriever, "fetch", return_value="/path/to/downloaded/file1.txt"):
        download_file(retriever, file_name, downloader=None)


def test_download_single_file_ignore_registry(data_registry_file: str):
    retriever = make_retriever(registry_file=data_registry_file)
    file_name = "file1.txt"
    downloader = pooch.HTTPDownloader(headers={"User-Agent": "LePHARE"})
    with patch("pooch.retrieve", return_value="/path/to/downloaded/file1.txt") as mock_retrieve:
        download_file(retriever, file_name, ignore_registry=True, downloader=downloader)
        mock_retrieve.assert_any_call(
            url=f"{retriever.base_url}{file_name}",
            known_hash=None,
            fname=file_name,
            path=retriever.path,
            downloader=downloader,
        )


@patch("lephare.data_retrieval.download_file")
def test_download_all_files(mock_download_file, data_registry_file: str):
    """This test checks to make sure that we're calling download_file for each file in the list.
    and also accounts for the fact that we retry downloading each file MAX_RETRY_ATTEMPTS times."""
    retriever = make_retriever(registry_file=data_registry_file)

    # Test with an empty list of files
    download_all_files(retriever, [])
    assert mock_download_file.call_count == 0

    # Test with list of 2 files
    file_names = ["file1.txt", "file2.txt"]
    download_all_files(retriever, file_names)
    mock_download_file.assert_any_call(retriever, "file1.txt", ignore_registry=False)
    mock_download_file.assert_any_call(retriever, "file2.txt", ignore_registry=False)
    # `1 + MAX_RETRY_ATTEMPTS` because we try once, then retry MAX_RETRY_ATTEMPTS more times
    assert mock_download_file.call_count == len(file_names) * (1 + MAX_RETRY_ATTEMPTS)


@patch("lephare.data_retrieval.download_file")
def test_download_all_files_non_default_retry(mock_download_file, data_registry_file: str):
    """This test checks to make sure that we're calling download_file for each file in the list.
    and sets the retries to 0. This should only try to download each file once."""
    retriever = make_retriever(registry_file=data_registry_file)
    file_names = ["file1.txt", "file2.txt"]
    download_all_files(retriever, file_names, retry=0)
    mock_download_file.assert_any_call(retriever, "file1.txt", ignore_registry=False)
    mock_download_file.assert_any_call(retriever, "file2.txt", ignore_registry=False)
    assert mock_download_file.call_count == len(file_names)


def test_get_auxiliary_data(test_data_dir: str):
    """Check the file downloader works for a new filter file."""
    test_dir = os.path.abspath(os.path.dirname(__file__))
    os.environ["LEPHAREDIR"] = os.path.join(test_dir, "../data")
    os.environ["LEPHAREWORK"] = os.path.join(test_dir, "../tmp")
    config = lp.default_cosmos_config.copy()
    # Set to a filter that is not in the test directory
    config["FILTER_LIST"] = "lsst/total_g.dat"
    # Check the downloader runs
    lp.data_retrieval.get_auxiliary_data(keymap=config)
    filter_file_path = os.path.join(test_dir, "../data/filt/lsst/total_g.dat")
    # Check it was downloaded
    assert os.path.exists(filter_file_path)
    os.remove(filter_file_path)
    # Check it was deleted
    assert not os.path.exists(filter_file_path)
