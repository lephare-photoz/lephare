import os
import tempfile
from pathlib import Path
from unittest.mock import mock_open, patch

import pytest
import requests
from lephare.data_retrieval import (
    _check_downloaded_files,
    _create_directories_from_files,
    download_all_files,
    download_registry_from_github,
    filter_files_by_prefix,
    make_default_retriever,
    make_retriever,
    read_list_file,
)

# TODO: this will be bundled into a module in the future,
# so replace the hardcoding when that happens
DEFAULT_BASE_DATA_URL = "https://raw.githubusercontent.com/lephare-photoz/lephare-data/main/"
DEFAULT_REGISTRY_FILE = "data_registry.txt"
DEFAULT_LOCAL_DATA_PATH = "./data"


def test_filter_file_by_prefix(test_data_dir):
    file_path = os.path.join(test_data_dir, "test_file_names.list")
    target_prefixes = ["prefix1", "prefix2"]
    expected_lines = ["prefix1_file1", "prefix2_file2"]
    assert filter_files_by_prefix(file_path, target_prefixes) == expected_lines


@patch("requests.get")
def test_download_registry_from_github_success(mock_get):
    mock_get.return_value.status_code = 200
    mock_get.return_value.text = "file1\nfile2\nfile3"

    with tempfile.TemporaryDirectory() as tmpdir:
        download_registry_from_github(outfile=os.path.join(tmpdir, "registry.txt"))

        with open(os.path.join(tmpdir, "registry.txt"), "r") as file:
            assert file.read() == "file1\nfile2\nfile3"


@patch("requests.get")
def test_download_registry_from_github_failure(mock_get):
    mock_get.return_value.status_code = 404
    with pytest.raises(requests.exceptions.HTTPError):
        download_registry_from_github()


def test_read_list_file(test_data_dir):
    file_path = os.path.join(test_data_dir, "test_file_names.list")
    expected_files = ["prefix1_file1", "prefix2_file2"]
    assert read_list_file(file_path) == expected_files


@patch("requests.get")
def test_read_list_file_remote(mock_get):
    mock_get.return_value.status_code = 200
    mock_get.return_value.text = "file1\tA\nfile2\tA\nfile3\tA"

    expected_files = ["file1", "file2", "file3"]
    assert read_list_file("http://example.com") == expected_files


def test_make_default_retriever(data_registry_file):
    mock_registry_content = "file1.txt hash1\nfile2.txt hash2"
    with patch("builtins.open", mock_open(read_data=mock_registry_content)):
        retriever = make_default_retriever()
        assert retriever.base_url == DEFAULT_BASE_DATA_URL
        assert retriever.path == Path(DEFAULT_LOCAL_DATA_PATH)
        assert retriever.registry["file1.txt"] == "hash1"
        assert retriever.registry["file2.txt"] == "hash2"


def test_make_retriever(data_registry_file):
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


@patch("lephare.data_retrieval.download_file")
def test_download_all_files(mock_download_file, data_registry_file):
    retriever = make_retriever(registry_file=data_registry_file)
    file_names = ["file1.txt", "file2.txt"]
    download_all_files(retriever, file_names)
    mock_download_file.assert_any_call(retriever, "file1.txt", ignore_registry=False)
    mock_download_file.assert_any_call(retriever, "file2.txt", ignore_registry=False)
    assert mock_download_file.call_count == len(file_names)
    # TODO could stand to expand this test
    # Additionally, would be nice to explicitly test single file download
