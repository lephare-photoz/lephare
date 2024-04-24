"""There's a lot of branching logic in download_registry_from_github.

I did my best to keep it simple, but it's still a little tricky.

Here's how it works:

Call download_registry_from_github
1. Local registry file does not exist → Download remote registry file
    1. Successfully downloaded registry file → Exit (now have local registry)
    2. Fail to download registry file → raise Exception
2. Local registry file exists
    → Need to check if local registry is up to date
    → Call _check_registry_is_latest_version
        1. True → Exit (confirmed local registry is up to date)
        2. False → Download updated version
            1. Successfully downloaded registry file → Exit (local registry updated)
            2. Fail to download registry file → raise Exception
        3. Failed to download hash file → raise Exception
"""

import os
import tempfile
from unittest.mock import patch

import pytest
import requests
from lephare.data_retrieval import (
    download_registry_from_github,
)


@patch("requests.get")
def test_download_registry_success(mock_get):
    # 1. Local registry does not exist
    #     1. Successfully downloaded registry file
    mock_get.return_value.status_code = 200
    mock_get.return_value.text = "file1\nfile2\nfile3"

    with tempfile.TemporaryDirectory() as tmpdir:
        download_registry_from_github(outfile=os.path.join(tmpdir, "registry.txt"))

        with open(os.path.join(tmpdir, "registry.txt"), "r") as file:
            assert file.read() == "file1\nfile2\nfile3"


@patch("requests.get")
def test_download_registry_failure(mock_get):
    # 1. Local registry does not exist
    #     2. Fail to download registry file
    mock_get.return_value.status_code = 404
    with pytest.raises(requests.exceptions.HTTPError):
        download_registry_from_github()


@patch("os.path.isfile")
@patch("pooch.file_hash")
@patch("requests.get")
def test_update_registry_hash_matches(mock_get, mock_file_hash, mock_isfile):
    # 2. Local registry exists
    #     1. _check_registry_is_latest_version == True
    mock_isfile.return_value = True
    mock_file_hash.return_value = "registryhash123"
    mock_get.return_value.status_code = 200
    mock_get.return_value.text = "registryhash123"

    with tempfile.TemporaryDirectory() as tmpdir:
        outfile = os.path.join(tmpdir, "registry.txt")
        download_registry_from_github(url="http://example.com/data_registry.txt", outfile=outfile)

        assert mock_get.call_count == 1
        assert mock_get.call_args[0][0] == "http://example.com/data_registry_hash.sha256"


@patch("os.path.isfile")
@patch("pooch.file_hash")
@patch("requests.get")
def test_update_registry_hash_mismatches(mock_get, mock_file_hash, mock_isfile):
    # 2. Local registry exists
    #     2. _check_registry_is_latest_version == False
    #       2. Successfully downloaded registry file
    mock_isfile.return_value = True
    mock_file_hash.return_value = "registryhash123"
    mock_get.return_value.status_code = 200
    mock_get.return_value.text = "hash_doesn't_match123"

    with tempfile.TemporaryDirectory() as tmpdir:
        outfile = os.path.join(tmpdir, "registry.txt")
        download_registry_from_github(url="http://example.com/data_registry.txt", outfile=outfile)

        # One call to download the hash file, one call to download the full registry file:
        assert mock_get.call_count == 2
        # The following set of [0][0][0] and such is because the call args list is
        #     [call('http://example.com/data_registry_hash.sha256', timeout=60),
        #      call('http://example.com/data_registry.txt', timeout=60)]
        # and we're only interested in checking the urls:
        assert mock_get.call_args_list[0][0][0] == "http://example.com/data_registry_hash.sha256"
        assert mock_get.call_args_list[1][0][0] == "http://example.com/data_registry.txt"


@patch("os.path.isfile")
@patch("requests.get")
def test_update_registry_hash_mismatches_and_download_fails(mock_get, mock_isfile):
    # 2. Local registry exists
    #     2. _check_registry_is_latest_version == False
    #       2. Fail to download registry file
    mock_isfile.return_value = True
    mock_get.return_value.text = "file1\nfile2\nfile3"

    mock_get.return_value.status_code = 404
    with pytest.raises(requests.exceptions.HTTPError):
        download_registry_from_github()
