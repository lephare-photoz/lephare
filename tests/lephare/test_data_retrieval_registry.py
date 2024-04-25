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
from unittest.mock import Mock, patch

import pytest
import requests
from lephare.data_retrieval import (
    download_registry_from_github,
)


def test_download_registry_success():
    # 1. Local registry does not exist (so no mocking needed here)
    #     1. Successfully downloaded registry file

    # Mock remote registry file
    mock_response = Mock()
    mock_response.status_code = 200
    mock_response.text = "file1\nfile2\nfile3"
    with patch("requests.get", return_value=mock_response) as mock_get_remote_registry:  # noqa: F841
        with tempfile.TemporaryDirectory() as tmp_dir:
            registry_outfile = os.path.join(tmp_dir, "registry.txt")
            download_registry_from_github(outfile=registry_outfile)
            # Check that we can open it (and it contains expected content)
            with open(registry_outfile, "r") as file:
                assert file.read() == "file1\nfile2\nfile3"


def test_download_registry_failure():
    # 1. Local registry does not exist (no mocking needed)
    #     2. Fail to download registry file

    # Mock failed registry file download
    mock_response = Mock()
    mock_response.raise_for_status.side_effect = requests.exceptions.HTTPError(
        "404 Client Error: Not Found for url"
    )
    with patch("requests.get", return_value=mock_response) as mock_get_remote_registry:  # noqa: F841
        with pytest.raises(requests.exceptions.HTTPError):
            download_registry_from_github()


def test_update_registry_hash_matches():
    # 2. Local registry exists
    #     1. _check_registry_is_latest_version == True

    # Mock the local registry file existing
    with patch("os.path.isfile", return_value=True) as mock_local_registry_existing:  # noqa: F841
        # Mock local registry having a certain pooch hash
        with patch("pooch.file_hash", return_value="registryhash123") as mock_local_registry_hash:  # noqa: F841
            # Mock getting the remote registry hash file
            mock_response = Mock()
            mock_response.status_code = 200
            mock_response.text = "registryhash123"
            with patch("requests.get", return_value=mock_response) as mock_get_remote_hash_file:
                # Call the function
                with tempfile.TemporaryDirectory() as tmp_dir:
                    registry_outfile = os.path.join(tmp_dir, "registry.txt")
                    download_registry_from_github(
                        url="http://example.com/data_registry.txt", outfile=registry_outfile
                    )
                    # Assert that we only make 1 request for the hash file, and that we used the right url:
                    assert mock_get_remote_hash_file.call_count == 1
                    assert (
                        mock_get_remote_hash_file.call_args[0][0]
                        == "http://example.com/data_registry_hash.sha256"
                    )


def test_update_registry_hash_mismatches():
    # 2. Local registry exists
    #     2. _check_registry_is_latest_version == False
    #         1. Successfully downloaded registry file

    # Mock the local registry file existing
    with patch("os.path.isfile", return_value=True) as mock_local_registry_existing:  # noqa: F841
        # Mock local registry having a certain pooch hash
        with patch("pooch.file_hash", return_value="registryhash123") as mock_local_registry_hash:  # noqa: F841
            # Mock getting the remote hash/registry files
            mock_hash_response = Mock()
            mock_hash_response.status_code = 200
            mock_hash_response.text = "hash_doesn't_match123"

            mock_registry_response = Mock()
            mock_registry_response.status_code = 200
            mock_registry_response.text = "file1\nfile2\nfile3"

            def which_mock_get(*args, **kwargs):
                url = args[0]
                if "hash.sha256" in url:
                    return mock_hash_response
                else:
                    return mock_registry_response

            with patch("requests.get", side_effect=which_mock_get) as mock_remote_files_get:
                # Call the function
                with tempfile.TemporaryDirectory() as tmp_dir:
                    registry_outfile = os.path.join(tmp_dir, "registry.txt")
                    download_registry_from_github(
                        url="http://example.com/data_registry.txt", outfile=registry_outfile
                    )
                    # Checks:
                    # One call to download the hash file, one call to download the full registry file:
                    assert mock_remote_files_get.call_count == 2
                    # The following set of [0][0][0] and such is because the call args list is
                    #     [call('http://example.com/data_registry_hash.sha256', timeout=60),
                    #      call('http://example.com/data_registry.txt', timeout=60)]
                    # and we're only interested in checking the urls:
                    assert (
                        mock_remote_files_get.call_args_list[0][0][0]
                        == "http://example.com/data_registry_hash.sha256"
                    )
                    assert (
                        mock_remote_files_get.call_args_list[1][0][0]
                        == "http://example.com/data_registry.txt"
                    )


def test_update_registry_hash_mismatches_and_download_fails():
    # 2. Local registry exists
    #     2. _check_registry_is_latest_version == False
    #         2. Fail to download registry file

    # Mock the local registry file existing
    with patch("os.path.isfile", return_value=True) as mock_local_registry_existing:  # noqa: F841
        # Mock local registry having a certain pooch hash
        with patch("pooch.file_hash", return_value="registryhash123") as mock_local_registry_hash:  # noqa: F841
            # Mock getting the remote hash/registry files
            mock_hash_response = Mock()
            mock_hash_response.status_code = 200
            mock_hash_response.text = "hash_doesn't_match123"

            mock_registry_response = Mock()
            mock_registry_response.raise_for_status.side_effect = requests.exceptions.HTTPError(
                "404 Client Error: Not Found for url"
            )

            def which_mock_get(*args, **kwargs):
                url = args[0]
                if "hash.sha256" in url:
                    return mock_hash_response
                else:
                    return mock_registry_response

            with patch("requests.get", side_effect=which_mock_get) as mock_remote_files_get:  # noqa: F841
                # Check that we raise HTTPError as expected
                with pytest.raises(requests.exceptions.HTTPError):
                    download_registry_from_github()


def test_update_registry_hash_download_fails():
    # 2. Local registry exists
    #     3. Fail to download registry hash file

    # Mock the local registry file existing
    with patch("os.path.isfile", return_value=True) as mock_local_registry_existing:  # noqa: F841
        # Mock local registry having a certain pooch hash
        with patch("pooch.file_hash", return_value="registryhash123") as mock_local_registry_hash:  # noqa: F841
            # Mock getting the remote hash/registry files
            mock_hash_response = Mock()
            mock_hash_response.raise_for_status.side_effect = requests.exceptions.HTTPError(
                "404 Client Error: Not Found for url"
            )

            with patch("requests.get", return_value=mock_hash_response) as mock_get_remote_hash_file:  # noqa: F841
                # Check that we get the expected exception
                with pytest.raises(requests.exceptions.HTTPError):
                    download_registry_from_github()
