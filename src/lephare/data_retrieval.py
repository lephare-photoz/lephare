"""This module provides functionality for downloading and managing data files using pooch."""

import concurrent.futures
import os
import warnings
from functools import partial
from urllib.parse import urljoin, urlparse

import numpy as np
import pooch
import requests

from lephare import LEPHAREDIR
from lephare.prepare import all_types_to_keymap

DEFAULT_BASE_DATA_URL = "https://raw.githubusercontent.com/lephare-photoz/lephare-data/main/"
DEFAULT_REGISTRY_FILE = "data_registry.txt"
DEFAULT_LOCAL_DATA_PATH = LEPHAREDIR

# If a file is not downloaded the first time, retry this many times
MAX_RETRY_ATTEMPTS = 2

__all__ = [
    "download_all_files",
    "download_file",
    "download_registry_from_github",
    "filter_files_by_prefix",
    "make_default_retriever",
    "make_retriever",
    "read_list_file",
]


def filter_files_by_prefix(file_path, target_prefixes):
    """Returns all lines in a file that contain any of the target prefixes.

    Parameters
    ----------
    file_path : str
        The path to the file.
    target_prefixes : list
        A list of target prefixes to check for in each line.

    Returns
    -------
    list
        A list of lines that contain one of the target prefixes.
    """
    matching_lines = []
    with open(file_path, "r", encoding="utf-8") as file:
        for line in file:
            if any(line.startswith(prefix) for prefix in target_prefixes):
                matching_lines.append(line.split(" ")[0].strip())
    return matching_lines


def _check_registry_is_latest_version(remote_registry_url, local_registry_file):
    """Checks whether the local registry file is the latest version compared to a remote registry.

    Parameters
    ----------
    remote_registry_url : str
        The URL to the remote registry file, used to construct the URL to fetch the remote hash.
    local_registry_file : str
        The path to the local registry file whose up-to-date status is to be checked.

    Returns
    -------
    bool
        Returns True if the local registry file is up to date, otherwise False.

    Notes
    -----
    We make the assumption that the hash file for the registry will be stored in
    the same directory as the registry file, with the same name (sans extension)
    plus "_hash.sha256".

    Raises
    ------
    Exception
        If there is any problem fetching the registry hash file, including network issues,
        server errors, or other HTTP errors.
    """
    local_registry_hash = pooch.file_hash(local_registry_file, alg="sha256")
    remote_hash_url = os.path.splitext(remote_registry_url)[0] + "_hash.sha256"

    remote_hash_response = requests.get(remote_hash_url, timeout=60)
    remote_hash_response.raise_for_status()  # Raise exceptions for non-200 status codes

    return remote_hash_response.text.strip() == local_registry_hash


def download_registry_from_github(url="", outfile=""):
    """Fetch the contents of a file from a GitHub repository.

    Parameters
    ----------
    url : str
        The URL of the registry file. Defaults to a "data_registry.txt" file at
        DEFAULT_BASE_DATA_URL.
    outfile : str
        The path where the file will be saved. Defaults to DEFAULT_REGISTRY_FILE.

    Raises
    ------
    Exception
        If there is any problem fetching the registry hash file or full registry file,
        including network issues, server errors, or other HTTP errors.
    """
    remote_registry_name = "data_registry.txt"

    # Assign defaults if keywords left blank
    if url == "":
        url = urljoin(DEFAULT_BASE_DATA_URL, remote_registry_name)
    if outfile == "":
        outfile = DEFAULT_REGISTRY_FILE

    # If local registry hash matches remote hash, our registry is already up-to-date:
    if os.path.isfile(outfile) and _check_registry_is_latest_version(url, outfile):
        print(f"Local registry file is up to date: {outfile}")
        return

    # Download the registry file
    response = requests.get(url, timeout=120)
    response.raise_for_status()  # Raise exceptions for non-200 status codes

    with open(outfile, "w", encoding="utf-8") as file:
        file.write(response.text)

    print(f"Registry file downloaded and saved as {outfile}.")


def read_list_file(list_file, prefix=""):
    """Reads file names from a list file and returns a list of file paths.

    Parameters
    ----------
    list_file : str
        The name of the file containing the list of filenames. Can be local or a URL.

    prefix : str
        Optional prefix to add to all file names. When downloaded, file paths
        must be relative to the "base url," which is the top-level directory.

        Prefixes will be inferred from list_file paths or urls that contain
        "sed" or "filt"; otherwise; they should be manually specified.


    Returns
    -------
    list of str
        A list of file paths read from the list file.
    """
    file_names = []

    # Check if the list_file is a URL
    if urlparse(list_file).scheme in ("http", "https"):
        response = requests.get(list_file, timeout=60)
        response.raise_for_status()
        content = response.text
    else:
        with open(list_file, "r", encoding="utf-8") as file:
            content = file.read()

    # Infer the prefix if not provided
    # Note: pooch docs specify that registry files use Unix separators
    # Note as well: this may be phased out, if we decide to specify list
    #   files as containing paths relative to the root dir
    if prefix == "":
        if "sed" in list_file:
            start_index = list_file.find("sed/")
            end_index = list_file.rfind("/")
            prefix = list_file[start_index:end_index]
        elif "filt" in list_file:
            start_index = list_file.find("filt/")
            end_index = list_file.rfind("/")
            prefix = list_file[start_index:end_index]

    # Read in file
    for line in content.splitlines():
        file_name = line.split()[0].strip()
        if file_name[0] != "#":
            file_names.append(os.path.join(prefix, file_name))
    return file_names


def make_default_retriever():
    """Create a retriever with the default settings."""
    return make_retriever(
        base_url=DEFAULT_BASE_DATA_URL, registry_file=DEFAULT_REGISTRY_FILE, data_path=DEFAULT_LOCAL_DATA_PATH
    )


def make_retriever(
    base_url=DEFAULT_BASE_DATA_URL,
    registry_file=DEFAULT_REGISTRY_FILE,
    data_path=DEFAULT_LOCAL_DATA_PATH,
):
    """Create a retriever for downloading files.

    Parameters
    ----------
    base_url : str, optional
        The base URL for the data files.
    registry_file : str, optional
        The path to the registry file that lists the files and their hashes.
    data_path : str, optional
        The local path where the files will be downloaded.

    Returns
    -------
    pooch.Pooch
        The retriever object for downloading files.
    """
    retriever = pooch.create(
        base_url=base_url,
        path=data_path,
        registry=None,  # We're using a registry file instead (set below)
    )
    retriever.load_registry(registry_file)
    return retriever


def _create_directories_from_files(file_names):
    """Create directories for the given file names if they do not already exist.

    This function is for thread safety when downloading files in parallel.

    Parameters
    ----------
    file_names : list of str
        List of file names with absolute paths.
    """
    unique_directories = set(
        os.path.dirname(file_name) for file_name in file_names if os.path.dirname(file_name)
    )
    for directory in unique_directories:
        if not os.path.exists(directory):
            os.makedirs(directory)
            print(f"Created directory: {directory}")


def download_file(retriever, file_name, ignore_registry=False):
    """Download a file using the retriever, optionally ignoring the registry.

    Parameters
    ----------
    retriever : pooch.Pooch
        The retriever object for downloading files.
    file_name : str
        The name of the file to download.
    ignore_registry : bool
        If True, download the file without checking its hash against the registry.

    Returns
    -------
    str
        The path to the downloaded file.
    """
    if ignore_registry:
        print(f"Downloading without registry: {file_name}...")
        return pooch.retrieve(
            url=urljoin(retriever.base_url, file_name),
            known_hash=None,
            fname=file_name,
            path=retriever.path,
        )
    else:
        return retriever.fetch(file_name)


def download_all_files(retriever, file_names, ignore_registry=False, retry=MAX_RETRY_ATTEMPTS):
    """Download all files in the given list using the retriever.

    Parameters
    ----------
    retriever : pooch.Pooch
        The retriever object for downloading files.
    file_names : list of str
        List of file names to download.
    ignore_registry : bool
        If True, download the files without checking their hashes against the registry.
    retry : int
        Number of times to retry downloading a file if first attempt fails.
    """
    if len(file_names) == 0:
        print("Download all files called for list of 0 files; done.")
        return

    # First make directories, for thread safety
    absolute_file_names = [os.path.join(retriever.path, file_name) for file_name in file_names]
    _create_directories_from_files(absolute_file_names)

    # Now the downloading
    print(f"Checking/downloading {len(file_names)} files...")
    with concurrent.futures.ThreadPoolExecutor() as executor:
        download_fn = partial(download_file, retriever, ignore_registry=ignore_registry)
        futures = [executor.submit(download_fn, file_name) for file_name in file_names]

        # We're gathering the completed futures here to make sure we aren't skipping any files,
        # which seemed to be happening earlier, when using an executor mapping function instead
        completed_futures = []
        for future in concurrent.futures.as_completed(futures):
            try:
                completed_futures.append(future.result(timeout=60))  # timeout is in seconds
            except TimeoutError as e:
                print(f"Future completed with a timeout exception: {e}")
            except Exception as e:
                print(f"Future completed with an exception: {e}")

        print(f"{len(completed_futures)} completed.")

    # Finish with some checks on our downloaded files
    all_files_present = _check_downloaded_files(absolute_file_names, completed_futures)

    if not all_files_present and retry > 0:
        print("Retrying download for missing files...")
        download_all_files(retriever, file_names, ignore_registry=ignore_registry, retry=retry - 1)


def _check_downloaded_files(file_names, completed_futures):
    """Check if all files have been downloaded successfully and are not empty.

    Parameters
    ----------
    file_names : list of str
        List of expected file names.
    completed_futures : list of str
        List of file names that were downloaded.

    Returns
    -------
    bool
        True if all files are downloaded and non-empty, False otherwise.
    """
    # Check if all files were downloaded
    missing_files = False
    potentially_missing_files = set(file_names) - set(completed_futures)
    for file in potentially_missing_files:
        if not os.path.exists(file):
            print("The following file was not downloaded:", file)
            missing_files = True
    if missing_files:
        return False

    # Check if any downloaded file is empty
    for file_name in completed_futures:
        if os.path.getsize(file_name) == 0:
            print(f"The file {file_name} is empty.")
            return False

    print("All files downloaded successfully and are non-empty.")
    return True


def config_to_required_files(keymap, base_url=None):
    """Take a lephare config and return list of auxiliary files required for run.

    For the sed lists these must be present in the auxiliary files directory. If
    local full paths set the code will only retrieve opa, vega, and filters.

    In addition to the specified files we also add opa and vega files. These
    are always required.

    We use the tau opacities by default.

    Parameters
    ==========
    keymap : dict of lephare.keyvalue
        The dictionary of config keys containing filters etc required.
    base_url : str
        Url to overwrite default base.
    """
    if base_url is None:
        base_url = DEFAULT_BASE_DATA_URL
    required_files = []
    # We always need alloutputkeys.txt
    required_files += ["alloutputkeys.txt"]
    # Opacity always required
    opa_list = ["opa/OPACITY.dat"] + [f"opa/tau{i:02d}.out" for i in np.arange(81)]
    required_files += opa_list
    # vega always required
    vega_list = [
        "vega/BD+17.sed",
        "vega/BD+17o4708.sed",
        "vega/SunLCB.sed",
        "vega/VegaLCB.sed",
        "vega/a0v.sed",
        "vega/a0v_n.sed",
    ]
    required_files += vega_list
    required_files += [f"filt/{f}" for f in keymap["FILTER_LIST"].value.split(",")]
    # Get user specified sed lists
    sed_keys = ["STAR_SED", "GAL_SED", "QSO_SED"]
    for key in sed_keys:
        try:
            list_file = base_url + keymap[key].value
            required_files += [keymap[key].value]
            file_names = read_list_file(list_file, prefix=f"sed/{key.split('_')[0]}/")
            required_files += file_names
        except KeyError:
            warnings.warn(f"{key} keyword not set or not present in auxiliary files directory.")
    # Bethermin12 always required
    bet_list = "sed/GAL/BETHERMIN12/BETHERMIN12_MOD.list"
    required_files += [bet_list]
    required_files += read_list_file(base_url + bet_list, prefix="sed/GAL/")
    # Get extinction law files
    ext_list = [f"ext/{f}" for f in keymap["EXTINC_LAW"].value.split(",")]
    ext_list += ["ext/MW_seaton.dat"]  # Appears to be always required
    required_files += ext_list
    return required_files


def get_auxiliary_data(lephare_dir=LEPHAREDIR, keymap=None, additional_files=None):
    """Get all auxiliary data required to run lephare.

    This gets all the filters, seds, and other data files.

    If no keymap is set this will git clone the full repository.

    Parameters
    ==========
    lephare_dir : str
        The path to the lephare directory for auxiliary files.
    keymap : dict
        The config dictionary.
    additional_files : list
        Any additional files to be downloaded from the auxiliary file repo.
    """

    # ensure that all values in the keymap are keyword objects
    keymap = all_types_to_keymap(keymap)

    # Get the registry file
    download_registry_from_github()
    base_url = DEFAULT_BASE_DATA_URL
    repo_name = "lephare-data"
    repo_url = f"https://github.com/lephare-photoz/{repo_name}"
    registry_file = DEFAULT_REGISTRY_FILE
    data_path = lephare_dir
    if keymap is None:
        # Assume if filt is present assume everything is.
        if os.path.isdir(f"{lephare_dir}/filt"):
            warnings.warn(
                "Some data appears present. Not downloading."
                "Consider setting a keymap to download a subset."
            )
        else:
            # Get the full repository
            print(f"Downloading all auxiliary data (~1.5Gb) to {lephare_dir}.")
            print(f"Getting data from {repo_url}.")
            os.system(f"git clone {repo_url} {lephare_dir}")
    else:
        retriever = make_retriever(base_url=base_url, registry_file=registry_file, data_path=data_path)
        file_list = config_to_required_files(keymap)
        download_all_files(retriever, file_list, ignore_registry=False)
    if additional_files is not None:
        download_all_files(retriever, additional_files, ignore_registry=False)
