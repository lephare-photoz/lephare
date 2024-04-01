""" TODO write this docstring"""

"""
thoughts after talking to Sandro - loading data-registry as data frame
"""

# TODO isort at the end
import concurrent.futures
import os
from urllib.parse import urlparse

import pooch
import requests

DEFAULT_BASE_DATA_URL = "https://raw.githubusercontent.com/OliviaLynn/LEPHARE-data/main/"
DEFAULT_REGISTRY_FILE = "data_registry.txt"

#! Replace DEFAULT_LOCAL_DATA_PATH with the following:
# from lephare import data_marshaller
# DEFAULT_LOCAL_DATA_PATH = data_marshaller.get_data_path() # likely something like: ~/Library/Caches/lephare/data/
DEFAULT_LOCAL_DATA_PATH = "./data"


__all__ = [  # TODO alphabetical maybe. also at the end, prune this to be sure it's all used
    "read_list_file",
    "filter_files_by_prefix",
    "download_registry_from_github",
    "make_retriever",
    "download_file",
    "download_all_files",
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
    with open(file_path, "r") as file:
        for line in file:
            if any(line.startswith(prefix) for prefix in target_prefixes):
                matching_lines.append(line.split(" ")[0].strip())
    return matching_lines


def download_registry_from_github(url, outfile):
    """Fetch the contents of a file from a GitHub repository.

    Parameters
    ----------
    url : str
        The URL of the file in the GitHub repository.

    Returns
    -------
    str
        The contents of the file.
    """
    response = requests.get(url)
    if response.status_code == 200:
        with open(outfile, "w") as file:
            file.write(response.text)
        print(f"File downloaded and saved as {outfile}")
    else:
        raise Exception(f"Failed to fetch file: {response.status_code}")


def read_list_file(list_file, prefix=""):
    """Reads file names from a list file, returns list of file paths.

    Parameters
    ----------
    list_file : str
        The name of the file containing the list of filenames.
        Can be local or a url

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
        response = requests.get(list_file)
        response.raise_for_status()
        content = response.text
    else:
        with open(list_file, "r") as file:
            content = file.read()

    # Infer the prefix if not provided
    if prefix == "":
        if "sed" in list_file:
            start_index = list_file.find("sed/")
            end_index = list_file.rfind("/")
            prefix = list_file[start_index:end_index]
        elif "filt" in list_file:
            start_index = list_file.find("filt/")
            end_index = list_file.rfind("/")
            prefix = list_file[start_index:end_index]

    # Avoid ending up with double seperators after prefix
    # TODO : add consideration about using '/'. Pooch uses this separator in
    # registries, and we've discussed not targetting windows for our builds,
    # but it could be nice to say something about this somewhere.
    # Also, pooch says it will handle conversions when checking registry,
    # so we actually could go ahead and add this
    if len(prefix) > 0 and prefix[-1] == "/":
        prefix = prefix[:-1]

    # Read in file
    for line in content.splitlines():
        file_name = line.split()[0].strip()
        file_names.append(f"{prefix}/{file_name}")
    return file_names


def make_default_retriever():
    """Create a retriever with the default settings."""
    return make_retriever(
        base_url=DEFAULT_BASE_DATA_URL, registry_file=DEFAULT_REGISTRY_FILE, data_path=DEFAULT_LOCAL_DATA_PATH
    )


def make_retriever(base_url, registry_file, data_path=None):
    """TODO docstring"""
    if not data_path:
        data_path = pooch.os_cache("lephare")

    retriever = pooch.create(
        path=data_path,  # TODO does the default work
        base_url=base_url,
        registry=None,  # We're using a registry file instead
    )
    retriever.load_registry(registry_file)
    return retriever


def download_file(retriever, file):
    """TODO docstring"""
    # TODO would be nice to have a part here that will let someone proceed
    # with download even if hashes don't match up
    return retriever.fetch(file)


def create_directories_from_files(file_names):
    """TODO docstring. Basically, we need this for thread safety it seems."""
    unique_directories = set(
        os.path.dirname(file_name) for file_name in file_names if os.path.dirname(file_name)
    )
    for directory in unique_directories:
        if not os.path.exists(directory):
            os.makedirs(directory)
            print(f"Created directory: {directory}")


def download_all_files(retriever, file_names, ignore_registry=False):
    """TODO docstring"""
    #! talk to Drew
    # maybe TODO, check if we need to first? but doesn't pooch do this for us already?
    # based on sorcha's: determine if we should attempt to download or create any files.
    # perhaps it's mostly there because of the force-refresh functionality they have
    # should we have that too? I'm leaning probably not, as I don't see us changing template
    #     files all that often. when/if we add versioning, that should handle the rare case we do

    # First make directories, for thread safety
    create_directories_from_files(file_names)

    # Create an unregistered retriever if needed
    if ignore_registry:
        """
        unregistered_retriever = pooch.create(
            path=retriever.abspath,
            base_url=retriever.base_url,
            registry=None,  # No registry
            urls={file_name: urljoin(retriever.base_url, file_name) for file_name in file_names}  # Construct URLs manually
        )
        """

        def unregistered_retrieve(file_name):
            print(f"unregistered_retrieve({file_name})...")
            return pooch.retrieve(
                url=retriever.base_url + file_name,
                known_hash=None,
                fname=file_name,
                path=retriever.abspath,
            )

    # Now the downloading
    print(f"Checking/downloading {len(file_names)} files...")
    with concurrent.futures.ThreadPoolExecutor() as executor:
        #! talk to Drew
        # executor.map(retriever.fetch, file_names, timeout=0.1) # timeout in seconds. Doesn't seem to make a difference?
        # executor.map(fetch_partial, file_names) # doesn't seem to work at all? is there something wrong with the partial? is this blocking somehow?

        # TODO might be better to move this branching into our single file download_file. actually yes I like that a lot.
        if ignore_registry:
            # futures = [executor.submit(unregistered_retriever.fetch, file_name) for file_name in file_names]
            futures = [executor.submit(unregistered_retrieve, file_name) for file_name in file_names]
            results = [future.result() for future in concurrent.futures.as_completed(futures)]
            print(f"{len(results)} completed.")
            for result in results:
                print(result)
        else:
            futures = [executor.submit(retriever.fetch, file_name) for file_name in file_names]

            results = [future.result() for future in concurrent.futures.as_completed(futures)]
            print(f"{len(results)} completed.")

    #! talk to Drew
    # TODO perhaps programmatically check that they all downloaded correctly?
    # or just compare len file_names and len results? raise exception if not same?
