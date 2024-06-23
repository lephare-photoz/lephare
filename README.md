# LePHARE

[![Template](https://img.shields.io/badge/Template-LINCC%20Frameworks%20Python%20Project%20Template-brightgreen)](https://lincc-ppt.readthedocs.io/en/latest/)

[![PyPI](https://img.shields.io/pypi/v/lephare?color=blue&logo=pypi&logoColor=white)](https://pypi.org/project/lephare/)
[![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/lincc-frameworks/lephare/smoke-test.yml)](https://github.com/lephare-photoz/lephare/actions/workflows/smoke-test.yml)
[![Codecov](https://codecov.io/gh/lephare-photoz/lephare/branch/main/graph/badge.svg)](https://codecov.io/gh/lephare-photoz/lephare)
[![Read The Docs](https://img.shields.io/readthedocs/lephare)](https://lephare.readthedocs.io/)

IMPORTANT! This project is in an early development stage. If you wish to use and run LePHARE please download it from the [official repository](https://gitlab.lam.fr/Galaxies/LEPHARE/).

LePHARE (PHotometric Analysis for Redshift Estimation) is a Python package built on a complete rewrite in C++ of the [Fortran code](https://www.cfht.hawaii.edu/~arnouts/LEPHARE/acknowledgement.html) LePhare.
LePHARE computes photometric redshifts and physical parameters by fitting spectral energy distributions (SED) to a dataset of photometric fluxes or apparent magnitudes.

## Installation

LePHARE has been compiled and deployed on pypi for several systems; a list can be found at https://pypi.org/project/lephare/#files.
The simplest way to install LePHARE is using pip:

```
pip install lephare
```

If you have any problems using pip install please consider creating an issue and informing us what system you are using in order to help us improve robustness.

## Ancillary data files and environment variables
The legacy code LePHARE has been split into two repositorie: (lephare)(https://github.com/lephare-photoz/lephare?tab=readme-ov-file)
and (lephare-data)(https://github.com/lephare-photoz/lephare-data). The module `data_retrieval` helps the user with downloading the latter, based on
the setting of the environment variable `LEPHAREDIR`.
A second environment variable, `LEPHAREWORK`, controls where intermediary files are saved to and retrieved from, during LePHARE execution.

## Example usage

We provide a number of [Jupyter notebooks](docs/notebooks/) demonstrating various aspects of the Python code.

## Dev Guide - Getting Started

Before installing any dependencies or writing code, it's a great idea to create a
virtual environment. LINCC-Frameworks engineers primarily use `conda` to manage virtual
environments. If you have conda installed locally, you can run the following to
create and activate a new environment.

```
conda create -n <env_name>
conda activate <env_name>
```

Once you have created a new environment, you can install this project for local
development using the following commands:

```
git clone git@github.com:lephare-photoz/lephare.git
git submodule update --init --recursive
conda install -c conda-forge cxx-compiler #needed for MACOSX mostly
pip install -e .'[dev]'
pre-commit install
conda install pandoc
```
The installation can then be tested by running `python -m pytest tests` from the root directory of the repository.

Developers can also build the documentation in the following way:

```
cd docs/
pip install -r requirements.txt #install sphinx dependencies
make html
```

The doc entry will then be located at `../_readthedocs/html/index.html`. The documentation includes a rendering of
the notebooks, which thus need to be executed. You can bypass this stage by replacing `make html` above by `make no-notebooks`.
Executing `make` will list further options.

Notes:

1. The single quotes around `'[dev]'` may not be required for your operating system.
2. `pre-commit install` will initialize pre-commit for this local repository, so
   that a set of tests will be run prior to completing a local commit. For more
   information, see the Python Project Template documentation on 
   [pre-commit](https://lincc-ppt.readthedocs.io/en/latest/practices/precommit.html)
3. Install `pandoc` allows you to verify that automatic rendering of Jupyter notebooks
   into documentation for ReadTheDocs works as expected. For more information, see
   the Python Project Template documentation on
   [Sphinx and Python Notebooks](https://lincc-ppt.readthedocs.io/en/latest/practices/sphinx.html#python-notebooks)


This project was automatically generated using the LINCC-Frameworks 
[python-project-template](https://github.com/lincc-frameworks/python-project-template).
For more information about the project template see the 
[documentation](https://lincc-ppt.readthedocs.io/en/latest/).
