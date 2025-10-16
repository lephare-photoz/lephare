  | <img src="https://avatars.githubusercontent.com/u/165841626?s=400&u=ff86bd4c19a9d36958cf1b47d84849dbe25c274a&v=4" width="140"/> | LePHARE <br> PHotometric Analysis for Redshift Estimation <br> <br> [![PyPI](https://img.shields.io/pypi/v/lephare?color=blue&logo=pypi&logoColor=white)](https://pypi.org/project/lephare/) [![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/lincc-frameworks/lephare/smoke-test.yml)](https://github.com/lephare-photoz/lephare/actions/workflows/smoke-test.yml) [![Codecov](https://codecov.io/gh/lephare-photoz/lephare/branch/main/graph/badge.svg)](https://codecov.io/gh/lephare-photoz/lephare) [![Read The Docs](https://img.shields.io/readthedocs/lephare)](https://lephare.readthedocs.io/) [![DOI](https://zenodo.org/badge/307380211.svg)](https://zenodo.org/records/14162574)|
|---|---|

LePHARE is a code for estimating galaxy redshifts and physical parameters using template fitting. It is a complete rewrite in C++ of the [Fortran code](https://www.cfht.hawaii.edu/~arnouts/LEPHARE), with python binding and extension using [pybind11](https://github.com/pybind/pybind11).

Please find the code [documentation](https://lephare.readthedocs.io/) here, including how to install it and to get started.

# Requests and help

If you need help with the code, or if you have feature requests, please use the github issue system to let us know.

# Citing LePHARE

If you use LePHARE, please acknowledge it with the following references:

- Arnouts, S.; Cristiani, S.; Moscardini, L., Matarrese, S., Lucchin, F.  et al., 1999, MNRAS,  310, 540

- Ilbert, O.; Arnouts, S.; McCracken, H. J.; Bolzonella, M.; Bertin, E et al., 2006, A&A, 457, 841

# Contributors

LePHARE was originally developped in Fortran by [Stéphane Arnouts](https://people.lam.fr/arnouts.stephane/) and [Olivier Ilbert](https://people.lam.fr/ilbert.olivier/).

The C++ and python rewriting of the code is the work of Olivier Ilbert, [Johann Cohen-Tanugi](https://github.com/johannct), and [Raphael Shirley](http://raphaelshirley.co.uk/).

Other contributors include:
Iary Davidzon, Mara Salvato (MPE), Cédric Dubois (LAM), and Maria Petkova.

We acknowledge fruitful discussions with
Emeric Le Floc'h (CEA), Léo Michel-Dansac (LAM), Jean-Charles Lambert (LAM).


# Acknowledgements
[![Template](https://img.shields.io/badge/Template-LINCC%20Frameworks%20Python%20Project%20Template-brightgreen)](https://lincc-ppt.readthedocs.io/en/latest/)

This project was automatically generated using the LINCC-Frameworks [python-project-template](https://github.com/lincc-frameworks/python-project-template).

The authors gratefully acknowledge the important contribution of LINCC Framework members (notably [Olivia Lynn](https://github.com/OliviaLynn) and [Drew Oldag](https://github.com/drewoldag)) to the
construction of the LePHARE github infratructure.

The LePHARE logo is the work of [Eve Barlier](https://www.instagram.com/eve.barlier/). Thank you!
