Downloads
=========

GPSP Atlas Releases
===================

The complete GPSP Atlas v1.0 release is archived on Zenodo and can be cited using

**DOI:** `10.5281/zenodo.21403524 <https://doi.org/10.5281/zenodo.21403524>`_

The release contains four atlas products and associated helper modules,
comprising approximately 33 GB of compressed data products.

.. note::

   All atlas products are distributed as compressed ``.fits.gz`` archives.
   For optimal performance, users are encouraged to extract the files prior
   to analysis. The helper modules and all examples throughout this
   documentation assume the use of uncompressed ``.fits`` files.

Standard GPSP Atlases
---------------------

The standard GPSP atlas provides Galactic gamma-ray survival probabilities
computed using two independent GALPROP-based interstellar radiation field
(ISRF) models.

============================= ======== ==========
Atlas Product                 Size     Download
============================= ======== ==========
GPSP Atlas R12                11.2 GB  `Zenodo <https://zenodo.org/records/21403524/files/GPSP_Atlas_R12.fits.gz?download=1>`_
GPSP Atlas F98                11.2 GB  `Zenodo <https://zenodo.org/records/21403524/files/GPSP_Atlas_F98.fits.gz?download=1>`_
============================= ======== ==========

The R12 atlas uses the ISRF model of Robitaille et al. (2012), while
the F98 atlas uses the ISRF model of Freudenreich (1998).

GPSP-LIV Atlases
----------------

The GPSP-LIV atlases incorporate the modified subluminal Lorentz-Invariance
Violation (LIV) pair-production cross section introduced by Carmona et al.
(2024) directly into the gamma-ray propagation calculations.

================================= ======== ==========
Atlas Product                     Size     Download
================================= ======== ==========
GPSP-LIV Atlas R12                5.5 GB   `Zenodo <https://zenodo.org/records/21403524/files/GPSP_Atlas_LIV_R12.fits.gz?download=1>`_
GPSP-LIV Atlas F98                5.5 GB   `Zenodo <https://zenodo.org/records/21403524/files/GPSP_Atlas_LIV_F98.fits.gz?download=1>`_
================================= ======== ==========

The LIV parameter space covers

.. math::

   -4.5 \leq \log_{10}(\lambda/E_{\rm Pl}) \leq -2.0

with steps of 0.05.

Python Helper Modules
---------------------

The following helper modules provide a high-level interface for loading,
interpolating, and analyzing GPSP Atlas products.

================================= ======== ==========
Module                            Size     Download
================================= ======== ==========
gpsp_atlas_helper.py              30.6 kB  `Zenodo <https://zenodo.org/records/21403524/files/gpsp_atlas_helper.py?download=1>`_
gpsp_liv_atlas_helper.py          20.9 kB  `Zenodo <https://zenodo.org/records/21403524/files/gpsp_liv_atlas_helper.py?download=1>`_
================================= ======== ==========

Example Datasets & Scripts
--------------------------

Example datasets and scripts used throughout the tutorials can be obtained
from the corresponding documentation pages and GitHub repository. These files
are intended solely for demonstration purposes and are not required for the
use of the full atlas products.
