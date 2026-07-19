Quick Start
===========

.. note::

   The examples presented throughout this documentation assume that the GPSP
   Atlas FITS files have been extracted from the distributed ``.fits.gz``
   archives prior to use. While the atlas products are distributed in compressed
   form to reduce download size and storage requirements, the helper modules are
   designed to operate on the uncompressed ``.fits`` files.

   Using uncompressed FITS files allows direct memory-mapped access to the data
   through Astropy, enabling efficient random access to the multi-dimensional
   atlas without loading the entire dataset into memory. In contrast, compressed
   ``.fits.gz`` files must first be decompressed during access, which can
   significantly increase I/O overhead and may require substantially more
   temporary memory resources for large atlas products.

   For optimal performance, users are therefore encouraged to extract the
   downloaded ``.fits.gz`` archives before performing analyses.

Loading a GPSP Atlas
--------------------

.. code-block:: python

   from gpsp_atlas_helper import GPSPAtlasHelper

   atlas = GPSPAtlasHelper(
       "GPSP_Atlas_R12.fits"
   )

Computing a Survival Profile
----------------------------

.. code-block:: python

   energies, p_surv = atlas.get_survival_vs_energy(
       l=0.0,
       b=0.0,
       d=8.5
   )

Computing a Gamma-Ray Horizon
-----------------------------

.. code-block:: python

   horizons = atlas.get_horizon_distances(
       l=0.0,
       b=0.0,
       energies_ev=[1e15],
       target_p=0.37
   )

See the Tutorials section for more advanced examples.
