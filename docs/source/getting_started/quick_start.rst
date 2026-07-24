Quick Start
===========

.. note::

   The GPSP Atlas products are distributed as compressed ``.fits.gz`` archives
   to reduce download size and storage requirements. Both compressed
   ``.fits.gz`` files and uncompressed ``.fits`` files can be opened directly
   by the helper modules.

   However, for large atlas products, the use of uncompressed ``.fits`` files
   is recommended. Astropy can efficiently memory-map uncompressed FITS files,
   allowing random access to the multi-dimensional atlas without repeatedly
   decompressing the underlying data. In contrast, access to compressed
   ``.fits.gz`` files generally involves additional decompression overhead,
   which may result in slower I/O performance for large datasets.

   For optimal performance during scientific analyses, users are therefore
   encouraged to extract the downloaded ``.fits.gz`` archives before use.

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
