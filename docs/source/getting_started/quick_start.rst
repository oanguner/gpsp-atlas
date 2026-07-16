Quick Start
===========

Loading a GPSP Atlas
--------------------

.. code-block:: python

   from gpsp_atlas_helper import GPSPAtlasHelper

   atlas = GPSPAtlasHelper(
       "GPSP_Atlas_4D_R12_20bins.fits"
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
