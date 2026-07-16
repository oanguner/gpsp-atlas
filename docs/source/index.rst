GPSP Atlas Documentation
========================

**Galactic Photon Survival Probability (GPSP) Atlas
for Very-High-Energy and Ultra-High-Energy
Gamma-Ray Astronomy.**

The GPSP Atlas provides precalculated Galactic gamma-ray survival probabilities
together with Python helper modules for scientific analyses of VHE and UHE
gamma-ray sources.

The atlas describes the Galactic photon survival probability over the full
parameter space

.. math::

   P(l,b,d,E_{\gamma})

and covers gamma-ray energies from 1 TeV to 10 PeV.

Two independent GALPROP-based interstellar radiation field (ISRF) models are
provided:

* R12 (Robitaille et al. 2012) --> `R12 Paper <https://doi.org/10.1051/0004-6361/201219073>`_
* F98 (Freudenreich 1998) --> `F98 Paper <https://doi.org/10.1086/305065>`_

These models are further discussed in Porter et al. 2017 (`See their paper <https://doi.org/10.3847/1538-4357/aa844d>`_).

A dedicated Lorentz-Invariance Violation version of the atlas (GPSP-LIV)
extends the framework by incorporating modified subluminal pair-production
cross sections introduced in `Carmona et al. 2024 <https://link.aps.org/doi/10.1103/PhysRevD.110.063035>`_ 
into the gamma-ray propagation calculations.

GPSP Atlas Features
-------------------

* Gamma-ray Survival probability calculations
* Gamma-ray horizon determination
* Galactic bird’s-eye view transparency maps
* Galactic Plane transparency maps
* Extended source analyses
* Application of photon absorption to intrinsic spectrum 
* Correcting observed spectrum for absorption effects
* Subluminal Lorentz-Invariance Violation studies

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   getting_started/installation
   getting_started/atlas_structure
   getting_started/quick_start

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/survival_probabilities
   tutorials/galactic_transparancy_horizons
   tutorials/gc_pevatron_region
   tutorials/cygnus_x3_region

.. toctree::
   :maxdepth: 2
   :caption: Example Datasets & Scripts

   example_datasets
   example_scripts

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/gpsp_atlas_helper
   api/gpsp_liv_atlas_helper

.. toctree::
   :maxdepth: 1
   :caption: Additional Information

   downloads
   citations
