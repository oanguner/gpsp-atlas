# GPSP Atlas

The **GPSP Atlas** is a computational framework designed to map the Galactic Plane Survival Probability (GPSP) for high-energy gamma-ray sources. 

This repository provides helper modules to facilitate the analysis of FITS-based survival probability maps, implementing modified cross-section models—including Lorentz Invariance Violation (LIV) scenarios—to interpret observations of violent astrophysical accelerators.

## Documentation

Full documentation, including API references and tutorial notebooks, is hosted on Read the Docs:

**[https://gpsp-atlas.readthedocs.io](https://gpsp-atlas.readthedocs.io)**

## Features

* **Data Handling:** Utilities to parse and manipulate GPSP atlas FITS files.
* **LIV Modeling:** Explicit calculation of photon absorption cross-sections (e.g., Model A) as detailed in [arXiv:2404.07842].
* **Visualization:** Helper modules for plotting survival probabilities as a function of energy and distance.

## Installation

```bash
pip install .
