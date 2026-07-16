import os
import sys

# Be sure that the GPSPAtlasHelper module is in your system path!
sys.path.insert(0, os.path.abspath("../../../gpsp-atlas"))

import numpy as np
from astropy.io import fits
from gpsp_atlas_helper import GPSPAtlasHelper
import matplotlib.pyplot as plt

helper = GPSPAtlasHelper("GPSP_Atlas_4D_R12_20bins.fits")

# Desired energy slice in eV. Currently set to 150 TeV
target_energy_ev = 150e12

# Desired distance in kpc. Currently set to 10 kpc
target_distance_kpc = 10.0

# bmin and bmax. Currently set to produce Galactic plane map between -5 deg and 5 deg in latitude
b_lower = -5
b_upper = 5

# Output FITS file. 
output_fits="GPSP_2D_Map_10kpc_150TeV.fits"

map_2d, l_grid, b_grid = helper.get_galactic_plane_map(target_energy_ev=target_energy_ev, target_distance_kpc=target_distance_kpc, b_lower=b_lower, b_upper=b_upper,
                                                       output_fits=output_fits)
                                                       
print("Done...")
