import os
import sys

# Be sure that the GPSPAtlasHelper is in the system path!
sys.path.insert(0, os.path.abspath("../../../gpsp-atlas"))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.ndimage import gaussian_filter
from gpsp_atlas_helper import GPSPAtlasHelper

plt.rcParams['font.weight'] = 'bold'
plt.rcParams['axes.labelweight'] = 'bold'
plt.rcParams['axes.titleweight'] = 'bold'
plt.rcParams['xtick.labelsize'] = 8  # Changes the 90, 180, 270 labels

# Force X11/Xcb for stable rendering
os.environ["QT_QPA_PLATFORM"] = "xcb"

# --- 1. CONFIGURATION ---
FILE_PATH = "GPSP_Atlas_4D_R12_20bins.fits"
ENERGIES_TEV = [100, 150, 200, 250]
B_LOWER, B_UPPER = -1.0, 1.0
vmin = 0.7
vmax = 1.0
# Contour Settings
CONTOUR_LEVELS = [0.75, 0.80, 0.85]
CONTOUR_COLORS = ['white']
CONTOUR_LINESTYLES = ['dotted', 'dashed', 'solid']

# Tangential Arms
ARMS = {
    "Norma-3 kpc": [329, 3],
    "Scutum–Crux": [27, 3],
}

# --- 2. DATA EXTRACTION & PLOTTING ---
helper = GPSPAtlasHelper(FILE_PATH)
num_maps = len(ENERGIES_TEV)

fig, axes = plt.subplots(1, num_maps, subplot_kw={'projection': 'polar'}, 
                         figsize=(18, 5), constrained_layout=True)

for i, energy_tev in enumerate(ENERGIES_TEV):
    ax = axes[i]
    target_energy_ev = energy_tev * 1e12
    
    # 1. Extraction
    l_grid, d_axis, map_raw = helper.get_bird_view_map(
        target_energy_ev=target_energy_ev,
        b_lower=B_LOWER,
        b_upper=B_UPPER,
        method='percentile',
        percentile=50
    )

    L_mesh, D_mesh = np.meshgrid(np.radians(l_grid), d_axis, indexing='ij')
    map_smooth = gaussian_filter(map_raw, sigma=(12.0, 0.0), mode='wrap')

    # 2. Polar Setup
    ax.set_facecolor('black')
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(1)
    ax.set_rasterization_zorder(1)

    # Base Map
    pc = ax.pcolormesh(L_mesh, D_mesh, map_raw, cmap='gist_heat', 
                       norm=LogNorm(vmin=vmin, vmax=vmax), shading='auto', zorder=0, edgecolor='face')

    # --- 3. CONTOURS WITH BLACK SHADOWS ---
    # Solid black "trench" underneath for contrast
    ax.contour(L_mesh, D_mesh, map_smooth, levels=CONTOUR_LEVELS, 
               colors='black', linewidths=2.2, zorder=0.4, linestyles='solid')

    # Main white contours
    cs = ax.contour(L_mesh, D_mesh, map_smooth, levels=CONTOUR_LEVELS, 
                    colors=CONTOUR_COLORS, linestyles=CONTOUR_LINESTYLES, 
                    linewidths=1.4, zorder=0.5)

    # 4. Arm Tangents
    for name, (l_mean, l_rms) in ARMS.items():
        theta_mid = np.radians(l_mean)
        
        ax.vlines(theta_mid, 0, 20, color='black', linestyle='-', 
                      linewidth=1.8, alpha=0.7, zorder=4.9)
                      
        ax.vlines(theta_mid, 0, 20, color='lime', linestyle='--', 
                      linewidth=1.2, alpha=0.7, zorder=5.0)

    # 5. Grid & Labels Formatting
    ax.set_xticks(np.radians([0, 90, 180 ,270]))
    #ax.set_xticklabels(['$0^\circ$', '$90^\circ$', '$180^\circ$', '$270^\circ$'], 
    #                   color='none', fontsize=8, fontweight='bold')
                       
    labels = ['$0^\circ$', '', '$180^\circ$', '']
    if i == 0: # Far left map
        labels = ['$0^\circ$', '$90^\circ$', '$180^\circ$', ''] 
    if i == num_maps - 1: # Far right map
        labels = ['$0^\circ$', '', '$180^\circ$', '$270^\circ$']
        
    ax.set_xticklabels(labels, fontweight='bold', fontsize=11, color='none')
    
    ax.set_rlim(0, 20)
    ax.set_yticks([5, 10, 15, 20])
    ax.set_yticklabels([]) 
    ax.grid(True, color='black', linestyle='--', alpha=0.3)

    # Markers
    ax.scatter(0, 0, marker='*', s=80, color='lime', zorder=5, edgecolor='black', linewidths=0.5)
    ax.scatter(0, 8.5, marker='o', s=30, color='cyan', zorder=5, edgecolor='black', linewidths=0.5)
    ax.scatter(1.393, 7.4, marker='*', s=80, color='purple', zorder=5, edgecolor='black', linewidths=0.5)
    
    #ax.set_title(f"{energy_tev} TeV", fontsize=15, pad=25, fontweight='bold')

# --- 6. COLORBAR BOLDING ---
cbar = fig.colorbar(pc, ax=axes, location='right', pad=0.02, shrink=0.7, aspect=25)
cbar.set_label('Survival Probability ($P_{surv}$)', fontsize=14, fontweight='bold')

# 3. Optional: If you want exactly "0.7, 0.8..." and they aren't showing up
# because it's a log scale, manually set the ticks
cbar.set_ticks([0.7, 0.8, 0.9, 1.0])
cbar.ax.set_yticklabels(['0.7', '0.8', '0.9', '1.0'], fontweight='bold', fontsize=12)

plt.savefig("bird_view_map.png", dpi=300)

plt.show()
