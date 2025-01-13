"""
Script Information
-------------------
Writer:      Tillmann Stralka
Date:        2020-05-05
Owner:       HLP University of Leipzig
Python Ver.: 2.7

Description:
This script automates the processing of AFM files, including finding, listing, 
applying fit functions, generating images, extracting statistics, and plotting. 
It also correlates current and topographic information and convolves dynamic scans.

"""

# -------------------------------
# System and Gwyddion-specific Configuration
# -------------------------------
# Some Python libraries have to be added via PATH since they are near Gwyddion
import sys
# Adding paths for Gwyddion and Python 2.7 specific modules
sys.path.append('/usr/local/opt/python@2/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages')
sys.path.append('/usr/share/gwyddion/pygwy')
# Importing pygtk for GUI support (requires GTK-2.0)
import pygtk
pygtk.require20()  # Ensures GTK-2.0 compatibility
# Adding Gwyddion module paths
sys.path.append('/usr/local/Cellar/gwyddion/2.52/share/gwyddion/pygwy')
# Importing Gwyddion utilities and main modules
import gwyutils
import gwy

# -------------------------------
# General Python Standard Libraries
# -------------------------------
# System and OS-related utilities
import os
import time
import shutil

# -------------------------------
# Numerical Computation and Data Analysis Libraries
# -------------------------------
# Core numerical library
import numpy as np
# Data manipulation and analysis
import pandas as pd

# -------------------------------
# Scientific Computing and Signal Processing
# -------------------------------
# Statistical distributions
from scipy.stats import norm, halfnorm
# Curve fitting
from scipy.optimize import curve_fit
# Signal processing utilities
from scipy.signal import find_peaks
import scipy.signal  # Importing the full module for extended utilities

# -------------------------------
# Visualization Libraries
# -------------------------------
# Comprehensive plotting library
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl

# -------------------------------
# Image Processing and File Handling
# -------------------------------
# Image I/O and video handling
import imageio

# -------------------------------
# PDF Handling Libraries
# -------------------------------
# Exception handling for PDF-to-image conversions
from pdf2image.exceptions import (
    PDFInfoNotInstalledError,
    PDFPageCountError,
    PDFSyntaxError,
)
# PDF-to-image conversion utilities
from pdf2image import convert_from_path, convert_from_bytes

# -------------------------------
# # Local Imports
# -------------------------------
# Import utils_afm from the src directory
src_path = os.path.abspath(os.path.join('src'))
print("src path:", src_path)
# Add the src directory to the system path
sys.path.append(src_path)
from utils_afm import *
# Define the path to the 'data' folder in the local repository
data_path = os.path.abspath(os.path.join(os.getcwd(), 'data'))
# Define the general working directory / path 
path = os.path.abspath(os.path.join(os.getcwd()))


###############################################################################
#                               Declare Empty Lines and DataFrames             #
###############################################################################

# Initialize empty lists to store statistical values and data
statistics_current = []  # For storing current statistical values
statistics_topo = []     # For storing topological statistical values
statistics_error = []    # For storing error-related statistical values
statistics_current_gb = []  # For storing current statistical values grouped by something (e.g., geography)
statistics_current_grain = []  # For storing current statistical values by grain

# Initialize empty lists for data frames to be populated later
dataframes_current = []   # For storing current data frames
dataframes_topo = []      # For storing topological data frames
dataframes_error = []     # For storing error data frames

# Empty lists for data lines (e.g., raw input lines or processed data)
datalines_current = []    # For storing lines of current data
datalines_topo = []       # For storing lines of topological data
datalines_distance = []   # For storing lines of distance data

###############################################################################
#                              Convolution Drift Initialization             #
###############################################################################

# Drift initialization with zeros
array_old = 0
array_old2 = 0

# Offsets used for drift correction (initialize as tuples of zeros)
offset_by_drift = (0, 0)
offset_by_drift2 = (0, 0)

###############################################################################
#                             Line Extraction Settings                       #
###############################################################################

# Define starting and ending positions for the line extraction process
line_x_start = 46  # X-coordinate of the line start
line_y_start = 47  # Y-coordinate of the line start
line_x_end = 58    # X-coordinate of the line end
line_y_end = 34    # Y-coordinate of the line end

# Resolution for line extraction (how finely data is sampled)
line_res = 1000

# Initial peak position (set to 0 initially, will be overwritten later)
peak_x_position = 0

###############################################################################
#                           Presets: File Sorting, Info Retrieval     #
###############################################################################

# Sort and list files based on the working path
N, fnames_topo, fnames_current, fnames_amp, fnames_phase, fnames_error = sortandlist(path)

# Retrieve information from 'Infos.csv' file (e.g., voltage list and cutoff lines)
voltage_list, lines_cutoff = get_info_sheet(path, name='Infos.csv')

# Create necessary folders for output files (e.g., PDFs, statistics, images)
path_pdfs, path_histo, path_statistics, path_gifs, path_jpgs, path_lines, path_fitted, path_stableframe = make_folders(path)

###############################################################################
#                          Interaction Settings                               #
###############################################################################

# Define interaction modes for different processes
image_mode = gwy.RUN_NONINTERACTIVE  # Set to non-interactive for image handling
fit_mode = gwy.RUN_NONINTERACTIVE    # Set to non-interactive for fitting process
mask_mode = gwy.RUN_NONINTERACTIVE   # Set to non-interactive for masking process

###############################################################################
#               Main Evaluation Run: Data Processing and Analysis             #
###############################################################################

# First run: Process topological data, fit data, and calculate statistics
df_stat_topo, range_topo, dataframes_topo = first_run_topo(voltage_list, lines_cutoff, fnames_topo, factor=9, daten=statistics_topo)

# Second run: Refine topological data, apply drift correction, and store results
df_drift, dataframes_topo, lines_topo = second_run_topo(dataframes_topo, range_topo, path_fitted, array_old, offset_by_drift, factor=9)

# Uncomment and implement third and fourth runs if needed for more processing steps
# third_run_topo(dataframes_topo, range_topo, path_stableframe, df_drift, factor=9)
# df_drift2, dataframes_topo = fourth_run_topo(dataframes_topo, range_topo, path_stableframe, array_old2, offset_by_drift, factor=9)

# First run: Process current data, fit data, and calculate statistics
df_stat_current, range_current, dataframes_current = first_run_current(voltage_list, lines_cutoff, fnames_current, factor=9, daten=statistics_current)

# Second run: Refine current data and apply corrections
dataframes_current = second_run_current(dataframes_current, range_current, path, factor=9)

# Uncomment to perform further processing for current data (e.g., third run)
# third_run_current(dataframes_current, range_current, path_fitted, path_stableframe, df_drift, factor=9)

# Output check (for debugging purposes)
print('CHECKPOINT')

# Optional: Plot multi-line graph for topography and current (parameters adjusted for visualization)
# plot_line_multi(df_line_topo_small, df_line_current_small, df_line_x, name, path_lines, 
#                 plot_min_y1=topo_min + 50, plot_max_y1=topo_max, plot_min_y2=current_min + 5000, 
#                 plot_max_y2=current_max - 5000, xlabel='Distance ($\it{nm}$)', y1label='Topography ($\it{nm}$)', y2label='Current ($\it{nA}$)')

# Generate GIFs from the topological data for visualization
make_videos(dataframes_topo, range_topo, path_stableframe, path_gifs)


