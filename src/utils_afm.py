"""Legacy AFM pipelines and compatibility exports.

Independent I/O, processing, analysis, and visualization helpers live in their
respective modules and are re-exported here for existing callers. Pipeline
functions that rely on this module's global state remain in this module.
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

# Keep the legacy utils_afm import surface while implementations live in afm_io.
from afm_io import (
    SaveStatisticsToFile,
    data_save,
    extract_time,
    extract_voltage,
    get_info_sheet,
    get_time,
    load_data,
    load_statistics_data,
    make_drift_list,
    make_folders,
    remove,
    sortandlist,
)
from afm_processing import (
    area_extract,
    average_check,
    cut_edges,
    fit_functions,
    get_line,
    print_df_NumpyArray,
    select_dataframe,
    zero_check,
)
from afm_analysis import (
    Gauss,
    get_drift,
    get_maximum_drift,
    get_offset_from_first_image,
    get_statistic,
    peak_extraction,
)
from afm_visualization import (
    get_colormap,
    image_save,
    make_histogram,
    make_histogram_all,
    make_statistics_plot,
    make_statistics_plot_multi,
    select_color,
    select_range,
)

# -------------------------------
# Path Handling
# -------------------------------
# Define the path to the 'data' folder in the local repository
data_path = os.path.abspath(os.path.join(os.getcwd(), 'data'))



###############################################################################
#   Global Variables    
###############################################################################
# Statistics lists for storing measurement values
statistics_current = []
statistics_topo = []
statistics_error = []
statistics_current_gb = []
statistics_current_grain = []

# Data frame lists for storing processed data
dataframes_current = []
dataframes_topo = []
dataframes_error = []

# Line scan data lists
datalines_current = []
datalines_topo = []
datalines_distance = []

# Drift tracking variables
array_old = 0  # Previous array for drift calculation
array_old2 = 0  # Secondary previous array for drift calculation
offset_by_drift = (0, 0)  # Cumulative drift offset
offset_by_drift2 = (0, 0)  # Secondary drift offset


##################################################################################################
# Legacy orchestration and helpers that depend on utils_afm module globals
##################################################################################################

def make_gif(path, time, name):
    #Making a video of the evaluated scans with opencv 
    #Change to image folder for getting pics and saving video
    os.chdir(path) 
    #Collect Images and put them in list
    #print("The unsorted list:")
    filenames = os.listdir(path)
    #print(filenames)
    jpgs_list = [k for k in filenames if name in k]
    #print(jpgs_list)
    # Sort them with key definition for voltage extraction and order them with int Voltage
    jpgs_list = sorted(jpgs_list,key=extract_time)
    #print("The sorted list:")
    #print(jpgs_list    )
    #jpgs_list.reverse()
    #print(jpgs_list)
    cwd = os.getcwd()
    #print(cwd)
    if cwd == path_pdfs:
        #print("We are making AFM scan Gifs")
        resolution = 1000
    else:
        #print("We are making histograms or linescan Gifs")
        resolution = 100
        
    print("Making Video from image list:"    )
    with imageio.get_writer(name + ".gif", mode="I", duration = time) as writer:
        for image in jpgs_list:
            #print("We are in the image loop, image type:")
            #print(type(image))
            #print("We are in the image loop, image name:")
            #print(image)
            image = make_pdf_to_jpg(image, resolution)
            #print("Image after conversion:")
            #print(image)
            #Change path into jpg folder where the converdet jpg filed (made from the pdf are stored)
            os.chdir(path_jpgs)
            img = imageio.imread(image)
            #print("Type with the image reading function")
            #print(type(img))
            os.chdir(path)
            writer.append_data(img)
    writer.close()            
    return

def make_pdf_to_jpg(image, resolution):
    #A little function to convert pdfs to pixel images, in that case jpgs to make gifs or
    #Videos, since the video makers can not handle vector graphics 
    #The function returns a string for a file name, where .pdf is replaced by .jpg 
    #print("We are in the image converter function, type and name:"
    name_long = image
    name = name_long[:-4]
    #print(type(name)
    #print(name
    images = convert_from_path(image, resolution)
    for image in images:
        os.chdir(path_jpgs)
        image.save(name + ".jpg", "JPEG") 
    #print("Conversion worked!"    
    return (name + ".jpg")

def gauss_fit(x_data_peak, y_data, x_data, range_idx):
    """
    Fit a Gaussian curve to data around a peak position.
    
    Args:
        x_data_peak (float): X position of peak to fit around
        y_data (array): Full Y-axis data
        x_data (array): Full X-axis data
        range_idx (int): Number of points to include on each side of peak
        
    Returns:
        None: Currently just plots and saves the fit
        
    Notes:
        - Extracts window of data centered on peak position
        - Calculates initial guess parameters from data
        - Fits Gaussian using scipy.optimize.curve_fit
        - Plots and saves the fit results
    """
    # Convert data to numpy arrays
    y_all = np.asarray(y_data)
    x_all = np.asarray(x_data)
    
    # Find array index closest to peak position
    peak_idx = (np.abs(x_all - x_data_peak)).argmin()
    
    # Extract window of data around peak
    y = y_all[(peak_idx-range_idx):(peak_idx+range_idx)]
    x = x_all[(peak_idx-range_idx):(peak_idx+range_idx)]

    # Calculate initial guess parameters
    mean = np.sum(x * y) * (1. / np.sum(y))
    sigma = np.sqrt(np.sum(y * (x - mean)**2) / np.sum(y))
    
    # Fit Gaussian curve
    popt, pcov = curve_fit(Gauss, x, y, p0=[np.max(y), mean, sigma])
    perr = np.sqrt(np.diag(pcov))
    
    # Extract fit parameters
    gauss_peak_y = popt[0]  # Amplitude
    gauss_peak_x = popt[1]  # Center position
    gauss_sigma = popt[2]   # Standard deviation
    
    # Plot results
    fig, ax = plt.subplots(figsize=(10,5))
    plt.plot(x, y, "b+:", label="data")
    plt.plot(x, Gauss(x, *popt), "r-", label="fit")
    plt.legend()
    plt.title("Gaussian Fit to Peak")
    plt.xlabel("Position (nm)")
    plt.ylabel("Signal")
    ax.set_ylim(-3000, 3000)
    
    # Save plot
    os.chdir(path_lines)
    plt.savefig("{}_Gaussfit.pdf".format(name), format="pdf")
    os.chdir(path)
    plt.show()

    return  "The fit parameters for the gauss function: popt, pcov" , popt ,pcov

def first_run_topo(voltage_list, lines_cutoff, fnames_topo, daten, factor=9):
    """
    Performs the first run of the topography pipeline.

    Steps:
    1. Extract data from ActCon files and preprocess it.
    2. Cut edges of the data and apply color selection.
    3. Fit functions and calculate statistics.
    4. Save processed DataFrames and compute an overview evaluation.

    Args:
        voltage_list (list): List of voltages corresponding to each file.
        lines_cutoff (int): Number of lines to cut off from the data.
        fnames_topo (list): List of file names to process.
        factor (int): Scaling factor for data adjustment.
        daten (list): Metadata associated with the dataset.

    Returns:
        tuple: Processed statistical DataFrame, range of topography, and processed DataFrames.
    """
    for i, fname in enumerate(fnames_topo):
        print("First Run Iteration: i={}".format(i))
        name = "{} V_Topo".format(voltage_list[i])
        print(name)

        ids, actcon = load_data(path, fname)
        df = select_dataframe(actcon, parameter_name=name, key=0)
        df = cut_edges(df, lines_cutoff)
        actcon = select_color(actcon, color="Gwyddion.net", key=0)
        actcon = fit_functions(actcon, df, i, mode=fit_mode)
        statistics, max_val, min_val, Sme, zwischendaten, df = get_statistic(
            df, actcon, fname, name, factor, daten, i, noise=0, key=0
        )
        df = df_save(df, actcon, dataframes_topo, fname, name, voltage_list[i])
        df = average_check(df, actcon, key=0)
        actcon = data_save(actcon, fname, path_saving=path_fitted, path_working=path)
        remove(actcon)

    # Perform overview evaluation and save statistical data
    SaveStatisticsToFile(daten, path_statistics, voltage_list, name="Statistics_Topo.csv", unit="nm")
    df_stat_topo = load_statistics_data(path_statistics, name="Statistics_Topo.csv")
    range_topo = get_opimum_range_topo(df_stat_topo, column_name="Maximum [nm]")

    # Create histograms for topography data
    make_histogram_all(
        N,
        dataframes_topo,
        path_histo,
        path,
        factor,
        plot_max=range_topo,
        plot_min=-range_topo,
        name_ylabel="Topography ($\it{nm}$)",
        name_file="Topo_Histogramms.pdf",
        list_names=voltage_list,
    )
    return df_stat_topo, range_topo, dataframes_topo

def second_run_topo(dataframes_topo, range_topo, path_fitted, array_old, offset_by_drift, factor = 9):
    """
    Performs the second run of the topography pipeline.

    Steps:
    1. Reprocess DataFrames with drift extraction.
    2. Generate and save images, videos, and histograms.
    3. Convolve the data to extract drift.

    Args:
        dataframes_topo (list): List of processed DataFrames from the first run.
        range_topo (float): Topography range for processing.
        path_fitted (str): Path to the fitted data.
        array_old (numpy.array): Initial array for drift extraction.
        offset_by_drift (tuple): Offset values to adjust drift.
        factor (int): Scaling factor for data adjustment.

    Returns:
        tuple: Updated drift DataFrame, processed DataFrames, and line scan data.
    """
    drift_list = []

    for i, data in enumerate(dataframes_topo):
        print("Second Run Topo Iteration: i={}".format(i))
        name, voltage, df, actcon, fname = data
        ids, actcon = load_data(path_fitted, fname + ".gwy")

        df = average_check(df, actcon, key=0)
        actcon = image_save(actcon, i, path_pdfs, path, mode=image_mode, dataname=fname)

        df_line_topo, df_line_x, step_size = get_line(
            df, factor, x_start=line_x_start, y_start=line_y_start, 
            x_end=line_x_end, y_end=line_y_end, res=line_res
        )
        plot_line(
            df_line_topo, df_line_x, name, path_lines, 
            plot_min=-range_topo, plot_max=range_topo, 
            xlabel="Distance ($\it{nm}$)", ylabel="Topo ($\it{nm}$)"
        )
        append_to_linescan_list(df_line_topo, name, daten=datalines_topo)

        array_old, drift = get_drift(df, array_old, i)
        offset_by_drift = get_offset_from_first_image(drift, offset_by_drift, i)
        drift_list = make_drift_list(
            drift, offset_by_drift, i, drift_list, path, path_statistics, 
            dataframes=dataframes_topo, name="Drift_List.csv"
        )
        df = df_save(df, actcon, dataframes_topo, fname, name, voltage_list[i])
        remove(actcon)

    df_drift = load_statistics_data(path_statistics, name="Drift_List.csv")
    df_stat_topo = load_statistics_data(path_statistics, name="Statistics_Topo.csv")
    make_statistics_plot_multi(
        df_stat_topo, 1, 4, 6, 5, path_statistics, 
        plot_max=(2 * range_topo), plot_min=0, 
        label="Height ($\it{nm}$)", name="Topo_MultiPlot"
    )
    make_statistics_plot(
        df_stat_topo, path_statistics, 13, [2, 3, 4, 5, 6], label="Topo"
    )
    return df_drift, dataframes_topo[:i + 1], datalines_topo

def third_run_topo(dataframes_topo, Range, path_fitted, path_stableframe, df_drift, factor = 9):
    #Here we extract from the fittted files a stable frame with the drift list
    maximum_drift_x, maximum_drift_y = get_maximum_drift(df_drift, column_name1="X Offset by drift from first image [px]", column_name2="Y Offset by drift from first image [px]")
    for i in range(len(dataframes_topo)):
        print("Third Run Topo Iteration: i=")
        print(i )
        name, voltage, df, actcon, fname = dataframes_topo[i]
        ids, actcon = load_data(path_fitted, fname + ".gwy")
 
        df_stable, key = area_extract(df, actcon, name, df_drift, i, maximum_drift_x, maximum_drift_y, lines_cutoff)
       
        df_stable= average_check(df_stable, actcon, key)
        actcon = select_color(actcon, color="Gwyddion.net", key=key)
        actcon = select_range(actcon, Range*(10**-factor), key)
        actcon = image_save(actcon, i, path_stableframe, path, mode = image_mode, dataname= fname)
        actcon = data_save(actcon, fname, path_saving = path_stableframe, path_working = path)
        #df_stable = df_save(df_stable, actcon, dataframes_topo, fname, name, voltage_list[i]) 
        remove(actcon)
        print(len(dataframes_topo))
    return

def fourth_run_topo(dataframes_topo, range_topo, path_stableframe, array_old2, offset_by_drift, factor):
    #Foruth run, make convolution again, since it often is inacurate at first
    drift_list2 = []
    for i in range(len(dataframes_topo)):
        print("Fourth Run Topo Iteration: i=")
        print(i )
        name, voltage, df, actcon, fname = dataframes_topo[i]
        ids, actcon = load_data(path_stableframe, fname + ".gwy")
        df = select_dataframe(actcon, parameter_name=voltage, key = 0)

        df = average_check(df, actcon, key = 0)
        #actcon = select_range(actcon, range_topo*(10**-factor), key=0)
        actcon = image_save(actcon, i, path_pdfs, path, mode = image_mode, dataname= fname)

        #make_histogram(df, name, fname, path_histo, path, factor, plot_max = range_topo, plot_min = -range_topo, label="Topo distribution", xlabel="Topography ($\it{nm}$)")
        ##########  Convolve    ##########
        array_old2, drift = get_drift(df, array_old2, i)
        offset_by_drift = get_offset_from_first_image(drift, offset_by_drift, i)
        drift_list2 = make_drift_list(drift,offset_by_drift,i,drift_list2, path, path_statistics, dataframes = dataframes_topo, name="Drift_List2.csv")
        df = df_save(df, actcon, dataframes_topo, fname, name, voltage_list[i]) 
        remove(actcon)

    df_drift2 = load_statistics_data(path_statistics, name = "Drift_List2.csv")    
    dataframes_topo = dataframes_topo[:(i+1)]
    
    #cut out for another convolution run (basically what we did in thrid topo run)
    maximum_drift_x, maximum_drift_y = get_maximum_drift(df_drift2, column_name1="X Offset by drift from first image [px]", column_name2="Y Offset by drift from first image [px]")
    for i in range(len(dataframes_topo)):
        print("Fourth2 Run Topo Iteration: i=")
        print(i )
        name, voltage, df, actcon, fname = dataframes_topo[i]
        ids, actcon = load_data(path_stableframe, fname + ".gwy")
        df = select_dataframe(actcon, parameter_name=voltage, key = 0)
        print(df)
        
        df_stable, key = area_extract(df, actcon, name, df_drift2, i, maximum_drift_x, maximum_drift_y, lines_cutoff)
       
        #df_stable= average_check(df_stable, actcon, key)
        #actcon = select_color(actcon, color="Gwyddion.net", key=key)
        #actcon = select_range(actcon, Range*(10**-factor), key)
        #actcon = image_save(actcon, i, path_stableframe, path, mode = image_mode, dataname= fname)
        actcon = data_save(actcon, fname, path_saving = path_stableframe, path_working = path)
        #df_stable = df_save(df_stable, actcon, dataframes_topo, fname, name, voltage_list[i]) 
        remove(actcon)
        print(len(dataframes_topo))


    #iteration loop again, for even better convolution of images 
    drift_list3 = []
    array_old3 = 0
    offset_by_drift3 = (0,0)
    for i in range(len(dataframes_topo)):
        print("Fourth Run Topo Iteration: i=")
        print(i )
        name, voltage, df, actcon, fname = dataframes_topo[i]
        ids, actcon = load_data(path_stableframe, fname + ".gwy")
        df = select_dataframe(actcon, parameter_name=voltage, key = 0)

        df = average_check(df, actcon, key = 0)
        #actcon = select_range(actcon, range_topo*(10**-factor), key=1)
        actcon = image_save(actcon, i, path_pdfs, path, mode = image_mode, dataname= fname)

        #make_histogram(df, name, fname, path_histo, path, factor, plot_max = range_topo, plot_min = -range_topo, label="Topo distribution", xlabel="Topography ($\it{nm}$)")
        ##########  Convolve    ##########
        array_old3, drift = get_drift(df, array_old3, i)
        offset_by_drift3 = get_offset_from_first_image(drift, offset_by_drift3, i)
        drift_list3 = make_drift_list(drift,offset_by_drift3,i,drift_list3, path, path_statistics, dataframes = dataframes_topo, name="Drift_List3.csv")
        df = df_save(df, actcon, dataframes_topo, fname, name, voltage_list[i]) 
        remove(actcon)

    df_drift3 = load_statistics_data(path_statistics, name = "Drift_List3.csv")    
    dataframes_topo = dataframes_topo[:(i+1)]    
    return df_drift, dataframes_topo[:(i+1)]

def first_run_current(voltage_list, lines_cutoff, fnames_current, factor, daten):
    #First run for current files makes: Extraction of DF from Actcon, Edges, Color, ZeroCheck, Statistics, Save after Fitting  
    for i in range(len(fnames_current)):
        print("Iteration: i=")
        print(i )
        fname = fnames_current[i]
        name =  str(voltage_list[i]) + " V" + "_Current"
        ids, actcon = load_data(path, fname)
        df = select_dataframe(actcon, parameter_name=voltage_list[i], key = 0)
        df = cut_edges(df, lines_cutoff)
        #df, offset = zero_check(df, actcon, key = 0)
        actcon = select_color(actcon, color="Red-Cyan", key = 0)  
        statistics, max, min, Sme, zwischendaten, df =  get_statistic(df, actcon, fname, name, factor, daten, i, noise = 0, key=0)        
        df = df_save(df, actcon, dataframes_current, fname, name, voltage_list[i])    
        remove(actcon)
    #After first run overview evaluation and collection of values for image making and so on            
    SaveStatisticsToFile(daten, path_statistics, voltage_list, name="Statistics_Current.csv", unit="nA")
    df_stat_current = load_statistics_data(path_statistics, name = "Statistics_Current.csv")    
    range_current = get_opimum_range_current(df_stat_current, column_name1 = "Maximum [nA]", column_name2 = "Minimum [nA]")        
    #the df lists consit of name and gwy objects 
    return df_stat_current, range_current, dataframes_current

def second_run_current(dataframes, Range, path, factor):
    #Second run makes: average on 0, make and save: images-videos-histograms, convolution to get drift  
    #Set array_old for beginning, since its used by drift extraction function, will be overwritten in loop and reused
    print(type(dataframes))
    print(dataframes)
    for i in range(len(dataframes)):
        print("Second Run Current Iteration: i=")
        print(i )
        name, voltage, df, actcon, fname = dataframes[i]
        ids, actcon = load_data(path, fname)
        df = cut_edges(df, lines_cutoff)
        df, offset = zero_check(df, actcon, key=0)
        actcon = select_color(actcon, color="Red-Cyan", key = 0)
        df = select_dataframe(actcon, parameter_name=voltage_list[i], key = 0)
        actcon = select_range(actcon, Range*(10**-factor), key=0)
        actcon = image_save(actcon, i, path_pdfs, path, mode = image_mode, dataname= fname)
        
        df_line, df_line_x, step_size = get_line(df, factor, x_start=line_x_start, y_start=line_y_start, x_end=line_x_end, y_end=line_y_end, res=line_res) 
        save_line(df_line, df_line_x, path_lines, path, name, parameter=voltage, name_x = "Topography", name_y = "Current")
        
        append_to_linescan_list(df_line, name, daten=datalines_current)
        #peak_x_position, x_value_Peak = peak_extraction(i, df_line_x, df_line_current, step_size, peak_x_position, startingpoint = 0.75, width = 200)
        #gauss_fit(x_value_Peak, df_line_current, df_line_x, range_idx=150)
        plot_line(df_line, df_line_x, name, path_lines, plot_min=-Range, plot_max=Range, xlabel="Distance ($\it{nm}$)", ylabel="Current ($\it{unit_current}$)")      

        
        #make_histogram(df, name, fname, path_histo, path, factor, plot_max = Range, plot_min = -Range, label="Current distribution", xlabel="Current ($\it{nA}$)")
        actcon = data_save(actcon, fname, path_saving = path_fitted, path_working = path)
        remove(actcon)

    make_histogram_all(N, dataframes, path_histo, path, factor, plot_max = Range, plot_min = -Range, name_ylabel="Current ($\it{nA}$)", name_file="Current_Histogramms.pdf", list_names = voltage_list)
    df_stat_current = load_statistics_data(path_statistics, name = "Statistics_Current.csv")
    make_statistics_plot_multi(df_stat_current, 1, 4 ,6 ,5 , path_statistics, plot_max = Range, plot_min = -Range, label="Current ($\it{nA}$)", name = "Current_MultiPlot")
    make_statistics_plot(df_stat_current, path_statistics, 13, [2, 3, 4,5,6], label="Current")

    return dataframes

def third_run_current(dataframes, Range, path_fitted, path_stableframe, df_drift, factor = 9):
    #Here we extract from the fittted files a stable frame with the drift list
    maximum_drift_x, maximum_drift_y = get_maximum_drift(df_drift, column_name1="X Offset by drift from first image [px]", column_name2="Y Offset by drift from first image [px]")
    for i in range(len(dataframes)):
        print("Third Run Current Iteration: i=")
        print(i )
        name, voltage, df, actcon, fname = dataframes[i]
        ids, actcon = load_data(path_fitted, fname + ".gwy")
 
        df_stable, key = area_extract(df, actcon, name, df_drift, i, maximum_drift_x, maximum_drift_y, lines_cutoff)
       
        #df_stable= average_check(df_stable, actcon, key)
        actcon = select_color(actcon, color="Red-Cyan", key=0)
        actcon = select_range(actcon, Range*(10**-factor), key=0)
        actcon = image_save(actcon, i, path_stableframe, path, mode = image_mode, dataname= fname)
        df_stable = df_save(df_stable, actcon, dataframes_topo, fname, name, voltage_list[i]) 
        actcon = data_save(actcon, fname, path_saving = path_stableframe, path_working = path)
    return

def first_run_error(voltage_list, lines_cutoff, fnames_error, factor = 0, daten = statistics_error):
    #Evaluation of all error files 
    for i in range(len(fnames_error)):
        print("Iteration of first run error: i=")
        print(i )
        fname = fnames_error[i]
        name =  str(voltage_list[i]) + " V" + "_Error"
        ids, actcon = load_data(path, fname)
        df = select_dataframe(actcon, parameter_name=voltage_list[i], key = 0)
        df = cut_edges(df, lines_cutoff)
        df = average_check(df, actcon, key = 0)
        actcon = select_color(actcon, color="Gray", key = 0)  
        actcon = image_save(actcon, i, path_pdfs, path, mode = image_mode, dataname= fname)        
        statistics, max, min, Sme, zwischendaten, df =  get_statistic(df, actcon, fname, name, factor, daten, i, noise = 0, key=1)        
        df = df_save(df, actcon, dataframes_error, fname, name, voltage_list[i]) 
        actcon = data_save(actcon, fname, path_saving = path_fitted, path_working = path)         
        remove(actcon)
    #After first run overview evaluation and collection of values for image making and so on            
    SaveStatisticsToFile(daten, path_statistics, voltage_list, name="Statistics_Error.csv", unit="V")
    df_stat_error = load_statistics_data(path_statistics, name = "Statistics_Error.csv")    
    range_error = get_opimum_range_error(df_stat_error, column_name1 = "Maximum [V]", column_name2 = "Minimum [V]")    
    #the df lists consit of name and gwy objects 
    return df_stat_error, int(range_error), dataframes_error

def second_run_error(dataframes_error, range_error, path, factor, key=0):
    #Second run makes: average on 0, make and save: images-videos-histograms, convolution to get drift  
    #Set array_old for beginning, since its used by drift extraction function, will be overwritten in loop and reused
    print(type(dataframes_error))
    print(dataframes_error)
    for i in range(len(dataframes_error)):
        print("Second Run Error Iteration: i=")
        print(i )
        name, voltage, df, actcon, fname = dataframes_error[i]
        ids, actcon = load_data(path, fname)
        
        df = select_dataframe(actcon, parameter_name=voltage_list[i], key=0)
        print(df)
        actcon = select_range(actcon, range_error*(10**-factor), key=0)
        actcon = image_save(actcon, i, path_pdfs, path, mode = image_mode, dataname= fname)
        #make_histogram(df, name, fname, path_histo, path, factor, plot_max = range_error, plot_min = -range_error, label="Error distribution", xlabel="Error ($\it{V}$)")

        remove(actcon)
    make_histogram_all(N, dataframes_error, path_histo, path, factor, plot_max = range_error, plot_min = -range_error, name_ylabel="Error ($\it{V}$)", name_file="Error_Histogramms.pdf", list_names = voltage_list)
    df_stat_error = load_statistics_data(path_statistics, name = "Statistics_Error.csv")
    make_statistics_plot_multi(df_stat_error, 1, 4 ,6 ,5 , path_statistics, plot_max = (range_error/2), plot_min = -range_error/2, label="Error ($\it{V}$)", name = "Error_MultiPlot")
    return

def assemble(n, topo, error, current, amp, phase, time_list, path_new):
    """
    Assemble multiple scan types into combined .gwy files.
    
    Args:
        n (int): Number of scans
        topo (list): Topography scan files
        error (list): Error scan files  
        current (list): Current scan files
        amp (list): Amplitude scan files
        phase (list): Phase scan files
        time_list (list): Timestamps for each scan
        path_new (str): Output directory path
    """
    # Sort files by timestamp
    topo = sorted(topo, key=extract_time)
    error = sorted(error, key=extract_time)
    current = sorted(current, key=extract_time)
    amp = sorted(amp, key=extract_time)
    phase = sorted(phase, key=extract_time)

    for i in range(len(topo)):
        print('Iteration run: i=')
        print(i)
        
        # Load topography data
        TOPO = topo[i]
        ids_topo, actcon_topo = load_data(path, TOPO)
        df_topo, name_topo = select_dataframe(actcon_topo, ids_topo[0])
        
        # Get output filename from timestamp
        name = time_list[i]
        print(name)
        
        # Add error data if available
        if len(error) == 0:
            print('The Error list is empty')
        else:
            ERROR = error[i]
            ids_error, actcon_error = load_data(path, ERROR)
            df_error, name_error = select_dataframe(actcon_error, ids_error[0])
            id_error = gwy.gwy_app_data_browser_add_data_field(df_error, actcon_topo, 1)
            actcon_topo['/' + str(id_error) +'/data/title'] = 'Error'
            remove(actcon_error)
       
        # Add current data if available
        if len(current) == 0:
            print('The current list is empty')
        else:
            CURRENT = current[i]
            ids_current, actcon_current = load_data(path, CURRENT)
            df_current, name_current = select_dataframe(actcon_current, ids_current[0])
            id_current = gwy.gwy_app_data_browser_add_data_field(df_current, actcon_topo, 2)
            actcon_topo['/' + str(id_current) +'/data/title'] = 'Current'
            remove(actcon_current)
                        
        # Add amplitude data if available
        if len(amp) == 0:
            print('The Amplitude list is empty')
        else:
            AMP = amp[i]
            ids_amp, actcon_amp = load_data(path, AMP)
            df_amp, name_amp = select_dataframe(actcon_amp, ids_amp[0])
            id_amp = gwy.gwy_app_data_browser_add_data_field(df_amp, actcon_topo, 3)
            actcon_topo['/' + str(id_amp) +'/data/title'] = 'Amplitude'
            remove(actcon_amp)    
            
        # Add phase data if available
        if len(phase) == 0:
            print('The Phase list is empty')
        else:
            PHASE = phase[i]
            ids_phase, actcon_phase = load_data(path, PHASE)
            df_phase, name_phase = select_dataframe(actcon_phase, ids_phase[0])
            id_phase = gwy.gwy_app_data_browser_add_data_field(df_phase, actcon_topo, 4)
            actcon_topo['/' + str(id_phase) +'/data/title'] = 'Phase'
            remove(actcon_phase) 
            
        # Save combined file
        actcon_topo = data_save(actcon_topo, name, path_new, path)
        remove(actcon_topo)

    return
