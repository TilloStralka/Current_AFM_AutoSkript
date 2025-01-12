"""
AFM Data Processing Utility Functions

This module contains functions for processing and analyzing AFM (Atomic Force Microscopy) data.
Functions are organized into the following categories:

1. Loading and Saving
   - File I/O operations
   - Path handling
   - Data import/export

2. Data Processing
   - Data frame selection and manipulation
   - Edge cutting
   - Area extraction
   - Fitting and averaging

3. Visualization
   - Color selection
   - Range setting
   - Image saving
   - GIF creation
   - Histogram generation
   - Statistical plotting

4. Data Analysis
   - Statistical calculations
   - Line profile extraction
   - Peak detection
   - Gaussian fitting

5. Drift Analysis
   - Drift calculation
   - Offset determination
   - Drift list management
   - Maximum drift detection

6. Helper Functions
   - String parsing
   - Time extraction
   - Color mapping
   - Array conversion

7. Pipeline Functions
   - Topography processing pipeline
   - Current processing pipeline
   - Error processing pipeline
   - Data assembly

Each function is documented with its specific purpose and parameters.
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
########################## Functions for Loading and Saving ######################################
##################################################################################################

def sortandlist(path):
    files_topo = [k for k in names if "Topography" in k]
    files_amp = [k for k in names if "Amplitude" in k]
    files_phase = [k for k in names if "Phase" in k]  
    files_error = [k for k in names if "Error" in k]
    #Sort files by name, so the program goes through them chronologically
    files_current = sorted(files_current)
    files_topo = sorted(files_topo)
    files_amp = sorted(files_amp)
    files_phase = sorted(files_phase) 
    files_error = sorted(files_error)       
    #All files are listed here as an overview before the process starts
    print( "Folgende Current Dateien werden bearbeitet:" )
    print(files_current)
    print( "Folgende Topography Dateien werden bearbeitet:") 
    print(files_topo)
    print( "Folgende Amplituden Dateien werden bearbeitet:" )
    print(files_amp)
    print("Folgende Phasen Dateien werden bearbeitet:" )
    print(files_phase)
    print("Folgende Error Dateien werden bearbeitet:" )
    print(files_error)
    #Get Number of Elements in List for overall Histogram Plot 
    N = len(files_topo)
    print("Number of treated elements:")
    print(N)
    return N, files_topo, files_current, files_amp, files_phase, files_error

def get_info_sheet(path, name): 
    #Function to load the info data for all those scans, which must include: 
    #voltasge List, Number of lines which should be cut away, 
    #Read function 
    #Change to statistics folder for loading csv data
    os.chdir(path)
    df_info = pd.read_csv(name, sep=";", header=[0])
    print(df_info)
    #Call the column Voltage and extract it directly as list, and other infos  
    print("Cheeeeeeeeeeeeeeeeeeeeeeeeeeeck")
    voltage = df_info["Voltage [V]"].tolist()
    print("Cheeeeeeeeeeeeeeeeeeeeeeeeeeeck")
    lines_cutoff = int(df_info.iloc[0]["Lines Cutoff"])
    #Go back to working directory 
    os.chdir(path)
    return voltage, lines_cutoff

def make_folders(path):
    #Working path in which the pdfs will be saved 
    try:
        os.makedirs(path + "PDFs")        
    except OSError:
        print("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(path + "PDFs")
        os.makedirs(path + "PDFs")        
    else:
        print("Successfully made new folder") 
    #Working path in which the jpgs will be saved 
    try:
        os.makedirs(path + "gifs")        
    except OSError:
        print("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(path + "gifs")
        os.makedirs(path + "gifs")        
    else:
        print("Successfully made new folder") 
    #Working path in which the Histogramms will be saved 
    try:
        os.makedirs(path + "Histograms")        
    except OSError:
        print("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(path + "Histograms")
        os.makedirs(path + "Histograms")
    else:
        print("Successfully made new folder") 
    #Working path in which the Linescans will be saved 
    try:
        os.makedirs(path + "Linescans")        
    except OSError:
        print("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(path + "Linescans")
        os.makedirs(path + "Linescans")
    else:
        print("Successfully made new folder")         
    #Working path in which the Statistics csv plus Plots will be saved 
    try:
        os.makedirs(path + "Statistics")        
    except OSError:
        print("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(path + "Statistics")
        os.makedirs(path + "Statistics")
    else:
        print("Successfully deleted and made new")     
        #Working path in which the fitted data will be stored  
    try:
        os.makedirs(path + "Fitted")        
    except OSError:
        print("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(path + "Fitted")
        os.makedirs(path + "Fitted")
    else:
        print("Successfully deleted and made new")     
            #Working path in which the stable frame data will be stored  
    try:
        os.makedirs(path + "StableFrame")        
    except OSError:
        print("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(path + "StableFrame")
        os.makedirs(path + "StableFrame")
    else:
        print("Successfully deleted and made new")      
                #Working path in which the jpgs will be stored  
    try:
        os.makedirs(path + "JPGs")        
    except OSError:
        print("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(path + "JPGs")
        os.makedirs(path + "JPGs")
    else:
        print("Successfully deleted and made new")      
    path_pdfs = path + "PDFs"
    path_histo = path + "Histograms"    
    path_statistics = path + "Statistics"  
    path_gifs = path + "gifs"
    path_lines = path + "Linescans"
    path_fitted = path + "Fitted"
    path_stableframe = path + "StableFrame"
    path_jpgs = path + "JPGs"
    return path_pdfs, path_histo, path_statistics, path_gifs, path_jpgs, path_lines, path_fitted, path_stableframe


def load_data(path, name):
    print("Es wird gerade folgende Datei bearbeitet:")
    print(name)
    
    
    os.chdir(path)
    #Load data from path
    actcon = gwy.gwy_file_load(name, gwy.RUN_NONINTERACTIVE)
    #Load container into browser
    gwy.gwy_app_data_browser_add(actcon)
    #Ids is identifikation code or key of the active container 
    ids = gwy.gwy_app_data_browser_get_data_ids(actcon)
    print("With the following ids:")
    print(ids)
    return ids, actcon

def data_save(actcon, fname, path_saving, path_working):
    #Save files as a new .gwy file with a _fitted ending
    os.chdir(path_saving)
    #gwy.gwy_file_save(actcon, fname)
    
    gwy.gwy_file_save(actcon, fname + ".gwy")
    os.chdir(path_working)
    return actcon

def remove(actcon):
    #Remove file from browser to prevent program crashes due to overload
    gwy.gwy_app_data_browser_remove(actcon)	

def SaveStatisticsToFile(daten, path_statistics, voltage_list, name, unit):
    #Saving of all statistic data which is in the data list into a csv data 
    #Also adds the voltage list to the csv data 
    print("Save statistics data....")    
    #before saving change into statistics folder
    os.chdir(path_statistics)
    try:
        d = open(name, "w")
    except:
        print("Dateizugriff nicht erfolgreich")
        sys.exit(0)
    #Saving the data in a csv file
    d.write("Name" + ";"  +"Voltage [V]" + ";" + "R_q square roughness ["+unit+"]" + ";" + "R_a mean roughness ["+unit+"]" + ";" + "Maximum ["+unit+"]" + ";" + "Minimum ["+unit+"]" + ";" + "Average ["+unit+"]" + ";" + "SSK Skew" + ";" + "Kurtosis" + ";" + "Median ["+unit+"]" + ";" + "Sum ["+unit+"]" + ";" + "Area [mym^2]" + ";" + "Datei" + ";" + "Number of Scan[#]" + "\n" )
    i = 0
    for zwischendaten in daten:
        # zwischendaten = [str(element) for element in zwischendaten]
        # d.write(";".join(zwischendaten) + "\n")
        d.write(zwischendaten[0] + ";" + str(voltage_list[i]) + ";" + str(zwischendaten[1]) + ";" + str(zwischendaten[2]) + ";" + str(zwischendaten[3]) + ";" +   	str(zwischendaten[4]) + ";" + str(zwischendaten[5]) + ";" + str(zwischendaten[6]) + ";" + str(zwischendaten[7]) + ";" + str(zwischendaten[8]) + ";" + str(zwischendaten[9]) + ";" + str(zwischendaten[10]) + ";" + str(zwischendaten[11]) + ";" + str(zwischendaten[12])  + "\n")
        i = i + 1
    #close the just written data an return         
    d.close()
	#Change back to working directory
    os.chdir(path)    
    return 

def load_statistics_data(path_statistics, name):
    #Function to load statistics data csv data into the browser and make statistics data frame 
    #Change to statistics folder for loading csv data
    os.chdir(path_statistics)
    #Read function 
    df_stat = pd.read_csv(name, sep=";", header=[0])
    #Go back to working directory 
    os.chdir(path)
    print("Loading of csv data worked!")
    return df_stat

##################################################################################################
########################## Functions for Data Processing ##########################################
##################################################################################################

def select_dataframe(actcon, parameter_name, key):
    #Select df with key id 0 since the tiff only has one 
    df = actcon[gwy.gwy_app_get_data_key_for_id(key)] 
    df_name = actcon["/" + str(key) +"/data/title"]
    print("The name of the datafield as represented in the gwy container:")
    print(df_name)
    title = df_name + " " + str(parameter_name) + " V"
    actcon["/" + str(key) +"/data/title"] = title
    #print(actcon["/"key"/data/title"])
    #Select datafield so that the fit functions (gwy_process_func_run) are not confused
    gwy.gwy_app_data_browser_select_data_field(actcon, key)
    return df

def cut_edges(df, number_of_lines):
    #For longtime measurements sometimes the voltage is not set right for the
    #the first 2 or 4 lines (given here as number_of_lines), those will be cut away 
    x_size = df.get_xres()
    y_size = df.get_yres()
    #since counting starts with 0 we make ther staring point the number of lines
    #minus 1 since 0 is 1 
    df.resize(number_of_lines-1, number_of_lines-1, (x_size - number_of_lines -1), (y_size - number_of_lines -1))
    return df

def area_extract(dataframe, actcon, name, drift, ite, x_total_offset, y_total_offset, lines_cutoff):
    #Area extraction for a stable frame
    #The gwyddion dataframe.area_extract function takes starting point at upper eft corner and wants to know the witdth and height of rectangular image
    #Duplicate dataframe for extraction
    stable = dataframe.duplicate()
    #Get the offset for the i image 
    x_offset_fromfirst = drift.iloc[ite, 1]
    y_offset_fromfirst = drift.iloc[ite, 2]
    print("Offset from first image")
    print(x_offset_fromfirst)
    print(y_offset_fromfirst)
    #get total resolution of whole image 
    xres = dataframe.get_xres() 
    yres = dataframe.get_yres() 
    #get maximum image size of stable frame 
    x_total = xres - abs(x_total_offset) 
    y_total = yres - abs(y_total_offset) 
    print("Total image size possible:")
    print(x_total)
    print(y_total)
    #Case differentiation for drift directions 
    #x drift positive, sample moves right     
    if x_total_offset>=0:
        #y drift positive, sample moves downwards 
        if y_total_offset >= 0: 
            print("Case 1")
            a = 0 +(x_offset_fromfirst)
            b = 0 - y_offset_fromfirst
            c = x_total
            d = y_total
            stable2 = stable.area_extract(a, b, c, d)        
        #y drift positive, sample moves upwards 
        else: 
            print("Case 2")
            print("The offsets:")
            print(x_total_offset)
            print(x_offset_fromfirst)
            print(y_total_offset)
            print(y_offset_fromfirst)
            a = x_offset_fromfirst
            b = abs(y_total_offset) - abs(y_offset_fromfirst) 
            c = x_total
            d = y_total
            print("a,b,c,d is")
            print(a, b, c, d    ) 
            print("Size x and y")
            print(c - a)
            print(d - b)
            stable2 = stable.area_extract(a, b, c, d)    
            print(stable)
    #x drift negative, sample moves left 
    else:    
        #y drift positive, sample moves downwards 
        if y_total_offset >= 0: 
            print("Case 3")
            print("The offsets:")
            print(x_total_offset)
            print(x_offset_fromfirst)
            print(y_total_offset)
            print(y_offset_fromfirst)
            a = abs(x_total_offset) + x_offset_fromfirst
            b = 0 + (y_offset_fromfirst)
            c = x_total
            d = y_total
            print("a,b,c,d is")
            print(a, b, c, d   )  
            print("Size x and y")
            print(c - a)
            print(d - b)
            stable2 = stable.area_extract(a, b, c, d)
      
        #y drift positive, sample moves upwards
        else: 
            print("Case 4")
            a = 0 -(x_offset_fromfirst)
            b = 0 + y_offset_fromfirst
            c = x_total
            d = y_total
            stable2 = stable.area_extract(a, b, c, d)
    
    #Select dataframe and add to browser, so it can be treated 
    gwy.gwy_app_data_browser_add_data_field(stable2, actcon, True)

def fit_functions(actcon, df, i, mode):
    #Classical Gwyddion Fitfunctions, leveling, align rows, scar removement, fix zero 
    if i>0: 
        gwy.gwy_process_func_run("level", actcon, gwy.RUN_IMMEDIATE)
        gwy.gwy_process_func_run("align_rows", actcon, gwy.RUN_IMMEDIATE)
        gwy.gwy_process_func_run("scars_remove", actcon, gwy.RUN_IMMEDIATE)
        gwy.gwy_process_func_run("fix_zero", actcon, gwy.RUN_IMMEDIATE)
    else: 
        gwy.gwy_process_func_run("level", actcon, gwy.RUN_IMMEDIATE)
        gwy.gwy_process_func_run("align_rows", actcon, mode)
        gwy.gwy_process_func_run("scars_remove", actcon, gwy.RUN_IMMEDIATE)
        gwy.gwy_process_func_run("fix_zero", actcon, gwy.RUN_IMMEDIATE)
    df.data_changed()
    return actcon

def average_check(x, actcon, key):
    #Find the average value of the dataframe and set it as new zero to clamp 
    #around it and hopefully compensate drift in z-direction pizeo 
    statistics = x.get_stats()
    avg = statistics[0]
    x.add(-avg)
    #Select datafield again, due to duplocate there are now 2 df inside the container 
    x = actcon[gwy.gwy_app_get_data_key_for_id(key)]
    gwy.gwy_app_data_browser_select_data_field(actcon, key)    
    return x

def zero_check(df, actcon, key):
    #Get the zero offset and set it zero, as well as the width of the gaussian noise around 0
    #First make histogram and find maximum y value (here is the 0 current)        
    y, x = np.histogram(df, bins = 40000)
    index_of_maximum = np.where(y == y.max())
    array_y_max = y, x[index_of_maximum]
    offset = array_y_max[1]
    print(offset)
    count = offset.size
    if count > 1: 
        offset = offset[0]
        print("We have double maximum count value, therefore we are chosing the first one:")
        print(offset )
    print("Point of Maximum:")
    print(offset)
    #Add or Substrate to whole df 
    df.add(-offset)
    df.data_changed()
    #Make second similar df and make df absolut and fit a half gaussian distribution around it 
    df_fit = df.duplicate()
    #print(type(df_fit))
    #df_fit = df_fit.multiply_fields(df_fit, df_fit)
    #df_fit = 
    y_fit, x_fit = np.histogram(df_fit, bins = 4000)  
    #print(y_fit, x_fit)
    
    mean, var = halfnorm.stats(moments="mv")
    #print(mean, var   )
    #Select datafield again, due to duplocate there are now 2 df inside the container 
    df = actcon[gwy.gwy_app_get_data_key_for_id(key)]
    gwy.gwy_app_data_browser_select_data_field(actcon, key)
    return df, offset

##################################################################################################
########################## Functions for Visualization ###########################################
##################################################################################################

def select_color(actcon, color, key):
    """Set color palette for data visualization.
    
    Args:
        actcon: Active container
        color (str): Color palette name (e.g. 'Gwyddion.net', 'Red-Cyan')
        key (int): Data field key
        
    Returns:
        actcon: Updated container with new color palette
    """   
    actcon.set_string_by_name("/" + str(key) + "/base/palette", color)
    return actcon

def select_range(actcon, range_topo, key):
    """Set fixed scale range for visualization.
    
    Args:
        actcon: Active container
        range_val (float): Maximum absolute value for scale
        key (int): Data field key
        
    Returns:
        actcon: Updated container with fixed range
    """
    actcon.set_double_by_name("/{}/base/max".format(key), range_val)
    actcon.set_double_by_name("/{}/base/min".format(key), -range_val)
    actcon.set_int32_by_name("/{}/base/range-type".format(key), 1)
    #set range type: 0 = Full, 1 = Fixed, 2 = Automatic, 3 = Adaptive     
    #actcon.set_logarithmic(is_logarithmic)
    return actcon

def image_save(actcon, i, path_pdfs, path, mode, dataname): 
    #Save files as pdf and .tiff file 
    #First change into the pdf folder for saving and afterwards go back to working path
    os.chdir(path_pdfs)
    if i>0: 
        gwy.gwy_file_save(actcon, dataname + ".pdf", gwy.RUN_NONINTERACTIVE)
    else: 
        gwy.gwy_file_save(actcon, dataname + ".pdf", mode)
    os.chdir(path)       
    return actcon

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

def make_histogram(df, name, fname, path_histo, path, factor, plot_max,  plot_min, label, xlabel):
    """Create histogram of data field values.
    
    Args:
        df: Data field to analyze
        name (str): Display name for legend
        fname (str): Output filename
        path_histo (str): Output directory for histograms
        path (str): Working directory
        factor (int): Scaling factor for values
        plot_max/min (float): Plot axis limits
        label (str): Plot title
        xlabel (str): X-axis label
    """
    # Calculate histogram
    resolution = df.get_xres() * df.get_yres()
    y, x = np.histogram(df, bins=resolution)
    x = x * (10**factor)  # Scale x values
    
    # Create plot
    fig, ax = plt.subplots(figsize=(10,5))
    ax.plot(x[:-1], y, "b", label=name)
    ax.legend(loc=1, fancybox=True, framealpha=1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    ax.set_xlim(plot_min, plot_max)
    plt.semilogy()
    plt.tight_layout()
    
    # Save plot
    os.chdir(path_histo)
    plt.savefig(fname.replace(".tiff", "_Histo.pdf"), format="pdf")
    os.chdir(path)
    return 

def make_histogram_all(N, dataframes, path_histo, path, factor, plot_max, plot_min, name_ylabel, name_file, list_names):
    #Histogram-Plot von allen Histogrammen in einem Plot dazu durch alle eintraege in der Dataframes Liste gehen
    #Settings for the plot  
    fig, ax = plt.subplots(figsize=(10,4))
    #Function is looping through all dfs to make multi histo plot, its reversed so that first one is most non transparent 
    for i in range(N):
        #print("Index i:")
        zeile = dataframes[i]
        #Make Histogram bins with numpy         
        y, x = np.histogram(zeile[2], bins = 500)
        #Normierung on and off
        #y = (y/float(np.max(y)))
        #multiply the x value with the factor, important for topo scans in nm range 
        #print("We are in the histogram function, here comes the x value hopefully in nm:")
        x = x*(10**factor) 
        #Give histograms name of samples, they will be renamed in the plot function to get rid of data type ending   
        ax.plot(x[:-1],y,color=get_colormap(i,N), label = str(list_names[i]) + " V",  alpha=0.7) 
        #for legend label = name.replace(".tiff", " V"),
        #Fill area under histo curve with rising transparancy over N dfs 
        ax.fill_between(x[:-1], y, 0, color=get_colormap(i,N), alpha=0.2)
        #for transparency alpha=(0.6 -(i*0.01))           
    #Settings for the plot        
    ax.legend(loc = 1, fancybox = True, framealpha = 1, ncol=3)
    ax.set_xlabel(name_ylabel)
    ax.set_ylabel("Count")
    ax.set_xlim(plot_min, plot_max)
    #ax.set_ylim(10**0,4*10**2)
    plt.semilogy()
    #plt.semilogx()
    plt.tight_layout()
    #If required, save funktion of plot
    #before saving change into histo folder
    os.chdir(path_histo)
    #Save function
    plt.savefig(name_file, format="pdf")
    #Change back to working directory
    os.chdir(path)
    #plt.show()
    return 

def make_statistics_plot(df_stat, path_save, x_column, y_column, label):
    print("We are in the statistics plot function!")
    #Make a string list of the header of the statistics df to call them in plot
    names = list(df_stat.columns)
    #Define the number of plots with the number of elements in y_column
    n = len(y_column) 
    print(n)
    for i in y_column:
        print(i)
        #indexing the list to get a consistent spread of color change
        N = y_column.index(i)
        print(N)
        #Making a plot of the statistic values of the evaluated scans
        fig, ax = plt.subplots(figsize=(5,5))
        # determine x-data
        #print(df_stat.iloc[:, i])
        y = df_stat.iloc[:, i]
        #print(x.max)
        #print(df_stat.iloc[:, 10])
        x = df_stat.iloc[:, x_column]
        #more settings for the plot, color scale over n number of plots with index N in the list
        ax.plot(x,y, marker = "p", color = get_colormap(N ,n))        
        #Settings for the plot        
        ax.legend(loc = "best", fontsize = 14, fancybox = True, framealpha = 1)
        ax.set_xlabel(names[x_column], fontsize = 14)
        ax.set_ylabel(names[i], fontsize = 14)
        #Set axis labeling
        ax.tick_params(axis="x", labelsize=14)
        ax.tick_params(axis="y", labelsize=14)
        #tight layout prevents cutting of during pdf making 
        plt.tight_layout()
        os.chdir(path_save)
        print(os.getcwd())
        plt.savefig(names[i] + label + ".pdf", format="pdf")
        plt.show() 
    return  

def make_statistics_plot_multi(df, x_column, y_column1, y_column2, y_column3, 
                             path_statistics, plot_min, plot_max, label, name):
    """Create multi-line plot of statistics data.
    
    Args:
        df: Statistics dataframe
        x_column (int): Column index for x-axis data
        y_column1/2/3 (int): Column indices for y-axis data
        path_statistics (str): Output directory
        plot_min/max (float): Y-axis limits
        label (str): Y-axis label
        name (str): Output filename base
    """
    fig, ax = plt.subplots(figsize=(10,5))
    
    # Plot data series
    x = df.iloc[:,x_column]
    for y_col, color in [(y_column1, "red"), 
                        (y_column2, "blue"),
                        (y_column3, "black")]:
        y = df.iloc[:,y_col]
        ax.plot(x, y, marker="p", color=color, linestyle="None")
    
    # Configure plot
    ax.legend(loc=4, fontsize=14, fancybox=True, framealpha=1)
    ax.set_xlabel("Voltage ($\it{V}$)", fontsize=14)
    ax.set_ylabel(label, fontsize=14)
    ax.tick_params(axis="both", labelsize=14)
    ax.set_ylim(plot_min, plot_max)
    plt.tight_layout()
    
    # Save plot
    os.chdir(path_statistics)
    plt.savefig("{}.pdf".format(name), format="pdf")
    os.chdir(path)
    plt.show()
    return 

##################################################################################################
########################## Functions for Data Analysis ##########################################
##################################################################################################

def get_statistic(x, actcon, filename, Name, factor, daten, i, noise, key):
    """
    Extract statistical values from a data field and append to statistics list.
    
    Args:
        x: Input data field
        actcon: Active container
        filename: Name of file
        Name: Display name
        factor: Scaling factor for values
        daten: List to store statistics
        i: Current iteration
        noise: Noise threshold
        key: Data field key
        
    Returns:
        tuple: (statistics, max, min, sum, intermediate_data, data_field)
    """
    # Get basic statistics
    max_val = x.get_max()
    stats = x.get_stats()
    
    # Extract individual statistics
    avg, ra, rms, skew, kurtosis = stats
    min_val = x.get_min()
    median = x.get_median()
    sum_val = x.get_sum()
    surface_area = x.get_surface_area()

    # Select data field
    df = actcon[gwy.gwy_app_get_data_key_for_id(key)]
    gwy.gwy_app_data_browser_select_data_field(actcon, key)

    # Scale values by factor (surface area uses different scaling)
    scaled_vals = [
        val * (10**factor) for val in 
        [avg, ra, rms, max_val, min_val, median, sum_val]
    ]
    surface_area *= 10**12  # Convert to square micrometers
    
    # Round all values to 2 decimal places
    avg, ra, rms, max_val, min_val, median, sum_val = [
        round(val, 2) for val in scaled_vals
    ]
    skew = round(skew, 2)
    kurtosis = round(kurtosis, 2)
    surface_area = round(surface_area, 2)

    # Compile intermediate data
    intermediate_data = [
        Name, ra, rms, max_val, min_val, avg, 
        skew, kurtosis, median, sum_val, surface_area,
        filename, i
    ]
    
    print("Statistical values for {}:".format(filename))
    print(intermediate_data)
    
    daten.append(intermediate_data)
    
    return stats, max_val, min_val, sum_val, intermediate_data, df

def get_line(df, factor, x_start, y_start, x_end, y_end, res):
    """
    Extract line profile from data field.
    
    Args:
        df: Input data field
        factor: Scaling factor
        x_start, y_start: Starting coordinates
        x_end, y_end: Ending coordinates
        res: Resolution of extracted line
        
    Returns:
        tuple: (line_data, x_positions, point_size)
            - line_data: Extracted line profile
            - x_positions: Array of x coordinates
            - point_size: Physical size of each point
    """
    # Extract line profile (thickness=5, interpolation=2)
    df_line = df.get_profile(x_start, y_start, x_end, y_end, res, 5, 2)
    
    # Scale values
    df_line.multiply(10**factor)
    
    # Calculate physical dimensions
    real_size = df_line.get_real()
    line_res = df_line.get_res()
    point_size = (10**factor) * (real_size/line_res)
    
    print("Physical size per point: {} nm".format(point_size))
    
    # Generate x coordinates
    x_positions = [l * point_size for l in range(line_res)]
    
    return df_line, x_positions, point_size

def peak_extraction(i, x_data, y_data, step_size, peak_x_position, startingpoint, width):
    """
    Extract peak from line scan data within a specified window.
    
    Args:
        i (int): Iteration index (0 for first scan, >0 for subsequent)
        x_data (array): X-axis data points
        y_data (array): Y-axis data points (e.g. current/height values)
        step_size (float): Physical distance between data points
        peak_x_position (float): Previous peak position (used if i>0)
        startingpoint (float): Initial x position to start search (used if i=0)
        width (float): Width of window to search for peak in physical units
        
    Returns:
        tuple: (peak_x_position, x_value)
            - peak_x_position (int): Index position of peak
            - x_value (float): X coordinate of peak
            
    Notes:
        - For first scan (i=0), searches around startingpoint
        - For subsequent scans, searches around previous peak position
        - Window width is converted from physical units to array indices
        - Returns both array index and physical position of peak
    """
    # Convert y_data to numpy array and take absolute values
    y_all = y_data.get_data()    
    y = np.asarray(y_all)
    y = np.absolute(y)
    x = np.asarray(x_data)

    # Convert window width from physical units to number of data points
    n_points = int(round(width/step_size))
    print("Number of points in search window: {}".format(n_points))
    
    # Extract window of data to search for peak
    if i == 0:
        # First scan: center window on startingpoint
        idx = (np.abs(x-startingpoint)).argmin()
        x = x[idx-n_points:idx+n_points]
        y = y[idx-n_points:idx+n_points]
    else:
        # Subsequent scans: center window on previous peak
        startingpoint = x[peak_x_position] 
        idx = (np.abs(x-startingpoint)).argmin()
        x = x[idx-n_points:idx+n_points]
        y = y[idx-n_points:idx+n_points]

    # Find peaks in the windowed data
    peaks, _ = find_peaks(y, distance=1000)  
    peaks = int(round(peaks))
    
    # Convert peak position to full array index and get physical x value
    peak_x_position = int(round(x[peaks]))
    x_value = x[peaks]

    return peak_x_position, x_value

def Gauss(x, a, x0, sigma):
    """
    Gaussian function for curve fitting.
    
    Args:
        x (array): X values to evaluate function at
        a (float): Amplitude
        x0 (float): Center position
        sigma (float): Standard deviation
        
    Returns:
        array: Gaussian values evaluated at x positions
    """
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

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


##################################################################################################
########################## Functions for Drift Analysis ##########################################
##################################################################################################

def get_drift(x, old_array, i):
    """
    Calculate drift between consecutive AFM images using FFT convolution.
    
    Args:
        x: Current image datafield
        old_array: Previous image array (or None for first image)
        i: Current iteration index
        
    Returns:
        tuple: (array_for_next_iteration, drift_vector)
            - array_for_next_iteration: Current image as numpy array
            - drift_vector: (x,y) tuple of drift in pixels
            
    Notes:
        Uses FFT convolution to find the offset between consecutive images.
        For first image (i=0), establishes baseline and returns zero drift.
        For subsequent images, calculates drift relative to previous image.
    """
    #print("The start position is the half of the resolution (middle of image)"
    start_position = ((x.get_xres() / 2), (x.get_yres() / 2))
    
    if i == 0:
        # First image - establish baseline
        old_array = print_df_NumpyArray(x)
        # Self-convolve to get reference correlation
        corr_img = scipy.signal.fftconvolve(old_array, old_array[::-1,::-1], mode="same")
        # Find correlation peak
        start_position = np.unravel_index(np.argmax(corr_img), corr_img.shape)
        return old_array, (0, 0)
    
    else:
        # Convert current image to array
        array_new = print_df_NumpyArray(x)
        
        # Cross-correlate with previous image using FFT convolution
        corr_img = scipy.signal.fftconvolve(old_array, array_new[::-1,::-1], mode="same")
        
        # Find correlation peak position
        position = np.unravel_index(np.argmax(corr_img), corr_img.shape)
        
        # Calculate drift as difference from start position
        drift = (
            start_position[0] - position[0],
            start_position[1] - position[1]
        )
        
        return array_new, drift
    
def get_offset_from_first_image(drift, offset, i):
    """
    Calculate cumulative offset from first image in series.
    
    Args:
        drift: Current (x,y) drift vector
        offset: Previous cumulative offset
        i: Current iteration index
        
    Returns:
        tuple: Updated cumulative (x,y) offset from first image
        
    Notes:
        For i=0, returns initial offset.
        For i>0, adds current drift to previous offset.
    """
    if i == 0:
        return offset
        
    else:
        new_offset = (
            offset[0] + drift[0],
            offset[1] + drift[1]
        )
        print("For iteration > 0 the offset from the image before:", new_offset)
        return new_offset

def make_drift_list(drift, offset_by_drift, i, dlist, path, path_statistics, dataframes, name):
    """
    Record drift measurements and save to CSV file.
    
    Args:
        drift: Current (x,y) drift vector
        offset_by_drift: Cumulative (x,y) offset from first image
        i: Current iteration index
        dlist: List to store drift measurements
        path: Working directory path
        path_statistics: Statistics output directory path
        dataframes: List of image dataframes
        name: Output filename
        
    Returns:
        list: Updated drift measurements list
        
    Notes:
        Saves drift data as CSV with columns:
        - Iteration
        - X/Y offset from first image
        - X/Y drift from previous image
    """
    # Add current measurements to list
    current_data = (i, offset_by_drift[0], offset_by_drift[1], drift[0], drift[1])
    dlist.append(current_data)
    
    # Save data on all but last iteration
    if i != len(dataframes) - 1:
        os.chdir(path_statistics)
        print("Saving drift data...")
        
        try:
            with open(name, "w") as d:
                # Write header
                d.write("Iteration;X Offset by drift from first image [px];" +
                       "Y Offset by drift from first image [px];" +
                       "X Drift from image before [px];" +
                       "Y Drift from image before [px]\n")
                
                # Write data rows
                for data in dlist:
                    d.write(";".join(str(x) for x in data) + "\n")
                    
        except IOError:
            print("Failed to access file")
            sys.exit(0)
            
        os.chdir(path)
        
    else:
        print("Current drift list:", dlist)
        
    return dlist

def get_maximum_drift(df, column_name1, column_name2):
    """
    Calculate maximum drift extent in X and Y directions.
    
    Args:
        df: Drift measurements dataframe
        column_name1: Column name for X drift values
        column_name2: Column name for Y drift values
        
    Returns:
        tuple: (max_x_drift, max_y_drift)
            Maximum absolute drift values in each direction
            
    Notes:
        Compares positive and negative drift extremes to find
        the largest magnitude in each direction.
    """
    # Get X drift values and find maximum magnitude
    x_drifts = df[column_name1]
    max_x = x_drifts.max()
    min_x = abs(x_drifts.min())
    total_drift_x = max_x if max_x >= min_x else -min_x
    
    # Get Y drift values and find maximum magnitude  
    y_drifts = df[column_name2]
    max_y = y_drifts.max()
    min_y = abs(y_drifts.min())
    total_drift_y = max_y if max_y >= min_y else -min_y
    
    return total_drift_x, total_drift_y

##################################################################################################
########################## Helper Functions ####################################################
##################################################################################################

def extract_voltage(string):
    """
    Extract voltage value from a filename string that contains '_V' pattern.
    
    Args:
        string (str): Filename containing voltage information (e.g. 'scan_2.5V_001.dat')
        
    Returns:
        float: Extracted voltage value, or None if no voltage found
        
    Example:
        >>> extract_voltage('scan_2.5V_001.dat')
        2.5
    """
    for element in string.split("_"):
        # Look for elements containing both 'V' and '.' (e.g. '2.5V')
        if "V" in element and "." in element:
            return float(element.replace("V",""))
    return None

def extract_time(string):
    """
    Extract timestamp from a filename string that contains timestamp pattern.
    
    Args:
        string (str): Filename containing timestamp (e.g. 'scan_2.5V_20230401_001234.dat')
        
    Returns:
        str: Extracted timestamp string, or None if no timestamp found
        
    Example:
        >>> extract_time('scan_2.5V_20230401_001234.dat') 
        '20230401_001234'
    """
    for element in string.split("_"):
        # Look for elements starting with '00' and containing '.' 
        # This pattern matches typical timestamp formats
        if "00" in element and "." in element:
            return str(element)
    return None

def get_time(list_of_filenames):
    """
    Extract timestamps from a list of filenames and return sorted list.
    
    Args:
        list_of_filenames (list): List of filenames containing timestamps
        
    Returns:
        list: List of extracted timestamps (first 9 chars only)
        
    Example:
        >>> get_time(['scan_001_20230401.dat', 'scan_002_20230402.dat'])
        ['20230401', '20230402']
    """
    # Sort filenames by timestamp
    sorted(list_of_filenames, key=extract_time)
    list_TIME = []
    
    for filename in list_of_filenames:
        ending = extract_time(filename)
        if ending:
            # Take first 9 chars of timestamp
            timestamp = ending[0:9]
            list_TIME.append(timestamp)
    
    return list_TIME

def get_colormap(index, N):
    """
    Generate a color from blue->red colormap based on index position.
    
    Args:
        index (int): Current index position
        N (int): Total number of colors needed
        
    Returns:
        tuple: RGBA color values
        
    Notes:
        Creates a continuous color gradient from blue (negative/electrons) 
        to red (positive) with N segments.
    """
    # Create custom colormap from blue to red
    colormap = mpl.colors.LinearSegmentedColormap.from_list(
        "custom", 
        [(0,"blue"), (1,"red")], 
        N=N
    )
    
    # Set up normalization and mapping
    color_norm = colors.Normalize(vmin=0, vmax=N)
    scalar_map = cmx.ScalarMappable(norm=color_norm, cmap=colormap)
    
    return scalar_map.to_rgba(index)

def print_df_NumpyArray(df):
    """
    Convert a Gwyddion DataField object to a numpy array.
    
    Args:
        df: Gwyddion DataField object
        
    Returns:
        numpy.ndarray: Array containing the DataField values
    """
    return gwyutils.data_field_data_as_array(df)

##################################################################################################
########################## Pipeline Functions #################################################
##################################################################################################

def first_run_topo(voltage_list, lines_cutoff, fnames_topo, factor = 9, daten):
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
        print(f"First Run Iteration: i={i}")
        name = f"{voltage_list[i]} V_Topo"
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
        print(f"Second Run Topo Iteration: i={i}")
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
    sorted(topo, key=extract_time)
    sorted(error, key=extract_time)
    sorted(current, key=extract_time)
    sorted(amp, key=extract_time)
    sorted(phase, key=extract_time)

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
        actcon_topo = data_save(actcon_topo, name, path_newData, path)
        remove(actcon_topo)

    return
