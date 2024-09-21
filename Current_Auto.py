"""
Writer: 
Tillmann Stralka 2020.Mai.05
Owner: 
HLP University of Leipzig
About:
Python version 2.7
Automatic procession of AFM files.
Finding, listing, fitfunctions, image making, statistic extraction and plotting 
Correlation of Current and topographic information, convolution of dynamicscans
"""


#Some python libs have to be added via PATH since they are near gwyddion
import sys
sys.path.append("/usr/local/opt/python@2/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages")
sys.path.append("/usr/share/gwyddion/pygwy")
import pygtk
pygtk.require20() # adds gtk-2.0 folder to sys.path
sys.path.append("/usr/local/Cellar/gwyddion/2.52/share/gwyddion/pygwy")
import gwyutils
#Importing all other necessary libraries from the python surounding 
import gwy
import os, time
import numpy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl
import numpy as np
from scipy.stats import norm
from scipy.stats import halfnorm
from scipy.optimize import curve_fit
from scipy.signal import find_peaks
import scipy.signal
import sys 
import pandas as pd  
import shutil  
import cv2 
import imageio
from pdf2image.exceptions import (PDFInfoNotInstalledError, PDFPageCountError, PDFSyntaxError)
from pdf2image import convert_from_path, convert_from_bytes
import glob


###############################################################################
        #   Declare empty lines and dataframes   #
###############################################################################
#Create an empty line in which all data is stored before it gets saved to another file 
#statistics_current and so are empty lines where the statistic values will be saved, and the lists will be made to csv datas 
statistics_current = []
statistics_topo = []
statistics_error = []
statistics_current_gb = []
statistics_current_grain = []
dataframes_current = []
dataframes_topo = []
dataframes_error = []
datalines_current = []
datalines_topo = []
datalines_distance = []
#For convolution declare drift at the beginning, zeros     
array_old = 0
array_old2 = 0
offset_by_drift = (0,0)
offset_by_drift2 = (0,0)
#For Line extraction
line_x_start = 46
line_y_start = 47
line_x_end = 58
line_y_end = 34 
line_res = 1000
peak_x_position = 0 #only for the beginning set 0, then it will be overwritten 





###############################################################################
        #   Subfunctions    #
###############################################################################
def working_path():
    #Working path in which the files to be processed are located, either way Linux or Mac
    path_linux = "/home/tillmann/Desktop/DataToEvaluate/"
    path_mac = "/Users/tillo/DataToEvaluate/"
    if os.path.isdir(path_linux):
        print("Seems to be a linux")
        path = path_linux
    else:
        print( "Seems to be a mac")
        path = path_mac
    return path
        
def sortandlist(path): 
    #Get file names from the working path, only .tiff names 
    names = os.listdir(path)
    #The data is in the folder sorted after name (which is chronologically)
    files_current = [k for k in names if "Current" in k]
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
    print "Folgende Phasen Dateien werden bearbeitet:" 
    print(files_phase)
    print "Folgende Error Dateien werden bearbeitet:" 
    print(files_error)
    #Get Number of Elements in List for overall Histogram Plot 
    N = len(files_topo)
    print "Number of treated elements:"
    print(N)
    return N, files_topo, files_current, files_amp, files_phase, files_error

def get_info_sheet(path, name):
    #Function to load the info data for all those scans, which must include: 
    #voltasge List, Number of lines which should be cut away, 
    #Read function 
    #Change to statistics folder for loading csv data
    os.chdir(path)
    df_info = pd.read_csv(name, sep=";", header=[0])
    print df_info
    #Call the column Voltage and extract it directly as list, and other infos  
    print "Cheeeeeeeeeeeeeeeeeeeeeeeeeeeck"
    voltage = df_info["Voltage [V]"].tolist()
    print "Cheeeeeeeeeeeeeeeeeeeeeeeeeeeck"
    lines_cutoff = int(df_info.iloc[0]["Lines Cutoff"])
    #Go back to working directory 
    os.chdir(path)
    return voltage, lines_cutoff

def make_folders(path):
    #Working path in which the pdfs will be saved 
    try:
        os.makedirs(path + "PDFs")        
    except OSError:
        print ("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(path + "PDFs")
        os.makedirs(path + "PDFs")        
    else:
        print ("Successfully made new folder") 
    #Working path in which the jpgs will be saved 
    try:
        os.makedirs(path + "gifs")        
    except OSError:
        print ("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(path + "gifs")
        os.makedirs(path + "gifs")        
    else:
        print ("Successfully made new folder") 
    #Working path in which the Histogramms will be saved 
    try:
        os.makedirs(path + "Histograms")        
    except OSError:
        print ("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(path + "Histograms")
        os.makedirs(path + "Histograms")
    else:
        print ("Successfully made new folder") 
    #Working path in which the Linescans will be saved 
    try:
        os.makedirs(path + "Linescans")        
    except OSError:
        print ("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(path + "Linescans")
        os.makedirs(path + "Linescans")
    else:
        print ("Successfully made new folder")         
    #Working path in which the Statistics csv plus Plots will be saved 
    try:
        os.makedirs(path + "Statistics")        
    except OSError:
        print ("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(path + "Statistics")
        os.makedirs(path + "Statistics")
    else:
        print ("Successfully deleted and made new")     
        #Working path in which the fitted data will be stored  
    try:
        os.makedirs(path + "Fitted")        
    except OSError:
        print ("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(path + "Fitted")
        os.makedirs(path + "Fitted")
    else:
        print ("Successfully deleted and made new")     
            #Working path in which the stable frame data will be stored  
    try:
        os.makedirs(path + "StableFrame")        
    except OSError:
        print ("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(path + "StableFrame")
        os.makedirs(path + "StableFrame")
    else:
        print ("Successfully deleted and made new")      
                #Working path in which the jpgs will be stored  
    try:
        os.makedirs(path + "JPGs")        
    except OSError:
        print ("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(path + "JPGs")
        os.makedirs(path + "JPGs")
    else:
        print ("Successfully deleted and made new")      
    path_pdfs = path + "PDFs"
    path_histo = path + "Histograms"    
    path_statistics = path + "Statistics"  
    path_gifs = path + "gifs"
    path_lines = path + "Linescans"
    path_fitted = path + "Fitted"
    path_stableframe = path + "StableFrame"
    path_jpgs = path + "JPGs"
    return path_pdfs, path_histo, path_statistics, path_gifs, path_jpgs, path_lines, path_fitted, path_stableframe

def extract_voltage(string):
    for element in string.split("_"):
        if "V" in element and "." in element:
            #print "First element in list only the Voltage number extraction:"
            #print (float(element.replace("V","")))
            return float(element.replace("V",""))
        
def extract_time(string):        
    for element in string.split("_"):
        if "00" in element and "." in element:
            #print "First element in list only the Voltage number extraction:"
            #print (float(element.replace("V","")))
            return str(element) 
    
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
    print "With the following ids:"
    print ids
    return ids, actcon

def select_dataframe(actcon, parameter_name, key):
    #Select df with key id 0 since the tiff only has one 
    df = actcon[gwy.gwy_app_get_data_key_for_id(key)] 
    df_name = actcon["/" + str(key) +"/data/title"]
    print "The name of the datafield as represented in the gwy container:"
    print df_name
    title = df_name + " " + str(parameter_name) + " V"
    actcon["/" + str(key) +"/data/title"] = title
    #print actcon["/"key"/data/title"]
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
    print "Offset from first image"
    print x_offset_fromfirst
    print y_offset_fromfirst
    #get total resolution of whole image 
    xres = dataframe.get_xres() 
    yres = dataframe.get_yres() 
    #get maximum image size of stable frame 
    x_total = xres - abs(x_total_offset) 
    y_total = yres - abs(y_total_offset) 
    print "Total image size possible:"
    print x_total
    print y_total
    #Case differentiation for drift directions 
    #x drift positive, sample moves right     
    if x_total_offset>=0:
        #y drift positive, sample moves downwards 
        if y_total_offset >= 0: 
            print "Case 1"
            a = 0 +(x_offset_fromfirst)
            b = 0 - y_offset_fromfirst
            c = x_total
            d = y_total
            stable2 = stable.area_extract(a, b, c, d)        
        #y drift positive, sample moves upwards 
        else: 
            print "Case 2"
            print "The offsets:"
            print x_total_offset
            print x_offset_fromfirst
            print y_total_offset
            print y_offset_fromfirst
            a = x_offset_fromfirst
            b = abs(y_total_offset) - abs(y_offset_fromfirst) 
            c = x_total
            d = y_total
            print "a,b,c,d is"
            print a, b, c, d     
            print "Size x and y"
            print c - a
            print d - b
            stable2 = stable.area_extract(a, b, c, d)    
            print stable
    #x drift negative, sample moves left 
    else:    
        #y drift positive, sample moves downwards 
        if y_total_offset >= 0: 
            print "Case 3"
            print "The offsets:"
            print x_total_offset
            print x_offset_fromfirst
            print y_total_offset
            print y_offset_fromfirst
            a = abs(x_total_offset) + x_offset_fromfirst
            b = 0 + (y_offset_fromfirst)
            c = x_total
            d = y_total
            print "a,b,c,d is"
            print a, b, c, d     
            print "Size x and y"
            print c - a
            print d - b
            stable2 = stable.area_extract(a, b, c, d)
      
        #y drift positive, sample moves upwards
        else: 
            print "Case 4"
            a = 0 -(x_offset_fromfirst)
            b = 0 + y_offset_fromfirst
            c = x_total
            d = y_total
            stable2 = stable.area_extract(a, b, c, d)
    
    #Select dataframe and add to browser, so it can be treated 
    gwy.gwy_app_data_browser_add_data_field(stable2, actcon, True)


























    
    print "The ids of all df in the active container:"
    ids_list =  gwy.gwy_app_data_browser_get_data_ids(actcon)
    print ids_list
    print type(ids_list)
    print "Size of new df:"
    print stable2.get_xres()
    print stable2.get_yres()
    print "ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ"
    print ids_list[-1]
    print name 
    
    key = ids_list[-1]
    
    
    gwy.gwy_app_data_browser_select_data_field(actcon, ids_list[-1])
    #df_small = select_dataframe(actcon, parameter_name=voltage_list[i], key = ids_list[-1])
    #rename titel in the new created channel to get a better subscription within the gifs     
    actcon["/" + str(key) +"/data/title"] = name
    return stable2, key
  
    
def select_color(actcon, color, key):
    #Datafield output with setting for color palette    
    actcon.set_string_by_name("/" + str(key) + "/base/palette", color)
    return actcon

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

def select_range(actcon, range_topo, key):
    # Setting limits for the scale in the resulting image (Fixed imaging)
    actcon.set_double_by_name("/" + str(key) + "/base/max", range_topo)
    actcon.set_double_by_name("/" + str(key) + "/base/min", -range_topo)
    actcon.set_int32_by_name("/" + str(key) + "/base/range-type", 1)
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
    #print "The unsorted list:"
    filenames = os.listdir(path)
    #print filenames
    jpgs_list = [k for k in filenames if name in k]
    #print jpgs_list
    # Sort them with key definition for voltage extraction and order them with int Voltage
    jpgs_list = sorted(jpgs_list,key=extract_time)
    #print "The sorted list:"
    #print jpgs_list    
    #jpgs_list.reverse()
    #print jpgs_list
    cwd = os.getcwd()
    #print cwd
    if cwd == path_pdfs:
        #print "We are making AFM scan Gifs"
        resolution = 1000
    else:
        #print "We are making histograms or linescan Gifs"
        resolution = 100
        
    print "Making Video from image list:"    
    with imageio.get_writer(name + ".gif", mode="I", duration = time) as writer:
        for image in jpgs_list:
            #print "We are in the image loop, image type:"
            #print type(image)
            #print "We are in the image loop, image name:"
            #print image
            image = make_pdf_to_jpg(image, resolution)
            #print "Image after conversion:"
            #print image
            #Change path into jpg folder where the converdet jpg filed (made from the pdf are stored)
            os.chdir(path_jpgs)
            img = imageio.imread(image)
            #print "Type with the image reading function"
            #print type(img)
            os.chdir(path)
            writer.append_data(img)
    writer.close()            
    return

#################


def make_pdf_to_jpg(image, resolution):
    #A little function to convert pdfs to pixel images, in that case jpgs to make gifs or
    #Videos, since the video makers can not handle vector graphics 
    #The function returns a string for a file name, where .pdf is replaced by .jpg 
    #print "We are in the image converter function, type and name:"
    name_long = image
    name = name_long[:-4]
    #print type(name)
    #print name
    images = convert_from_path(image, resolution)
    for image in images:
        os.chdir(path_jpgs)
        image.save(name + ".jpg", "JPEG") 
    #print "Conversion worked!"    
    return (name + ".jpg")

def get_statistic(x, actcon, filename, Name, factor, daten, i, noise, key):
    #Extraktion of the statistic values avg, ra, rms, skwe, kurosis and append to the empty statistics list 
    #First the noise will be removed from the dataframe (a mask is beeing applied)
    #Steps are: duplicate, get max, set range (cutoff), get stat, select normal df again 
    print "We are in the statistics function" 
    #Get maximum for the clamps 
    max = x.get_max()
    #Extraktion of the statistic values avg, ra, rms, skwe, kurosis and append to the empty statistics list 
    statistics = x.get_stats()
    avg = statistics[0]
    ra = statistics[1]
    rms = statistics[2]
    skwe = statistics[3]
    kurosis = statistics[4]
    min = x.get_min()
    median = x.get_median()
    Sme = x.get_sum()
    surfc = x.get_surface_area()
    #Select datafield again, due to duplocate there are now 2 df inside the container 
    df = actcon[gwy.gwy_app_get_data_key_for_id(key)]
    gwy.gwy_app_data_browser_select_data_field(actcon, key)    
    #recalculte all statistic values to the order of magnitude that makes sense 
    #the area SURFC is given in square mycrometer, therefore here no recalculation
    avg, ra, rms, max, min, median, Sme, surfc = avg*(10**factor), ra*(10**factor), rms*(10**factor), max*(10**factor), min*(10**factor), median*(10**factor), Sme*(10**factor), surfc*(10**12)
    #round them to 2 digits after point 
    avg = round(avg, 2)    
    ra = round(ra, 2)
    rms = round(rms, 2)
    skwe = round(skwe, 2)
    kurosis = round(kurosis, 2)
    max = round(max, 2)
    min = round(min, 2)
    Sme = round(Sme, 2)
    median = round(median,2)    
    surfc = round(surfc,2)    
    zwischendaten = [Name, ra, rms, max, min, avg, skwe, kurosis, median, Sme, surfc, filename, i]
    print "Here the statistic values"
    print zwischendaten
    print filename
    daten.append(zwischendaten)
    return statistics, max, min, Sme, zwischendaten, df 

def df_save(df, actcon, dataframes, fname, name, voltage):
    #saves df AFTER! fitting into the list of dataframes for the overview histogram  
    zwischendaten2 = [name, voltage, df, actcon, fname]
    dataframes.append(zwischendaten2)
    print "Dataframe has been added to the list of dataframes"
    print df
    return df

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
    print "Loading of csv data worked!"
    return df_stat

def make_statistics_plot(df_stat, x_column, y_column, path_statistics):
    print "We are in the statistics plot function!"
    #Make a string list of the header of the statistics df to call them in plot
    names = list(df_stat.columns)
    #Define the number of plots with the number of elements in y_column
    n = len(y_column)    
    #make single plots 
    for i in y_column:
        #indexing the list to get a consistent spread of color change
        N = y_column.index(i)
        #Making a plot of the statistic values of the evaluated scans
        fig, ax = plt.subplots(figsize=(5,5))
        # determine x-data
        #print df_stat.iloc[:, i]
        y = df_stat.iloc[:, i]
        #print x.max
        #print df_stat.iloc[:, 10]
        x = df_stat.iloc[:, x_column]
        #more settings for the plot, color scale over n number of plots with index N in the list
        ax.plot(x,y, marker = "p", color = get_colormap(N ,n))        
        #Settings for the plot        
        ax.legend(loc = 2, fontsize = 14, fancybox = True, framealpha = 1)
        ax.set_xlabel(names[x_column], fontsize = 14)
        ax.set_ylabel(names[i], fontsize = 14)
        #Set axis labeling
        ax.tick_params(axis="x", labelsize=14)
        ax.tick_params(axis="y", labelsize=14)
        #tight layout prevents cutting of during pdf making 
        plt.tight_layout()
        
        #Change to statistics folder for saving statistic plots
        os.chdir(path_statistics)
        #Save function
        plt.savefig(names[i] +".pdf", format="pdf")
        #Go back to working directory 
        os.chdir(path)
        plt.show() 
    return  

def make_statistics_plot(df_stat, path_save, x_column, y_column, label):
    print "We are in the statistics plot function!"
    #Make a string list of the header of the statistics df to call them in plot
    names = list(df_stat.columns)
    #Define the number of plots with the number of elements in y_column
    n = len(y_column) 
    print n
    for i in y_column:
        print i
        #indexing the list to get a consistent spread of color change
        N = y_column.index(i)
        print N
        #Making a plot of the statistic values of the evaluated scans
        fig, ax = plt.subplots(figsize=(5,5))
        # determine x-data
        #print df_stat.iloc[:, i]
        y = df_stat.iloc[:, i]
        #print x.max
        #print df_stat.iloc[:, 10]
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
        print os.getcwd()
        plt.savefig(names[i] + label + ".pdf", format="pdf")
        plt.show() 
    return  

def make_statistics_plot_multi(df, x_column, y_column1, y_column2, y_column3, path_statistics, plot_min, plot_max, label, name):
    #make multiplots with one x-axis and several y data
    print "We are in the statistics multi plot function!"
    #Making a plot of the statistic values of the evaluated scans
    fig, ax = plt.subplots(figsize=(10,5))        
    #determine x-data
    #print df_stat.iloc[:, i]
    x = df.iloc[:,x_column]
    #print x.max
    #determin y datas 
    y1 = df.iloc[:,y_column1]
    y2 = df.iloc[:,y_column2]
    y3 = df.iloc[:,y_column3]
    #setting color and mark style for graphs
    ax.plot(x,y1, marker = "p", color = "red", linestyle = "None")
    ax.plot(x,y2, marker = "p", color = "blue", linestyle = "None")
    ax.plot(x,y3, marker = "p", color = "black", linestyle = "None")    
    #more settings for the plot, color scale over n number of plots with index N in the list
    #plt.plot(x, y1, marker = "p", color = "r", x,y3, marker = "p", color = "b")
    #, color = "b" ,x,y3, marker = "p", color = "b"        
    #Settings for the plot  
    #Legend      
    ax.legend(loc = 4, fontsize = 14, fancybox = True, framealpha = 1)
    #axis, name, size of numbers
    ax.set_xlabel("Voltage ($\it{V}$)", fontsize = 14)
    ax.set_ylabel(label, fontsize = 14)
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    ax.set_ylim(plot_min, plot_max)
    
    #ax.set_ylim(10**-2,10**5)
    
    #y1 = abs(y1)
    #y2 = abs(y2)
    #y3 = abs(y3)
    #plt.semilogy()
    #tight layout prevents cutting of during pdf making 
    plt.tight_layout()
    #Change to statistics folder for saving statistic plots
    os.chdir(path_statistics)
    #Save function
    plt.savefig(name + ".pdf", format="pdf")
    #Go back to working directory 
    os.chdir(path)
    plt.show() 
    return 

def get_opimum_range_topo(dframe, column_name):
    #get optimum range out of statistic values df for imaging 
    #make only half of value, since average will be set to zero    
    listofvalues = dframe[column_name]
    print listofvalues
    print type(listofvalues)
    maximum = listofvalues.max()
    maximum = maximum/2
    maximum = round(maximum, -1)
    maximum = int(maximum) 
    print "Optimum +- range for all Images:"
    print maximum
    return maximum

def get_opimum_range_current(dframe, column_name1, column_name2):
    #get optimum range out of statistic values df for imaging   
    listofvalues1 = dframe[column_name1]
    listofvalues2 = dframe[column_name2]    
    maximum1 = int(listofvalues1.max())
    maximum2 = int(abs(listofvalues2.min()))    
    if maximum1 >= maximum2:
        m = maximum1
    else:
        m = maximum2
    m = round(m, -1)
    print "Optimum +- range for all Current Images:"
    print m
    return m

def get_opimum_range_error(dframe, column_name1, column_name2):
    #get optimum range out of statistic values df for imaging  
    listofvalues1 = dframe[column_name1]
    listofvalues2 = dframe[column_name2]    
    print listofvalues1
    print listofvalues2
    maximum1 = listofvalues1.max()
    maximum2 = abs(listofvalues2.min())    
    print maximum1
    print maximum2
    if maximum1 >= maximum2:
        m = maximum1
    else:
        m = maximum2
    m = round(m, 1)
    m = m
    print "Optimum +- range for all Error Images:"
    print m
    return m

def get_maximum_drift(df, column_name1, column_name2):
    #get the optimum range for the cut out of stable image 
    listofvaluesX = df[column_name1]
    listofvaluesY = df[column_name2]
    maxX = listofvaluesX.max()
    minX = abs(listofvaluesX.min())
    if maxX >= minX:
        total_drift_x = maxX 
    else:
        total_drift_x = -minX
    maxY = listofvaluesY.max()
    minY = abs(listofvaluesY.min())
    if maxY >= minY:
        total_drift_y = maxY 
    else:
        total_drift_y = -minY            
    return total_drift_x, total_drift_y
    
def get_line(df, factor, x_start, y_start, x_end, y_end, res):
    #Function takes the following values: 
    #scol, srow, ecol, erow, res, thickness, interpolation
    # Gives back a list of values, which are hight or current a.s.f. and the size of one datapoint in lateral resolution (nano meter)
    print "check"
    df_line = df.get_profile(x_start, y_start, x_end, y_end, res,5,2)
    #print "The real physical size of the extracted dataline"    
    df_line.multiply(10**factor)
    real_size = df_line.get_real()
    #print real_size
    line_res = df_line.get_res()
    #print "Line size"
    #print line_res
    #print (real_size/line_res)
    #line_point_size = (10**factor) * df_line.get_dx()
    line_point_size = (10**factor) *(real_size/line_res)
    print "The size of one idx in the numpy array /which is the line (in nm)"
    print line_point_size
    df_line_x = [] 
    for l in range(line_res):
        df_line_x.append(l*line_point_size)
    return df_line, df_line_x, line_point_size

def peak_extraction(i, x_data, y_data, step_size, peak_x_position, startingpoint, width): 
    #Function finds peak of segment in the dataline, set starting point bevor the peak starts 
    print "We are in the peak find function" 
    #Declare the y data from the gwy-Liste file and make it an array, then make all values positive
    y_all = y_data.get_data()    
    y = np.asarray(y_all)
    y = np.absolute(y)
    x = np.asarray(x_data)

    n_points = int(round(width/step_size))
    print "number of points of region we like to extract" 
    print n_points
    
 
    if i==0:
        idx = (np.abs(x-startingpoint)).argmin()
        print idx
        print x[idx]
        x = x[idx-n_points:idx+n_points]
        #print x 
        #x.flat[np.abs(a - 600).argmin()]
        y = y[idx-n_points:idx+n_points]
    else:
        startingpoint = x[peak_x_position] 
        idx = (np.abs(x-startingpoint)).argmin()
        print idx
        print x[idx]
        x = x[idx-n_points:idx+n_points]
    #print x 
    #x.flat[np.abs(a - 600).argmin()]
        y = y[idx-n_points:idx+n_points]

    peaks, _ = find_peaks(y,  distance = 1000)  
    #print peaks
    #print type(y)
    #print type(y)
    #plt.plot(y)
    #plt.plot(peaks, y[peaks], "x")
    #plt.show()
    print "Hellllooooo here comes the peak" 
    print peaks
    peaks = int(round(peaks))
    
    #print "X position of peak in nm"
    peak_x_position = int(round(x[peaks]))
    #print peak_x_position
    #print "x position of peak in int of small cut out area"
    #print peaks
    
    #print "Number of elements in array"
    #print len(y)
    #print len(y_all)

    print "Make it easy, just give backe x value of the peak and search for it later" 
    x_value = x[peaks]
    print x_value

    #y_out = y[(peaks-n_points):(peaks+n_points)]
    #x_out = x[(peaks-n_points):(peaks+n_points)]
    #plt.plot(y_out)
    #plt.plot(x_out)
    #plt.show()
    return peak_x_position, x_value

def Gauss(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def gauss_fit(x_data_peak, y_data, x_data, range_idx):
    #Fits gauss distribution to the segment an gives back maximum value, full half witdh 
    print "We are in the gauss-fit function" 
    #print type(y_data)    
    #print type(x_data)    
        
    y_all = np.asarray(y_data)
    x_all = np.asarray(x_data)
    peak_idx = (np.abs(x_all - x_data_peak)).argmin()
    #print peak_idx
    
    y = y_all[(peak_idx-range_idx):(peak_idx+range_idx)]
    x = x_all[(peak_idx-range_idx):(peak_idx+range_idx)]

    
    mean = np.sum(x * y)* (1. / np.sum(y))
    print "Mean and sigma of data"
    print mean
    sigma = np.sqrt(np.sum(y * (x - mean)**2) / np.sum(y))
    print sigma
    
    popt,pcov = curve_fit(Gauss, x, y, p0=[np.max(y), mean, sigma])
    print "The fit parameters for the gauss function" 
    print popt,pcov
    perr = np.sqrt(np.diag(pcov))
    print perr
    gauss_peak_y = popt[0]
    gauss_peak_x = popt[1]
    gauss_sigma = popt[2]
    
    print gauss_peak_y, gauss_peak_x, gauss_sigma
    
    # popt,pcov = curve_fit(Gauss, x, y, p0=[max(y), mean, sigma])
    fig, ax = plt.subplots(figsize=(10,5))
    plt.plot(x, y, "b+:", label="data")
    plt.plot(x, Gauss(x, *popt), "r-", label="fit")
    plt.legend()
    plt.title("Fig. 3 - Fit for Time Constant")
    plt.xlabel("Time (s)")
    plt.ylabel("Voltage (V)")
    
    
    ax.set_ylim(-3000,3000)
    
    
    #Change to linescan folder for saving statistic plots
    os.chdir(path_lines)
    #Save function
    plt.savefig(name + " " + "Gaussfit" + " " + "_.pdf", format="pdf")
    #Go back to working directory 
    os.chdir(path)
    plt.show() 
    #####Not ready yet, but somehow like this it should work 
    
    #curve_fit(gaus,x,y,p0=[1,mean,sigma])
    #return gauss_ax, gauss_mean, gauss_sigma
    return 


def plot_line(y_column, x_column, name, path_, plot_min, plot_max, xlabel, ylabel):
    #make multiplots with one x-axis and several y data
    print "We are in the statistics multi plot function!"
    #Making a plot of the statistic values of the evaluated scans
    fig, ax = plt.subplots(figsize=(10,5))        
    #determine x-data
    y = y_column.get_data()  
    x = x_column
    #setting color and mark style for graphs
    ax.plot(x, y, marker = "p", color = "red", label =str(extract_voltage(name)) + " V")
    #more settings for the plot, color scale over n number of plots with index N in the list
    #plt.plot(x, y1, marker = "p", color = "r", x,y3, marker = "p", color = "b")
    #, color = "b" ,x,y3, marker = "p", color = "b"        
    #Settings for the plot  
    #Legend      
    ax.legend(loc = 2, fontsize = 14, fancybox = True, framealpha = 1) 
    #axis, name, size of numbers
    ax.set_xlabel(xlabel, fontsize = 14)
    ax.set_ylabel(ylabel, fontsize = 14)
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    ax.set_ylim(plot_min, plot_max)
    #plt.semilogy()
    #tight layout prevents cutting of during pdf making 
    plt.tight_layout()
    #Change to folder for saving statistic plots
    os.chdir(path_)
    #Save function
    plt.savefig(name + "_Linescan.pdf", format="pdf")
    #Go back to working directory 
    os.chdir(path)
    plt.show() 
    return 

def append_to_linescan_list(df_line, name, daten):
    #Get gwy object df_line and make it a numpy array to add it to a list, so we can make multiple line plots later on 
    array = df_line.get_data() 
    #print "We are in the append line to linescan list function:"
    daten.append(array)
    print (type(daten))
    return daten

def save_line(line, line_x, path_line, path, name, parameter, name_x, name_y):
    #save line as csv data in the lines folder 
    print "cheeeeeeeeeeeck"
    print type(line)
    print line
    
    #determine x-data
    y = line.get_data()  
    x = line_x
    
    DF = pd.DataFrame()
    DF[name_x] = x
    DF[name_y] = y
    
    # save the dataframe as a csv file 
    os.chdir(path_lines)
    DF.to_csv(name + ".csv") 
    os.chdir(path)
    return 

def get_colormap(index, N):
    #Subfunction to get continuous color rgb-values for a multidata plot  
    #Define the scale for color from blue (neg for electrones) to red (pos) with N segmentations
    map = mpl.colors.LinearSegmentedColormap.from_list("custom", [(0,"blue"),(1,"red")], N=N)
    colorMap = plt.get_cmap(map)
    colorNorm  = colors.Normalize(vmin=0, vmax=N)
    scalarMap = cmx.ScalarMappable(norm=colorNorm, cmap=colorMap)
    #print scalarMap
    return scalarMap.to_rgba(index)

def zero_check(df, actcon, key):
    #Get the zero offset and set it zero, as well as the width of the gaussian noise around 0
    #First make histogram and find maximum y value (here is the 0 current)        
    y, x = np.histogram(df, bins = 40000)
    index_of_maximum = np.where(y == y.max())
    array_y_max = y, x[index_of_maximum]
    offset = array_y_max[1]
    print offset
    count = offset.size
    if count > 1: 
        offset = offset[0]
        print "We have double maximum count value, therefore we are chosing the first one:"
        print offset 
    print "Point of Maximum:"
    print offset
    #Add or Substrate to whole df 
    df.add(-offset)
    df.data_changed()
    #Make second similar df and make df absolut and fit a half gaussian distribution around it 
    df_fit = df.duplicate()
    #print type(df_fit)
    #df_fit = df_fit.multiply_fields(df_fit, df_fit)
    #df_fit = 
    y_fit, x_fit = np.histogram(df_fit, bins = 4000)  
    #print y_fit, x_fit
    
    mean, var = halfnorm.stats(moments="mv")
    #print mean, var   
    #Select datafield again, due to duplocate there are now 2 df inside the container 
    df = actcon[gwy.gwy_app_get_data_key_for_id(key)]
    gwy.gwy_app_data_browser_select_data_field(actcon, key)
    return df, offset

def make_histogram(df, name, fname, path_histo, path, factor, plot_max,  plot_min, label, xlabel):
    #Histogram zwischenfunktionen         
    xres = df.get_xres()
    print xres
    yres = df.get_yres()
    resolution = yres*xres
    print "The resolution of the scan is:"
    print resolution
    y, x = np.histogram(df, bins = resolution)
    #multiply the x value with the factor, important for topo scans in nm range 
    #print "We are in the histogram function, here comes the x value hopefully in nm:"
    x = x*(10**factor)
    #Settings for the plot    
    fig, ax = plt.subplots(figsize=(10,5))    
    #ax.scatter(x[:-1], y, "-b", label)
    ax.plot(x[:-1], y, "b", label = name)
    ax.legend(loc = 1, fancybox = True, framealpha = 1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Count")
    ax.set_xlim(plot_min, plot_max)
    plt.semilogy()
    plt.tight_layout()
    #before saving change into histo folder
    os.chdir(path_histo)
    #Save function
    plt.savefig(fname.replace(".tiff", "_Histo.pdf"), format="pdf")
    #Change back to working directory
    os.chdir(path)
    #plt.show()
    return 

def make_histogram_all(N, dataframes, path_histo, path, factor, plot_max, plot_min, name_ylabel, name_file, list_names):
    #Histogram-Plot von allen Histogrammen in einem Plot dazu durch alle eintraege in der Dataframes Liste gehen
    #Settings for the plot  
    fig, ax = plt.subplots(figsize=(10,4))
    #Function is looping through all dfs to make multi histo plot, its reversed so that first one is most non transparent 
    for i in range(N):
        #print "Index i:"
        zeile = dataframes[i]
        #Make Histogram bins with numpy         
        y, x = np.histogram(zeile[2], bins = 500)
        #Normierung on and off
        #y = (y/float(np.max(y)))
        #multiply the x value with the factor, important for topo scans in nm range 
        #print "We are in the histogram function, here comes the x value hopefully in nm:"
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

def print_df_NumpyArray(df):
    array = gwyutils.data_field_data_as_array(df)
    #gwy.gwy_app_data_browser_get_current(APP_DATA_VIEW)
    #print array
    return array

def get_drift(x, old_array, i):
    #Subfunktion to convolute two consecutive dataframes to find the drift 
    #function makes dataframes an array and saves df from the loop bevore as
    #old-array and compares them therefore
    #print "We are in the get drift function"
    #print "The start position is the half of the resolution (middle of image)"
    start_position = ((x.get_xres() / 2), (x.get_yres() / 2))
    
    if i == 0:
        #print "Firsttime, therefore print out array of first image:"
        old_array = print_df_NumpyArray(x) 
        corr_img = scipy.signal.fftconvolve(old_array, old_array[::-1,::-1], mode="same")
        #print corr_img
        start_position = np.unravel_index(np.argmax(corr_img), corr_img.shape)
        #print start_position
        #print type(start_position)
        drift = (0, 0)
        #print old_array
        return old_array, drift
    #have to give start_position and star-array back to make them global variables 
    
    else:
        #Also irgendwie muss man das so machen, das man das erste bild mit sich selbst convoluted und der unterschied zwischen dem        
        #punkt und dem punkt wenn ich das naechste mal zwei unterschiedliche mit sich convolute
        #dieser unterschied gibt mir den drift zwischen image 1 und 2 an 
        array_new = print_df_NumpyArray(x)
        #print "Here you get the convolution image as arrays:"
        #print type(array_new)
        #here the magic happens, we use the fastfouriertransformation to convolute
        #the two images and get an convolution image which we search for the 
        #highest/brightest point 
        corr_img = scipy.signal.fftconvolve(old_array, array_new[::-1,::-1], mode="same")
        #print corr_img
        position = np.unravel_index(np.argmax(corr_img), corr_img.shape)
        #print "Position of overlay center of first image with the following images:"
        #print position
        drift = start_position[0] - (position[0]), start_position[1] - (position[1])
        #print "Offset/drif of image from the image before (x/y):"
        #print drift
        #Here we return the data which is beeing treated right now as array_old (even though its been declared as array_new) 
        #so we can use it in the next round as array_old
        return array_new, drift
    
def get_offset_from_first_image(drift, offset, i):
    #Function to find the total offset from the first image in stack, normal 
    #iteration addition over i-loop 
    
    if i == 0: 
        #print "Type of offset and offset:"
        #print offset 
        #print type(offset) 
        return offset 
        
    else: 
        new_offset = ((offset[0] + drift[0]), (offset[1] + drift[1]))
        print "For iteration > 0 the offset from the image before"
        print new_offset
        return new_offset   
    
def make_drift_list(drift,offset_by_drift,i,dlist, path, path_statistics, dataframes, name):
    #Listing these calculated values into a list with header
    zwischendaten = i, offset_by_drift[0],offset_by_drift[1], drift[0],drift[1]
    dlist.append(zwischendaten)
    #For last run, make header and safe it as csv in statistics folgder 
    if i != len(dataframes) - 1:
        os.chdir(path_statistics)
        print("Save drift data....")    
        try:
            d = open(name, "w")
        except:
            print("Dateizugriff nicht erfolgreich")
            sys.exit(0)
        #Saving the data in a csv file
        d.write("Iteration" + ";" + "X Offset by drift from first image [px]" + ";" + "Y Offset by drift from first image [px]" + ";" + "X Drift from image before [px]" + ";" + "Y Drift from image before [px]" + "\n" )
        for zwischendaten in dlist:
            d.write(str(zwischendaten[0]) + ";" + str(zwischendaten[1]) + ";" + str(zwischendaten[2]) + ";" + str(zwischendaten[3]) + ";" + str(zwischendaten[4])+ "\n")
            #close the just written data an return         
        d.close()	    
        os.chdir(path)
        return dlist
    else:
        print "The drift list so far:"
        print dlist
        print 
        return dlist


def first_run_topo(voltage_list, lines_cutoff, fnames_topo, factor, daten):
    #First run makes: Extraction of DF from Actcon, Edges, Color, Fit, Statistics  
    for i in range(len(fnames_topo)):
        print "First Run Iteration: i="
        print i 
        fname = fnames_topo[i]
        name =  str(voltage_list[i]) + " V" + "_Topo"
        print name
        ids, actcon = load_data(path, fname)
        df = select_dataframe(actcon, parameter_name=name, key = 0)
        df = cut_edges(df, lines_cutoff)
        actcon = select_color(actcon, color="Gwyddion.net", key = 0)  
        actcon = fit_functions(actcon, df, i, mode = fit_mode)  
        statistics, max, min, Sme, zwischendaten, df =  get_statistic(df, actcon, fname, name, factor, daten, i, noise = 0, key = 0)        
        df = df_save(df, actcon, dataframes_topo, fname, name, voltage_list[i]) 
        df = average_check(df, actcon, key = 0)
        actcon = data_save(actcon, fname, path_saving = path_fitted, path_working = path)         
        remove(actcon)
    #After first run overview evaluation and collection of values for image making and so on            
    SaveStatisticsToFile(daten, path_statistics, voltage_list, name="Statistics_Topo.csv", unit="nm")
    df_stat_topo = load_statistics_data(path_statistics, name = "Statistics_Topo.csv")
    range_topo = get_opimum_range_topo(df_stat_topo, column_name = "Maximum [nm]")
    make_histogram_all(N, dataframes_topo, path_histo, path, factor, plot_max = range_topo, plot_min = -range_topo, name_ylabel="Topography ($\it{nm}$)", name_file="Topo_Histogramms.pdf", list_names = voltage_list)
    
    #the df lists consit of name and gwy objects 
    return df_stat_topo, range_topo, dataframes_topo

def second_run_topo(dataframes_topo, range_topo, path_fitted, array_old, offset_by_drift, factor):
    #Second run makes: average on 0, make and save: images-videos-histograms, convolution to get drift  
    #Set array_old for beginning, since its used by drift extraction function, will be overwritten in loop and reused
    #Make empty drift list (local) to fill in iteration loop
    drift_list = []
    for i in range(len(dataframes_topo)):
        print "Second Run Topo Iteration: i="
        print i 
        name, voltage, df, actcon, fname = dataframes_topo[i]
        ids, actcon = load_data(path_fitted, fname + ".gwy")

        df = average_check(df, actcon, key = 0)
        #actcon = select_range(actcon, range_topo*(10**-factor), key=0)
        actcon = image_save(actcon, i, path_pdfs, path, mode = image_mode, dataname= fname)
        
        df_line_topo, df_line_x, step_size = get_line(df, factor, x_start=line_x_start, y_start=line_y_start, x_end=line_x_end, y_end=line_y_end, res=line_res) 
        plot_line(df_line_topo, df_line_x, name, path_lines, plot_min=-range_topo, plot_max=range_topo, xlabel="Distance ($\it{nm}$)", ylabel="Topo ($\it{nm}$)")  
        data_lines_topo = append_to_linescan_list(df_line_topo, name, daten=datalines_topo)

        
        #make_histogram(df, name, fname, path_histo, path, factor, plot_max = range_topo, plot_min = -range_topo, label="Topo distribution", xlabel="Topography ($\it{nm}$)")
        ##########  Convolve    ##########
        array_old, drift = get_drift(df, array_old, i)
        offset_by_drift = get_offset_from_first_image(drift, offset_by_drift, i)
        drift_list = make_drift_list(drift,offset_by_drift,i,drift_list, path, path_statistics, dataframes = dataframes_topo, name="Drift_List.csv")
        df = df_save(df, actcon, dataframes_topo, fname, name, voltage_list[i]) 
        remove(actcon)

    df_drift = load_statistics_data(path_statistics, name = "Drift_List.csv")
    df_stat_topo = load_statistics_data(path_statistics, name = "Statistics_Topo.csv")
    make_statistics_plot_multi(df_stat_topo, 1, 4 ,6 ,5 , path_statistics, plot_max = (2*range_topo), plot_min = 0, label="Height ($\it{nm}$)", name = "Topo_MultiPlot")    
    make_statistics_plot(df_stat_topo, path_statistics, 13, [2, 3, 4,5,6], label="Topo")
    return df_drift, dataframes_topo[:(i+1)], data_lines_topo 

def third_run_topo(dataframes_topo, Range, path_fitted, path_stableframe, df_drift, factor = 9):
    #Here we extract from the fittted files a stable frame with the drift list
    maximum_drift_x, maximum_drift_y = get_maximum_drift(df_drift, column_name1="X Offset by drift from first image [px]", column_name2="Y Offset by drift from first image [px]")
    for i in range(len(dataframes_topo)):
        print "Third Run Topo Iteration: i="
        print i 
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
        print len(dataframes_topo)
    return 

def fourth_run_topo(dataframes_topo, range_topo, path_stableframe, array_old2, offset_by_drift, factor):
    #Foruth run, make convolution again, since it often is inacurate at first
    drift_list2 = []
    for i in range(len(dataframes_topo)):
        print "Fourth Run Topo Iteration: i="
        print i 
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
        print "Fourth2 Run Topo Iteration: i="
        print i 
        name, voltage, df, actcon, fname = dataframes_topo[i]
        ids, actcon = load_data(path_stableframe, fname + ".gwy")
        df = select_dataframe(actcon, parameter_name=voltage, key = 0)
        print df
        
        df_stable, key = area_extract(df, actcon, name, df_drift2, i, maximum_drift_x, maximum_drift_y, lines_cutoff)
       
        #df_stable= average_check(df_stable, actcon, key)
        #actcon = select_color(actcon, color="Gwyddion.net", key=key)
        #actcon = select_range(actcon, Range*(10**-factor), key)
        #actcon = image_save(actcon, i, path_stableframe, path, mode = image_mode, dataname= fname)
        actcon = data_save(actcon, fname, path_saving = path_stableframe, path_working = path)
        #df_stable = df_save(df_stable, actcon, dataframes_topo, fname, name, voltage_list[i]) 
        remove(actcon)
        print len(dataframes_topo)


    #iteration loop again, for even better convolution of images 
    drift_list3 = []
    array_old3 = 0
    offset_by_drift3 = (0,0)
    for i in range(len(dataframes_topo)):
        print "Fourth Run Topo Iteration: i="
        print i 
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


def make_videos(dataframes_topo, Range, path_stableframe, path_gifs):
    #make gifs, cut out stable frame 
    #make_gif(path_stableframe, time = 2, name = "Topography")
    make_gif(path_pdfs, time = 2, name = "Current")
    make_gif(path_pdfs, time = 2, name = "Topography")    
    make_gif(path_histo, time = 2, name = "Topography")
    make_gif(path_histo, time = 2, name = "Current")

    return

def first_run_current(voltage_list, lines_cutoff, fnames_current, factor, daten):
    #First run for current files makes: Extraction of DF from Actcon, Edges, Color, ZeroCheck, Statistics, Save after Fitting  
    for i in range(len(fnames_current)):
        print "Iteration: i="
        print i 
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
    print type(dataframes)
    print dataframes
    for i in range(len(dataframes)):
        print "Second Run Current Iteration: i="
        print i 
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
        print "Third Run Current Iteration: i="
        print i 
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
        print "Iteration of first run error: i="
        print i 
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
    print type(dataframes_error)
    print dataframes_error
    for i in range(len(dataframes_error)):
        print "Second Run Error Iteration: i="
        print i 
        name, voltage, df, actcon, fname = dataframes_error[i]
        ids, actcon = load_data(path, fname)
        
        df = select_dataframe(actcon, parameter_name=voltage_list[i], key=0)
        print df
        actcon = select_range(actcon, range_error*(10**-factor), key=0)
        actcon = image_save(actcon, i, path_pdfs, path, mode = image_mode, dataname= fname)
        #make_histogram(df, name, fname, path_histo, path, factor, plot_max = range_error, plot_min = -range_error, label="Error distribution", xlabel="Error ($\it{V}$)")

        remove(actcon)
    make_histogram_all(N, dataframes_error, path_histo, path, factor, plot_max = range_error, plot_min = -range_error, name_ylabel="Error ($\it{V}$)", name_file="Error_Histogramms.pdf", list_names = voltage_list)
    df_stat_error = load_statistics_data(path_statistics, name = "Statistics_Error.csv")
    make_statistics_plot_multi(df_stat_error, 1, 4 ,6 ,5 , path_statistics, plot_max = (range_error/2), plot_min = -range_error/2, label="Error ($\it{V}$)", name = "Error_MultiPlot")
    return

###############################################################################
    ### Presets: Make path, Sort files, Get Infos ... 
###############################################################################
  
path = working_path()
N, fnames_topo, fnames_current, fnames_amp, fnames_phase, fnames_error = sortandlist(path) 
voltage_list, lines_cutoff = get_info_sheet(path, name = "Infos.csv")
path_pdfs, path_histo, path_statistics, path_gifs, path_jpgs, path_lines, path_fitted, path_stableframe = make_folders(path)

###############################################################################
    ### Interatcion settings ...  
###############################################################################

image_mode = gwy.RUN_NONINTERACTIVE
fit_mode = gwy.RUN_NONINTERACTIVE 
mask_mode = gwy.RUN_NONINTERACTIVE


###############################################################################
    ### First evaluation run: open data, putting in lists, fitting, cutting edges, get statistics, get drift 
###############################################################################

df_stat_topo, range_topo, dataframes_topo = first_run_topo(voltage_list, lines_cutoff, fnames_topo, factor = 9, daten = statistics_topo)
df_drift, dataframes_topo, lines_topo = second_run_topo(dataframes_topo, range_topo, path_fitted, array_old, offset_by_drift, factor = 9)
#third_run_topo(dataframes_topo, range_topo, path_fitted, path_stableframe, df_drift, factor = 9)
#df_drift2, dataframes_topo = fourth_run_topo(dataframes_topo, range_topo, path_stableframe, array_old2, offset_by_drift, factor = 9)

#print df_drift2.iloc[:,1]
#print df_drift2.iloc[:,2]

#third_run_topo(dataframes_topo, range_topo, path_stableframe, path_stableframe, df_drift2, factor = 9)

#df_stat_error, range_error, dataframes_error = first_run_error(voltage_list, lines_cutoff, fnames_error, factor = 0, daten = statistics_error)
#second_run_error(dataframes_error, range_error, path, factor=0)

df_stat_current, range_current, dataframes_current = first_run_current(voltage_list, lines_cutoff, fnames_current, factor = 9, daten = statistics_current)
dataframes_current = second_run_current(dataframes_current, range_current, path, factor = 9)
#third_run_current(dataframes_current, range_current, path_fitted, path_stableframe, df_drift, factor = 9)


print "cheeeeeeck"


#plot_line_multi(df_line_topo_small, df_line_current_small, df_line_x, name, path_lines, plot_min_y1=topo_min + 50, plot_max_y1=topo_max ,plot_min_y2=current_min+5000, plot_max_y2=current_max-5000, xlabel="Distance ($\it{nm}$)", y1label="Topography ($\it{nm}$)", y2label="Current ($\it{nA}$)")


make_videos(dataframes_topo, range_topo, path_stableframe, path_gifs)


###############################################################################
    ### First run: fitting, cutting edges, get statistics, get drift 
###############################################################################

