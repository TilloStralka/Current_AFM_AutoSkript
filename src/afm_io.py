"""File and metadata helpers for AFM analysis."""

import os
import shutil
import sys

import pandas as pd


def get_info_sheet(path, name):
    #Function to load the info data for all those scans, which must include:
    #voltage list, number of lines which should be cut away.
    os.chdir(path)
    df_info = pd.read_csv(name, sep=";", header=[0])
    print(df_info)
    print("Cheeeeeeeeeeeeeeeeeeeeeeeeeeeck")
    voltage = df_info["Voltage [V]"].tolist()
    print("Cheeeeeeeeeeeeeeeeeeeeeeeeeeeck")
    lines_cutoff = int(df_info.iloc[0]["Lines Cutoff"])
    os.chdir(path)
    return voltage, lines_cutoff


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
            return float(element.replace("V", ""))
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
    list_of_filenames = sorted(list_of_filenames, key=extract_time)
    list_TIME = []

    for filename in list_of_filenames:
        ending = extract_time(filename)
        if ending:
            # Take first 9 chars of timestamp
            timestamp = ending[0:9]
            list_TIME.append(timestamp)

    return list_TIME


def sortandlist(path):
    names = os.listdir(path)
    files_topo = [k for k in names if "Topography" in k]
    files_current = [k for k in names if "Current" in k]
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

def make_folders(path):
    #Working path in which the pdfs will be saved 
    try:
        os.makedirs(os.path.join(path, "PDFs"))
    except OSError:
        print("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(os.path.join(path, "PDFs"))
        os.makedirs(os.path.join(path, "PDFs"))
    else:
        print("Successfully made new folder") 
    #Working path in which the jpgs will be saved 
    try:
        os.makedirs(os.path.join(path, "gifs"))
    except OSError:
        print("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(os.path.join(path, "gifs"))
        os.makedirs(os.path.join(path, "gifs"))
    else:
        print("Successfully made new folder") 
    #Working path in which the Histogramms will be saved 
    try:
        os.makedirs(os.path.join(path, "Histograms"))
    except OSError:
        print("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(os.path.join(path, "Histograms"))
        os.makedirs(os.path.join(path, "Histograms"))
    else:
        print("Successfully made new folder") 
    #Working path in which the Linescans will be saved 
    try:
        os.makedirs(os.path.join(path, "Linescans"))
    except OSError:
        print("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(os.path.join(path, "Linescans"))
        os.makedirs(os.path.join(path, "Linescans"))
    else:
        print("Successfully made new folder")         
    #Working path in which the Statistics csv plus Plots will be saved 
    try:
        os.makedirs(os.path.join(path, "Statistics"))
    except OSError:
        print("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(os.path.join(path, "Statistics"))
        os.makedirs(os.path.join(path, "Statistics"))
    else:
        print("Successfully deleted and made new")     
        #Working path in which the fitted data will be stored  
    try:
        os.makedirs(os.path.join(path, "Fitted"))
    except OSError:
        print("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(os.path.join(path, "Fitted"))
        os.makedirs(os.path.join(path, "Fitted"))
    else:
        print("Successfully deleted and made new")     
            #Working path in which the stable frame data will be stored  
    try:
        os.makedirs(os.path.join(path, "StableFrame"))
    except OSError:
        print("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(os.path.join(path, "StableFrame"))
        os.makedirs(os.path.join(path, "StableFrame"))
    else:
        print("Successfully deleted and made new")      
                #Working path in which the jpgs will be stored  
    try:
        os.makedirs(os.path.join(path, "JPGs"))
    except OSError:
        print("Folder exists already, will be deleted and replaced by new one")
        shutil.rmtree(os.path.join(path, "JPGs"))
        os.makedirs(os.path.join(path, "JPGs"))
    else:
        print("Successfully deleted and made new")      
    path_pdfs = os.path.join(path, "PDFs")
    path_histo = os.path.join(path, "Histograms")
    path_statistics = os.path.join(path, "Statistics")
    path_gifs = os.path.join(path, "gifs")
    path_lines = os.path.join(path, "Linescans")
    path_fitted = os.path.join(path, "Fitted")
    path_stableframe = os.path.join(path, "StableFrame")
    path_jpgs = os.path.join(path, "JPGs")
    return path_pdfs, path_histo, path_statistics, path_gifs, path_jpgs, path_lines, path_fitted, path_stableframe

def load_data(path, name):
    import gwy
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
    import gwy
    #Save files as a new .gwy file with a _fitted ending
    os.chdir(path_saving)
    #gwy.gwy_file_save(actcon, fname)
    
    gwy.gwy_file_save(actcon, fname + ".gwy")
    os.chdir(path_working)
    return actcon

def remove(actcon):
    import gwy
    #Remove file from browser to prevent program crashes due to overload
    gwy.gwy_app_data_browser_remove(actcon)

def SaveStatisticsToFile(daten, path_statistics, voltage_list, name, unit):
    #Saving of all statistic data which is in the data list into a csv data 
    #Also adds the voltage list to the csv data 
    print("Save statistics data....")    
    # Write the statistics file to the output folder.
    try:
        d = open(os.path.join(path_statistics, name), "w")
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
    return

def load_statistics_data(path_statistics, name):
    #Load statistics CSV data into a DataFrame.
    #Read function 
    df_stat = pd.read_csv(os.path.join(path_statistics, name), sep=";", header=[0])
    print("Loading of csv data worked!")
    return df_stat

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
    
    # Save the complete drift list on the final iteration.
    if i == len(dataframes) - 1:
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
