"""AFM statistics, peak analysis, and drift calculations."""

import numpy as np
import gwy
import scipy.signal
from scipy.signal import find_peaks

from afm_processing import print_df_NumpyArray



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
