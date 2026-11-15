"""AFM data selection, preprocessing, and transformations."""

import numpy as np
import gwy
import gwyutils
from scipy.stats import halfnorm



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

def print_df_NumpyArray(df):
    """
    Convert a Gwyddion DataField object to a numpy array.
    
    Args:
        df: Gwyddion DataField object
        
    Returns:
        numpy.ndarray: Array containing the DataField values
    """
    return gwyutils.data_field_data_as_array(df)
