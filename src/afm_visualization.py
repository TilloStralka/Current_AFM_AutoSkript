"""AFM image, histogram, and statistical plotting helpers."""

import os

import numpy as np
import gwy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib as mpl



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
        range_topo (float): Maximum absolute value for scale
        key (int): Data field key
        
    Returns:
        actcon: Updated container with fixed range
    """
    actcon.set_double_by_name("/{}/base/max".format(key), range_topo)
    actcon.set_double_by_name("/{}/base/min".format(key), -range_topo)
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
    plt.savefig(os.path.join(path_statistics, "{}.pdf".format(name)), format="pdf")
    plt.show()
    return

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
