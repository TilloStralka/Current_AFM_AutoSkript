# Current AFM automatic evaluationd and correlation skript

## Overview

Python script for autonomous Current AFM scans evaluation. It will treat a number of scans simultaneously and correlate topographic and current data and make a statistic evaluation of the whole stack of images.

""" Writer: Tillmann Stralka 2020.Mai.05 Owner: HLP University of Leipzig About: Python version 2.7 Automatic procession of AFM files. Finding, listing, fitfunctions, image making, statistic extraction and plotting Correlation of Current and topographic information, convolution of dynamicscans """

The script was created as part of Tillmann Stralka's doctoral thesis. Please refer to it for further questions and information: "Electrical Atomic Force Microscopy Applied on Copper Iodide Thin Films and Crystals - 17.06.2024 - Tillmann Stralka" [Dissertation](https://nbn-resolving.org/urn:nbn:de:bsz:15-qucosa2-927428)


## **1. How to Run the Project**
It is recommended to work with a venv, since this project works with the outdated Python 2.7 version. 
==============================

## Installation
1. Clone the repository:
    ```bash
    git clone https://github.com/TilloStralka/Current_AFM_AutoSkript.git
    cd CO2_Emission_Predictor
    ```

2. Create a virtual environment and activate it:
    ```bash
    python -m venv env
    source env/bin/activate  # On Windows: env\Scripts\activate
    ```

3. Install dependencies:
    ```bash
    pip install -r requirements.txt
    ```
4. Run the main analysis program:
    ```bash
    python Current_Auto.py
    ```

5. Download and Install Gwyddion
    The Gwyddion software API will be used to perform specific functions in the script. This script is written for **Gwyddion version 3.52**. 
    You can download the software and find additional information here: [Gwyddion Download Page](https://gwyddion.net/download/).

    **Steps to Set Up Gwyddion:**
    a) Locate the Installation Path:
    After installing Gwyddion, identify the location of the `pygwy` library on your system. Typically, this is located in a directory such as:  `/usr/share/gwyddion/pygwy`


    b) Bind Gwyddion Libraries:
    To enable your script to use the Gwyddion API, append the library path to Python’s system path in the script. Add the following line to your main script:  
    ```python
    sys.path.append('/usr/share/gwyddion/pygwy')
    ```

    c) Import Required Libraries:
    After binding the library path, you can import the required Gwyddion modules in your script:
    ```python
    import gwy
    import pygtk
    import gwyutils
    ```
     After that the required libraries: gwy, pygtk, gwyutils can be used. 

==============================

## **2. Project Organization**

==============================

    ├── LICENSE
    ├── README.md          <- The top-level README for developers using this project.
    ├── data               <- Should be in your computer but not on Github (only in .gitignore)
    │                         The original, immutable data dump.
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's name, and a short `-` delimited description, e.g.
    │                         `1.0-alban-data-exploration`. (Non yet, because this runs on python 2.7)
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`           
    │
    ├── src                <- Source code for this project.
    │   ├── __init__.py    <- Makes src a Python module.
    │   └── utils_afm.py   <- This contains all the separate defined functions which are deployed 
    │                         in the main file: Current_Auto.py.
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, 
    │                         e.g.generated with `pip freeze > requirements.txt`                          
    │
    ├── Fitted             <- The fitted files (adapted topography and current) are saved here as .gwy files.
    │
    ├── Histograms         <- Here the histograms of the distributions of the scans will be stored.
    │
    ├── JPGs               <- Here the scans will be saved as JPGs, which are required for the gif 
    │                         making process. 
    │
    ├── Linescans          <- Here the Linescans are saved (the plots as well as the txt files containing
    │                         distance over measurement, which could be height or current intensity) 
    │
    ├── PDFs               <- Here the final images of the scan are saved as PDF formats.
    │
    ├── gifs               <- Here the created videos/gifs of the scan sequences are saved. 
    │
    ├── StableFrame        <- Here do lie the cut out images of the unmoved are, which is extracted via
    │                         convolution. 
    │
    ├── Statistics         <- Here lies the csv files which contain all statistic values for topography, and
    │                         current and furthermore. 

--------

==============================

## **2. Data**

==============================

The data is was collected by the author in the lab and is  **not included in the repository** due to size and privacy constraints. 
The files should be AFM scans or current AFM scans (I-AFM, cp-AFM, cAFM etc.) of the format gwy containing differnet channels (topography, error signal, laterals force, current, etc.). These channels will be treated separately, correlated and evaluated. Within the data file there should also be a .txt file. Which contains a list of voltages applied, or paramter changes during the scan sequence recording/measuring. For further informations see def get_info_sheet in the utils_afm.py file.

==============================

## **2. Programm Steps**

==============================


1. **Data Collection**: 
   - Collect data from the path
   - Unpack and 
   - bla 

2. **Data Organization**:
   - 
   
3. **Data Visualization**:
   - Create visualizations to show dependencies between brands, specifications, and fuel types.
   
## Technologies Used

- Python
- gwyddion /pygwy
- Pandas
- Matplotlib / Seaborn for visualization



Citing from the dissertation: 

Chapter 3.2.4 Sequential AFM and correlation with masks
#### fill in text #### 
<img width="858" alt="ImageProcessing" src="https://github.com/user-attachments/assets/1584d244-4594-42fc-a660-4563bc2de305">

Chapter 5.5 The deployed code and evaluation approach

#### fill in text #### 
<img width="489" alt="DataProcessing" src="https://github.com/user-attachments/assets/03cb276f-c5f9-4c4b-a7db-3a7f71e7abc1">

