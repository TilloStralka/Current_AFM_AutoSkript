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
    │   ├── __init__.py    <- Makes `src` a Python module.
    │   └── utils_afm.py   <- Contains all the separate defined functions deployed 
    │                         in the main file: `Current_Auto.py`.
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, 
    │                         e.g., generated with `pip freeze > requirements.txt`.
    │
    ├── Fitted             <- The fitted files (adapted topography and current) are saved here as .gwy files. 
    │
    ├── Histograms         <- Histograms of the distributions of the scans are stored here.
    │
    ├── JPGs               <- Scans are saved as JPGs, which are required for the GIF-making process.
    │
    ├── Linescans          <- Contains:
    │                         - Plots of linescans.
    │                         - TXT files with distance over measurement data (e.g., height or current intensity).
    │
    ├── PDFs               <- Final images of the scans are saved here in PDF format.
    │
    ├── gifs               <- Created videos/GIFs of the scan sequences are saved here.
    │
    ├── StableFrame        <- Contains the cut-out images of the unmoved areas, extracted via convolution.
    │
    ├── Statistics         <- CSV files containing all statistical values for topography, current, and more.

--------

==============================

## **3. Data**

==============================

The data is was collected by the author in the lab and is  **not included in the repository** due to size and privacy constraints. 
The files should be AFM scans or current AFM scans (I-AFM, cp-AFM, cAFM etc.) of the format gwy containing differnet channels (topography, error signal, laterals force, current, etc.). These channels will be treated separately, correlated and evaluated. Within the data file there should also be a .txt file. Which contains a list of voltages applied, or paramter changes during the scan sequence recording/measuring. For further informations see def get_info_sheet in the utils_afm.py file.


==============================

## 4. Sequential Electrical AFM Measurements

==============================

### Measurement Procedure
- Conduct cp-AFM, KPFM, and SCM measurements while varying parameters:
  - Voltage
  - Temperature
  - Time
  - Light
  - Position
  - Tip-sample distance
  - Tip pressure
- Record two stacks of 3D data:
  - Topography
  - Electrical measurements

### Data Region Classification
- Separate topography images into regions:
  - Crystal facets
  - Grain sides
  - Grain boundaries (GBs)
  - Defect-dense regions
- Create a 2D mask to map these regions onto electrical signals for correlation.

### Behavioral Analysis
- Observe physical property changes across varying parameters.
- Generate statistical insights surpassing point-probe spectroscopy.

### Crystallographic Classification
- Identify crystallographic facets using angular relationships.
- Confirm orientation via XRD (e.g., (111) out-of-plane γ-CuI).

### Region Segmentation
- For small areas:
  - Identify GB regions with slopes exceeding 25%.
- For larger areas:
  - Use watershed or height cut-off methods for grain classification and thin-film characterization.

> **Citing from the dissertation: Chapter 3.2.4 Sequential AFM and correlation with masks**  
![Image Processing](https://github.com/user-attachments/assets/1584d244-4594-42fc-a660-4563bc2de305)

---

==============================

## 5. Execution Flow of Program and Data Evaluation

==============================

### Data Preprocessing
- Locate and load measurement files and corresponding information sheets.
- Combine separate signal files (e.g., topography, current, slope) into container files for correlation.

### First Run
- Process data channels (e.g., topography, slope) for:
  - Distortion correction
  - Plane adjustment
  - Artifact removal
- Generate outputs:
  - Images and histograms
  - Basic statistics (e.g., max, min, roughness)

### Second Run
- Create masks:
  - Use slope thresholds to identify GB regions.
- Segregate areas into GB and grain regions.
- Save masked data and additional statistics (e.g., GB/grain proportions).

### Third Run
- Apply masks to electrical channels for grain/GB-specific analysis.
- Calculate:
  - Region-specific average currents
  - Current proportions
  - Noise thresholds
- Save outputs:
  - Separated and processed data files
  - Final statistics
  - Visual summaries (plots and videos)

### Parameter-Dependent Analysis
- Plot and analyze statistics (e.g., currents, region proportions) as functions of varying parameters.
- Combine results into holistic visualizations for clear presentation.

> **Citing from the dissertation: Chapter 5.5 The deployed code and evaluation approach**  
> _Fill in text as needed._  

![Data Processing](https://github.com/user-attachments/assets/03cb276f-c5f9-4c4b-a7db-3a7f71e7abc1)

---

==============================

## Technologies Used
- Python
- gwyddion / pygwy
- Pandas
- Matplotlib

---