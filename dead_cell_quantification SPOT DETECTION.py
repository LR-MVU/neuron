# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 09:57:27 2022

@author: Administrator
"""

import numpy as np
import os
import bigfish.stack as stack
import bigfish.detection as detection
import bigfish.multistack as multistack
import bigfish.plot as plot
from bigfishUI import *
from matplotlib import pyplot as plt
import pandas as pd
import glob

# Ask user for the path to directory containing MAX projection images
path_input = getImagesDir()
path_output = os.path.join(path_input, "output")

# Get channel for spot detection (ex. CH3 or CH4)
_, max_files, chan = get_TIF_and_MAX_Images(path_input, "Enter name of channel for which spot detection is being done: ", show=True, return_chan=True)

# Folder path for detection results
detection_output = os.path.join(path_output, 'Spot detection results '+chan+" automatic")
# Check to see if folder path exists, if it does not create it 
if not os.path.exists(detection_output):
    os.makedirs(detection_output)

# Create file path for spreadsheet of thresholds, if it does not exist we will create this file
threshold_file = os.path.join(path_output, 'Spot detection results '+chan+" automatic" + "/thresholds.xlsx")
# Ask if user wants to use pre-set threshold values from spreadsheet or automatically detect threshold values
man_thresh = int(input("Type 0 if you want automatic spot detection. Otherwise, type 1.\n"))
elbow = getYesorNo("Do you want to see the elbow plot (number of mRNAs detected vs selected thresholds? (Y/N):")


# We will save the results of this into a file labelled for Automatic spot detection for this channel
img_filenames = []
automatic_thresholds = []
num_detected_spots = []

# Pre-set spot radius based on the estimated minimum size of a spot and voxel size based on the microscope used
spot_radius = 250
voxel_size= 322.5 

# Use values inputted in spreadsheet (manuel detection) if instructed by user
if man_thresh != 0:
    
    # Create spot detection folder if it does not already exist
    detection_output = os.path.join(path_output, 'Spot detection results '+chan+" manual")
    if not os.path.exists(detection_output):
        os.makedirs(detection_output)
        
    
    # Access spreadsheet with updated thresholds by user in the automatic detection folder
    thresholds_df = pd.read_excel(threshold_file)
    # The user must input the manual threshold values in the last coloumn of the spreadsheet
    man_thresh = list(thresholds_df[thresholds_df.columns[-1]])
    print(man_thresh)
    manually_detected_spots = []


for j in range(len(max_files)):
    
    # Get full path of the max file and corresponding TIF file
    filename = max_files[j].group(0)[:-4]
    img_filenames.append(filename)
    
    max_filename = max_files[j].group(0)
    rnaMAX_path = glob.glob(path_input + "/" + max_filename)[0]
    rna_mip = np.squeeze(plt.imread(rnaMAX_path))[:,:,-1]
    
    # If automatic thresholding is selected by the user
    if man_thresh == 0:
        # Show elbow plots only if user asked for it
        if elbow: 
            plot.plot_elbow(
                images=rna_mip, 
                voxel_size=voxel_size, 
                spot_radius=spot_radius,
                path_output=os.path.join(detection_output, filename+"_elbow"))
              
        # Detect spots. Bigfish will choose threshold based on voxel size and spot radius
        spots,threshold = detection.detect_spots(images=rna_mip, return_threshold=True, 
            voxel_size=voxel_size, 
            spot_radius=spot_radius)
        # Add thresholds to list of automatic thresholds
        automatic_thresholds.append(threshold)
        # Keep track of the number of spots per image
        num_detected_spots.append(len(spots))
    
    # Use manual detection values if instructed by users
    else:
        # Extract threshold value from list and give to Bigfish for spot detection
        threshold = man_thresh[j]
        spots = detection.detect_spots(images=rna_mip, threshold=threshold, 
            voxel_size = voxel_size,
            spot_radius=spot_radius)
        # Keep track of the number of spots per image
        manually_detected_spots.append(spots.shape)
        
        # Show elbow plots if instructed by user
        if elbow: 
            plot.plot_elbow(
                images=rna_mip, 
                voxel_size= voxel_size, 
                spot_radius=spot_radius,
                path_output=os.path.join(detection_output, filename+"_elbow"))
        
    # Print Spot Detection Information for the user
    print(f"Image {j}: {filename}: detected spots")
    print("\r shape: {0}".format(spots.shape))
    print("\r dtype: {0}".format(spots.dtype))
    print("\r threshold: {0} \n".format(threshold))

    
    # Plot detected spots
    path = os.path.join(detection_output, filename+"_spots")
    plot.plot_detection(rna_mip, spots, title=filename+" Detected spots", contrast=True, path_output=path)
    
    # Plot spots with decomposed dense regions
    path = os.path.join(detection_output, filename+"_detection")

    # Save in csv files    
    path = os.path.join(detection_output, filename+ "_spots.csv")
    stack.save_data_to_csv(spots, path)
    
    # Save mRNA mask 
    path2 = os.path.join(detection_output,'spot detection masks')
    if not os.path.exists(path2):
        os.makedirs(path2)
    mask = np.zeros(rna_mip.shape)
    
    coords = []
    for p in spots:
        x, y = int(p[0]), int(p[1])
        coords.append([x,y])
            
    for spot in coords:
        mask[spot[0],spot[1]] = 255
        
    plt.imsave(os.path.join(path2, filename+ "_mRNA_mask.gif"), mask.astype("uint8"))

# Save results of spot detection to threshold spreadsheet which can be edited to adjust manual values in the future 
if man_thresh == 0:
    df = pd.DataFrame({"Index": list(range(len(max_files))), 
                       "Filename": img_filenames,
                       "Automatic Threshold": automatic_thresholds,
                       "Num Spots detected": num_detected_spots})
    saveDfDictToExcel({"Thresholds by img name": df}, threshold_file)
    
else:
    df = thresholds_df
    new_threshold_file = os.path.join(detection_output, "_new_thresholds.xlsx")
    saveDfDictToExcel({"Thresholds by img name": df}, new_threshold_file)
    
    


    