# -*- coding: utf-8 -*-
"""
Created on Fri May  3 13:09:28 2024

@author: yorub
"""
import numpy as np
import matplotlib.pyplot as plt
import os 
from skimage.filters import threshold_minimum
from scipy.ndimage import binary_dilation
import pandas as pd
from skimage.segmentation import flood

import bigfish.stack as stack
import bigfish.detection as detection
import bigfish.plot as plot

# Function converts a dictionary into an excel spreadsheet
def saveDfDictToExcel(data_dict, name):
  with pd.ExcelWriter(name) as writer:
    for div_treat, df in data_dict.items():
      cols = df.columns.values.tolist()
      df.to_excel(writer, sheet_name=div_treat, columns=cols, index=False)

# Performs automatic spot detection, is not necessary if spot detection has been performed with other code
def spot_detection(twoD_image, filename, detection_output):
    # Minimum radius of a spot, can be detected based on image
    spot_radius = 300
    
    # Generate an elbow plot
    plot.plot_elbow(
        images=twoD_image, 
        voxel_size=322.5, 
        spot_radius=spot_radius,
        path_output=os.path.join(detection_output, filename+"_elbow"))
    
    # Perform automatic spot detection
    spots,threshold = detection.detect_spots(images=twoD_image, return_threshold=True, 
        voxel_size=322.5,
        spot_radius=spot_radius)
    
    # Plot detected spots
    path = os.path.join(detection_output, filename+"_spots")
    plot.plot_detection(twoD_image, spots, title=filename+" Detected spots", contrast=True, path_output=path)
    
    # Save detected spots in CSV
    path = os.path.join(detection_output, filename+ "_spots.csv")
    stack.save_data_to_csv(spots, path)
    
    # Create and save mRNA mask
    path2 = os.path.join(detection_output,'spot detection masks')
    if not os.path.exists(path2):
        os.makedirs(path2)
    mask = np.zeros(twoD_image.shape)
    
    coords = []
    for p in spots:
        x, y = int(p[0]), int(p[1])
        coords.append([x,y])
    for spot in coords:
        mask[spot[0],spot[1]] = 255
        
    plt.imsave(os.path.join(path2, filename+ "_mRNA_mask.gif"), mask.astype("uint8"))

    # Return the path to the spot detection results
    return path
    
    
def create_mask(spot_detection_file, img):
    # Read in the spot detection by BigFish
    spots = []
    with open(spot_detection_file, 'r') as f:
        for line in f:
            p = line.split(';')
            x, y = int(p[0]), int(p[1])
            spots.append((x,y))
    
    # Create array of zeros for binary mask
    mask = np.zeros(img.shape)
    # Create a copy of the image to be used in the flood fill algorithem
    image_alt = np.copy(img)

    # Loop through all spots detected for image during spot detection
    for s in spots:
        # Extract the spot coordinates
        spot = (s[0], s[1])
        
        # Check to see if the spot has already been accounted for in mask
        not_done = mask[spot] == 0
        if not_done:
            # Extract the amplitude of the spot in the image
            image_value = img[spot]
            # Set tolerance level for flood fill based on the brightness of the spot
            if image_value >= 210:
                # The flood will fill any pixel with at least 85% intensity of the spot
                tolerance = round(255 - image_value*0.65)
                
            elif image_value >= 100:
                tolerance = 255 - image_value*0.60
            
            elif image_value > 2:
                # This pixel intensity range is found in smaller mit in less bright areas.
                tolerance = 255 - image_value*0.65
            
            else:
                continue
            
            # Set the spot detected to full brightness for the binary mask
            image_alt[spot] = 255
            # Flood fill the area around the spot detected based on the the tolerance set above
            m = flood(image_alt, seed_point=spot, tolerance=tolerance)
            
            # If no flood occurs continue and skip the rest of the loop
            if np.sum(m) < 2:
                continue
            
            # Update flood fill image to original value to prevent errors in furture flood fill operations
            image_alt[spot] = image_value
            # Add flooded binary mask to previous binary mask updating new flooded values to positive integers
            mask += m
    # Returns a True/False mask where if the pixel has been flooded the value will be true
    return mask > 0 # binary mask

# Declare empty arrays to be filled and run the main program which will call functions above
image_name, pix_in_filtered_dapis, pix_in_filtered_mags, pix_in_overlap, frac_pix_in_overlap = [], [], [], [], []
def main_program(path_mag, path_dapi):
    # Get folder path of current directory
    top_directory = os.getcwd()
    # Extract the filename from the CH4 pathway
    fname = os.path.basename(path_mag).split("_CH4.tif")[0]
    # Create new folder where results can be placed for this specfic file run allowing for validation of masking
    folder_path = os.path.join(os.getcwd(), fname)
    os.makedirs(folder_path, exist_ok=True)
    # Move to newly created folder to place files
    os.chdir(folder_path)
    
    
    ## Read in images from inputed file paths
    mag_img = np.squeeze(plt.imread(path_mag))[:,:,0]
    dapi_img =  np.squeeze(plt.imread(path_dapi))[:,:,-1]
    # Save images in folder for verfication of good masking
    plt.imsave(fname+"_Dapi_8bit.gif", dapi_img)
    plt.imsave(fname+"_Magenta_8bit.gif", mag_img)
    
    # Change back to the top directory, exit sub folder for image
    os.chdir(top_directory)

    # Extract the spot detection results for the dapi and the mag channels (previously been calculated using the prior program)
    spot_detection_file_dapi = os.path.join(detection_output_dapi, os.path.basename(path_dapi)[:-4]+ "_spots.csv")
    spot_detection_file_mag = os.path.join(detection_output_mag, os.path.basename(path_mag)[:-4]+ "_spots.csv")
    
    # Get masks for dapi and magenta channels using bigfish spot detection and the create_mask function described above    
    dapi_binary_stain = create_mask(spot_detection_file_dapi, dapi_img)
    mag_binary_stain = create_mask(spot_detection_file_mag, mag_img)
    
    # Change back to the sub folder to save the initial binary masks for verification
    os.chdir(folder_path)
    plt.imsave(fname+"_Dapi_Binarized_Mask.gif", dapi_binary_stain)
    plt.imsave(fname+"_Magenta_Binarized_Mask.gif", mag_binary_stain)

    # Identify the brightest spots in the magenta image as progenitors which are not neurons so the can be removed 
    try:
        # The threshold has an adjustment factor of 0.85 of what is calculated using the function based on visual analysis of our images
        thresh = threshold_minimum(mag_img) * 0.85
    except RuntimeError:
        thresh = 255
        print(f"For {fname}, did not identify any proginators")
        
    # Find pixels in magenta image that are above the calculated threshold and dialate them 6 pixels to delete the progenitors
    binary_mag = mag_img > thresh
    binary_to_delete = binary_dilation(binary_mag, iterations=6)
    plt.imsave(fname+"_progenitors.gif", binary_mag)

    # Delete the spots by inverting the image and multiplying it with the dapi image to remove progenitors
    binary = binary_to_delete < 1 
    filtered_dapi = dapi_img*binary
    plt.imsave(fname+"_filetered_dapi.gif", filtered_dapi)


    # Create the binary masks excluding the progenitors for the dapi and magenta channels and determine the overlap
    dapi_binary_stain_filtered = dapi_binary_stain*binary
    mag_binary_stain_filtered = mag_binary_stain*binary
    overlap_binary = dapi_binary_stain_filtered*mag_binary_stain_filtered
    # Save the filtered masks for visual validation
    plt.imsave(fname+"_Dapi_Binarized_Mask_filtered.gif", dapi_binary_stain_filtered)
    plt.imsave(fname+"_Magenta_Binarized_Mask_filtered.gif", mag_binary_stain_filtered)
    plt.imsave(fname+"_Overlap_binary.gif", overlap_binary)

    # Add information to lists about the total area of nerve cells (excluding progenitors) in the dapi and magenta channels to list
    # Also include the percent of the dapi binary mask that was also present in magenta (what percent of neuron area died)
    image_name.append(os.path.basename(path_mag).split("_CH4.tif")[0])
    pix_in_filtered_dapis.append(np.sum(dapi_binary_stain_filtered))
    pix_in_filtered_mags.append(np.sum(mag_binary_stain_filtered))
    pix_in_overlap.append(np.sum(overlap_binary))
    frac_pix_in_overlap.append(pix_in_overlap[-1]/pix_in_filtered_dapis[-1])
    # Return to the top directory
    os.chdir(top_directory)



# This part of the code deals with processing many images at one time

# Get folder where user is storing there image files and enter that folder
folder_name = input(r"What is the path to the folder containing your images: ")
os.chdir(folder_name)
# Create a new folder for the results, we called ours "Results no dialation"
folder_results = "Results no dilation"
if not os.path.exists(folder_results):
    os.makedirs(folder_results)
# Chaneg directory to the newly created folder
os.chdir(folder_results)

# Define path for spot detection results, this may have to be adjusted based on where spot detection results are
detection_output_mag = os.path.join(folder_name, "output", "Spot detection results CH4 manual")
detection_output_dapi = os.path.join(folder_name, "output", "Spot detection results CH3 manual")

# Create a new folder for the detection dapi output
if not os.path.exists(detection_output_dapi):
    os.makedirs(detection_output_dapi)
    
# Create a list of all images file names a subdirectories
image_files = os.listdir(folder_name)
# Create empty lists for dapi and magenta file names
mag_file_names, dapi_file_names = [], []
# Loop through all files in the folder and sort those that are dapi and magenta (should match pairwise) and ignore other files in folder
for file in image_files:
    if "CH4.tif" in file:
        mag_file_names += [file]
    elif "CH3.tif" in file:
        dapi_file_names += [file]

# Loop through the file names and run the program for pairwise dapi and magenta images by giving the main program there paths (described above)
for file_mag in mag_file_names:
    for file_dapi in dapi_file_names:
        if file_mag.split("CH4.tif")[0] == file_dapi.split("CH3.tif")[0]:
            path_mag = os.path.join(folder_name, file_mag)
            path_dapi = os.path.join(folder_name, file_dapi)
            main_program(path_mag,path_dapi)

# Create quantification dataframe based on lists created from running the main program
df = pd.DataFrame({"Image Name": image_name,
                    "Num. Pixels in Filtered Dapi image": pix_in_filtered_dapis,
                    "Num. Pixels in Filtered Magenta image": pix_in_filtered_mags, 
                    "Num. Pixels in Overlap of Dapi and Magenta": pix_in_overlap, 
                    "Frac. Pixels in Filtered Dapi in Overlap": frac_pix_in_overlap})

# Convert dataframe created above to a spreadsheet 
saveDfDictToExcel({"Data": df}, "CellDeath.xlsx")


