__author__ = "Zoe Wefers"
__credits__ = ["Ryan Huang"]
__version__ = "1.0.0"
__maintainer__ = "Zoe Wefers"
__email__ = "zoe.wefers@mail.mcgill.ca"

"""Purpose: Calculates four kinds of statistics about the localization of mRNA in dendrites/somas. For mRNA colocalization 
and synaptic localization performs 100 simulations of placing mRNAs and reports simulated statistics as a computational control."""

from re import L
import numpy as np
import pandas as pd
import glob
import os
from math import sqrt, inf, ceil
import random
from collections import Counter
from collections import defaultdict
import matplotlib.pyplot as plt

zones = ['0-25', '25-50','50-75','75-100','100-125','125-150','150-175','175-200','200-225','225-250', '250-inf']

#Find pixels in skel with only 1 neighboor
def findEndpoints(skel, points):
  endpoints = []
  for point in points:
    if np.sum((skel[(point[0] - 1):(point[0] + 2), (point[1] - 1):(point[1] + 2)])) == 510: #255*2 (pixel and its 1 neighboor)
      endpoints.append(point)
  return endpoints


def sortDistance(dist):
    '''    
    Parameters
    ----------
    dist : num
        Input a distance (um) and the function will sort it into a zone of 25 um each from 0 to 250.

    Returns
    -------
    zone : str
        a description of which zone the given distance falls into.

    '''
    zones = ['0-25', '25-50','50-75','75-100','100-125','125-150','150-175','175-200','200-225','225-250', '250-inf']
    
    if 0 <= dist < 25:
        zone = zones[0]   
    elif 25 <= dist < 50:
        zone = zones[1]       
    elif 50 <= dist < 75:
        zone = zones[2]       
    elif 75 <= dist < 100:
        zone = zones[3]       
    elif 100 <= dist < 125:
        zone = zones[4]       
    elif 125 <= dist < 150:
        zone = zones[5] 
    elif 150 <= dist < 175:
        zone = zones[6]   
    elif 175 <= dist < 200:
        zone = zones[7]      
    elif 200 <= dist < 225:
        zone = zones[8]    
    elif 225 <= dist < 250:
        zone = zones[9]    
    elif 250 <= dist:
        zone = zones[10]
    else:
        zone = "neg" #for negative numbers
        print("Negative!!! why?")
        print(dist)
    return zone


def calculate_densities(segs_byDivTreat, syn_chan, nanometers_per_pixel):
  densities_byDivTreat = {}
  
  for div_treat in segs_byDivTreat.keys():                                                                                                                                                                                                                                                                                                                                                                                                                             
    numRNA_area_density = [] #array of tuples (segmentation #, Channel, num of mRNA, area, density)
    
    for segmentation in segs_byDivTreat[div_treat]:
      print_points = np.argwhere(segmentation.print == 255)
      area = len(print_points) * nanometers_per_pixel**2
      
      if area != 0:
        for channel in segmentation.spot_coords.keys():
          if channel != syn_chan:
            num_mRNAs = len(segmentation.spot_coords[channel])
            numRNA_area_density.append([segmentation.num, channel, num_mRNAs, area, num_mRNAs/area])
      else:
        print("Print for " + segmentation.name() + " is empty")
        
    df = pd.DataFrame(numRNA_area_density, columns = ["Number", "Channel", "Number of mRNA", "Area (sq. nanometer)", "Density"])
    densities_byDivTreat[div_treat] = df
    
  return densities_byDivTreat



def dist(pt, sig):
  return sqrt((pt[0] - sig[0])**2 + (pt[1] - sig[1])**2)


#map a mRNA coordinate to the closest point in a list of pixel coordinates
def mapSig(sig, points):
  k = np.argmin(np.array([dist(pt, sig) for pt in points]))
  return points[k]



#Order coordinates of skeleton in order from one end to the other
def orderSkelPoints(skel, skel_points, endpoint, nm_per_pixel=107.5):
  dist_array = np.zeros(len(skel_points))
  skel_array = np.zeros(skel_points.shape) #an array of all the skel points, but in order from the endpoint
  prev_coord = endpoint
  cur_coord = endpoint
  dist_list_um = []
  
  #iterate through the number of pixels in the skeleton
  for i in range(len(skel_points)):
    skel_array[i] = cur_coord
    dist_array[i] = dist(prev_coord, cur_coord)
    dist_list_um.append(np.sum(dist_array)*nm_per_pixel/1000)
    
    #look for white pixels adjacent/diagonal to current pixel
    pixels = np.argwhere(skel[cur_coord[0] - 1:cur_coord[0] + 2, cur_coord[1] - 1:cur_coord[1] + 2] == 255)
    pixels = [(cur_coord - (1,1) + pix) for pix in pixels]

    cur_coord_copy = cur_coord
    prev_coord_copy = prev_coord
    
    for pix in pixels:
      if np.any(pix != cur_coord_copy) and np.any(pix != prev_coord_copy):
        prev_coord = cur_coord
        cur_coord = pix
  
  
  
  
  return skel_array, dist_array, dist_list_um #returns skel_points in order and an array of distance between adjacent points in skeleton




#Assume soma is the end of skeletons with more mRNA coordinates
def findSoma(distances, skel_len, nm_to_pixels):
    '''
    Reorders a list of distances if needed, to make sure the distances are measured from the soma end
    Parameters
    ----------
    distances : list of num
        a list of distances (in nm) of some spots. 
    skel_len : num
        length of skeleton in pixels.    
        can be found by getting sum of the list of distances between each point in the skeleton (in pixels)
        (the list would have each value be either 1 or sqrt(2))
    nm_to_pixels : float
        conversion factor: 1 pixel = nm_to_pixels nanometer. this is not a good name, should be nm_per_pixel

    Returns
    -------
    distances : list of num
        a list of distances (nm) of the points, reordered if necessary.

    '''
    skel_len_um = skel_len * nm_to_pixels / 1000 #convert skel to nm, then from nm to um. 1 nm = 1/1000 um
    
    sum_mrna_dists = np.sum(distances)
    #what is this criteria? what does it mean?
    
    if 2*sum_mrna_dists > skel_len_um * len(distances): #if initial guess was not soma end
        
        ret_distances = [] #take distance starting from oposite end of skel
        for dist in distances:
            new_dist =skel_len_um-dist 
            
            if -0.01 <= new_dist < 0: #due to computational or memory limits, points that are mapped to an endpoint of a skeleton will sometimes be off by a tiny amount
                #e.g. if total skel_length is 48.972831283, the distance of the mapped point might be 48.972831288 or something
                #This causes new_dist to be a negative number of order of magnitude e-13 or e-14
                #if it is a tiny negative value, we should consider it to be essentially zero
                print("Neg number:", new_dist)
                ret_distances.append(0)
                
            else:
                ret_distances.append(new_dist)
    else: #if initial guess was correct, change nothing
        ret_distances = distances
        
    return ret_distances

def calculateSimSpotsDists(dendrite, nm_to_pixels):
    '''
    Takes in a single Dendrite object, goes through the simulated spots.
    For each of the 100 simulations, sorts the sim spots into zones and saves in dendrite.simulated_spot_distances attribute
    Needed for mRNA perspective of synAnalysis, as simulated and real mRNA need to be 
    Parameters
    ----------
    dendrite : TYPE
        DESCRIPTION.
    nm_to_pixels : TYPE
        DESCRIPTION.

    Raises
    ------
    ValueError
        If dendrite.start_end is None. This means that distrAnalysis has not been done before.

    Returns
    -------
    None.

    '''
    for chan, simulations in dendrite.simulated_spots.items(): #simulations is the list of 100 tuples (sim_coords, sim_coords_px)
        total_dists = []
        for simulation in simulations:
            spots = simulation[1] # spots in pixels NOT nm
            #spots = [(x/nm_to_pixels, y/nm_to_pixels) for x, y in sim_coords]
            
            
            skeleton = dendrite.skeleton
            skel_points = np.argwhere(skeleton == 255)
            
            endpoints = findEndpoints(skeleton, skel_points) #find endpoints of skeleton
            if len(endpoints) == 2:
              if dendrite.start_end is not None:
                  startEnd = dendrite.start_end
              else:
                  raise ValueError("No endpoint known for this dendrite: Dend. "+ str(dendrite.num)+ " of " +str(dendrite.img_name))
              
              skel_points, skel_distances, skel_distances_um = orderSkelPoints(skeleton, skel_points, startEnd)
              
              #map each simulated spot to the skeleton
              mapped_spots = [mapSig(spot, skel_points) for spot in spots] 
              
              spot_distances = [] #calculate all the distances of the mrna along the skeleton
              #calculate distance of mappen mrna coord along skeleton
              for spot in mapped_spots:
                k = np.argwhere(np.all(skel_points == spot, axis=1))[0][0]
                dist = skel_distances_um[k]
                spot_distances.append(dist)
            
              #since this is for simulated spots, no need to check if the endpoint was correct.
              total_dists.append(spot_distances)
              
        dendrite.simulated_spot_distances[chan] = total_dists
              

def distrAnalysis(dendrites, syn_chan, nm_to_pixels):
  '''
  Goews through each dendrite and calculates the number of spots in each zone. 
  The zones are 25 um each from 0 to 250 um from the soma. 
  For each Dendrite, also initializes attributes dendrite.spot_distribution and dendrite.spot_distances
  
  Parameters
  ---------
  dendrites: list of Dendrite objects
  the distribution analysis will be performed for each dendrite individually
  
  syn_chan: str
  name of synapse channel. this is specified to avoid using synapse as an indicator of where the start end is. 
  
  nm_to_pixels : float
      conversion factor: 1 pixel = nm_to_pixels nanometer. this is not a good name, should be nm_per_pixel

  Returns
  -------
  df: Pandas DataFrame
  a breakdown of each dendrite, with each channel, and how many mrnas are in 0-25 um from the soma, 25-50 um, etc. until 150 um
  '''


  data = []
  
  for dendrite in dendrites:
    
    channels = list(dendrite.spot_coords.keys())
    if channels[0] == syn_chan:
        channels = channels[1:] + [syn_chan] #we do not want syn chan to be the first channel processed
        #bc it is not a good indicator of where the soma is 
        print("moving synchan")
        
    for channel in channels:
        
        spots = [(x/nm_to_pixels, y/nm_to_pixels) for x, y in dendrite.spot_coords[channel]]
        
        skeleton = dendrite.skeleton
        skel_points = np.argwhere(skeleton == 255)
        
        endpoints = findEndpoints(skeleton, skel_points) #find endpoints of skeleton
        
        if len(endpoints) == 2:
            
          if dendrite.start_end is not None:
              startEnd = dendrite.start_end
          else:
              startEnd = endpoints[0]
          
          
          skel_points, skel_distances, skel_distances_um = orderSkelPoints(skeleton, skel_points, startEnd)
          dendrite.skel_points = skel_points #list of skeleton coordinates, in order
          dendrite.skel_distances = skel_distances #distance bw skeleton coordinates, in order
          
          
          
          mapped_spots = [mapSig(spot, skel_points) for spot in spots] #map mrna coords to skeleton
          
          spot_distances_first = [] #first calculate all the distances of the mrna along the skeleton
    
          #calculate distance of mapped mrna coord along skeleton
          for spot in mapped_spots:
              
            k = np.argwhere(np.all(skel_points == spot, axis=1))[0][0]
            #dist1 = np.sum(skel_distances[0:k+1]*nm_to_pixels/1000)
            
            dist = skel_distances_um[k]
            
                
            spot_distances_first.append(dist)
        
          #then check if endpoint[0] was the correct starting endpoint. if not, update the attribute          
          skel_length = np.sum(skel_distances)
          if channel == channels[0]:
              spot_distances = findSoma(spot_distances_first, skel_length, nm_to_pixels)
              
              if not np.all(spot_distances == spot_distances_first):
                  dendrite.start_end = endpoints[1]
                  
              else:
                  dendrite.start_end = startEnd
                  
          else: 
              spot_distances = spot_distances_first
          
          
          dendrite.spot_distances[channel] = spot_distances
              
          #then do total distribution analysis to be saved in Excel sheet
          binned_distances, bins = np.histogram(spot_distances, bins=[0,25,50,75,100,125,150,175,200,225,250, inf])
          
          data.append([dendrite.num, channel] + list(binned_distances))
          
          dendrite.spot_distribution[channel] = list(binned_distances)
        else:
          print("\nSkeleton of " + dendrite.name() + " does not have 2 endpoints, so skipped it in distribution analysis")
  
  cols = ["Dendrite Num", "Channel", "0-25 (um)", "25-50 (um)", "50-75 (um)", "75-100 (um)", "100-125 (um)", "125-150 (um)","150-175 (um)","175-200 (um)","200-225 (um)", "225-250 (um)", ">= 250 (um)"]
  df = pd.DataFrame(data, columns=cols)
  
  return df


#Pick fake mRNA coordinate randomly from area of dendrite print
def simulateSpots(spot_distribution, print_distribution, nm_per_pix, return_pix=False):
  ''' Returns spots simulated in the entire dendrite, matching the distribution zones of mRNAs.
  If there is an error in spot detection, the dendrite will have mRNA detected in a zone where there are no print points
  In that case, PrintPointProblems.txt will show the spot distribution of the dendrite and the problematic zone, then select
  from the closest lesser zone with mRNA
  
  
  Parameters
  ----------      
  spot_distribution : list of num (1D)
        a list representing the distribution of mrna in the dendrite for a specific channel.
        e.g. [251,12,2,0,0] means from 0-25 um, there are 251 mrna, 12 mrna from 25-50 um, etc
        
  print_distribution : dictionary of form {'0-25':[list of coordinates], '25-50': list of coordinates,etc}
        a dictionary that has sorted all the spots in the print into the appropiate zones.
        Used to select simulated spots proportionally to the number of mrnas in a zone
        the x and y coordinates should represent pixels.
        
  nm_per_pix : num
        Conversion factor, nanometers per pixel
        
   

  Returns
  -------
  simNanomt : list of coordinates [x,y]
        a list of points, randomly chosen from the print, matching the distribution of the mrnas in each zone.
        the x and y coordinates represent nanometers.
        
  OR, if return_pix
  simPixels: a list of coordinates [x, y]
      a list of points, randomly chosen from the print, matching the distribution of the mrnas in each zone. 
      The x and y coordinates represent pixels.

  '''
  
  
  
  simPixels = [] #alist of [x,y] points that are randomly chosen within one zone only.
  
  #print("spot_distribution:", spot_distribution)
  
  
  for j in range(len(spot_distribution)):
      
      print_points_chosen = print_distribution[zones[j]]
      
      #check that it's not asking us to simulate mrnas where there are no print points
      ind = j
      while print_points_chosen == [] and spot_distribution[j] != 0:
          
          print("zone: ", zones[ind])
          print("no print points here")
          print("num mnra: ", spot_distribution[ind])
          print_points_chosen = print_distribution[zones[ind-1]]
          
          with open("PrintPointProblems.txt",'a') as f:
              f.write("\nSpot distribution: " + str(spot_distribution))
              f.write("\nZone: "+str(zones[ind]))
              f.write("\nnum mrna: "+str(spot_distribution[ind]))
              f.write("\nNo print points here\n")
          
          ind -= 1
          
          
          
      
      #edit: append method gives generator object which can not be indexed in the altered for loop
      simPixels_for_this_bin  = []
      
      for i in range(spot_distribution[j]):
          new_pix = list(random.choice(print_points_chosen))
          
          while new_pix in list(simPixels_for_this_bin):
              #print("choosing again...")
              new_pix = list(random.choice(print_points_chosen))
              
          simPixels_for_this_bin.append(new_pix)
      
      simPixels += simPixels_for_this_bin
  
  simNanomt = []
  for pix in simPixels:
      try:
          x_min = int(ceil(pix[0]*nm_per_pix))
          x_max = int(ceil(pix[0]*nm_per_pix + nm_per_pix))
          y_min = int(ceil(pix[1]*nm_per_pix))
          y_max = int(ceil(pix[1]*nm_per_pix + nm_per_pix))
          x_coord = random.choice(range(x_min, x_max))
          y_coord = random.choice(range(y_min, y_max))
          simNanomt.append((x_coord, y_coord))
      except:
          print("pix:", pix)
          
          global simp
          print(type(simPixels))
          
          simp = simPixels
          print(type(simp))
          print("index: ", simPixels.index(pix))
  
  if return_pix:
      return simNanomt, simPixels
  
  return simNanomt


#find the minimum distance of each spot of type 1 to any spot of type 2
def spots_to_spots(spots1, spots2, same=False):
  
  dists = []
  
  
  if len(spots1) != 0 and len(spots2) != 0:
      for spot1 in spots1:
          if same:
            if len(spots1) == 1: #added this if statement b/c taking np.min of empty list raises error
                return [-1] # self colocalization can not be calculated with a single spot
            
            #calculate the distance between this spot and every spot in spots2 (excluded itself), then choose the smallest one
            min_dist = np.min(np.array([dist(spot1, spot2) for spot2 in spots2 if np.all(spot1!=spot2)])) #need to find closest point other than itself
          
          else:  
            min_dist = np.min(np.array([dist(spot1, spot2) for spot2 in spots2]))
          
          
          dists.append(min_dist)
      
  if len(dists) == 0: #i.e. there may be some spots in in spots1 but there are zero spots in spots2 (colocalization can not be calculated)
      dists = [-1]*len(spots1)
  
  return dists




      
      
      

def dendriteCompartment(dendrites, nm_per_pixel):
  '''Sorts dendrite print points into zones. Initializes attribute dendrite.printDistribution
  Parameters
  ----------
  dendrites : list of Dendrite objects
  
  nm_to_pixels : num
      Conversion factor. Nanometers per pixel.
  Returns
  -------
  None.
  '''
  
  for dendrite in dendrites:
    for channel in dendrite.spot_coords.keys():
        
        
        spots = np.argwhere(dendrite.print)
        skeleton = dendrite.skeleton
        skel_points = np.argwhere(skeleton == 255)
        endpoints = findEndpoints(skeleton, skel_points) #find endpoints of skeleton
        
        if len(endpoints) == 2:
          skel_points, skel_distances, skel_distances_um = orderSkelPoints(skeleton, skel_points, dendrite.start_end)
          mapped_spots = [mapSig(spot, skel_points) for spot in spots] #map mrna coords to skeleton
          spot_distances = []

          #calculate distance of mapped mrna coord along skeleton
          for spot in mapped_spots:
            k = np.argwhere(np.all(skel_points == spot, axis=1))[0][0] 
            spot_distances.append(skel_distances_um[k])
          
          
          spot_distances_sorted = {} #{'0-25': [list of [x,y] points in print], '25-50': list of points, etc}
          zones = ['0-25', '25-50','50-75','75-100','100-125','125-150','150-175','175-200','200-225','225-250', '250-inf']
          for z in zones:
              spot_distances_sorted[z] = []
          
          
          for j, distance in enumerate(spot_distances):
              
              sorted_zone = sortDistance(distance)
              if sorted_zone == "neg" :
                  spot_distances_sorted['0-25'].append(spots[j])
                  print("Neg!")
              else:
                  spot_distances_sorted[sorted_zone].append(spots[j])
              
              
              
          dendrite.printDistribution = spot_distances_sorted
          
  return


def findDistAlongSkel(spots, dendrite, nm_per_pixel=107.5):
    '''
    Returns list of distances along skeleton for a list of spots (used in coloc analysis to find distances of simulated spots)

    Parameters
    ----------
    spots_nm : TYPE
        DESCRIPTION.
    dendrite : TYPE
        DESCRIPTION.
    nm_per_pixel : TYPE, optional
        DESCRIPTION. The default is 107.5.

    Returns
    -------
    spot_distances : TYPE
        DESCRIPTION.

    '''
    #spots = [(s[0]/nm_per_pixel, s[1]/nm_per_pixel) for s in spots_nm]
    skeleton = dendrite.skeleton
    skel_points = np.argwhere(skeleton == 255)
    
    skel_points, skel_distances, skel_distances_um = orderSkelPoints(skeleton, skel_points, dendrite.start_end)
    
    mapped_spots = [mapSig(spot, skel_points) for spot in spots] #map mrna coords to skeleton
    spot_distances = []

    #calculate distance of mapped mrna coord along skeleton
    for spot in mapped_spots:
      k = np.argwhere(np.all(skel_points == spot, axis=1))[0][0] 
      spot_distances.append(skel_distances_um[k])
      
    return spot_distances



def checkColocIssues(dendrite, mrna1, chan1, mrna2, chan2, mrna1_to_mrna2, mrna2_to_mrna1, check_self=False, mrna1_to_mrna1=None, mrna2_to_mrna2=None):
    ''' Will go through the spots of both channels. Consider following possibilities:
    1) There are some spots in channel1 but zero spots in channel2. Colocalization between these two can not be calculated. 
    Will replace mrna1_to_mrna2 with a list of the same length of mrna1, with every entry as -1. So, when np.histogram is run, it will not include this data and this dendrite will have 0s everywhere in the excel file. 
    Then will print a message explaining this. Same if mrna1 and mrna2 are reversed.
    2) There are some spots in channel1 and only one spot in channel2. Colocalization between them can be calculated and nothing is changed. Self colocalization between chan2 and itself can not be calculated.
    Will replace mrna2_to_mrna2 with a list of -1 entries so np.histogram does not count it. Print a message explaining this.
    3) There are some spots in both channel1 and channel2. Nothing is changed. 
    4) There are no spots in either channel 1 or channel2.. Raise error?
    5) There is only one spot in both channel 1 and channel 2...?
    Parameters
    ----------
    mrna1 : TYPE
        DESCRIPTION.
    chan1 : TYPE
        DESCRIPTION.
    mrna2 : TYPE
        DESCRIPTION.
    chan2 : TYPE
        DESCRIPTION.
    mrna1_to_mrna2 : TYPE
        DESCRIPTION.
    mrna2_to_mrna1 : TYPE
        DESCRIPTION.
    mrna1_to_mrna1 : TYPE
        DESCRIPTION.
    mrna2_to_mrna2 : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    len_mrna1, len_mrna2 = len(mrna1), len(mrna2)
    
    if len_mrna1 > 1 and len_mrna2 == 0:
        message = f"For image {dendrite.img_name} in dendrite number {dendrite.num}\n There are {len_mrna1} spots in {chan1} but zero spots in {chan2} \n Colocalization between {chan1} and {chan2} can not be calculated. This data will be missing in the excel sheet (zeros in every row)"
        with open("ColocIssues.txt", 'a') as f:
            f.write(message)
            f.write("\n \n")
        print(f"For image {dendrite.img_name} in dendrite number {dendrite.num}")
        print(f"There are {len_mrna1} spots in {chan1} but zero spots in {chan2}")
        print(f"Colocalization between {chan1} and {chan2} can not be calculated.")
        print("This data will be missing in the excel sheet (zeros in every row)")
        mrna1_to_mrna2 = [-1]*len_mrna1 #a list of -1 entries to hold the place of real data
    
    if len_mrna2 > 1 and len_mrna1 == 0:
        message = f"For image {dendrite.img_name} in dendrite number {dendrite.num}\nThere are {len_mrna2} spots in {chan2} but zero spots in {chan1}\nColocalization between {chan2} and {chan1} can not be calculated. This data will be missing in the excel sheet (zeros in every row)"
        with open("ColocIssues.txt", 'a') as f:
            f.write(message)
            f.write("\n \n")
        print(f"For image {dendrite.img_name} in dendrite number {dendrite.num}")
        print(f"There are {len_mrna2} spots in {chan2} but zero spots in {chan1}")
        print(f"Colocalization between {chan2} and {chan1} can not be calculated.")
        print("This data will be missing in the excel sheet (zeros in every row)")
        mrna2_to_mrna1 = [-1]*len_mrna2 #a list of -1 entries to hold the place of real data
    
    if not check_self:
        return  mrna1_to_mrna2, mrna2_to_mrna1
    
    if len_mrna1 > 1 and len_mrna2 == 1:
        message = f"For image {dendrite.img_name} in dendrite number {dendrite.num}\nThere are {len_mrna1} spots in {chan1} but only one spot in {chan2}\nColocalization between {chan1} and {chan2} can be calculated. Self colocalization of {chan2} can not be calculated and will be excluded from the excel sheet"
        with open("ColocIssues.txt", 'a') as f:
            f.write(message)
            f.write("\n \n")
        
        print(f"For image {dendrite.img_name} in dendrite number {dendrite.num}")
        print(f"There are {len_mrna1} spots in {chan1} but only one spot in {chan2}")
        print(f"Colocalization between {chan1} and {chan2} can be calculated.")
        print(f"Self colocalization of {chan2} can not be calculated and will be excluded from the excel sheet")
        mrna2_to_mrna2 = [-1] #a list of -1 entries to hold the place of real data
    
    if len_mrna2 > 1 and len_mrna1 == 1:
        message = f"For image {dendrite.img_name} in dendrite number {dendrite.num}\nThere are {len_mrna2} spots in {chan2} but only one spot in {chan1}\nColocalization between {chan2} and {chan1} can be calculated. Self colocalization of {chan1} can not be calculated and will be excluded from the excel sheet"
        with open("ColocIssues.txt", 'a') as f:
            f.write(message)
            f.write("\n \n")
            
        print(f"For image {dendrite.img_name} in dendrite number {dendrite.num}")
        print(f"There are {len_mrna2} spots in {chan2} but only one spot in {chan1}")
        print(f"Colocalization between {chan2} and {chan1} can be calculated.")
        print(f"Self colocalization of {chan1} can not be calculated and will be excluded from the excel sheet")
        mrna1_to_mrna1 = [-1] #a list of -1 entries to hold the place of real data
        
    return mrna1_to_mrna2, mrna2_to_mrna1, mrna1_to_mrna1, mrna2_to_mrna2

    
def colocAnalysis(dendrites, chan1, chan2, max_dist=250, incre=25, nanometers_per_pixel=107.5):
  '''
  Method: find colocalization of one whole dendrite alone. Simulate spots for the whole dendrite.
  Break up the counts into zones by using dendrite.spot_distances (indices correspond to dendrite.spot_coords[chan1])
  Do the simulation 1000 times. Once you save all the spots from all 100 sims, find distances along skeleton for each simulated spot.
  Use mrnaSim1_to_mrnaSim2 and the corresponding indices to mrnaSim1_dists and mrnaSim2_dists to break it up into zones
  Save as one dendrite, broken up into zones. 
  
    Parameters
    ----------
    dendrites : TYPE
        DESCRIPTION.
    chan1 : TYPE
        DESCRIPTION.
    chan2 : TYPE
        DESCRIPTION.
    max_dist : TYPE, optional
        DESCRIPTION. The default is 250.
    incre : TYPE, optional
        DESCRIPTION. The default is 25.
    nanometers_per_pixel : TYPE, optional
        DESCRIPTION. The default is 107.5.

    Returns
    -------
    data_df : TYPE
        DESCRIPTION.

  '''
  total_zone_list = []
  total_dend_num_list = []
  total_bin_labels = []
  
  total_chan1_to_chan2 = []
  total_chan2_to_chan1 = []
  total_chan1_to_chan1 = []
  total_chan2_to_chan2 = []
  total_Simchan1_to_Simchan2 = []
  total_Simchan2_to_Simchan1 = []
  total_Simchan1_to_Simchan1 = []
  total_Simchan2_to_Simchan2 = []
  
  for dendrite in dendrites:
    if chan1 in dendrite.spot_coords.keys() and chan2 in dendrite.spot_coords.keys():
      mrnas1 = dendrite.spot_coords[chan1] #len equal to mrnas1
      mrnas2 = dendrite.spot_coords[chan2]
      mrnas1_dists = dendrite.spot_distances[chan1] #len equal to mrnas1
      mrnas2_dists = dendrite.spot_distances[chan2]
      
      
      
      #get unsimulated coloc data
      mrna1_to_mrna2 = spots_to_spots(mrnas1, mrnas2) #len equal to mrnas1
      mrna2_to_mrna1 = spots_to_spots(mrnas2, mrnas1)
      mrna1_to_mrna1 = spots_to_spots(mrnas1, mrnas1) #len equal to mrnas1
      mrna2_to_mrna2 = spots_to_spots(mrnas2, mrnas2)
      
      #print(f"mRNAs1: {len(mrnas1)}, mRNAs1to2: {len(mrna1_to_mrna2)}, mRNAS1to1: {len(mrna1_to_mrna1)}")
      #print(f"mRNAs2: {len(mrnas2)}, mRNAs2to1: {len(mrna2_to_mrna1)}, mRNAS1to1: {len(mrna2_to_mrna2)}")
      
      #fill in any empty lists if the dend happens to have zero spots for some channel (1to2 or 2to1) or just one (self coloc)
      #mrna1_to_mrna2, mrna2_to_mrna1, mrna1_to_mrna1, mrna2_to_mrna2 = checkColocIssues(dendrite, mrnas1, chan1, mrnas2, chan2, mrna1_to_mrna2, mrna2_to_mrna1, check_self=True, mrna1_to_mrna1=mrna1_to_mrna1, mrna2_to_mrna2=mrna2_to_mrna2)
      
      # sort above data into zone from dendrite (25 um each)
      mrnas1_sorted, mrnas1_dists_sorted, mrnas2_sorted, mrnas2_dists_sorted = {}, {}, {}, {}
      mrna1_to_mrna2_sorted, mrna1_to_mrna1_sorted = {}, {}
      mrna2_to_mrna1_sorted, mrna2_to_mrna2_sorted = {}, {}
      
      for sorted_zone in zones:
          mrnas1_sorted[sorted_zone] = []
          mrnas1_dists_sorted[sorted_zone] = []
          mrna1_to_mrna2_sorted[sorted_zone] = []
          mrna1_to_mrna1_sorted[sorted_zone] = []
          mrnas2_sorted[sorted_zone] = []
          mrnas2_dists_sorted[sorted_zone] = []
          mrna2_to_mrna2_sorted[sorted_zone] = []
          mrna2_to_mrna1_sorted[sorted_zone] = []
      
      # sort colocalization data for chan 1 (to chan2 and to self) into zones
      for i, mrna1 in enumerate(mrnas1):
          mrna1_dist = mrnas1_dists[i]
          closest_mrna2 = mrna1_to_mrna2[i]
          closest_mrna1 = mrna1_to_mrna1[i]
          
          #find distance of this spot from soma
          sorted_zone = sortDistance(mrna1_dist)
          
          if sorted_zone == "neg":
              print("negative value: ", mrna1_dist)
              sorted_zone = '0-25'
          
          mrnas1_sorted[sorted_zone].append(mrna1)
          mrnas1_dists_sorted[sorted_zone].append(mrna1_dist)
          mrna1_to_mrna2_sorted[sorted_zone].append(closest_mrna2)
          mrna1_to_mrna1_sorted[sorted_zone].append(closest_mrna1)
          
             
    
      
      # sort colocalization data for chan 2 (to chan1 and to self) into zones
      for i, mrna2 in enumerate(mrnas2):
          mrna2_dist = mrnas2_dists[i]
          closest_mrna1 = mrna2_to_mrna1[i]
          closest_mrna2 = mrna2_to_mrna2[i]
          
          #find distance of this spot from soma
          sorted_zone = sortDistance(mrna2_dist)
          
          if sorted_zone == "neg":
              print("negative value: ", mrna2_dist)
              sorted_zone = '0-25'
          
          mrnas2_sorted[sorted_zone].append(mrna2)
          mrnas2_dists_sorted[sorted_zone].append(mrna2_dist)
          mrna2_to_mrna2_sorted[sorted_zone].append(closest_mrna2)
          mrna2_to_mrna1_sorted[sorted_zone].append(closest_mrna1)
          
     
      
      #find binned distances
      bins = [i for i in range(0, max_dist+1, incre)] + [inf]
      bin_labels = [str(bins[i]) + "-" + str(bins[i+1]) for i in range(len(bins)-1)]
      
      # save real data into lists (columns for data frame)
      chan1_to_chan2 = {}
      chan2_to_chan1 = {}
      chan1_to_chan1 = {}
      chan2_to_chan2 = {}
      
      for z in zones:
          chan1_to_chan2[z] = np.histogram(mrna1_to_mrna2_sorted[z], bins=bins)[0]
          chan2_to_chan1[z] = np.histogram(mrna2_to_mrna1_sorted[z], bins=bins)[0]
          chan1_to_chan1[z] = np.histogram(mrna1_to_mrna1_sorted[z], bins=bins)[0]
          chan2_to_chan2[z] = np.histogram(mrna2_to_mrna2_sorted[z], bins=bins)[0]
          
          # first three columns 
          total_zone_list += [z]*len(bin_labels)
          total_dend_num_list +=  [dendrite.num]*len(bin_labels)
          total_bin_labels += bin_labels
          
          # real data columns
          total_chan1_to_chan2 += list(chan1_to_chan2[z])
          total_chan2_to_chan1 += list(chan2_to_chan1[z])
          total_chan1_to_chan1 += list(chan1_to_chan1[z])
          total_chan2_to_chan2 += list(chan2_to_chan2[z])
          
      
      ### get simulated data 
      
      # check if simulation was previously done (e.g. in synapseAnalysis)
      #previously_done_sim = chan1 in dendrite.simulated_spots.keys() and chan2 in dendrite.simulated_spots.keys()
      previously_done_sim = False
      if previously_done_sim:
          # make sure the correct number of simulations exist
          if  len(dendrite.simulated_spots[chan1]) != 100 or  len(dendrite.simulated_spots[chan1]) != 100:
              print("We have ", len(dendrite.simulated_spots[chan1]), "simulations for", chan1)
              print("We have ", len(dendrite.simulated_spots[chan2]), "simulations for", chan2)
      else:
          dendrite.simulated_spots[chan1] = []
          dendrite.simulated_spots[chan2] = []
        
        
      Simchan1_to_Simchan2 = {}
      Simchan2_to_Simchan1 = {}
      Simchan1_to_Simchan1 = {}
      Simchan2_to_Simchan2 = {}
      for z in zones:
          Simchan1_to_Simchan2[z] = []
          Simchan2_to_Simchan1[z] = []
          Simchan1_to_Simchan1[z] = []
          Simchan2_to_Simchan2[z] = []
          
      for j in range(100):
          
          if previously_done_sim:
              sim1, sim1_pix  = dendrite.simulated_spots[chan1][j][0]
              sim2, sim2_pix  = dendrite.simulated_spots[chan2][j][0]
          else:
              sim1, sim1_pix = simulateSpots(dendrite.spot_distribution[chan1], dendrite.printDistribution, nanometers_per_pixel, return_pix=True)
              sim2, sim2_pix  = simulateSpots(dendrite.spot_distribution[chan2], dendrite.printDistribution, nanometers_per_pixel, return_pix=True)
              dendrite.simulated_spots[chan1].append((sim1, sim1_pix))
              dendrite.simulated_spots[chan2].append((sim2, sim2_pix))
              
          mrnaSim1_to_mrnaSim2 = spots_to_spots(sim1, sim2)
          mrnaSim2_to_mrnaSim1 = spots_to_spots(sim2, sim1)
          mrnaSim1_to_mrnaSim1 = spots_to_spots(sim1, sim1)
          mrnaSim2_to_mrnaSim2 = spots_to_spots(sim2, sim2)
          
          mrnasSim1_dists = findDistAlongSkel(sim1_pix, dendrite, nm_per_pixel=nanometers_per_pixel)
          mrnasSim2_dists = findDistAlongSkel(sim2_pix, dendrite, nm_per_pixel=nanometers_per_pixel)
          
          # sort above data into zone from dendrite (25 um each)
          mrnasSim1_sorted, mrnasSim1_dists_sorted, mrnasSim2_sorted, mrnasSim2_dists_sorted = {}, {}, {}, {}
          mrnaSim1_to_mrnaSim2_sorted, mrnaSim1_to_mrnaSim1_sorted = {}, {}
          mrnaSim2_to_mrnaSim1_sorted, mrnaSim2_to_mrnaSim2_sorted = {}, {}
          
          for sorted_zone in zones:
              mrnasSim1_sorted[sorted_zone] = []
              mrnasSim1_dists_sorted[sorted_zone] = []
              mrnaSim1_to_mrnaSim2_sorted[sorted_zone] = []
              mrnaSim1_to_mrnaSim1_sorted[sorted_zone] = []
              mrnasSim2_sorted[sorted_zone] = []
              mrnasSim2_dists_sorted[sorted_zone] = []
              mrnaSim2_to_mrnaSim2_sorted[sorted_zone] = []
              mrnaSim2_to_mrnaSim1_sorted[sorted_zone] = []
          
          # sort colocalization data for chan 1 (to chan2 and to self) into zones
          for i, mrnaSim1 in enumerate(sim1):
              mrnaSim1_dist = mrnasSim1_dists[i]
              closest_mrnaSim2 = mrnaSim1_to_mrnaSim2[i]
              closest_mrnaSim1 = mrnaSim1_to_mrnaSim1[i]
              
              #find distance of this spot from soma
              sorted_zone = sortDistance(mrnaSim1_dist)
              
              if sorted_zone == "neg":
                  print("negative value: ", mrnaSim1_dist)
                  sorted_zone = '0-25'
              
              mrnasSim1_sorted[sorted_zone].append(mrnaSim1)
              mrnasSim1_dists_sorted[sorted_zone].append(mrnaSim1_dist)
              mrnaSim1_to_mrnaSim2_sorted[sorted_zone].append(closest_mrnaSim2)
              mrnaSim1_to_mrnaSim1_sorted[sorted_zone].append(closest_mrnaSim1)
              
        
          
          # sort colocalization data for chan 2 (to chan1 and to self) into zones
          for i, mrnaSim2 in enumerate(sim2):
              mrnaSim2_dist = mrnasSim2_dists[i]
              closest_mrnaSim1 = mrnaSim2_to_mrnaSim1[i]
              closest_mrnaSim2 = mrnaSim2_to_mrnaSim2[i]
              
              #find distance of this spot from soma
              sorted_zone = sortDistance(mrnaSim2_dist)
              
              if sorted_zone == "neg":
                  print("negative value: ", mrnaSim2_dist)
                  sorted_zone = '0-25'
              
              mrnasSim2_sorted[sorted_zone].append(mrnaSim2)
              mrnasSim2_dists_sorted[sorted_zone].append(mrnaSim2_dist)
              mrnaSim2_to_mrnaSim2_sorted[sorted_zone].append(closest_mrnaSim2)
              mrnaSim2_to_mrnaSim1_sorted[sorted_zone].append(closest_mrnaSim1)
        
          for z in zones:
              if  len(mrnas1_sorted[z]) != len(mrnasSim1_sorted[z]) or  len(mrnas2_sorted[z]) != len(mrnasSim2_sorted[z]):
                  print(z, "mRNA1: ", len(mrnas1_sorted[z]), "Sim mRNA1:", len(mrnasSim1_sorted[z]))
                  print("mRNA2: ", len(mrnas2_sorted[z]), "Sim mRNA2:", len(mrnasSim2_sorted[z]))
                  print("mRNA1to2: ", len(mrna1_to_mrna2_sorted[z]), "Sim mRNA1to2:", len(mrnaSim1_to_mrnaSim2_sorted[z]))
                  print("mRNA2to1: ", len(mrna2_to_mrna1_sorted[z]), "Sim mRNA1to2:", len(mrnaSim2_to_mrnaSim1_sorted[z]))
                  
              
              
          # save binned ddata for this simulation
          for z in zones:
              Simchan1_to_Simchan2[z].append(np.histogram(mrnaSim1_to_mrnaSim2_sorted[z], bins=bins)[0])
              Simchan2_to_Simchan1[z].append(np.histogram(mrnaSim2_to_mrnaSim1_sorted[z], bins=bins)[0])
              Simchan1_to_Simchan1[z].append(np.histogram(mrnaSim1_to_mrnaSim1_sorted[z], bins=bins)[0])
              Simchan2_to_Simchan2[z].append(np.histogram(mrnaSim2_to_mrnaSim2_sorted[z], bins=bins)[0])
        
              
      # save average of simulated data into lists (columns for data frame)
      for z in zones:
          sim1to2 =  np.array(Simchan1_to_Simchan2[z])
          total_Simchan1_to_Simchan2 += list(sim1to2.mean(axis=0))
          
          sim2to1 =  np.array(Simchan2_to_Simchan1[z])
          total_Simchan2_to_Simchan1 += list(sim2to1.mean(axis=0))
          
          sim1to1 =  np.array(Simchan1_to_Simchan1[z])
          total_Simchan1_to_Simchan1 += list(sim1to1.mean(axis=0))
          
          sim2to2 =  np.array(Simchan2_to_Simchan2[z])
          total_Simchan2_to_Simchan2 += list(sim2to2.mean(axis=0))
          
  
  print("ch1 to ch2:", len(total_chan1_to_chan2))
  print("zone list", len(total_zone_list))
  print("dend num list: ", len(total_dend_num_list))
  print("bin list: ", len(total_bin_labels))
  data = np.concatenate(([total_zone_list],
                         [total_dend_num_list],
                         [total_bin_labels],
                         [total_chan1_to_chan2],
                         [total_chan2_to_chan1],
                         [total_chan1_to_chan1],
                         [total_chan2_to_chan2],
                         [total_Simchan1_to_Simchan2],
                         [total_Simchan2_to_Simchan1],
                         [total_Simchan1_to_Simchan1],
                         [total_Simchan2_to_Simchan2]), axis=0)

  cols = ["Zone (um from Soma)",
          "Dendrite number",
          "Distance (nm)",
          chan1 + " to closest " + chan2,
          chan2 + " to closest " + chan1, 
          chan1 + " to closest " + chan1,
          chan2 + " to closest " + chan2,
          "Sim-"+chan1 + " to closest " + "Sim-"+chan2,
          "Sim-"+chan2 + " to closest " + "Sim-"+chan1,
          "Sim-"+chan1 + " to closest " + "Sim-"+chan1,
          "Sim-"+chan2 + " to closest " + "Sim-"+chan2]
  
  data_total = {cols[0]: data[0],
                cols[1]: data[1],
                cols[2]: data[2],
                cols[3]: data[3],
                cols[4]: data[4],
                cols[5]: data[5],
                cols[6]: data[6],
                cols[7]: data[7],
                cols[8]: data[8],
                cols[9]: data[9],
                cols[10]: data[10]}
  
  #data_df = pd.DataFrame(data, columns=cols) # add dataframe 
  data_df = pd.DataFrame(data_total)
  #print(data_df)
  return data_df
    
    
#numer of mRNA within threshold distance of each synapse
def synapMRNACount(synap_pts, mrna_pts, dist_thresh):
    counts = []
    for synap in synap_pts:
        count = np.count_nonzero([dist(synap, mrna) < dist_thresh for mrna in mrna_pts])
        counts.append(count)
    return counts


def mRNA_to_synAnalysis(dendrites, syn_chan, mRNA_chans, thresh, nanometers_per_pixel=107.5):
    master_dict = {} #{Dendrite1: {chan1: {'0-25': Counter dict, '25-50': counter dict}, chan2:{'0-25': Counter dict, '25-50': counter dict}}, Dendrite2:etc }
    #{Dendrite: {channel:{zone: {count:num of synapses}}}}
    #the zones in sub dictionary represents distances of the SYNAPSE from the soma
    #Counter dict is in form {num of synapses: number of mrna with this number of synapses} e.g. {1:3,2:1,3:0} 3 synapses have 1 mrna, 1 synapse has 2 mrna, and no synapse has 3 mrnas
    
    for dendrite in dendrites:
        synapses = dendrite.spot_coords[syn_chan]
        
        
        master_dict[dendrite] = {}
        
        for chan in mRNA_chans:
            spots = dendrite.spot_coords[chan]
            spot_dists = dendrite.spot_distances[chan]
            sorted_mrnas = {} #{zone: [list of mRNA coordinates in that zone]}
            sorted_dists = {}
            
            for z in zones:
                sorted_mrnas[z] = []
                sorted_dists[z] = []
                
                
            #first sort each mrna into its zone
            for j in range(len(spot_dists)):
                sorted_zone = sortDistance(spot_dists[j])
                if sorted_zone == "neg":
                    print("negative mRNA: ", spots[j])
                    print("distance: ", spot_dists[j])
                    sorted_zone = '0-25'
                
                sorted_dists[sorted_zone].append(spot_dists[j])
                
                sorted_mrnas[sorted_zone].append(spots[j])
            
            sorted_counters = {} #{zone: Counter object}
            #now iterate through each zone, and quantify # of synapses near each mrna
            for z in zones:
                sorted_counters[z] = Counter() #a counter for each zone for real mRNA
                
            
            for z in zones:
                spots_in_zone = sorted_mrnas[z]
                real_counts = synapMRNACount(spots_in_zone, synapses, thresh)
                sorted_counters[z].update(real_counts)
                if len(spots_in_zone) == 0:
                    sorted_counters[z] = Counter({0:0})
            
            master_dict[dendrite][chan] = sorted_counters
            
            #now simulate mRNAs and repeat the same thing. Assume simulation has already happened
            calculateSimSpotsDists(dendrite, nanometers_per_pixel)
            
            
            sorted_sim_counters = {} #{zone: Counter object}
            for z in zones:
                sorted_sim_counters[z] = Counter() #a counter for each zone for real mRNA
            
            
            for i in range(100):
                sim_spots = dendrite.simulated_spots[chan][i][0]
                sim_dists = dendrite.simulated_spot_distances[chan][i]
                
                sorted_sim_mrnas = {} #{zone: [list of mRNA coordinates in that zone]}
                sorted_sim_dists = {}
                
                for z in zones:
                    sorted_sim_mrnas[z] = []
                    sorted_sim_dists[z] = []
                
                #first sort each mrna into its zone
                for j in range(len(sim_dists)):
                    sorted_zone = sortDistance(sim_dists[j])
                    if sorted_zone == "neg":
                        print("negative mRNA: ", sim_spots[j])
                        print("distance: ", sim_dists[j])
                        sorted_zone = '0-25'
                    
                    sorted_sim_dists[sorted_zone].append(sim_dists[j])
                    
                    sorted_sim_mrnas[sorted_zone].append(sim_spots[j])
                    
                for z in zones:
                    spots_in_zone = sorted_sim_mrnas[z]
                    real_counts = synapMRNACount(spots_in_zone, synapses, thresh)
                    sorted_sim_counters[z].update(real_counts)
                
            #after all the simulations, take average
            for z in zones:
                if len(sorted_sim_counters[z]) == 0:
                    sorted_sim_counters[z] = Counter({0:0})
                    
                for k in sorted_sim_counters[z]: 
                    
                    sorted_sim_counters[z][k] /= 100
            
            master_dict[dendrite]["Sim-"+chan] = sorted_sim_counters
    
    #after iterating through every dendrite, every channel, and quantifying synapses near each mrna
    dfs = []
    total_df = pd.DataFrame()
    #create a single df for each dendrite with dend num, zone, num of synapses, and num of mrna with that synapse for each chan (real, sim)
    #then append all these dfs together
    
    for dendrite in dendrites:
        data = master_dict[dendrite]
        max_nums = {}
        for z in zones:
            max_nums[z] = 0
            for chan, sorted_counters in data.items():
                counter = sorted_counters[z]
                
                #check what the max # of synapses found around 1 mrna is; this determines how many rows we need for each column
                if max(list(counter.keys())) > max_nums[z]:
                       max_nums[z] = max(list(counter.keys()))
        
        #now we know how many rows each zone in this dend needs
        
        #create a df for each zone, then append them to get the full df for this dendrite.
        df_dend = pd.DataFrame()   
        
        for z in zones:
            max_num = max_nums[z]
            
            
            df = pd.DataFrame({"Dendrite number": [dendrite.num]*(max_num+1),
                              "Zone (um) - distance of mRNA from zoma": [z]*(max_num+1),
                              "Number of synapses": list(range(max_num+1))})
                        
            #df will contain all the info for a specific zone in a specific dend
            for chan, sorted_counters in data.items():
                counter = sorted_counters[z]
                chan_counts = [counter[i] if i in counter.keys() else 0 for i in range(max_num+1)]
                
                df.insert(len(df.columns), chan, chan_counts)
                
                
            df_dend = df_dend.append(df)
        
        #df_dend will contain info for a specific dend, with every zone
        total_df = total_df.append(df_dend)
        
    #total_df contains all the info for all dends, and all zones
    
    return total_df
        
        
                
                
                
            
    
def synAnalysis(dendrites, syn_chan, thresh, nanometers_per_pixel=107.5):
  master_dict = {} #{Dendrite1: {chan1: {'0-25': Counter dict, '25-50': counter dict}, chan2:{'0-25': Counter dict, '25-50': counter dict}}, Dendrite2:etc }
  #{Dendrite: {channel:{zone: {count:num of synapses}}}}
  #the zones in sub dictionary represents distances of the SYNAPSE from the soma
  #Counter dict is in form {num of mrnas: number of synapses with this number of mrnas} e.g. {1:3,2:1,3:0} 3 synapses have 1 mrna, 1 synapse has 2 mrna, and no synapse has 3 mrnas
  
  
  #an added layer of data to see if different mrna serve the same synapses
  serving_mrna_df = {} #{Dendrite: df} df with columns zone, synapse, mrna counts chan 1, mrna counts chan 2...
  
  for dendrite in dendrites:
    master_dict[dendrite] = {}
    
    
    
    channels = []
    if syn_chan in dendrite.spot_coords.keys():
      synapses = dendrite.spot_coords[syn_chan] #[list of synapse spots x,y]
      #we will sort synapses into their dendrite zones. each will get a dedicated Counter
      #then we will run synapseMRNACount for each zone-set of synapses individually and update their counters 
      #then we will run 100 sims. after each individual sim, we will run synapseMRNACount for each zone-set of synapses
      #and update their counters each time
      sorted_synapses = {}
      sorted_dists = {}
      
      for z in zones:
          sorted_synapses[z] = []
          sorted_dists[z] = []
                    
      synapses_dists = dendrite.spot_distances[syn_chan]
      
      for j in range(len(synapses_dists)):
          sorted_zone = sortDistance(synapses_dists[j])
          if sorted_zone == "neg":
              print("negative synapse: ", synapses[j])
              print("distance: ", synapses_dists[j])
              sorted_zone = '0-25'
          
          sorted_dists[sorted_zone].append(synapses_dists[j])
          
          sorted_synapses[sorted_zone].append(synapses[j])
      #now synapses_sorted gives us the synapses sorted into each zone. 
      
      #initialize dataframe for this dendrite
      DN_list, dist_list, zone_list, syn_list = [],[],[],[] #DN stands for dendrite Number
      for z in zones:
          
          DN_list += [dendrite.num]*len(sorted_synapses[z])
          dist_list += sorted_dists[z]
          zone_list += [z]*len(sorted_synapses[z])
          syn_list += sorted_synapses[z]
          
      df_serving_mrna = pd.DataFrame({"Dendrite number": DN_list,
                                      "Dist from soma (um)": dist_list, 
                                      "Zone from soma (um)": zone_list,
                                      "Synapse coordinate (yx)": syn_list})
      
      #iterate through each channel
      for chan, mrnas in dendrite.spot_coords.items():
          
          real_mrna_counters = {} #a counter for each zone for real mRNA
          ctrl_mrna_counters = {} #a counter for each zone for simulated mRNA
          for z in zones:
              real_mrna_counters[z] = Counter() #a counter for each zone for real mRNA
              ctrl_mrna_counters[z] = Counter() #{'0-25': Counter dict, '25-50': counter dict... }
              
          
          if chan!= syn_chan:
              channels.append(chan)
              #counts is the list of how many mrna are in each synapse. 
              #counter is the dictionary that counts the number of synapses with each number of mrna
              chan_served_mrna_counts = []
              
              for z in zones: #do the real mRNA counts for just the synapses in each zone
                  synapses_in_zone = sorted_synapses[z]
                  
                  if len(synapses_in_zone) == 0:
                      real_mrna_counters[z] = Counter({0:0})
                      
                  else:
                      mrna_counts = synapMRNACount(synapses_in_zone, mrnas, thresh)
                      chan_served_mrna_counts += mrna_counts
                      
                      real_mrna_counters[z].update(mrna_counts)
              shape = df_serving_mrna.shape
              df_serving_mrna.insert(shape[1],"Num of "+chan,chan_served_mrna_counts)   
              
              #run 100 simulations now and update the ctrl_mrna_counters for each zone every time. 
              if chan not in dendrite.simulated_spots.keys():
                  dendrite.simulated_spots[chan] = [] #save the simulated spots in this attribute to avoid redundancies when running colocAnalysis
                  
              
              for i in range(100):
                  
                  if chan in dendrite.simulated_spots.keys() and len(dendrite.simulated_spots[chan]) == 100:
                      print("Syn no sim needed, i=",i)
                      ctrl_mrnas, ctrl_mrnas_pix = dendrite.simulated_spots[chan][i]
                  else:
                      
                      ctrl_mrnas, ctrl_mrnas_pix = simulateSpots(dendrite.spot_distribution[chan], dendrite.printDistribution, nanometers_per_pixel, return_pix=True)
                      dendrite.simulated_spots[chan].append((ctrl_mrnas, ctrl_mrnas_pix))
                      
                  #update counter using these simulated mrnas for synapses in each zone of dendrites
                  for z in zones: #do the real mRNA counts for just the synapses in each zone
                      synapses_in_zone = sorted_synapses[z]
                      if len(synapses_in_zone) == 0:
                          ctrl_mrna_counters[z] = Counter({0:0})
                      else:
                          ctrl_counts = synapMRNACount(synapses_in_zone, ctrl_mrnas, thresh)
                          ctrl_mrna_counters[z].update(ctrl_counts)
              
              #now divide the control amounts by 100 to average over 100 simulations
              for z in zones:
                  for k in ctrl_mrna_counters[z]: 
                      ctrl_mrna_counters[z][k] /= 100
              
              master_dict[dendrite][chan] = [real_mrna_counters, ctrl_mrna_counters] #counter dicts for mrna and ctrl, sorted in dict by zone
  
    serving_mrna_df[dendrite] = df_serving_mrna
    
  mrnas_serving_df = pd.DataFrame()
  
  for dendrite, df in serving_mrna_df.items():
      
      mrnas_serving_df = mrnas_serving_df.append(df)
  
  
  d_nums = {}
  for dendrite in master_dict.keys():
      nums = {} #{zone: all possible mrna counts from all channels}, used to determine list length
      for z in zones:
          nums[z] = []
          
      for chan, counters in master_dict[dendrite].items():
          counter1, counter2 = counters #[sorted_mrna_counters, sorted_ctrl_counters]
          for z in zones:
              c1, c2 = counter1[z], counter2[z] #mrna_counters, ctrl_counters for a particular zone
              nums[z] += list(c1.keys()) + list(c2.keys())
              
      d_nums[dendrite] = nums
        
     
        
  dfs = {} #{chan1: list of dfs, chan2: list of dfs}. after appending all dfs we will concatenate 
  for c in channels:
      dfs[c] = pd.DataFrame()
       
  for dendrite in master_dict.keys():
      nums = d_nums[dendrite]
      max_nums = {}
      for z in zones:
          if len(nums[z]) != 0:
              max_nums[z] = max(nums[z])
          else: 
              max_nums[z] = 0
          
      for chan in channels:
          sorted_mrna_counters, sorted_ctrl_counters = master_dict[dendrite][chan]
          for zone, mrna_counter in sorted_mrna_counters.items():
              ctrl_counter = sorted_ctrl_counters[zone]
                               
              dendrite_nums = [dendrite.num]*(max_nums[zone]+1)
              zone_list = [zone]*(max_nums[zone]+1)
              counts_list = [i for i in range(max_nums[zone]+1)]
              
              
              real_data = [0 if i not in mrna_counter.keys() else mrna_counter[i] for i in range(max_nums[zone] + 1)]
              sim_data = [0 if i not in ctrl_counter.keys() else ctrl_counter[i] for i in range(max_nums[zone] + 1)]
              if sum(real_data) != 0:
                  real_freq = [i/sum(real_data) for i in real_data]
              else: 
                  real_freq = [0]*len(real_data)
              if sum(sim_data) != 0:
                  sim_freq = [i/sum(sim_data) for i in sim_data]
              else:
                  sim_freq = [0]*len(sim_data)
               
              #data = np.array([dendrite_nums, zone_list, counts_list, real_data, real_freq, sim_data, sim_freq]).T
              data = [dendrite_nums, zone_list, counts_list, real_data, real_freq, sim_data, sim_freq]
              columns=["Dendrite number", "Zone from soma (um)", "Num of mRNA", "Real " + chan,"Real "+chan+" Frequency", "Sim "+ chan, "Sim "+chan+" Frequency"]
              
              
              t_data = {columns[0]:data[0],
                        columns[1]:data[1],
                        columns[2]:data[2],
                        columns[3]:data[3],
                        columns[4]:data[4],
                        columns[5]:data[5],
                        columns[6]:data[6]}
              
              df = pd.DataFrame(t_data)
              
              dfs[chan] = dfs[chan].append(df, ignore_index=True)
            
  dfs = list(dfs.values())    
  
  #this is taken from zoes code but i dont realy get it. if !=0 but elif ==1 ? the second block would never run...
  
  if len(dfs) != 0:
      data = pd.concat(dfs, axis = 1)
  elif len(dfs) == 1:
      data = pd.concat([dfs[0], pd.DataFrame([], columns=[3,4])], axis = 1)
  else:
      data = pd.concat([pd.DataFrame([], columns = [1,2,3]), pd.DataFrame([], columns=[4,5,6])], axis = 1)
  return data, mrnas_serving_df

    

