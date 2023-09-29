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

#Find pixels in skel with only 1 neighboor
def findEndpoints(skel, points):
  endpoints = []
  for point in points:
    if np.sum((skel[(point[0] - 1):(point[0] + 2), (point[1] - 1):(point[1] + 2)])) == 510: #255*2 (pixel and its 1 neighboor)
      endpoints.append(point)
  return endpoints

#Order coordinates of skeleton in order from one end to the other
def orderSkelPoints(skel, skel_points, endpoint):
  dist_array = np.zeros(len(skel_points))
  skel_array = np.zeros(skel_points.shape)
  prev_coord = endpoint
  cur_coord = endpoint

  for i in range(len(skel_points)):
    skel_array[i] = cur_coord
    dist_array[i] = dist(prev_coord, cur_coord)
    #look for white pixels adjacent/diagonal to current pixel
    pixels = np.argwhere(skel[cur_coord[0] - 1:cur_coord[0] + 2, cur_coord[1] - 1:cur_coord[1] + 2] == 255)
    pixels = [(cur_coord - (1,1) + pix) for pix in pixels]

    cur_coord_copy = cur_coord
    prev_coord_copy = prev_coord
    
    for pix in pixels:
      if np.any(pix != cur_coord_copy) and np.any(pix != prev_coord_copy):
        prev_coord = cur_coord
        cur_coord = pix
        
  
  return skel_array, dist_array

#map a mRNA coordinate to the closest point in a list of pixel coordinates
def mapSig(sig, points):
  k = np.argmin(np.array([dist(pt, sig) for pt in points]))
  return points[k]



#skel_array is list of pixels in skeleton in order
#dist_array is list of distances from each pixel in skel_array to the next
def orderSkelPoints(skel, skel_points, endpoint):
    skel_len = len(skel_points)
    dist_array = np.zeros(skel_len)
    skel_array = np.zeros(skel_points.shape)
    prev_coord = endpoint
    cur_coord = endpoint

    for i in range(skel_len):
      skel_array[i] = cur_coord
      dist_array[i] = dist(prev_coord, cur_coord) #sqrt(2) if diagonal step, 1 o/w

      pixels = np.argwhere(skel[cur_coord[0] - 1:cur_coord[0] + 2, cur_coord[1] - 1:cur_coord[1] + 2] == 255)
      pixels = [(cur_coord - (1,1) + pix) for pix in pixels]

      cur_coord_copy = cur_coord
      prev_coord_copy = prev_coord

      for pix in pixels:
        if np.any(pix != cur_coord_copy) and np.any(pix != prev_coord_copy):
            prev_coord = cur_coord
            cur_coord = pix
        
    return skel_array, dist_array #returns skel_points in order and an array of distance between adjacent points in skeleton


def dist(pt, sig):
  return sqrt((pt[0] - sig[0])**2 + (pt[1] - sig[1])**2)

#Assume soma is the end of skeletons with more mRNA coordinates
def findSoma(distances, skel_len, nm_to_pixels):
    skel_len_um = skel_len * nm_to_pixels / 1000
    sum_mrna_dists = np.sum(distances)
    if 2*sum_mrna_dists > skel_len_um * len(distances): #if initial guess was not soma end
        distances = [skel_len_um-dist for dist in distances] #take distance starting from oposite end of skel
    return distances


#Pick fake mRNA coordinate uniformly randomly from area of dendrite print
def simulateSpots(print_points, num_spots, nm_per_pix):
  simPixels  = [random.choice(print_points) for i in range(num_spots)]
  simNanomt = []
  for pix in simPixels:
      x_min = int(ceil(pix[0]*nm_per_pix))
      x_max = int(ceil(pix[0]*nm_per_pix + nm_per_pix))
      y_min = int(ceil(pix[1]*nm_per_pix))
      y_max = int(ceil(pix[1]*nm_per_pix + nm_per_pix))
      x_coord = random.choice(range(x_min, x_max))
      y_coord = random.choice(range(y_min, y_max))
      simNanomt.append((x_coord, y_coord))
  return simNanomt


#find the minimum distance of each spot of type 1 to any spot of type 2
def spots_to_spots(spots1, spots2, same=False):
  dists = []
  if len(spots1) != 0 and len(spots2) != 0:
    for spot1 in spots1:
      if same:
        #if there is only one spot in the list, then self colocalization can not be found
        if len(spots1) == 1:
            return dists
        min_dist = np.min(np.array([dist(spot1, spot2) for spot2 in spots2 if np.all(spot1!=spot2)])) #need to find closest point other than itself
      else:  
        min_dist = np.min(np.array([dist(spot1, spot2) for spot2 in spots2]))
      dists.append(min_dist)
  return dists

#numer of mRNA within threshold distance of each synapse
def synapMRNACount(synap_pts, mrna_pts, dist_thresh):
    counts = []
    for synap in synap_pts:
        count = np.count_nonzero([dist(synap, mrna) < dist_thresh for mrna in mrna_pts])
        counts.append(count)
    return counts


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


def distrAnalysis(dendrites, syn_chan, nm_to_pixels): #list of dends --> 2d array of data with cols ["Num", "Channel", "0-25 (nm)", "25-50 (nm)", "50-75 (nm)", "75-100 (nm)", "100-125 (nm)", "<125-1(nm)"]
  data = []
  for dendrite in dendrites:
    for channel in dendrite.spot_coords.keys():
      if channel != syn_chan:
        spots = [(x/nm_to_pixels, y/nm_to_pixels) for x, y in dendrite.spot_coords[channel]]
        skeleton = dendrite.skeleton
        skel_points = np.argwhere(skeleton == 255)
        endpoints = findEndpoints(skeleton, skel_points) #find endpoints of skeleton
        if len(endpoints) == 2:
          skel_points, skel_distances = orderSkelPoints(skeleton, skel_points, endpoints[0])
          mapped_spots = [mapSig(spot, skel_points) for spot in spots] #map mrna coords to skeleton
          spot_distances = []

          #calculate distance of mappen mrna coord along skeleton
          for spot in mapped_spots: 
            k = np.argwhere(np.all(skel_points == spot, axis=1))[0][0]
            spot_distances.append((np.sum(skel_distances[0:k+1]*nm_to_pixels/1000)))

          skel_length= np.sum(skel_distances)
          spot_distances = findSoma(spot_distances, skel_length, nm_to_pixels)
          binned_distances, bins = np.histogram(spot_distances, bins=[0,25,50,75,100,125,150, inf])
          data.append([dendrite.num, channel] + list(binned_distances))
        else:
          print("\nSkeleton of " + dendrite.name() + " does not have 2 endpoints, so skipped it in distribution analysis")
  
  cols = ["Dendrite Num", "Channel", "0-25 (um)", "25-50 (um)", "50-75 (um)", "75-100 (um)", "100-125 (um)", "125-150 (um)", ">= 150 (um)"]
  df = pd.DataFrame(data, columns=cols)
  
  return df


def colocAnalysis(dendrites, chan1, chan2, max_dist=300, incre=20, nanometers_per_pixel=107.5):
  
  chan1_to_chan2 = []
  chan2_to_chan1 = []
  chan1_to_chan1 = []
  chan2_to_chan2 = []
  mrnasSim1_to_mrnasSim2 = []
  mrnasSim2_to_mrnasSim1 = []
  
  for dendrite in dendrites:
    if chan1 in dendrite.spot_coords.keys() and chan2 in dendrite.spot_coords.keys():
      #get spot_coords
      mrnas1 = dendrite.spot_coords[chan1]
      mrnas2 = dendrite.spot_coords[chan2]

      #get unsimulated coloc data
      chan1_to_chan2 += spots_to_spots(mrnas1, mrnas2)
      chan2_to_chan1 += spots_to_spots(mrnas2, mrnas1)
      chan1_to_chan1 += spots_to_spots(mrnas1, mrnas1, same=True)
      chan2_to_chan2 += spots_to_spots(mrnas2, mrnas2, same=True)
      
      print_points = np.argwhere(dendrite.print == 255)

      #simulate mrna coords and find colocaliztion data 50 times
      for i in range(100):
        mrnasSim1 = simulateSpots(print_points, len(mrnas1), nanometers_per_pixel)
        mrnasSim2 = simulateSpots(dendrite.print, len(mrnas2), nanometers_per_pixel)
        mrnasSim1_to_mrnasSim2 += spots_to_spots(mrnasSim1, mrnasSim2)
        mrnasSim2_to_mrnasSim1 += spots_to_spots(mrnasSim2, mrnasSim1)
      
  bins = [i for i in range(0, max_dist+1, incre)] + [inf]
  bin_labels = [str(bins[i]) + "-" + str(bins[i+1]) for i in range(len(bins)-1)]

  chan1_to_chan2 = np.histogram(chan1_to_chan2, bins=bins)[0]
  chan2_to_chan1 = np.histogram(chan2_to_chan1, bins=bins)[0]
  chan1_to_chan1 = np.histogram(chan1_to_chan1, bins=bins)[0]
  chan2_to_chan2 = np.histogram(chan2_to_chan2, bins=bins)[0]
  mrnasSim1_to_mrnasSim2 = np.histogram(mrnasSim1_to_mrnasSim2, bins=bins)[0]/100
  mrnasSim2_to_mrnasSim1 = np.histogram(mrnasSim2_to_mrnasSim1, bins=bins)[0]/100

  data = np.concatenate(([bin_labels],
                        [chan1_to_chan2],
                        [chan2_to_chan1],
                        [chan1_to_chan1],
                        [chan2_to_chan2],
                        [mrnasSim1_to_mrnasSim2],
                        [mrnasSim2_to_mrnasSim1]), axis=0, dtype=object).T

  cols = ["Distance (nm)",
          chan1 + " to closest " + chan2,
          chan2 + " to closest " + chan1, 
          chan1 + " to closest " + chan1,
          chan2 + " to closest " + chan2,
          "Sim-"+chan1 + " to closest " + "Sim-"+chan2,
          "Sim-"+chan2 + " to closest " + "Sim-"+chan1]
  
  
  data_df = pd.DataFrame(data, columns=cols) # add dataframe 

  return data_df
    
    
  
def synAnalysis(dendrites, syn_chan, thresh, nanometers_per_pixel=107.5):
  mrna_counter_dict = defaultdict(Counter) #{mRNA_chan: (mrna_Counter, ctrl_Counter)}
  ctrl_counter_dict = defaultdict(Counter)
  for dendrite in dendrites:
    if syn_chan in dendrite.spot_coords.keys():
      synapses = dendrite.spot_coords[syn_chan]
      for chan, mrnas in dendrite.spot_coords.items():
        if chan!= syn_chan:
          ctrl_runs = Counter()
          mrna_counts = synapMRNACount(synapses, mrnas, thresh)

          print_points = np.argwhere(dendrite.print == 255)
          for i in range(100):
            ctrl_mrnas = simulateSpots(print_points, len(mrnas), nanometers_per_pixel)
            ctrl_counts = synapMRNACount(synapses, ctrl_mrnas, thresh)
            ctrl_runs.update(ctrl_counts)
          
          for k in ctrl_runs: ctrl_runs[k] /= 100 #Divide ctrl_runs bc we did 100 simulations
        
          #Update counter dict
          mrna_counter = mrna_counter_dict[chan]
          ctrl_counter = ctrl_counter_dict[chan]
          mrna_counter.update(mrna_counts)
          ctrl_counter.update(ctrl_runs)
          mrna_counter_dict[chan] = mrna_counter
          ctrl_counter_dict[chan] = ctrl_counter

    else:
      print("No synapse coords for " + dendrite.name())
      print("Will skip this dendrite")
    
  dfs = []
  for mrna_chan, mrna_counter in mrna_counter_dict.items():
    ctrl_counter = ctrl_counter_dict[mrna_chan]
    max_num = max(list(mrna_counter.keys()) + list(ctrl_counter.keys()))
    mrna_data = [0 if i not in mrna_counter.keys() else mrna_counter[i] for i in range(max_num + 1)]
    ctrl_data = [0 if i not in ctrl_counter.keys() else ctrl_counter[i] for i in range(max_num + 1)]
    counts = [i for i in range(max_num+1)]
    data = np.array([counts, mrna_data, ctrl_data]).T
    df = pd.DataFrame(data, columns=["Num of mRNA", "Real " + mrna_chan, "Sim " + mrna_chan])
    dfs.append(df)
  
  if len(dfs) != 0:
    data = pd.concat(dfs, axis = 1)
  elif len(dfs) == 1:
    data = pd.concat([dfs[0], pd.DataFrame([], columns=[3,4])], axis = 1)
  else:
    data = pd.concat([pd.DataFrame([], columns = [1,2,3]), pd.DataFrame([], columns=[4,5,6])], axis = 1)
  return data

          

    

