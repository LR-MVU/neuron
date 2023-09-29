__author__ = "Zoe Wefers"
__credits__ = ["Ryan Huang"]
__version__ = "1.0.0"
__maintainer__ = "Zoe Wefers"
__email__ = "zoe.wefers@mail.mcgill.ca"

"""Purpose: Performs all the image cleaning, smoothing, processing etc to generate prints for every annoated cellular compartment
and skeletons for every dendrite."""

import numpy as np
from matplotlib.pyplot import imread
from skimage.morphology import skeletonize, disk, binary_dilation, binary_erosion
from skimage.filters import threshold_otsu
from skimage.measure import label
from skimage import measure
from skimage import io
from scipy import ndimage as ndi
from datetime import date
from fil_finder import FilFinder2D
import astropy.units as u



#Use color mask and otsu threshold to extract single dendrite/soma
def colorDetect(max_file, color_file, color, segmentation):

  #reading in colored image and regular image
  color_img = imread(color_file)[:,:,0:3]
  img = imread(max_file)

  io.imshow(img)
 
  mask = np.where(np.all(color_img == color, axis=-1), 1, 0)

  wrongColor = False
  if np.count_nonzero(mask) == 0:
    wrongColor = True


  if "soma" in segmentation.name(): #for soma just use annotation
    ndi.binary_fill_holes(mask).astype(int)
    mask[mask==1] = 255
    mask = mask.astype(np.uint8)
    segmentation.print = mask
    
  else:
    p = np.multiply(img, mask)

    #initial thresholding using whole image
    thresh = threshold_otsu(img)
    binary = p > thresh
    
    if np.sum(binary) < np.sum(mask)*0.5:
        binary = p
    
    binary[binary==1] = 255
    binary = binary.astype(np.uint8)

    segmentation.print = binary

  return wrongColor

#Dilates outline (for making synapse outlines)
def dilateOutline(p, r):
  radius = r
  selem = disk(radius)
  dilated_print = binary_dilation(p, selem=selem)
  return dilated_print

#Gets nuxcleus print by using dapi signal behind soma annotation
def getNucleus(soma_print, dapi_gif, nucleus):
  missing_nuc = False
  mask = np.where(soma_print == 255, 1, 0)
  dapi_img = imread(dapi_gif)
  masked_dapi = np.multiply(dapi_img, mask)
  thresh = threshold_otsu(dapi_img) * 0.65
  nuc_print = np.where(masked_dapi > thresh, 1, 0)
  if np.count_nonzero(nuc_print) == 0:
    missing_nuc = True
  nuc_print[nuc_print==1] = 255
  nucleus.print = nuc_print.astype(np.uint8)
  return missing_nuc


def refineSomaOrNucPrint(segmentation):
  p = segmentation.print
  if np.count_nonzero(p) == 0:
    problem = True
  else:
    p, problem = cleanBinaryImg(p)
    if not problem:
      p = ndi.binary_fill_holes(p).astype(int)
      p[p==1] = 255
      segmentation.print = p.astype(np.uint8)
  return problem

#Sets skeleton and print attribute for a dendrite
def getPrintAndSkeleton(dendrite):

  #initial dilation of binary image of selected dendrite
  radius = 5
  selem = disk(radius)
  p = binary_dilation(dendrite.print, selem=selem)

  #Repeatedly dilate ouline until binary image of dendrite is completely connected
  print_copy = np.copy(p)
  labeled, num_components = label(print_copy, background=0, return_num=True)

  while num_components > 1:
    for i in range(10):
      print_copy = binary_dilation(print_copy, selem=disk(4))
      print_copy = binary_erosion(print_copy, selem=disk(3))
    labeled, num_components = label(print_copy, background=0, return_num=True)

  #Skeletonize dilate outline
  skeleton = skeletonize(print_copy)
  skeleton = trimSkeleton(skeleton)
  if np.count_nonzero(skeleton) == 0:
    skel_issue = True
  else:
    skeleton, skel_issue = cleanBinaryImg(skeleton)

  #This step only changes image if gaps in binary dendrite
  #Extract pixels of dilated skeleton not in outline
  adding = binary_dilation(skeleton, selem=disk(5))
  adding = adding * (adding != p)
  p = adding + p #Adding to print to fill gaps

  #Repeatedly dilate and erode to smooth outline
  for i in range(10):
    p = binary_dilation(p, selem=disk(3))
    p = binary_erosion(p, selem=disk(3))

  if np.count_nonzero(p) == 0:
    print_issue = True
  else:
    p, print_issue = cleanBinaryImg(p) #MAYBE MOVE THIS ABOVE THE FOR LOOP!!!
  
  p[p==1] = 255
  p = p.astype(np.uint8)
  skeleton[skeleton==1] = 255
  skeleton = skeleton.astype(np.uint8)

  #update dendrite object
  dendrite.print = p
  dendrite.skeleton = skeleton

  return skel_issue or print_issue
    
#delete small circles aside print/skeleton and check if empty
def cleanBinaryImg(img):
  labeled_img, num_components = label(img, return_num=True)
  size_list = [np.count_nonzero(labeled_img == i) for i in range(1, num_components+1)]
  max_label, problem = secondLargest(size_list)
  clean_img = (np.where(labeled_img == (max_label + 1), 1, 0)) #add 1 because 0 label is background
  if not problem: return clean_img, problem #If problem save original "bad" image for visualization
  else: return img, problem

def getSynapPrint(dendrite):
  dendrite.synap_print = dilateOutline(dendrite.print, 8)


#trim off small branches on skeleton
def trimSkeleton(skeleton):
  fil = FilFinder2D(skeleton, distance=1 * u.pix, mask=skeleton)
  #fil.preprocess_image(flatten_percent=85)
  fil.create_mask(border_masking=True, verbose=False,
                      use_existing_mask=True)
  fil.medskel(verbose=False)
  fil.analyze_skeletons(branch_thresh=100* u.pix, skel_thresh=40 * u.pix, prune_criteria='length')
  #parameters can be adjusted
  skel = fil.skeleton_longpath
  skel[skel == 1] = 255  
  return skel



#Sets outline (or synapse_outline) attribute in dendrite obj
def getOutline(segmentation, syn=False):

  problem = False

  if syn: p = segmentation.synap_print
      #print error statement if trying to get synap print for soma
  else: p = segmentation.print

  contours = measure.find_contours(p, 0)
  if len(contours) == 0: #No contours found --> empty print
    print("Empty print for " + segmentation.name())
    problem = True
  elif len(contours) >= 2: #Multiple contours found
    #print("Deleting extra shapes in outline from " + segmentation.name())
    len_list = [len(contour) for contour in contours]
    arg_max, problem = secondLargest(len_list) #Pick longest contour, unless largest is relatively large
    contour = contours[arg_max]  
  else:
    contour = contours[0]

  contour = [contour[i] for i in range(len(contour)) if i % 10 == 0] #use every tenth coordinate
  if syn: segmentation.synap_outline = contour
  else: segmentation.outline = contour

  return problem


#Returns longest list (largest contour)
def find_max_list(list):
    list_len = [len(i) for i in list]
    return list[np.argmax(np.array(list_len))]


#returns label of second largest component and returns True if area is more than 20% the area of largest component
def secondLargest(list):
  arg_max = np.argmax(np.array(list))
  max_size = list[arg_max]
  count = 0
  for size in list:
    if size >= max_size * 0.2:
      count +=1
  return arg_max, (count > 1)