
"""Purpose: Function for obtaining folder paths, channel names and numerical values from user.Formats and saves txt and excel files."""

import re
import os
import glob
from socket import AddressFamily
import segmentation
from PIL import Image
from datetime import date
import pandas as pd


class MaxImg:

    def __init__(self, match):
      self.fullname = match.group(0)
      self.experiment = match.group(1)
      self.DIV = match.group(2)
      self.treatment = match.group(3)
      self.FOV = match.group(4, 5) #Maybe 4 is plate not FOV
      self.dendAnnot = [] #list accommodates multiple annotation images
      self.somAnnot = []

    def name(self):
      name =  "EXP: " + self.experiment + ", DIV: " + self.DIV + ", Treatment: " + self.treatment + ", FOV: " + self.FOV[0] +  "-" + self.FOV[1]
      return name

class AnnotImg:
  
  def __init__(self, match):
    self.fullname = match.group(0)
    self.num_colors = int(match.group(1))
    self.segmentations = []


def pt1AnalysisType(seg):
  analysis = False
  cont = True
  while cont:
    print("\n Are you analyzing " + seg + " (Y/N)?", end=" ")
    analysis = input()
    if analysis == "Y":
      cont = False
      analysis = True
    elif analysis == "N":
      cont = False
      analysis = False

  return analysis
    
  
def pt2AnalysisType(typ, cellular_comp):
  cont = True
  while cont:
    print("\n Are you doing " + typ + " analysis for " + cellular_comp + " (Y/N)?", end=" ")
    ans = input()
    if ans == "Y":
      cont = False
      ans = True
    elif ans == "N":
      cont = False
      ans = False
  return ans


def regexMatchFilter(list, pattern):
  matches = []
  for l in list:
    m = re.fullmatch(pattern, l)
    if m: matches.append(m)
  return matches


def getImagesDir():
  # Getting the image folder
  print("\n Full path to folder of (annotated) images: ", end = "")
  directory = input()
    
  while not os.path.isdir(directory):
    print("\n Incorrect path, please try again:")
    print("\n Full path to folder: ", end = "")
    directory = input()

  return directory

def getMaxDistAndIncre():
  print("\n Max distance (nm) recorded for mRNA colocalization data: ", end = "")
  distance = input()

  while not distance.isdigit() or (distance.isdigit() and int(distance) <= 0):
    print("\n Please enter an integer greater than 0: ", end="")
    distance = input()

  distance=int(distance)

  # Getting the image folder
  print("\n Name desired increment (nm) for mRNA colocalization data: ", end = "")
  incre = input()
    
  while not incre.isdigit()or (incre.isdigit() and int(incre) <= 0):
    print("\n Please enter an integer greater than 0: ", end="")
    incre = input()

  incre = int(incre)

  return distance, incre

def getThreshDist(coloc=False):
  if coloc:
      print("\n Threshold distance (nm) for mRNA to be considered colocalized to each other: ", end = "")
  else:
      print("\n Threshold distance (nm) for mRNA to be considered localized at synapse: ", end = "")
  distance = input()

  while not distance.isdigit() or (distance.isdigit() and float(distance) <= 0):
    print("\n Please enter an decimal greater than 0: ", end="")
    distance = input()

  distance=float(distance)

  return distance

def getmRNAChans(dapi, map2, dir=""):
  mRNA_channels = []
  stop = False
  while not stop:
    mRNA_channel = getChannelName("mrna", dir=dir)
    if mRNA_channel == "0":
      stop = True
    elif mRNA_channel in mRNA_channels + [dapi, map2]:
      print("\n Channel already inputted, please enter new channel name \n")
    else: 
      mRNA_channels.append(mRNA_channel)
  return mRNA_channels

def getSynapChan(dapi, map2, mrna_chans):
  prev_chans = [dapi, map2] + mrna_chans
  synap_channel = ""
  cont = True
  while cont:
    print("\n Is there a synapse channel (Y/N): ", end = " ")
    synap = input()
    if synap == "Y":
        synap_channel = getChannelName("synapse")
        cont = False
        while synap_channel in prev_chans:
          print("\n Channel already inputted, please enter new channel name \n")
          synap = getChannelName("synapse")
    elif synap == "N":
        cont = False
  return synap_channel

def getChannelName(channel_type, dir=""):
  if channel_type == "mrna":
      print("\n Name of " + channel_type + " channel (if no more mrnas, press 0): ", end = "")
  else:
      print("\n Name of " + channel_type + " channel: ", end = "")
  channel = input()
  ''' 
  if channel_type != "mrna" or channel != "0":
    while len(glob.glob(dir + "*_" + channel + ".tif")) == 0 and channel != "0":
      print("\n " + channel + " is an incorrect channel name. Try again: ", end = "")
      print("\n Name of " + channel_type + " channel: ", end = "")
      channel = input()
  '''
  return channel



def processResultsFile(match, dir, channel, divTreat_dict, soma=False):
  
  file = open(dir + match.group(0), "r")
  '''
  #extract DIV/TREAT/FOV info from file_name
  filename = match.group(0)
  div_Treat = "DIV" + match.group(2) + " " + match.group(3)
  FOV = match.group(4) + "-" + match.group(5)
  '''
  filename = match.group(0)
  #print("match.group(4)", match.group(4))
  
  date = match.group(1) # should be same for all images in an experiment, so not saved in dendrite object attributes
  
  Treat1_Treat2 = match.group(2) + "_" +match.group(3)
  
  FOV = match.group(4) + "-" +match.group(5)
  
  # Nested dictionaries of this format:   {DIV TREAT : {FOV : {num : soma/dend object}}}
  if Treat1_Treat2 not in divTreat_dict.keys():
    divTreat_dict[Treat1_Treat2] = {FOV : {}}
  elif FOV not in divTreat_dict[Treat1_Treat2].keys():
    divTreat_dict[Treat1_Treat2][FOV] = {}

  #parse result txt file
  lines = file.read()
  cells = lines.split("CELL_START\tCell_")[1:]

  for cell in cells:
    mrna_coords = []
    num = cell.split()[0] ##### may be double digits
    if "SPOTS_START" in cell:
      spot_info = cell.split("SPOTS_START")[1]
      spot_lines = spot_info.split("\n")[2:]
      for line in spot_lines:
        if line.startswith("SPOTS_END"): break
        else:
          coords = line.split("\t")
          try: 
            x_coord = float(coords[0])
            y_coord = float(coords[1])
            mrna_coords.append((x_coord, y_coord))
          except ValueError:
            print("\nError with line in results file: " + filename + " in " + dir)
            print("Problematic line: " + line)
    
    if num not in divTreat_dict[Treat1_Treat2][FOV].keys(): #add soma/dend number if not already in nested dicts
      if soma:
        seg = segmentation.Soma(filename)
      else:
        seg = segmentation.Dendrite(filename)
      seg.num = num
      seg.spot_coords[channel] = mrna_coords
      divTreat_dict[Treat1_Treat2][FOV][num] = seg
    else:
      divTreat_dict[Treat1_Treat2][FOV][num].spot_coords[channel] = mrna_coords #add coords to nested dict

  file.close()



#Records all outlines of dendrites in txt file in FishQuant v1 (Matlab) format
def getOutlinesTxts(channel_img, dapi_img, directory, segmentations, nuc=False, synap=False):
  
  if len(segmentations) != 0:

    outline_txt = open(directory + "/" + channel_img[:-4] + "__outline.txt", "w")

    outline_txt.write("FISH-QUANT\n")
    outline_txt.write("File-version   3D_v1\n")
    outline_txt.write("RESULTS OF SPOT DETECTION PERFORMED ON " + date.today().strftime("%d-%b-%Y") + "\n")
    outline_txt.write("COMMENT	Automated outline definition (batch or quick-save)\n")
    outline_txt.write("IMG_Raw\t" + channel_img + "\n")
    outline_txt.write("IMG_Filtered\n")
    outline_txt.write("IMG_DAPI\t" + dapi_img + "\n")
    outline_txt.write("IMG_TS_label\n")
    outline_txt.write("FILE_settings\n")
    outline_txt.write("PARAMETERS\n")
    outline_txt.write("Pix-XY	Pix-Z	RI	Ex	Em	NA	Type\n")
    outline_txt.write("107.5	300	1.35	547	583	1.4	widefield\n")

    for segmentation in segmentations:
      if synap:
        outline = segmentation.synap_outline
      else:
        outline = segmentation.outline
      i = segmentation.num
      if len(outline) != 0:
        outline_txt.write("CELL_START\tCell_" + str(i) + "\n")

        outline_txt.write("X_POS") # write x coordinate
        for i in range(len(outline)):
          outline_txt.write("\t" + str(int(outline[i][1])))
        outline_txt.write("\t")

        outline_txt.write("\nY_POS") # write x coordinate
        for i in range(len(outline)):
          outline_txt.write("\t" + str(int(outline[i][0])))
        outline_txt.write("\t")

      outline_txt.write("\nZ_Pos\t\n")
      outline_txt.write("CELL_END\n")

      if nuc:
        nuc_outline = segmentation.nucleus.outline
        if len(outline) != 0:
          outline_txt.write("Nucleus_START\tNucleus_manual\n")

          outline_txt.write("X_POS")
          for i in range(len(nuc_outline)):
            outline_txt.write("\t" + str(int(nuc_outline[i][1])))
          outline_txt.write("\t")

          outline_txt.write("\nY_POS")
          for i in range(len(nuc_outline)):
            outline_txt.write("\t" + str(int(nuc_outline[i][0])))
          outline_txt.write("\t")

          outline_txt.write("\nZ_Pos\t\n")
          outline_txt.write("Nucleus_END\n")

    outline_txt.close()

  return

def saveImg(img, location, tag):
  Image.fromarray(img).save(location + tag + ".gif")


def saveDfDictToExcel(data_dict, name):
  with pd.ExcelWriter(name) as writer:
    for div_treat, df in data_dict.items():
      cols = df.columns.values.tolist()
      df.to_excel(writer, sheet_name=div_treat, columns=cols, index=False)
