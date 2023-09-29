#!/usr/bin/env python3

__author__ = "Zoe Wefers"
__credits__ = ["Ryan Huang"]
__version__ = "1.0.0"
__maintainer__ = "Zoe Wefers"
__email__ = "zoe.wefers@mail.mcgill.ca"

""""Purpose: Executable script which reads annotated and non-annotated MAP2 max projections gif images and generates a print for every
annnotated dendrite, soma (and nucleus of each soma). Also generates a 1-pixel wide skeleton for every annotated dendrite."""

import imgProcessing as ip
import ui
import os
import glob
import re
import segmentation


def main():
    #####Set colors######
    colors = {"red": [255,0,0], 
            "green": [0,255,0], 
            "blue" : [0,0,255], 
            "orange" : [240, 134, 51], 
            "yellow" : [255, 255, 0], 
            "purple" : [143, 57, 182],
            "teal" : [130, 210, 204],
            "mint" : [214, 253, 208],
            "salmon": [255, 128, 102]
            }
    ######################################################################################
    
    
        
    ########################### Getting Channel Names #################################
    #Making new subfolder to store skeletons and prints
    
    dend_analysis = ui.pt1AnalysisType("dendrites")
    soma_analysis = ui.pt1AnalysisType("somas")
    nucleus_analysis = False
    if soma_analysis:
        nucleus_analysis = ui.pt1AnalysisType("nuclei")

    images_dir = ui.getImagesDir()
    os.chdir(images_dir)

    print_skel_dir = "SkeletonsAndPrints"
    if not os.path.exists(print_skel_dir):
        os.makedirs(print_skel_dir)

    map2 = ui.getChannelName("MAP2")
    dapi = ui.getChannelName("DAPI")

    mRNA_channels = ui.getmRNAChans(dapi, map2)
            
    synap_channel = ui.getSynapChan(dapi, map2, mRNA_channels)
    ###################################################################################



    ######################## Getting max images and annotation images #################
    maxfile_pattern = re.compile(".*MAX_(.*?)_(?:div|DIV)([0-9]?[0-9])_(.*?)_.*_(.*?)_(?:xy|XY)([0-9]?[0-9])_" + map2 + ".gif")
    maxfile_matches = ui.regexMatchFilter(os.listdir(), maxfile_pattern)

    missing_data = []
    bad_prints_skels = []

    if len(maxfile_matches) == 0:
        print ("Could not find any Max images")
    else:
        dend_count = 1
        soma_count = 1
        for file_match in maxfile_matches:
    
            maxImg = ui.MaxImg(file_match)
            max_filename = maxImg.fullname[:-4]
            print("\n" + max_filename)

            annot_files = glob.glob(max_filename + "_*.gif") # grab all possible annotation files
            base_loc = print_skel_dir + "/" + maxImg.fullname[:-4] #base name for prints/skels

            if dend_analysis:
                dendAnnot_pattern = re.compile(max_filename + "_([0-9]?[0-9])_dendAnnot.gif")
                dendAnnot_imgs = ui.regexMatchFilter(annot_files, dendAnnot_pattern)
                if len(dendAnnot_imgs) == 0:
                    print("No annotated dendrites for " + maxImg.name())
                elif len(dendAnnot_imgs) > 1:
                    print("More than one annotated dendrite file for " + maxImg.name())
                else:
                    maxImg.dendAnnot = ui.AnnotImg(dendAnnot_imgs[0])

                    #processing dendrites and saving skeletons/prints
                    for color_name, color in list(colors.items())[0:maxImg.dendAnnot.num_colors]:
                        dendrite = segmentation.Dendrite(maxImg.fullname)
                        dendrite.color = color_name
                        dendrite.num = dend_count
                        wrongColor = ip.colorDetect(maxImg.fullname, maxImg.dendAnnot.fullname, color, dendrite)
                        if wrongColor:
                            print("\nWrong color used for " + dendrite.name())
                            bad_prints_skels.append(dendrite.name())
                        else:
                            problem = ip.getPrintAndSkeleton(dendrite)
                            if len(synap_channel) != 0:
                                    ip.getSynapPrint(dendrite)
                            if problem:
                                print("Problem with skeleton/print of " + dendrite.name())
                                bad_prints_skels.append(dendrite.name())
                            else:
                                ip.getOutline(dendrite)
                                maxImg.dendAnnot.segmentations.append(dendrite)
                                if len(synap_channel) != 0:
                                    ip.getOutline(dendrite, True)

                            #save skeletons/prints regardless of problems for visualisation of problem
                            ui.saveImg(dendrite.skeleton, base_loc, "_skel_" + str(dendrite.num)) #add num to tag
                            ui.saveImg(dendrite.print, base_loc, "_print_"  + str(dendrite.num))#add num to tag
                            if len(synap_channel) != 0:
                                ui.saveImg(dendrite.synap_print, base_loc, "_synprint_" + str(dendrite.num))#add num to tag
                        
                        dend_count +=1
            if soma_analysis:
                somAnnot_pattern = re.compile(max_filename + "_([0-9]?[0-9])_somAnnot.gif")
                somAnnot_imgs = ui.regexMatchFilter(annot_files, somAnnot_pattern)

                if len(somAnnot_imgs) == 0:
                    print("No annotated somas for " + maxImg.name())
                elif len(somAnnot_imgs) > 1:
                    print("More than one annotated soma file for " + maxImg.name())
                else:
                    maxImg.somAnnot = ui.AnnotImg(somAnnot_imgs[0])

                    #processing somas and saving prints
                    for color_name, color in list(colors.items())[0:maxImg.somAnnot.num_colors]:
                        soma = segmentation.Soma(maxImg.fullname)
                        soma.color = color_name
                        soma.num = soma_count
                        wrongColor = ip.colorDetect(maxImg.fullname, maxImg.somAnnot.fullname, color, soma)
                        if wrongColor:
                            print("\nWrong color used for " + soma.name())
                            bad_prints_skels.append(soma.name())
                        else:
                            problem = ip.refineSomaOrNucPrint(soma)
                            if problem:
                                print("\n Problem with print for " + soma.name())
                                print("--> Will not add " + soma.name() + " to txt outlines")
                                bad_prints_skels.append(soma.name())
                            else:
                                if nucleus_analysis:
                                    end_index = maxImg.fullname.index(map2)
                                    dapi_gifs = glob.glob(maxImg.fullname[0:end_index] + dapi + ".gif")
                                    if len(dapi_gifs) != 1:
                                        print("\nMissing (or extra) DAPI gif for " + maxImg.name())
                                        print("--> Will not add " + soma.name() + " to txt outlines")
                                    else:
                                        dapi_gif = dapi_gifs[0]
                                        nucleus = segmentation.Nucleus(maxImg.fullname)
                                        nucleus.num = soma_count
                                        missing_nuc = ip.getNucleus(soma.print, dapi_gif, nucleus)
                                        if missing_nuc:
                                            print("\n"+ soma.name() + " has no visible nucleus")
                                            print("--> Will not add " + soma.name() + " to txt outlines")
                                            bad_prints_skels.append(nucleus.name())
                                        else:
                                            problem = ip.refineSomaOrNucPrint(nucleus)
                                            if problem:
                                                print("\nProblem with print of " + nucleus.name() + " from " + soma.name())
                                                print("--> Will not add " + soma.name() + " or " + nucleus.name() + " to txt outlines")
                                                bad_prints_skels.append(nucleus.name())
                                            else:
                                                ip.getOutline(soma, False)
                                                ip.getOutline(nucleus, False)
                                                soma.nucleus = nucleus
                                                maxImg.somAnnot.segmentations.append(soma)

                                        ui.saveImg(nucleus.print, base_loc, "_nucprint_" + str(nucleus.num))

                                else:
                                    ip.getOutline(soma, False)
                                    maxImg.somAnnot.segmentations.append(soma)

                            ui.saveImg(soma.print, base_loc, "_somaprint_" + str(soma.num))
                        soma_count +=1

                    #get outline txts
                    if len(synap_channel != 0):
                        channels = mRNA_channels + [synap_channel]
                    for channel in channels:

                        #Get DAPI and mRNA images
                        start_index = maxImg.fullname.index(maxImg.experiment)
                        end_index = maxImg.fullname.index(map2)
                        mrna_imgs = glob.glob(maxImg.fullname[start_index:end_index] + channel + ".tif")
                        dapi_imgs = glob.glob(maxImg.fullname[start_index:end_index] + dapi + ".tif")

                        if len(mrna_imgs) != 1:
                            print("Missing (or extra)" + channel + " image for " + maxImg.name())
                            missing_data.append((maxImg.name(), channel))
                        elif len(dapi_imgs) != 1:
                            print("Missing (or extra) DAPI image for " + maxImg.name())
                            missing_data.append((maxImg, "DAPI"))
                        else:
                            mrna_img = mrna_imgs[0]
                            dapi_img = dapi_imgs[0]


                        if dend_analysis:
                            #make directories to save outlines
                            dend_outline_dir = "Outlines Dendrite " + channel
                            if not os.path.exists(dend_outline_dir):
                                os.makedirs(dend_outline_dir) 
                            ui.getOutlinesTxts(mrna_img, dapi_img, dend_outline_dir, maxImg.dendAnnot.segmentations, synap=True) #write outline text files

                        if soma_analysis and channel != synap_channel:
                            soma_outline_dir = "Outlines Soma " + channel
                            if not os.path.exists(soma_outline_dir):
                                os.makedirs(soma_outline_dir)
                            ui.getOutlinesTxts(mrna_img, dapi_img, soma_outline_dir, maxImg.somAnnot.segmentations, nuc=nucleus_analysis) #write outline text files
    ##################################################################################################################################


    #TO DO: PRINT REPORTS ABOUT MISSING_DATA AND BAD_PRINTS_SKELS (put in ui.py)
     
     
     

main()