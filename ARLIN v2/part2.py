#!/usr/bin/env python3

__author__ = "Zoe Wefers"
__credits__ = ["Ryan Huang"]
__version__ = "1.0.0"
__maintainer__ = "Zoe Wefers"
__email__ = "zoe.wefers@mail.mcgill.ca"

"""Purpose: Executable file to collect localization statistics froms users result files and skeletons/prints."""

import localization_analysis as la
import ui
import os
import glob
import re
import segmentation
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle

ExpName = 'Exp25'
bad_dends = []
dends_byDIVTreat = {}

def main():
    

    nanometers_per_pixel = 107.5
    ########################### GET USER INFO ######################################################
    
    images_dir = ui.getImagesDir()
    
    map2 = ui.getChannelName("MAP2", dir=images_dir + "/")
    dapi = ui.getChannelName("DAPI", dir=images_dir + "/")

    mRNA_channels = ui.getmRNAChans(dapi, map2, dir=images_dir + "/")

    density_dend = ui.pt2AnalysisType("density", "dendrites")
    density_som = ui.pt2AnalysisType("density", "somas")
    
    distr_analysis = True #distribution analysis is mandatory to get mrnas sorted and to do any coloc or synpase analysis
    
    coloc_analysis = ui.pt2AnalysisType("colocalization", "dendrites")
    
    if coloc_analysis:
        coloc_mrna1_chan = ui.getChannelName("first colocalization mRNA", dir=images_dir+ "/")
        coloc_mrna2_chan = ui.getChannelName("second colocalizatoin mRNA", dir=images_dir+ "/")
        max_dist, incre = ui.getMaxDistAndIncre()
        
        
    syn_analysis = ui.pt2AnalysisType("synapse", "dendrites")
    if syn_analysis:
        syn_chan = ui.getChannelName("synapse", dir=images_dir + "/")
        
        threshold_dist = ui.getThreshDist()
        
        
    else:
        syn_chan = None
    ###################################################################################################


    skel_print_dir = images_dir + "/SkeletonsAndPrints"
    
    #initializing dict of somas with spot detection files
    somas_byDIVTreatFOV = {} #Nested dictionaries of this format:   {DIV TREAT : {FOV : {num : dend object}}}
    if density_som:
        for channel in mRNA_channels:
            result_dir = images_dir + "/Results Soma " + channel + "/"
            results_pattern = re.compile("(.*?)_(.*?)_(.*?)_(.*?)_(?:xy|XY)([0-9]?[0-9])_" + channel + "__outline_spots.*?.txt")
            if not os.path.isdir(result_dir):
                print("No folder called " + result_dir)
            else:
                results_matches = ui.regexMatchFilter(os.listdir(result_dir), results_pattern)
                for match in results_matches:
                    ui.processResultsFile(match, result_dir, channel, somas_byDIVTreatFOV, soma=True)
        
        # adding soma prints
        somas_byDIVTreat = {}
        for div_treat in somas_byDIVTreatFOV.keys():
            FOV_dict = somas_byDIVTreatFOV[div_treat]
            somas = []
            for FOV in FOV_dict.keys():
                for num in FOV_dict[FOV].keys():
                    soma = FOV_dict[FOV][num]
                    filename = FOV_dict[FOV][num].img_name
                    filename = filename.split("_")[0:-4] #delete "_chan___outline_spots.txt"
                    filename = "MAX_" + "_".join(filename) + "_" + map2 + "_somaprint_" + str(num) + ".gif"
                    print_files = glob.glob(skel_print_dir + "/" + filename)
                    if len(print_files) == 1:
                        print_file = print_files[0]
                        #soma.print = plt.imread(print_file)
                        soma.print = np.squeeze(plt.imread(print_file))[:,:,0]
                        somas.append(soma)
                    elif len(print_files) == 0:
                        print("\nMissing prints for soma #" + str(num) + " for " + div_treat + " for FOV " + FOV)
                    else:
                        print("\nMulitple prints for soma #" + str(num) + " for " + div_treat + " for FOV " + FOV)
            somas_byDIVTreat[div_treat] = somas


    #initializing dict of dendrites with spot detection files
    dends_byDIVTreatFOV = {}  #Nested dictionaries of this format:   {DIV TREAT : {FOV : {num : dend object}}}
    for channel in mRNA_channels:
        result_dir = images_dir + "/Results Dendrite " + channel + "/"
        results_pattern = re.compile("(.*?)_(.*?)_(.*?)_(.*?)_(?:xy|XY)([0-9]?[0-9])_" + channel + "__outline_spots.*?.txt")
        if not os.path.isdir(result_dir):
                print("No folder called " + result_dir)
        else:
            results_matches = ui.regexMatchFilter(os.listdir(result_dir), results_pattern)
            
            for match in results_matches:
                ui.processResultsFile(match, result_dir, channel, dends_byDIVTreatFOV)
    
    #Adding synapse coords if doing synapse analysis
    if syn_analysis:
        result_dir = images_dir + "/Results Dendrite " + syn_chan + "/"
        results_pattern = re.compile("(.*?)_(.*?)_(.*?)_(.*?)_(?:xy|XY)([0-9]?[0-9])_" + channel + "__outline_spots.*?.txt")
        if not os.path.isdir(result_dir):
                print("No folder called " + result_dir)
        else:
            results_matches = ui.regexMatchFilter(os.listdir(result_dir), results_pattern)
            for match in results_matches:
                ui.processResultsFile(match, result_dir, syn_chan, dends_byDIVTreatFOV)
    
    #Adding skels and prints to dendrite
    #dends_byDIVTreat = {}
    for div_treat in dends_byDIVTreatFOV.keys():
        FOV_dict = dends_byDIVTreatFOV[div_treat]
        dendrites = []
        for FOV in FOV_dict.keys():
            for num in FOV_dict[FOV].keys():
                
                dendrite = FOV_dict[FOV][num]
                filename = FOV_dict[FOV][num].img_name
                filename = filename.split("_")[0:-4] #delete "_chan___outline_spots.txt"
                
                printname = "MAX_" + "_".join(filename) + "_" + map2 + "_print_" + str(num) + ".gif"
                skelname = "MAX_" + "_".join(filename) + "_" + map2 + "_skel_" + str(num) + ".gif"
                print_files = glob.glob(skel_print_dir + "/" + printname)
                skel_files = glob.glob(skel_print_dir + "/" + skelname)
                
                if len(print_files) == 1 and len(skel_files)==1: #adding print/skel to dendrite
                    print_file = print_files[0]
                    skel_file = skel_files[0]
                    
                    file = plt.imread(print_file)
                    if file.shape == (2304, 2304):
                        dendrite.print = plt.imread(print_file)
                        dendrite.skeleton = plt.imread(skel_file)
                    else:
                        dendrite.print = np.squeeze(plt.imread(print_file))[:,:,0]
                        dendrite.skeleton = np.squeeze(plt.imread(skel_file))[:,:,0]
                    
                    
                    #print("print shape: ", dendrite.print.shape)
                    #print("skel shape: ", dendrite.skeleton.shape)
                    
                    dendrites.append(dendrite)
                    
                elif len(print_files) == 0:
                    print("\nMissing dendrite print #" + str(num) + " for " + div_treat + " for FOV " + FOV)
                
                elif len(skel_files) == 0:
                    print("\nMissing dendrite skeleton #" + str(num) + " for " + div_treat + " for FOV " + FOV)
                else:
                    print("\nMulitple prints/skeletons for dendrite #" + str(num) + " for " + div_treat + " for FOV " + FOV)
        dends_byDIVTreat[div_treat] = dendrites
        print("Done reading prints and skels")
        

    #density analysis
    if density_som:
        print("Calculating soma densitites...")
        densities_byDIVTreat = la.calculate_densities(somas_byDIVTreat, syn_chan, nanometers_per_pixel=nanometers_per_pixel)
        ui.saveDfDictToExcel(densities_byDIVTreat, "SomaDensities.xlsx")
    
    if density_dend:
        print("Calculating dendrite densities...")
        densities_byDIVTreat = la.calculate_densities(dends_byDIVTreat, syn_chan, nanometers_per_pixel=nanometers_per_pixel)
        ui.saveDfDictToExcel(densities_byDIVTreat, "DendriteDensities.xlsx")


    #distribution analysis
    if distr_analysis:
        #the distrAnalysis function also sorts all mrna into their zones, saved in dendrite.spot_coords_sorted
        #this is necessary for coloc and distr analysis
        
        print("Calculating distribution along dendrites...")
        data_dict = {} #{div_treat : dataframe of distance along skeleton}
        for sheet_name, dendrites in dends_byDIVTreat.items():
            
            
            df = la.distrAnalysis(dendrites, syn_chan, nanometers_per_pixel) #this function also sorts mrna into dendrite.spot_coords_sorted dictionary
            data_dict[sheet_name] = df
            
        ui.saveDfDictToExcel(data_dict, "DistrAnalysis.xlsx")

    for sheet_name, dendrites in dends_byDIVTreat.items():
        for d in dendrites:
            if d.start_end is None:
                print("Removing ", d.img_name," Num: ", d.num)
                bad_dends.append(d)
                dendrites.remove(d)
                
                
    #synaptic localization analysis
    if syn_analysis:
        #three synapse analysis done here.
        
        #Analyze 
        print("Calculating degree of localization at synapses...")
        data_dict = {}
        data_dict2 = {}
        data_dict_mRNAtoSyn = {}
          
        
        for sheet_name, dendrites in dends_byDIVTreat.items():
            
            print("compartmentalizing", sheet_name)
            la.dendriteCompartment(dendrites, nanometers_per_pixel)
            
            
            print("Synapse Analysis", sheet_name)
            print(sheet_name)
            
            df, df2 = la.synAnalysis(dendrites, syn_chan, threshold_dist, nanometers_per_pixel)
            data_dict[sheet_name] = df
            data_dict2[sheet_name] = df2
            
            df = la.mRNA_to_synAnalysis(dendrites, syn_chan, mRNA_channels, threshold_dist, nanometers_per_pixel=107.5)
            data_dict_mRNAtoSyn[sheet_name] = df
            
            
        ui.saveDfDictToExcel(data_dict, "SynAnalysis.xlsx")
        ui.saveDfDictToExcel(data_dict2, "SynAnalysis2.xlsx")
        ui.saveDfDictToExcel(data_dict_mRNAtoSyn, "mRNAtoSynapseAnalysis.xlsx")
        
    #mRNA colocalization analysis
    if coloc_analysis:
        print("Calculating degree of colocalization...")
                
        coloc_data_dict = {}
        
        for sheet_name, dendrites in dends_byDIVTreat.items():
            
            print(sheet_name)
            
            #If SynAnalysis was not done, dendriteCompartments so that the dendrite.print_distribution variable is initialized before la.colocAnalysis is called
            if not syn_analysis:
                print("compartmentalizing...")
                la.dendriteCompartment(dendrites, nanometers_per_pixel)
                 
            print("coloc...")
            print(sheet_name)
            
            df = la.colocAnalysis(dendrites, coloc_mrna1_chan, coloc_mrna2_chan, max_dist, incre, nanometers_per_pixel)
            
            coloc_data_dict[sheet_name] = df
            
        
        ui.saveDfDictToExcel(coloc_data_dict, "ColocAnalysisByZoneAndDend.xlsx")
    
    

        
    return dends_byDIVTreat

dic = main()


if __name__ == "__main__2":
    
    dic = main()
    
    for sheet_name, dendrites in dic.items():
            for d in dendrites:
                if d.start_end is None:
                    print("Removing ", d.img_name," Num: ", d.num)
                    bad_dends.append(d)
                    dendrites.remove(d)
                
    
    with open(ExpName+".pkl","wb") as f:
        pickle.dump(dic, f, protocol=pickle.HIGHEST_PROTOCOL)