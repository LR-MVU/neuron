#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 18 16:52:57 2023

@author: stellaxu
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 30 08:45:48 2023

@author: stellaxu
"""
# -*- coding: utf-8 -*-
"""
Created on Wed May 3 14:06:43 2023

@author: Lokha, Sethu
"""

import numpy as np
import math
import os
import re


channel1Path = input("what is the full path of your mRNA channel? ") #mRNA channel output file from fishquant - summary file all spots 
channel2Path = input("what is the full path of your peptide channel? ") #peptide channel output file from fishquant - before threshold cutoff

ch1=np.loadtxt(channel1Path, skiprows=14, dtype=str)   
ch2=np.loadtxt(channel2Path, skiprows=14, dtype=str)   

knowMean = input("What is the mean of peptide? If you don't know enter 0 ")


if (knowMean!="0"):
    mean = float(knowMean)
if (knowMean=="0"):
    thresholdedSpotFile = input("What is the full path of the summary thresholded spots? ")
    mean_estimation_file=np.loadtxt(thresholdedSpotFile, skiprows=14,dtype=str) #peptide channel output file from fishquant for SINGLES
    mean=np.mean(np.array(mean_estimation_file[:,5].astype(float)))

print("Mean peptide intensity is: ", mean)
#print(ch1.shape) #Ch1- mRNA, prints number of rows and columns with data entries in the mRNA fishquant results
#print(ch2.shape) #Ch2- nascent peptide, prints number of rows and columns with data entries in the peptide fishquant results

translating={}

thresholdedDistance = input("what is the distance threshold between mRNA and peptide?")   #if the distance b/w peptide and mRNA is less than or equal to 100
outputFileName = input("Save the output file as:")
file=open(outputFileName+"_fishQuantTransltion_Results.txt",'w+')
file.write('\t'.join(['image_name','cell_number','mRNA_entry','peptide_entry','pep_amp','normalized_pep_amp','distance','x-coordinate','y-coordinate','z-coordinate'])+"\n")

for entry in range(ch1.shape[0]):
    for peptide_entry in range(ch2.shape[0]):
        #print(peptide_entry)
        if (ch1[entry,1]==ch2[peptide_entry,1]) and (ch1[entry,0][0:-10]==ch2[peptide_entry,0][0:-10]):
            #get the x y z coordinates from the files (located in 2,3 and 4th column)
            p1= ch1[entry,2:5].astype(float)   #mRNA
            p2= ch2[peptide_entry,2:5].astype(float)    #peptide
            #compute the distance d (P1,P2)
            image_name=ch1[entry,0].astype(str) # getting the filename from the entry for each
            cell = ch1[entry,1].astype(str)
            distance=math.sqrt(((p1[0]-p2[0]))**2+((p1[1]-p2[1]))**2+((p1[2]-p1[2]))**2)
            if distance<=float(thresholdedDistance):      #if the distance b/w peptide and mRNA is less than or equal to 100
                summary=[]
                summary.append(str(image_name))
                summary.append(str(cell))
                summary.append(str(entry))
                summary.append(str(peptide_entry))
                summary.append(str(ch2[peptide_entry,5]))
                summary.append(str(ch2[peptide_entry,5].astype(float)/mean))
                summary.append(str(distance))
                summary.append(str(p1[1])) # x-coordinate
                summary.append(str(p1[0])) # y-coordinate
                summary.append(str(p1[2])) # z-coordinate
                translating[entry]=[summary]
                file.write("\t".join(summary)+"\n")
file.close()

################################remove repeated mRNA, based the brightness ################################
outputfilepath = str(os.getcwd()+"/"+outputFileName+"_fishQuantTransltion_Results.txt")

summary=np.loadtxt(outputfilepath, skiprows=1, dtype=str)
uniquemRNAs=np.unique(summary[:,2])
textAmp=open("tmp_amplitude.txt","w+")
textAmp.write("\t".join(['image_name','cell_number','mRNA_entry','peptide_entry','pep_amp','normalized_pep_amp','distance'])+"\n")

for mRNA in uniquemRNAs:
    repeatedEntries=np.where((summary[:,2]==mRNA))    
    if len(repeatedEntries[0])>1:
        pep_amp =[]
        for entry in repeatedEntries[0]:
            pep_amp.append(summary[entry,-2])
        entry_id=repeatedEntries[0][np.argmin(np.array(pep_amp))]
        textAmp.write("\t".join(summary[entry_id,:]))
        textAmp.write("\n")
    else:
        entry_id=repeatedEntries[0][0]
        textAmp.write("\t".join(summary[entry_id,:]))
        textAmp.write("\n")
textAmp.close()

###############################remove repeated peptides, based on distance ###############################

summaryAmpDist=np.loadtxt("tmp_amplitude.txt", skiprows=1, dtype=str)
uniquemPEPs=np.unique(summaryAmpDist[:,3])
textDistance=open("tmp_distance.txt","w+")
textDistance.write("\t".join(['image_name','cell_number', 'mRNA_entry','peptide_entry','pep_amp','normalized_pep_amp','distance'])+"\n")

for PEP in uniquemPEPs:
    repeatedEntries=np.where((summaryAmpDist[:,3]==PEP))    
    if len(repeatedEntries[0])>1:
        distances=[]
        for entry in repeatedEntries[0]:
            distances.append(summary[entry,-1])
        entry_id=repeatedEntries[0][np.argmax(np.array(distances))]
        textDistance.write("\t".join(summaryAmpDist[entry_id,:]))
        textDistance.write("\n")
    else:
        entry_id=repeatedEntries[0][0]
        textDistance.write("\t".join(summaryAmpDist[entry_id,:]))
        textDistance.write("\n")
textDistance.close()

############################### for prism output --> translation efficiency ###############################
#the cell_translation.txt shows the # of translating mRNA & non translating mRNA
summaryAmpDist=np.loadtxt("tmp_distance.txt", skiprows=1, dtype=str)

unified_translating_cellnames=np.char.add(summaryAmpDist[:,0].astype(np.str_),summaryAmpDist[:,1].astype(np.str_))
unified_no_translating_cellnames=np.char.add(ch1[:,0].astype(np.str_),ch1[:,1].astype(np.str_))

translatingFilepath = str(os.getcwd()+"/"+outputFileName+"_cell_translatingmRNA.txt")
text=open(translatingFilepath,"w+")
text.write("\t".join(['image_name','cell_number','translating','untranslating','total mRNA','% of translating mRNA'])+"\n")

for name in np.unique(unified_translating_cellnames):
    translating_count=len(np.where(unified_translating_cellnames==name)[0])
    untranslating_count=len(np.where(unified_no_translating_cellnames==name)[0])
    total = translating_count+untranslating_count
    translation_percentage = translating_count/total
    print(name,translating_count,untranslating_count,total,translation_percentage)
    itms=name.split('.tif')
    text.write('\t'.join([itms[0]+'.tif',itms[1]])+"\t"+str(translating_count)+"\t"+str(untranslating_count)+"\t"+str(total)+"\t"+str(translation_percentage)+"\n")
text.close()

############################### for prism output --> ribosome occupency ###############################
summaryAmpDist=np.loadtxt("tmp_distance.txt", skiprows=1, dtype=str)

unified_cellnames=np.char.add(summaryAmpDist[:,0].astype(np.str_),summaryAmpDist[:,1].astype(np.str_))

countfilepath = str(os.getcwd()+"/"+outputFileName+"_count.txt")
text=open(countfilepath,"w+")
text.write('\t'.join(['image_name','cell_number','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','20+'])+"\n")

for name in np.unique(unified_cellnames):
    count=len(np.where((unified_cellnames==name) * (summaryAmpDist[:,5].astype(float)>=0.2))[0])
    count1=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>=0.2) & (summaryAmpDist[:,5].astype(float)<=1.5))[0])
    count2=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>1.5) & (summaryAmpDist[:,5].astype(float)<=2.5))[0])
    count3=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>2.5) & (summaryAmpDist[:,5].astype(float)<=3.5))[0])
    count4=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>3.5) & (summaryAmpDist[:,5].astype(float)<=4.5))[0])
    count5=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>4.5) & (summaryAmpDist[:,5].astype(float)<=5.5))[0])
    count6=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>5.5) & (summaryAmpDist[:,5].astype(float)<=6.5))[0])
    count7=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>6.5) & (summaryAmpDist[:,5].astype(float)<=7.5))[0])
    count8=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>7.5) & (summaryAmpDist[:,5].astype(float)<=8.5))[0])
    count9=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>8.5) & (summaryAmpDist[:,5].astype(float)<=9.5))[0])
    count10=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>9.5) & (summaryAmpDist[:,5].astype(float)<=10.5))[0])
    count11=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>10.5) & (summaryAmpDist[:,5].astype(float)<=11.5))[0])
    count12=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>11.5) & (summaryAmpDist[:,5].astype(float)<=12.5))[0])
    count13=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>12.5) & (summaryAmpDist[:,5].astype(float)<=13.5))[0])
    count14=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>13.5) & (summaryAmpDist[:,5].astype(float)<=14.5))[0])
    count15=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>14.5) & (summaryAmpDist[:,5].astype(float)<=15.5))[0])
    count16=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>15.5) & (summaryAmpDist[:,5].astype(float)<=16.5))[0])
    count17=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>16.5) & (summaryAmpDist[:,5].astype(float)<=17.5))[0])
    count18=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>17.5) & (summaryAmpDist[:,5].astype(float)<=18.5))[0])
    count19=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>18.5) & (summaryAmpDist[:,5].astype(float)<=19.5))[0])
    count20=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>19.5) & (summaryAmpDist[:,5].astype(float)<=20.5))[0])
    countMore=len(np.where((unified_cellnames==name) & (summaryAmpDist[:,5].astype(float)>20.5) & (summaryAmpDist[:,5].astype(float)<=40.5))[0])
    itms=name.split('.tif')
  
    text.write("\t".join([itms[0]+'.tif',itms[1]])+"\t"+str(count1)+"\t"+str(count2)+"\t"+str(count3)+"\t"+str(count4)+"\t"+str(count5)+"\t"+str(count6)+"\t"+str(count7)+"\t"+str(count8)
               +"\t"+str(count9)+"\t"+str(count10)+"\t"+str(count11)+"\t"+str(count12)+"\t"+str(count13)
               +"\t"+str(count14)+"\t"+str(count15)+"\t"+str(count16)+"\t"+str(count17)+"\t"+str(count18)+"\t"+str(count19)+"\t"+str(count20)+"\t"+str(countMore)+"\n")
    
text.close()


############################### clean up ###############################
os.remove('tmp_amplitude.txt')
os.remove('tmp_distance.txt')



print("done")








