'''
Gets rid of all files in the varna directories that shouldn't be of use for the paper. I only want to keep the following: 
centroid, mfe, mea, locarna if present, paired constrained, consensus without dms and rnaalifold--all with paired coloring only
'''

import os


directory = "dlta_varna"
name = directory.split("_")[0] #e.g. if directory="dusc_varna", the split will produce the array ["dusc", "varna"]--the name of interest is "dusc"

files_to_remove = ["_consensus_paired.varna", "_dms_constrained_paired_colormap.varna", "_downstream_paired.varna", "_mfe_paired_dms_consensus.varna", "_mfe_paired_dms_consensus_paired.varna", "_paired_constrained_dms_colormap.varna", "_paired_constrained_downstream.varna", "_subtilis_consensus_paired.varna", "_without_dms_paired.varna", "_without_rnaalifold_consensus_paired.varna"]

for i in range(len(files_to_remove)): 
	to_remove = files_to_remove[i]
	
	command = "rm "+directory+"/"+name+to_remove
	os.system(command)


#show contents of directory
os.system("ls "+directory)
