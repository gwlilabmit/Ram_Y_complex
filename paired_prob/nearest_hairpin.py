#calculates the distribution of distances from the cleavage site to the nearest downstream hairpin

import glob
import numpy as np
import matplotlib.pyplot as plt
import pylab

#This dictionary specifies the cleavage site location for each site of interest--note that this is 1-indexed. 
cleavage_sites = {"atpi":100,"cggr":100,"cwlo":86,"daca":12,"ddl":43,"glna":100,"mena":23,"metn":100,"metq":100,"ydih":100,\
"rpsl":79,"rpob":100,"pyrg":93,"ruva":47,"taga":12,"ybfg":100,"ffh":100,"yonr":100,"yqhl":99,"dlta":100,"dusc":100}

infiles = glob.glob("*_paired.txt") #this is the naming convention I've used for my consensus files 

distances = {} #dictionary with distances from the cleavage site to the nearest downstream hairpin

unzipped_distances = {}

for structure in infiles: 
	infile = open(structure, "r")

	for line in infile: 
		if ">" in line: #i.e., if it specifies the name of the file as per fasta conventions
			name = line[1:-1].lower() #>atpI\n becomes "atpi"
			site = cleavage_sites[name]-1 #correct for one-indexing of the site 
		elif "G" in line and ">" not in line: 
			seq = line[:-1] #sequence corresponding to consensus structure--take off the \n at the end 

		elif "." in line: #i.e., it contains a structure
			line = line.split(" ")[0] #only consider the structure, not the delta G [they're separated by a space]
			post_cleavage = line[site:] #only consider the portion of the structure including/after the cleavage site
			post_seq = seq[site:] #slice the sequence in the same manner

			starting = post_cleavage.find(".") #start looking for the first downstream hairpin ater the first unpaired nt downstream of the cleavage site
			downstream = post_cleavage[starting:] #only consider regions downstream of the first unpaired nt	
			downstreamseq = post_seq[starting:]
			firstleft = downstream.find("(") #all portions of post-cleavage to the left of the first (, i.e. perceivably
			firstright = post_cleavage.find(")") #I'm using this to see if a ) appears before a (, which would mean that the site is within a hairpin
			
			
			if firstleft<firstright: #If the structure is either not in a hairpin or on the ( side of a stem
				distances[name] = firstleft+starting

			elif firstright<=firstleft:
				#the site's in a bulge--set the distance to be the distance to the nearest downstream paired region
				distances[name] = firstright
		
			###comment this out if you want to only look at downstream sites--this sets the
			if post_cleavage[0] != ".": #if the site is paired	
				distances[name]=0 #it's in a stem
				
	
		
####plotting						

#####distances######
distvals = distances.values()
distkeys = distances.keys()
meandistance = np.mean(distvals) #mean distance

####plot a histogram of the values in the dictionary distances and add a dashed line at the mean distance
plt.figure(1)
plt.hist(distvals,color="skyblue", bins=np.arange(min(distvals),max(distvals)+2,1.0), edgecolor="white", align='mid') #plot a histogram of distances 
plt.xticks(np.arange(min(distvals), max(distvals)+2, 1.0),rotation='vertical') #I redid the x axes--ticks in increments of +1 and rotated to be vertical
plt.autoscale() 

plt.xlabel("Distance")
plt.ylabel("Frequency")
plt.title("Distribution of distance from the cleavage site to the nearest\n downstream stem. Values of 0 indicate that the site is paired.",fontsize=11)
plt.draw()

###plot a bar plot of distances--one bar for each cleavage site
plt.figure(2)
plt.bar(np.arange(len(distvals)), distvals, color="salmon", edgecolor="white")
plt.autoscale()
plt.xlabel("Cleavage site")
plt.xticks(np.arange(len(distvals)), distkeys, rotation="vertical")
plt.ylabel("Distances")
plt.title("Distances from cleavage sites to the nearest downstream stem")
plt.draw()

plt.show()
