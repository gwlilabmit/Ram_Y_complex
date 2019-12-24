#calculates the distribution of distances from the cleavage site to the nearest downstream hairpin

import glob
import numpy as np
import matplotlib.pyplot as plt
import pylab

#This dictionary specifies the cleavage site location for each site of interest--note that this is 1-indexed. 
cleavage_sites = {"atpi":100,"cggr":100,"cwlo":86,"daca":12,"ddl":43,"glna":100,"mena":23,"metn":100,"metq":100,"ydih":100,\
"rpsl":79,"rpob":100,"pyrg":93,"ruva":47,"taga":12,"ybfg":100,"ffh":100,"yonr":100,"yqhl":99,"dlta":100,"dusc":100}

infiles = glob.glob("*_mfe.txt") #this is the naming convention I've used for my consensus files 

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
			
			#"unzipping" stems until we hit the first GC pair
			working_seq = downstreamseq[:] #only the portion after the cleavage site--same length as downstream
			working_struct = downstream[:]
		
			working_seqright = downstreamseq[:]
			working_structright = downstream[:]
			
			if firstleft<=len(working_seq): 
				startleft = working_seq[firstleft] #location of naive first left paren
						
				i = 0 #unzip left (--basically if the first ( is A/T I keep on shortening the sequence until the first ( is C/G 
				newstart = 0
				while startleft == "A" or startleft=="T" or startleft=="U": 
					working_seq = working_seq[1:] #I shorten the sequence and structure by one on the left--e.g. (()) becomes ())
					working_struct = working_struct[1:]
					newstart = working_struct.find("(")
					startleft = working_seq[newstart]
					i += 1
				unzippedleft = i+newstart #so the index of the first G/C left paren is 
										#the number of characters I chopped off the sequence+the location of the first left paren in that shortened sequence
				if i == 0: 
					unzippedleft = firstleft
			else: 
				firstleft = len(working_seq)+1
				unzippedleft = len(working_seq)+1
		
			#for debugging
			print name
			print "post_seq: "+post_seq
			print "post_cle: "+post_cleavage
			print "downstrseq: "+ downstreamseq
			print "downstream: "+ downstream
			print "starting: "+ str(starting)
			print "working_seq: "+ working_seq
			print "working_str: "+ working_struct
			print "i: "+ str(i)
			print "firstleft: "+ str(firstleft)
			print "unzippedleft: "+ str(unzippedleft)			
	
			if firstright<=len(working_seqright):
				startright = working_seqright[firstright] #location of naive first right paren 

				j = 0 #unzip right )--same idea as left
				newstart = 0
				while startright == "A" or startright=="T" or startright=="U": 
					working_seqright = working_seqright[1:]
					working_structright = working_structright[1:]
					newstart = working_structright.find(")")
					startright = working_seqright[newstart]
					j += 1
				unzippedright = j+newstart
				if j == 0: 
					unzippedright = firstright		
		
			else: 
				firstright = len(working_seqright)+1
				unzippedright = len(working_seqright)+1				

			print "firstright: "+str(firstright)
			print "unzippedright: "+str(unzippedright)+"\n"

			if firstleft<firstright: #If the structure is either not in a hairpin or on the ( side of a stem
				distances[name] = firstleft+starting

			elif firstright<=firstleft:
				#the site's in a bulge--set the distance to be the distance to the nearest downstream paired region
				distances[name] = firstright
			
			if unzippedleft<unzippedright: #If the structure is either not in a hairpin or on the ( side of a stem
				unzipped_distances[name] = unzippedleft+starting 

			elif unzippedright<=unzippedleft:
				 #the site's in a bulge--set the distance to be the distance to the nearest downstream paired region
				unzipped_distances[name] = unzippedright

			###comment this out if you want to only look at downstream sites--this sets the
			if post_cleavage[0] != ".": #if the site is paired	
				distances[name]=0 #it's in a stem
				if post_seq[0] == "G" or post_seq[0]=="C":
					unzipped_distances[name]=0 #set this to the distance from the site to the nearest downstream paired G/C--can be from the same stem as the site
				else: 
					for i in range(len(post_seq)):
						if post_seq[i]=="G" or post_seq[i]=="C":
							if post_cleavage[i]!=".": 
								unzipped_distances[name]=i
								break
	
	
		
####plotting						

print "Distances: "+ str(distances)

#####distances######
distvals = distances.values()
distkeys = distances.keys()
meandistance = np.mean(distvals) #mean distance

####plot a histogram of the values in the dictionary distances and add a dashed line at the mean distance
plt.figure(1)
plt.hist(distvals,color="skyblue", bins=np.arange(min(distvals),max(distvals)+2,1.0), edgecolor="white", align='mid') #plot a histogram of distances 
plt.axvline(x=meandistance, color="r", linestyle="--", label="Mean distance") #plot a dashed vertical line for mean distance
plt.xticks(np.arange(min(distvals), max(distvals)+2, 1.0),rotation='vertical') #I redid the x axes--ticks in increments of +1 and rotated to be vertical
pylab.legend(loc='upper right') #show a legend for the mean distance line 
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

####unzipped####
unzippedvals = unzipped_distances.values()
unzippedsites = unzipped_distances.keys()
meanzip = np.mean(unzippedvals) #mean unzipped distance
print "Unzipped:"+str(unzipped_distances)

###plot a histogram of the values in the dictionary unzipped_distances and add a dashed line at the mean distance
plt.figure(3)
plt.hist(unzippedvals,color="goldenrod", bins=np.arange(min(unzippedvals),max(unzippedvals)+2,1.0), edgecolor="white", align='mid') #plot a histogram of distances 
plt.axvline(x=meanzip, color="r", linestyle="--", label="Mean distance") #plot a dashed vertical line for mean distance
plt.xticks(np.arange(min(unzippedvals), max(unzippedvals)+2, 1.0),rotation='vertical') #I redid the x axes--ticks in increments of +1 and rotated to be vertical
pylab.legend(loc='upper right') #show a legend for the mean distance line 
plt.autoscale() 

plt.xlabel("Distance")
plt.ylabel("Frequency")
plt.title("Distribution of 'unzipped' distance from the cleavage site to the nearest G/C in a\n downstream stem. Values of 0 indicate that the site is paired.",fontsize=11)
plt.draw()

####plot a bar plot of unzipped
plt.figure(4)
plt.bar(np.arange(len(unzippedvals)), unzippedvals, color="forestgreen", edgecolor="white")
plt.autoscale()
plt.xlabel("Cleavage site")
plt.xticks(np.arange(len(unzippedvals)), unzipped_distances.keys(), rotation="vertical")
plt.ylabel("Unzipped distances")
plt.title("Unzipped distances from cleavage sites to the nearest downstream stem")
plt.draw()



plt.show()

