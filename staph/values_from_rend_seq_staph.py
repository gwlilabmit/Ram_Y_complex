import numpy as np
import pylab
import matplotlib.pyplot as plt
from scipy.stats import mstats

'''
This script plots data for the 5 S. aureus sites of interest. The 6th site corresponds to the sequence TACTTACTAAATTTTATTTAACCTAAAAATGAACCACCTGGATGTGTGGG and doesn't appear to be Y-dependent.
'''

####opening rend-seq files######

#open forward and reverse rend-seq files. Change the path as needed--my script resides in the same directory as these wig files but yours may not.
rend_forward_5 = open("staph_ylbf_5_f.wig", "r") #We can just look at 5' reads here, it doesn't really matter 
rend_reverse_5 = open("staph_ylbf_5_r.wig", "r")

rend_forward_3 = open("staph_ylbf_3_f.wig", "r") #We can just look at 5' reads here, it doesn't really matter 
rend_reverse_3 = open("staph_ylbf_3_r.wig", "r")


#####opening bedfile######

#open bedfile with indices pertaining to the window of interest as well as a text file for output. Again, change path as needed. 
bedfile=open("staph_sites.bed12", "r")

######initializing arrays######

#initialize arrays that you'll populate with data from the bedfile
starting = [] #starting window indices
ending = [] #ending window indices
direction = [] #strandedness of the window you're considering
name = [] #name of the region

######initializing dictionaries######

#For my use: this is where the cleavage sites are within the windows I've chosen
cleavage_sites = {"site1": 25, "site2": 25, "site3": 25, "site4": 25, "site5": 25}


######initializing normalization offset####
offset = 4 #so exclude 4nts on each side of the cleavage site when doing normalization 


#######get starting/ending indices of windows of interest from the bedfile#####

#Extract salient values from the bedfile and populate starting/ending/direction/name arrays
for line in bedfile:
	if '_' in line: #i.e. if the line isn't a header--all the genome names have a '_' in them. It's a dumb way of doing it but it works.
		splitline = line.split('\t') #bedfiles are tab-separated--I need to get rid of the tabs in between values of interest for each line
		for i in range(len(splitline)): 
			word = splitline[i] 
			if '\n' in word: #the last word in each line has an attached \n so I need to get rid of that
				splitline[i] = word[:-1]
		#the bedfile is structured as "genome window_start window_end window_name another_parameter strand_direction"
		#that's consistent throughout the file so I can just, say, add each second word in each line to my "starting" array
		starting.append(int(splitline[1])) 
		ending.append(int(splitline[2]))
		direction.append(splitline[5])
		name.append(splitline[3])
bedfile.close()


######grabbing dms values for each window, keeping strand directionality in mind#########

#Based on the strand directionality of each sequence I'm going to grab the corresponding dms values for the window
#I'm also going to write that to a text file along with the name of the window and the starting/ending indices of that window as specified by the bedfile
#FOR EACH WINDOW I parse through the wig file line by line and look for reads whose locations fall within the window
#So the wig file is structured as "location reads" so I need to read through the locations first to see if they're within the starting and ending values I'm currently considering
#Let's say that my starting location is 300 and I find 3 reads at location 302. In that case, this part will set readcounts[2]=3. 

#additionally, I'm exporting normalized dms-seq values into a .dat file [see: https://www.tbi.univie.ac.at/RNA/RNAfold.1.html for more info]
#I also export arcsinh-transformed normalized dms-seq values into a .arcsinh file--the transform just helps with visualization/prevents high-valued nts from saturating the colormap
#basically each line is [one-indexed] index \t readcount
#you can use RNAfold's SHAPE-seq support to constrain folding by DMS-seq values using the --shape=filename.dat option. 
#####3' reads!!!#######

j = 1 #to keep track of figures to plot
for i in range(len(starting)):

	window_start = starting[i]
	window_end = ending[i]
	window_direction = direction[i]
	rend_forward_3.seek(0) #"Seek" sets the pointer back to the top of the page
	rend_reverse_3.seek(0)
	
	if window_direction=="+":
		rend_file = rend_forward_3
	elif window_direction=="-": 
		rend_file = rend_reverse_3

	#Write the name of the current site 
	window_name = name[i]

	#Create an array for the readcounts of each index of interest of size [window_end-window_start]
	#I'm initializing it with a bunch of 0's because I'm not guaranteed a read at every location within my window
	#So, in that instance, there would be 0 reads. Best to initialize with that in mind.

	readcounts = [0]*(window_end-window_start+1)

	#Parse through the lines in the dms file and extract the readcounts for all locations in your window of interest. 
	#We'll add them to the readcounts array.

	for line in rend_file:
		#The lines are structured as 'location\treads\n' so I need to get rid of the \t and \n to access index and reads
		index = int(line.split('\t')[0]) #The wig file is structured by columns--first column is the location/index so we extract the first word of the line here
		reads = line.split('\t')[1] #Second column is the readcount, so we're extracting the second word of the line
		reads = float(reads[:-1]) #Get rid of the \n on the end of the reads--it's treated as a single character
		if window_start <= index <= window_end: #For locations within the specified window of interest
			readcounts[index-window_start]=reads #Add the readcounts and indices to their respective arrays

	#Write the array readcounts to your output file
	#If the window is on the minus strand you need to reverse the array

	if window_direction == "-":
		readcounts = readcounts[::-1] #if we're looking at the reverse strand, reverse the readcounts otherwise the window will be backwards
	else: 
		pass
	readcounts = readcounts[:-1] #I'm adding an extra index when defining readcounts--for some reason I get an error when I don't do that so I'm just removing it here

	#Here I'm normalizing readcounts by downstream mRNA expression. 
	#Specifically, I normalize the readcounts by the average reads downstream of the cleavage site
	#Recall: I specified start codon locations relative to my window start earlier on. 
	#I divide all reads prior to the start codon by the average read count after/including the cleavage site. 

	normalizedreadcounts = readcounts[:]
	
	#apply 98% winsorization 
	#normalizedreadcounts = mstats.winsorize(normalizedreadcounts, limits=[None, 0.02])

	cleavage = cleavage_sites[window_name] #-1 #subtract 1 to account for zero indexing--the 100th base in the window is the 99th index 
	
	downstream = normalizedreadcounts[cleavage+offset:-offset-1] #exclude the number of nts specified by offset from each side of the window to normalize	
	upstream = normalizedreadcounts[offset:cleavage-offset-1] #the -1 here accounts for python's indexing here and in downstream

	averagedownstream = np.mean(downstream) #Get average reads downstream of the cleavage site
	averageupstream = np.mean(upstream)
	
	
	#########plotting unnormalized rend-seq values#############
		
	plt.figure(j)
	
	color = ["#336699" for i in normalizedreadcounts]
	edgecolor = ["#336699" for i in normalizedreadcounts]

	plt.scatter(len(normalizedreadcounts)-1, 0, s=1, color=color, label="3' reads") #add a dot to the cleavage site bar so we can label it

	#color[cleavage]="b" #color the cleavage site bar red
	#edgecolor[cleavage]="b" #color the cleavage site bar's edge red
	
	plt.bar(np.arange(len(normalizedreadcounts)),normalizedreadcounts, color=color, edgecolor=edgecolor)
	#add lines corresponding to the average downstream and upstream read counts 
	#from https://stackoverflow.com/questions/36470343/how-to-draw-a-line-with-matplotlib
	
	plt.xticks(np.arange(0,len(normalizedreadcounts),10),rotation="vertical")
	plt.autoscale()	
	plt.xlabel("Location along window")
	plt.ylabel("Read counts")
	plt.ylim(top=4.5)
	plt.title(window_name+" Rend-seq reads")
	j += 1
	plt.draw()


####5' reads!!!!!#######

j = 1 #to keep track of figures to plot
for i in range(len(starting)):

	window_start = starting[i]
	window_end = ending[i]
	window_direction = direction[i]
	rend_forward_5.seek(0) #"Seek" sets the pointer back to the top of the page
	rend_reverse_5.seek(0)
	
	if window_direction=="+":
		rend_file = rend_forward_5
	elif window_direction=="-": 
		rend_file = rend_reverse_5

	#Write the name of the current site 
	window_name = name[i]

	#Create an array for the readcounts of each index of interest of size [window_end-window_start]
	#I'm initializing it with a bunch of 0's because I'm not guaranteed a read at every location within my window
	#So, in that instance, there would be 0 reads. Best to initialize with that in mind.

	readcounts = [0]*(window_end-window_start+1)

	#Parse through the lines in the dms file and extract the readcounts for all locations in your window of interest. 
	#We'll add them to the readcounts array.

	for line in rend_file:
		#The lines are structured as 'location\treads\n' so I need to get rid of the \t and \n to access index and reads
		index = int(line.split('\t')[0]) #The wig file is structured by columns--first column is the location/index so we extract the first word of the line here
		reads = line.split('\t')[1] #Second column is the readcount, so we're extracting the second word of the line
		reads = float(reads[:-1]) #Get rid of the \n on the end of the reads--it's treated as a single character
		if window_start <= index <= window_end: #For locations within the specified window of interest
			readcounts[index-window_start]=reads #Add the readcounts and indices to their respective arrays

	#Write the array readcounts to your output file
	#If the window is on the minus strand you need to reverse the array

	#if window_direction == "-":
	#	readcounts = readcounts[::-1] #if we're looking at the reverse strand, reverse the readcounts otherwise the window will be backwards
	#else: 
	#	pass
	readcounts = readcounts[:-1] #I'm adding an extra index when defining readcounts--for some reason I get an error when I don't do that so I'm just removing it here

	#Here I'm normalizing readcounts by downstream mRNA expression. 
	#Specifically, I normalize the readcounts by the average reads downstream of the cleavage site
	#Recall: I specified start codon locations relative to my window start earlier on. 
	#I divide all reads prior to the start codon by the average read count after/including the cleavage site. 

	normalizedreadcounts = readcounts[:]
	
	#apply 98% winsorization 
	#normalizedreadcounts = mstats.winsorize(normalizedreadcounts, limits=[None, 0.02])

	cleavage = cleavage_sites[window_name] #-1 #subtract 1 to account for zero indexing--the 100th base in the window is the 99th index 
	
	downstream = normalizedreadcounts[cleavage+offset:-offset-1] #exclude the number of nts specified by offset from each side of the window to normalize	
	upstream = normalizedreadcounts[offset:cleavage-offset-1] #the -1 here accounts for python's indexing here and in downstream

	averagedownstream = np.mean(downstream) #Get average reads downstream of the cleavage site
	averageupstream = np.mean(upstream)
	
			
	#########plotting unnormalized rend-seq values#############
		
	plt.figure(j)
	color = ["#cc6633" for i in normalizedreadcounts]
	edgecolor = ["#cc6633" for i in normalizedreadcounts]

	plt.scatter(cleavage, 0, s=1, color=color, label="5' reads") #add a dot to the cleavage site bar so we can label it	

	
	plt.bar(np.arange(len(normalizedreadcounts)),normalizedreadcounts, color=color, edgecolor=edgecolor)
		
	plt.xticks(np.arange(0,len(normalizedreadcounts),10),rotation="vertical")
	plt.autoscale()	
	plt.xlabel("Location along window")
	plt.ylabel("Read counts")
	plt.ylim(top=4.5)
	plt.title(window_name+" Rend-seq reads")
	j += 1
	plt.draw()
	
	
	#first we have the [one-indexed, hence the +1] index column and then we have the normalized DMS reads corresponding to that index. 
	#the columns need to be tab-separated. 




plt.show()

##close rend-seq files
rend_reverse_5.close()
rend_forward_5.close()
rend_reverse_3.close()
rend_forward_3.close()
