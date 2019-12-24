import numpy as np
import pylab
import matplotlib.pyplot as plt
from scipy.stats import mstats

'''
This runs a metagene analysis on Lydia's bedfile of >400 "simple" [1 isoform] operons. 
I take a window of user-defined size, find the mean read count at each position and plot it. 
This can take a window around both the 5' and 3'. Currently it's set to a window at the 5'. 
'''
import time
start_time = time.time()

####opening dms files######

#open forward and reverse dms files. Change the path as needed--my script resides in the same directory as these wig files but yours may not.
dms_forward = open("rend_seq_5_f.wig", "r") #We can just look at 5' reads here, it doesn't really matter 
dms_reverse = open("rend_seq_5_r.wig", "r")


#####opening bedfile######

#open bedfile with indices pertaining to the window of interest as well as a text file for output. Again, change path as needed. 
bedfile=open("monocistrons.bed", "r")

######initializing arrays######

#initialize arrays that you'll populate with data from the bedfile
starting = [] #starting window indices
ending = [] #ending window indices
direction = [] #strandedness of the window you're considering
name = [] #name of the region

######initializing window size####
window_size = 100 #take 50nt on either side of the start of the window [the 5']

######initializing 2d array of values at each index######
indices = []
normalizedindices = []

for i in range(window_size): 
	indices.append([])
	normalizedindices.append([])

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
		cds = int(splitline[6]) #cds in the window
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

j = 1 #to keep track of figures to plot
for i in range(len(starting)):

	fiveprime = starting[i]
	window_start = starting[i] #-window_size
	window_end = ending[i] #+window_size+1
	window_direction = direction[i]
	dms_forward.seek(0) #"Seek" sets the pointer back to the top of the page
	dms_reverse.seek(0)
	
	if window_direction=="+":
		dms_file = dms_forward
	elif window_direction=="-": 
		dms_file = dms_reverse

	#Write the name of the current site 
	window_name = name[i]
	
	#Create an array for the readcounts of each index of interest of size [window_end-window_start]
	#I'm initializing it with a bunch of 0's because I'm not guaranteed a read at every location within my window
	#So, in that instance, there would be 0 reads. Best to initialize with that in mind.

	readcounts = [0]*(window_end-window_start+1)

	#Parse through the lines in the dms file and extract the readcounts for all locations in your window of interest. 
	#We'll add them to the readcounts array.

	for line in dms_file:
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
	
	
	#So at this point we have a big transcript from the bedfile but we only want something of size 2*window_size 
	#Now I shave down readcounts to the first 2*window_size nts, the 1 is added because of python's string slicing

	readcounts = readcounts[:window_size] #get a window_size-sized window around the 5'
	#readcounts = readcounts[len(readcounts)-window_size:] #get a window_size-sized window around the 3'
	
	fiveprime = int(np.ceil(float(len(readcounts))/2))	

	#Here I'm normalizing readcounts by downstream mRNA expression. 
	#Specifically, I normalize the readcounts by the average reads downstream of the cleavage site
	#Recall: I specified start codon locations relative to my window start earlier on. 
	#I divide all reads prior to the start codon by the average read count after/including the cleavage site. 

	normalizedreadcounts = readcounts[:]
	averagereads = np.mean(normalizedreadcounts)

	######normalizing dms reads upstream/downstream of the cleavage site and writing to dat and arcsinh files#####

	#divide all values upstream of the cleavage site by averageupstream and all values downstream of cleavage by averagedownstream
	#the regions upstream and downstream of the site will have different expression levels so we should expect different average read counts
	#they need to be treated differently because they're part of two different transcripts whose expression levels will differ
	#normalizing by average expression helps us better determine what bases have a higher/lower signal than average 
	#this helps us make informed decisions about which regions are likely paired/unpaired

	for i in range(len(normalizedreadcounts)):
		normalizedreadcounts[i] = normalizedreadcounts[i]/averagereads

	#####adding values from each index in readcounts to the appropriate array#####

	for i in range(len(readcounts)): 
		read = readcounts[i]
		indices[i].append(read) #so the first array in indices has all the index 1 values

		normalizedread = normalizedreadcounts[i]
		normalizedindices[i].append(normalizedread)



##########finding mean/median values for each index########

medians = [0]*len(indices)
means = [0]*len(indices)

normalizedmedians = [0]*len(normalizedindices)
normalizedmeans = [0]*len(normalizedindices)


for i in range(len(indices)): 
	indexvals = indices[i]
	normalizedindexvals = normalizedindices[i]

	medians[i] = np.median(indexvals)
	means[i] = np.mean(indexvals)
	
	normalizedmedians[i] = np.median(normalizedindexvals)
	normalizedmeans[i] = np.mean(normalizedindexvals)

'''
#########plotting median values#############

fiveprime = 0 #the 5' is the first index 
#threeprime = int(np.ceil(float(len(medians))))

plt.figure(j)
plt.scatter(fiveprime, 0, s=1, color="r", label="5'") #add a dot to the cleavage site bar so we can label it

color = ["skyblue" for i in medians]
edgecolor = ["skyblue" for i in medians]

color[fiveprime]="r" #color the cleavage site bar red
edgecolor[fiveprime]="r" #color the cleavage site bar's edge red

plt.bar(np.arange(len(medians)),medians, color=color, edgecolor=edgecolor)
#add lines corresponding to the average downstream and upstream read counts 
#from https://stackoverflow.com/questions/36470343/how-to-draw-a-line-with-matplotlib

pylab.legend(loc='upper right') #show a legend for the cleavage site and start codon 
plt.xticks(np.arange(0,len(medians),10),rotation="vertical")
plt.autoscale()	
plt.xlabel("Location along window")
plt.ylabel("Median read counts for "+str(len(starting)) + " genes")
plt.title("Metagene analysis of median dms-seq reads near 5' end")
j += 1
plt.draw()

#########plotting mean values#############

fiveprime = 0 #the 5' is the first index 
#threeprime = int(np.ceil(float(len(medians))))

plt.figure(j)
plt.scatter(fiveprime, 0, s=1, color="r", label="5'") #add a dot to the cleavage site bar so we can label it

color = ["goldenrod" for i in means]
edgecolor = ["goldenrod" for i in means]

color[fiveprime]="r" #color the cleavage site bar red
edgecolor[fiveprime]="r" #color the cleavage site bar's edge red

plt.bar(np.arange(len(means)),means, color=color, edgecolor=edgecolor)
#add lines corresponding to the average downstream and upstream read counts 
#from https://stackoverflow.com/questions/36470343/how-to-draw-a-line-with-matplotlib

pylab.legend(loc='upper right') #show a legend for the cleavage site and start codon 
plt.xticks(np.arange(0,len(means),10),rotation="vertical")
plt.autoscale()	
plt.xlabel("Location along window")
plt.ylabel("Mean read counts for "+str(len(starting)) + " genes")
plt.title("Metagene analysis of mean dms-seq reads near 5' end")
j += 1
plt.draw()


#########plotting normalized median values#############

fiveprime = 0 #the 5' is the first index 
#threeprime = int(np.ceil(float(len(medians))))

plt.figure(j)
plt.scatter(fiveprime, 0, s=1, color="r", label="5'") #add a dot to the cleavage site bar so we can label it

color = ["thistle" for i in normalizedmedians]
edgecolor = ["thistle" for i in normalizedmedians]

color[fiveprime]="r" #color the cleavage site bar red
edgecolor[fiveprime]="r" #color the cleavage site bar's edge red

plt.bar(np.arange(len(normalizedmedians)),normalizedmedians, color=color, edgecolor=edgecolor)
#add lines corresponding to the average downstream and upstream read counts 
#from https://stackoverflow.com/questions/36470343/how-to-draw-a-line-with-matplotlib

pylab.legend(loc='upper right') #show a legend for the cleavage site and start codon 
plt.xticks(np.arange(0,len(normalizedmedians),10),rotation="vertical")
plt.autoscale()	
plt.xlabel("Location along window")
plt.ylabel("Median normalized read counts for "+str(len(starting)) + " genes")
plt.title("Metagene analysis of normalized median dms-seq reads near 5' end")
j += 1
plt.draw()

'''

#########plotting normalized mean values#############
fiveprime = 0 #the 5' is the first index 
#threeprime = int(np.ceil(float(len(medians))))

plt.figure(j)
plt.scatter(fiveprime, 0, s=1, color="r", label="5'") #add a dot to the cleavage site bar so we can label it

color = ["grey" for i in normalizedmeans]
edgecolor = ["grey" for i in normalizedmeans]

color[fiveprime]="r" #color the cleavage site bar red
edgecolor[fiveprime]="r" #color the cleavage site bar's edge red

plt.bar(np.arange(len(normalizedmeans)),normalizedmeans, color=color, edgecolor=edgecolor)
#add lines corresponding to the average downstream and upstream read counts 
#from https://stackoverflow.com/questions/36470343/how-to-draw-a-line-with-matplotlib

pylab.legend(loc='upper right') #show a legend for the cleavage site and start codon 
plt.xticks(np.arange(0,len(normalizedmeans),10),rotation="vertical")
plt.autoscale()	
plt.xlabel("Location along window")
plt.ylabel("Mean normalized read counts for "+str(len(starting)) + " genes")
plt.title("Metagene analysis of mean normalized dms-seq reads near 5' end")
j += 1
plt.draw()



#first we have the [one-indexed, hence the +1] index column and then we have the normalized DMS reads corresponding to that index. 
#the columns need to be tab-separated. 

plt.show()

print "Runtime: "+str((time.time()-start_time)/60)+ " minutes"

##close dms files
'''output.close()'''
dms_reverse.close()
dms_forward.close()
