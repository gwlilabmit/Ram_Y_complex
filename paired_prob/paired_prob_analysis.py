'''
Produces box plots of paired probabilities for each position in the windows of interest. Size of windows is standardized to be the size of the shortest window. 
'''

import glob
import re 
import numpy as np
import pylab
import matplotlib.pyplot as plt

#################
#opening files
#################

infiles = glob.glob("*.dat") #normalized dms sequences without -999's

#########################
#global variables/arrays
#########################

sequences = [] #2d array--sequences[i] gives you an array of pairing probabilities for that sequence i

shortenedsequences = [] #same information as sequences--this time with shortened sequences

probsbyposition = [] #probsbyposition[i] gives you pairing probabilities for all sites at position i

distancefromcleavage = 100 #number of nts from the cleavage site to the end of the sequence

##############################################
#read files and set them to be the same size
#######

for infile in infiles: 
	datfile = open(infile, 'r')

	pairingprobs = [] #array of pairing probabilities for that site

	for line in datfile: 
		line = line.split('\t') #each line is "position \t probability\n"
		prob = float(line[1][:-1]) #[:-1] gets rid of the \n at the end of the probability, convert to a float so we can plot 
		
		pairingprobs.append(prob)

	sequences.append(pairingprobs)

	datfile.close()


####################################################################
#standardize length of all sequences to the shortest window length
####################################################################

#note: I chop length from the beginning of the window, not the end! This is because all windows have 100nt after the cleavage site but
#have a variable length before 

#I add everything to shortenedsequences

shortestlength = len(min(sequences, key=len))

cleavage = shortestlength-distancefromcleavage-1 #-1 accounts for 0-indexing

for sequence in sequences:

	if len(sequence)>shortestlength: 
		lengthtocut = len(sequence)-shortestlength
		shortenedsequence = sequence[lengthtocut:]
		
	else: 
		shortenedsequence = sequence
	
	shortenedsequences.append(shortenedsequence)

##########################################################################
#find pairing probabilities for each position and add to probsbyposition
##########################################################################

#first populate probsbyposition with a bunch of empty arrays

for i in range(shortestlength): 
	probsbyposition.append([])

#then add the jth index of each sequence in shortenedsequences to probsbyposition[j]

for i in range(len(shortenedsequences)):
	sequence = shortenedsequences[i]
	for j in range(len(sequence)):
		probsbyposition[j].append(sequence[j])

#######################
########plotting
#######################

###box plot--sections lifted from http://blog.bharatbhole.com/creating-boxplots-with-matplotlib/ ###

plt.scatter(cleavage+1, 0, s=1, color="b", label="Cleavage site") #add a dot to the cleavage site bar so we can label it--for some reason the +1 is necessary else
																  #the dot is shifted over one to the left

#ax = fig.add_subplot(111)
bp = plt.boxplot(probsbyposition, 0, "", patch_artist=True) #the 0, '' options are to not show outliers--see https://matplotlib.org/3.1.1/gallery/statistics/boxplot_demo.html

for i in range(len(bp['boxes'])):
	box = bp['boxes'][i]
	box.set(color='goldenrod')
	if i == cleavage: #set the cleavage site to have a red bar
		box.set(color='blue')

for i in range(len(bp['whiskers'])):
	whisker = bp['whiskers'][i]
	whisker.set(color='goldenrod')
	if i == cleavage*2 or i==cleavage*2+1: #because whiskers has 2x as many entries as boxes--the+1 gets both whiskers 
		whisker.set(color='blue')

for i in range(len(bp['caps'])):
	cap = bp['caps'][i]
	cap.set(color='goldenrod')
	if i == cleavage*2 or i==cleavage*2+1: 
		cap.set(color='blue')

for median in bp['medians']:
	median.set(color='navy')


#pylab.legend(loc='upper right') #show a legend for the cleavage site and start codon 
plt.xticks(np.arange(0,shortestlength,10),np.arange(0, shortestlength,10),rotation="vertical")
plt.autoscale()	
plt.xlabel("Location along window")
plt.ylabel("p(unpaired) for all genes")
plt.title("Boxplot of pairing probabilities for the 21 sites")
plt.draw()

plt.show()

	



