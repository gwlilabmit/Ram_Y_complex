'''
Plots a histogram showing the distribution of staph rny sites on the + and - strand, from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4608709/
'''

import numpy as np
import matplotlib.pyplot as plt

infile = open("staph.txt", "r") #contains the location of all sites, grouped by + or - strand
#basically, everything above the "-" is on the + strand

num_bins = 10000 #number of bins for histogram

seen_minus = 0

plus = [] #all sites on + strand
minus = [] #all sites on - strand

for line in infile: 
	if "-" in line: 
		seen_minus = 1
	else: 
		if seen_minus == 0: #if the line contains a + strand cleavage site
			plus.append(int(line[:-1])) #the [:-1] gets rid of the \n at the end

		else: #if the line contains a - strand cleavage site
			minus.append(int(line[:-1]))

#draw a histogram for the + sites
plt.figure(0)
plt.hist(plus, color='skyblue', bins = num_bins, edgecolor='skyblue')
plt.xlabel("Location")
plt.ylabel("Frequency")
plt.title("Distribution of + strand RNase Y cleavage sites in Staph")
plt.draw()
			
plt.figure(1)
plt.hist(minus, color='goldenrod', bins = num_bins, edgecolor='goldenrod')
plt.xlabel("Location")
plt.ylabel("Frequency")
plt.title("Distribution of - strand RNase Y cleavage sites in Staph")
plt.draw()

plt.show()

infile.close()

