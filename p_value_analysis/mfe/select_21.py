'''
Randomly finds num_windows windows and sees how many of them have paired cleavage sites. Repeats num_repeats times and tells you how many times there were >= cleavage_paired paired sites. 
'''

import glob
import random

structures = glob.glob("*_random.txt") #get all dms-constrained folded sites 

num_windows = 21 #number of windows you want to get, same as what you have for cleavage sites
cleavage_paired = 13 #number of cleavage sites you actually witness being paired 

num_repeats = 100000 #number of times you want to repeat this analysis 

greater_than_cleavage_paired = 0 #number of times the number of paired sites in a sample is greater than cleavage_paired

for i in range(num_repeats): 

	random_sample = random.sample(structures, num_windows) #randomly sample num_windows files from structures

	num_paired = 0

	for sample in random_sample: 
		infile = open(sample, "r")
		
		for line in infile: 
			if "." in line: #if that line has the structure
				structure = line.split(" ")[0]
				cleavage = len(structure)/2 -1 #the cleavage site is in the middle of the window, -1 accounts for 1-indexing
				if structure[cleavage] != ".": 
					num_paired += 1

		infile.close()

	if num_paired >= cleavage_paired: 
		greater_than_cleavage_paired += 1


print "Out of " + str(num_repeats) + " iterations, there were "+ str(greater_than_cleavage_paired) + " instances in which there were an equal or greater number of 'cleavage sites' paired than what was observed for the actual cleavage sites. This amounts to a fraction of " + str(float(greater_than_cleavage_paired)/float(num_repeats)) + " of instances in which one would expect to see "+ str(cleavage_paired) + "/" + str(num_windows) + " sites paired by chance."			
			
			
