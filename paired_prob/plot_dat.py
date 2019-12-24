#Extracts second-column values from .dat files and prints them out, comma-separated, so they can be used as a colormap in VARNA
#It'll do this for all .dat files you have in your directory. If you don't want this feature just comment out everything with read_files in it
#and unindent as needed. 
#I also plot out the values for just A/C reads.  

#I'm also printing out Yeo-Johnson or arcsinh-transformed reads--this is useful if there's a wide range of values [0 included] and you don't want a high-read nt to affect your colormap visualization dramatically. 

#I also plot reads for a given sequence transformed both ways for the sake of comparison. 

#If you're curious about Yeo-Johnson--its main benefit is that it can transform exponentially distributed data into normally-distributed data, with the additional perk of being able to deal with negative/zero values [unlike a Boxcox transform]
#https://machinelearningmastery.com/how-to-transform-data-to-fit-the-normal-distribution/ does a nice job explaining what the Yeo-Johnson is/what it does. 

import re 
import numpy as np
import glob
import matplotlib.pyplot as plt
from sklearn.preprocessing import PowerTransformer

read_files = glob.glob("*.dat")
sequences = open("21_cleavage_windows_final.txt", "r")
all_sequences = {}

for line in sequences: 
	if ">" in line: 
		seqname = line[1:-1]
	else: 
		all_sequences[seqname]=line[:-1]
sequences.close()

j = 1
for datfile in read_files: 
	infile = open(datfile, "r")
	#comment out this regex stuff if your .dat file isn't named "gene.dat"--with my naming convention this extracts the gene name for me
	regex = r"^[a-zA-Z]+"
	matches = re.findall(regex, datfile) #say the filename is atpi.dat. This extracts "atpi"
	name = matches[0]

	values = [] #array of all second-column values, i.e. the values of interest for the colormap

	for line in infile: 
		reads = line.split("\t")[1] #Each line is tab-separated. We want the value in the second column. 
		reads = reads[:-1] #There's a \n at the end of the "reads" value, which counts as a single character. 
		values.append(reads)
	
	values = np.array(values[:]).astype(float)
	ac_values = []
	sequence = all_sequences[name]

	for i in range(len(sequence)): 
		if sequence[i]=="A" or sequence[i]=="C": 
			ac_values.append(values[i]) #only add dms reads corresponding to A/C nts to ac_values

	#########plotting reads for all nts###########
	'''
	plt.figure(j)
	plt.hist(values, color="lemonchiffon", bins=np.arange(0, max(values)+2,1.0), edgecolor="darkgoldenrod",align="mid")
	plt.xticks(np.arange(min(values), max(values)+2, 1.0),rotation="vertical")
	plt.autoscale()	
	plt.xlabel("Read count")
	plt.ylabel("Frequency")
	plt.title(name+" DMS untransformed reads")
	j += 1
	plt.draw()
	'''

	values_to_transform = values[:] #The dms values were strings earlier--we need to convert to floats to manipulate

	#log transform
	for i in range(len(values_to_transform)): 
		value = values_to_transform[i]
		if value == 0: 
			values_to_transform[i] = 1e-7 #add as a pseudocount 
	transformed_vals = np.log(values_to_transform)	

	#This gets a bit convoluted. Basically I find the second-smallest value in transformedvals [so, the smallest nonzero value], add that value to all values in
	#transformedvals and then set any negative values to 0	

	findmin = transformed_vals[:]
	minval = min(findmin)
	findmin = findmin[findmin!=minval] #from https://stackoverflow.com/questions/53541156/how-to-remove-all-occurrences-of-an-element-from-numpy-array
	smallestnonzero = min(findmin)
	offset = 1 #set the second-lowest values to 1
	transformed_vals = [i+np.abs(smallestnonzero)+offset for i in transformed_vals]
	for i in range(len(transformed_vals)): 
		value = transformed_vals[i]
		if value < offset: #if it's <offset it's smaller than smallestnonzero
			transformed_vals[i] = 0


	#arcsinh transform 
	#transformed_vals = np.arcsinh(values_to_transform)
	
	#implementing Yeo-Johnson as per https://stackoverflow.com/questions/53624804/how-to-normalize-a-non-normal-distribution
	#values_to_transform = values_to_transform.reshape(-1,1) #convert to a 2d array
	#pt = PowerTransformer(method='yeo-johnson')
	#calculate the right parameters to fit the data [this is lambda from the transform]
	#pt.fit(values_to_transform)
	#transformed_vals = pt.transform(values_to_transform) 

	plt.figure(j)
	plt.hist(transformed_vals, color="tomato", bins=np.arange(0, max(transformed_vals)+2,1.0), edgecolor="white",align="mid")
	plt.xticks(np.arange(min(transformed_vals), max(transformed_vals)+2, 1.0),rotation="vertical")
	plt.autoscale()	
	plt.xlabel("Read count")
	plt.ylabel("Frequency")
	plt.title(name+" DMS log-transformed reads")
	j += 1
	plt.draw()
	
	#######plotting reads for a/c only########
	'''
	plt.figure(j)
	plt.hist(ac_values, color="goldenrod", bins=np.arange(0, max(ac_values)+2,1.0), edgecolor="white",align="mid")
	plt.xticks(np.arange(min(ac_values), max(ac_values)+2, 1.0),rotation="vertical")
	plt.autoscale()	
	plt.xlabel("Read count")
	plt.ylabel("Frequency")
	plt.title(name+" DMS untransformed A/C reads")
	j += 1
	plt.draw()
	'''
	ac_values_to_transform = ac_values[:] #The dms values were strings earlier--we need to convert to floats to manipulate

	#log transform
	for i in range(len(ac_values_to_transform)): 
		value = ac_values_to_transform[i]
		if value == 0: 
			ac_values_to_transform[i] = 1e-7
	ac_transformed_vals = np.log(ac_values_to_transform)
	
	#This gets a bit convoluted. Basically I find the second-smallest value in transformedvals [so, the smallest nonzero value], add that value to all values in
	#transformedvals and then set any negative values to 0	

	findminac = ac_transformed_vals[:]
	minac = min(findminac)
	findminac = findminac[findminac!=minac] #findminac with all instances of the smallest value removed
	smallestnonzeroac = min(findminac)
	offset = 1 #the difference you want between the smallest [0] value and the second-smallest value 
	ac_transformed_vals = [i+np.abs(smallestnonzeroac)+offset for i in ac_transformed_vals]
	
	for i in range(len(ac_transformed_vals)): 
		value = ac_transformed_vals[i]
		if value < offset:  
			ac_transformed_vals[i] = 0



	#arcsinh transform 
	#ac_transformed_vals = np.arcsinh(ac_values_to_transform)
	
	'''
	#implementing Yeo-Johnson as per https://stackoverflow.com/questions/53624804/how-to-normalize-a-non-normal-distribution
	ac_values_to_transform = np.array(ac_values_to_transform).astype(float).reshape(-1,1) #convert to a 2d array
	pt = PowerTransformer(method='yeo-johnson')
	#calculate the right parameters to fit the data [this is lambda from the transform]
	pt.fit(ac_values_to_transform)
	ac_transformed_vals = pt.transform(ac_values_to_transform) 
	'''
	plt.figure(j)
	plt.hist(ac_transformed_vals, color="skyblue", bins=np.arange(0, max(ac_transformed_vals)+2,1.0), edgecolor="white",align="mid")
	plt.xticks(np.arange(min(ac_transformed_vals), max(ac_transformed_vals)+2, 1.0),rotation="vertical")
	plt.autoscale()	
	plt.xlabel("Read count")
	plt.ylabel("Frequency")
	plt.title(name+" DMS log-transformed A/C reads")
	j += 1
	plt.draw()
	

	#print name+" reads:\n" + ",".join(values.astype(str))+"\n" #i.e. print "atpI reads: \n" followed by the reads 
	#print "Arcsinh-transformed "+name+" reads:\n" + ",".join(transformed_vals.astype(str))+"\n" #i.e. print "arcsinh-transformed atpI reads: \n" followed by the transformed reads 
	
	
	infile.close()
plt.show()
