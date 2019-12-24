#Extracts second-column values from .dat files and prints them out, comma-separated, so they can be used as a colormap in VARNA
#It'll do this for all .dat files you have in your directory. If you don't want this feature just comment out everything with read_files in it
#and unindent as needed. 

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

i = 1
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

	plt.figure(i)
	plt.hist(values, color="lemonchiffon", bins=np.arange(0, max(values)+2,1.0), edgecolor="darkgoldenrod",align="mid")
	plt.xticks(np.arange(min(values), max(values)+2, 1.0),rotation="vertical")
	plt.autoscale()	
	plt.xlabel("Read count")
	plt.ylabel("Frequency")
	plt.title(name+" DMS untransformed reads")
	i += 1
	plt.draw()

	values_to_transform = values[:] #The dms values were strings earlier--we need to convert to floats to manipulate

	#arcsinh transform 
	transformed_vals = np.arcsinh(values_to_transform)
	
	#implementing Yeo-Johnson as per https://stackoverflow.com/questions/53624804/how-to-normalize-a-non-normal-distribution
	#values_to_transform = values_to_transform.reshape(-1,1) #convert to a 2d array
	#pt = PowerTransformer(method='yeo-johnson')
	#calculate the right parameters to fit the data [this is lambda from the transform]
	#pt.fit(values_to_transform)
	#transformed_vals = pt.transform(values_to_transform) 

	plt.figure(i)
	plt.hist(transformed_vals, color="lavender", bins=np.arange(0, max(values)+2,1.0), edgecolor="indigo",align="mid")
	plt.xticks(np.arange(min(transformed_vals), max(transformed_vals)+2, 1.0),rotation="vertical")
	plt.autoscale()	
	plt.xlabel("Read count")
	plt.ylabel("Frequency")
	plt.title(name+" DMS arcsinh-transformed reads")
	i += 1
	plt.draw()


	print name+" reads:\n" + ",".join(values.astype(str))+"\n" #i.e. print "atpI reads: \n" followed by the reads 
	print "Arcsinh-transformed "+name+" reads:\n" + ",".join(transformed_vals.astype(str))+"\n" #i.e. print "arcsinh-transformed atpI reads: \n" followed by the transformed reads 
	
	infile.close()
plt.show()
