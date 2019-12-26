'''
Plots the distribution of a/u/c/g for all sites with good alignments. 
The input file has 100nt windows of the msa for all sites [50nt on either side of the cleavage site]
From there I calculate the fraction of a/u/c/g at each position for all MSAs
This lets me see if there's any global conservation of a given nucleotide at some position relative to all cleavage sites. 
'''

import numpy as np
import matplotlib.pyplot as plt


infile = open("ydih.txt", "r")



#count the number of A/U/C/G at each position

for line in infile: 

	#initialize arrays to count the number of A/U/G/C at each position 

	if "adenine" not in globals():  #adenine[1] = number of a's at index 1 in all strings  
		adenine = np.zeros(len(line))
	if "guanine" not in globals(): 
		guanine = np.zeros(len(line))
	if "uracil" not in globals(): 
		uracil = np.zeros(len(line))
	if "cysteine" not in globals(): 
		cysteine = np.zeros(len(line))

	for i in range(len(line)): 
		if line[i] == "A": 
			adenine[i] += 1
		elif line[i] == "G": 
			guanine[i] += 1
		elif line[i] == "U": 
			uracil[i] += 1
		elif line[i] == "C": 
			cysteine[i] += 1

infile.close()

#get the fraction of A/U/C/G for each position

for i in range(len(adenine)): 
	if adenine[i]+guanine[i]+uracil[i]+cysteine[i] == 0: 
		adenine[i] = 0.0
		guanine[i] = 0.0
		cysteine[i] = 0.0
		uracil[i] = 0.0

	else: 
		denominator = float(adenine[i]+guanine[i]+uracil[i]+cysteine[i])
		adenine[i] = float(adenine[i])/denominator
		guanine[i] = float(guanine[i])/denominator
		cysteine[i] = float(cysteine[i])/denominator
		uracil[i] = float(uracil[i])/denominator


#plot distribution of A/U/C/G at each position

j = 0 #this just keeps track of the number of figues, so each plot opens in a separate figure

#A

plt.figure(j)

color = ["salmon" for i in adenine]
edgecolor = ["salmon" for i in adenine]

plt.bar(np.arange(len(adenine)),adenine, color=color, edgecolor=edgecolor)
#add lines corresponding to the average downstream and upstream read counts 
#from https://stackoverflow.com/questions/36470343/how-to-draw-a-line-with-matplotlib

plt.xticks(np.arange(0,len(adenine),5),rotation="vertical")
plt.autoscale()	
plt.xlabel("Location along window")
plt.ylabel("Fraction of A's")
plt.title("Fraction of A's in MSA")
j += 1
plt.draw()

#U

plt.figure(j)

color = ["salmon" for i in uracil]
edgecolor = ["salmon" for i in uracil]

plt.bar(np.arange(len(uracil)), uracil, color=color, edgecolor=edgecolor)
#add lines corresponding to the average downstream and upstream read counts 
#from https://stackoverflow.com/questions/36470343/how-to-draw-a-line-with-matplotlib

plt.xticks(np.arange(0,len(uracil),5),rotation="vertical")
plt.autoscale()	
plt.xlabel("Location along window")
plt.ylabel("Fraction of U's")
plt.title("Fraction of U's in MSA")
j += 1
plt.draw()

#C

plt.figure(j)

color = ["salmon" for i in cysteine]
edgecolor = ["salmon" for i in cysteine]

plt.bar(np.arange(len(cysteine)),cysteine, color=color, edgecolor=edgecolor)
#add lines corresponding to the average downstream and upstream read counts 
#from https://stackoverflow.com/questions/36470343/how-to-draw-a-line-with-matplotlib

plt.xticks(np.arange(0,len(cysteine),5),rotation="vertical")
plt.autoscale()	
plt.xlabel("Location along window")
plt.ylabel("Fraction of C's")
plt.title("Fraction of C's in MSA")
j += 1
plt.draw()

#G

plt.figure(j)

color = ["salmon" for i in guanine]
edgecolor = ["salmon" for i in guanine]

plt.bar(np.arange(len(guanine)),guanine, color=color, edgecolor=edgecolor)
#add lines corresponding to the average downstream and upstream read counts 
#from https://stackoverflow.com/questions/36470343/how-to-draw-a-line-with-matplotlib

plt.xticks(np.arange(0,len(guanine),5),rotation="vertical")
plt.autoscale()	
plt.xlabel("Location along window")
plt.ylabel("Fraction of G's")
plt.title("Fraction of G's in MSA")
plt.draw()

plt.show()
