#Takes RNAfold -p0 structures and outputs probabilities of a given nt being unpaired into a .dat file
#You can then use that .dat file to constrain structure formation using RNAfold --shape=___.dat
#unconstrained is a file containing the kTln(Z) value [from RNAfold -p0] for the sequence of interest
#constrained is a file with -p0 -C output for the same sequence--it has output for each nt individually constrained
#The assumption here is that the first sequence has the first nt constrained, the second has the second nt constrained and so on
#
#Basically I grab the kTln(Z) values for all sequences in constrained, from there get the value of Z and divide that by the Z value from unconstrained 
#This gives me the pairing probabilities for each nt
#I output that into a .dat file that I can then use to constrain structure output with RNAfold --shape
#Additionally, if just you want to use values from that .dat file as a colormap in VARNA run colormap_from_dat.py 

import re
import numpy as np


name = "metq_correct"

unconstrained = open(name+"_partition.txt", "r") #RNAfold -p0 output for sequence without constraints
constrained = open(name+"_constrained_structures.txt", "r") #RNAfold -p0 -C output for sequences with a single nt constrained to be unpaired
probabilities = open(name+"_paired_probs.dat", "w") #This is the .dat file you'll be outputting the paired probabilities to 

####please don't change these parameters!!

T = 310 #This is the folding temperature in Kelvin as you've specified in RNAfold 
k = 0.0019872041 #kcal/(mol*K), this is the Boltzmann constant 
regex = r"\d+\.\d+"	#captures anything of the form "int.int," e.g. 51.05
					#I'm purposely getting rid of the - sign because I would've gotten rid of it otherwise anyway



for line in unconstrained: 
	if "free" in line: #so if this line contains the kTln(Z) value--RNAfold will always put this on its own line as "free energy of ensemble = -51.05 kcal/mol"
		unconstrained_energy = re.findall(regex, line) #final all strings matching regex

energies = [] #array for all constrained kTln(Z) values

for line in constrained: 
	if "free" in line: #so if this line contains the kTln(Z) value--RNAfold will always put this on its own line as "free energy of ensemble = -51.05 kcal/mol"
		energy = re.findall(regex, line) #final all strings matching regex
		energies.extend(energy)

constrained.close()
unconstrained.close()

Z = float(unconstrained_energy[0]) #Z is the value of kTln(Z) for an unconstrained structure
#so Z is really k*T*ln(Z)
Z = Z/(k*T)
Z = np.exp(Z)

energies = np.array(energies) #convert all these energy values to floats [they're currently strings] so we can work with them
energies = energies.astype(np.float)

for i in range(len(energies)): 
	energy = energies[i]
	#again, each energy value is really kT ln(Z) and we just want Z  
	energy = energy/(k*T)
	energy = np.exp(energy)
	energy = energy/Z
	#So the structure of the datfile is index \t value--I'm just writing to the datfile here
	probabilities.write(str(i+1)+"\t"+str(energy)+"\n") #divide each constrained kTln(Z) by Z--I'm adding +1 to the i value because .dat files are 1-indexed
	energies[i] = energy

probabilities.close()

