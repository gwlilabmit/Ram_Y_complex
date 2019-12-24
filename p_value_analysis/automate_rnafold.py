import re
import os
import glob


###open files###

sequences = glob.glob("*_silvi_constraint.txt")

###regexes###

nameregex = r"(^[a-zA-Z]\d+)|(^[a-zA-Z]+)"
bsuregex = r"\d+"

for sequence in sequences: 
	nametuple = re.findall(nameregex, sequence)[0]
	for i in range(len(nametuple)): 
		if len(nametuple[i])> 0: 
			name = nametuple[i]
	if name == "BSU": 
		bsunumber = re.findall(bsuregex, sequence)[0]
		name = "BSU_misc_RNA_"+bsunumber

	command = "rnafold --noPS -C --infile="+name+"_silvi_constraint.txt --outfile="+name+"_silvi_only.txt" 

	os.system(command)
