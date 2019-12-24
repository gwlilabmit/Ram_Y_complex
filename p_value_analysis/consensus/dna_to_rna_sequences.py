'''
Literally takes a DNA sequence and converts it to RNA, i.e. literally converts T's to U's from a fasta-format input
'''

infile = open("yqhl.txt", "r")
outfile = open("rna_yqhl.txt", "w")

for line in infile: 
	if ">" in line: #i.e. if it's a fasta header
		outfile.write(line)
	else: 
		sequence = line[:-1] #:-1 gets rid of the \n at the end 
		sequence = list(sequence) #converting to a list so I can mutate it
		for i in range(len(sequence)): 
			char = sequence[i]
			char = char.upper() #just standardize everything as uppercase

			if char == "T": #convert T->U
				sequence[i] = "U" 
		sequence = "".join(sequence) #convert back into a string
		outfile.write(sequence+"\n")

outfile.close()
infile.close()
