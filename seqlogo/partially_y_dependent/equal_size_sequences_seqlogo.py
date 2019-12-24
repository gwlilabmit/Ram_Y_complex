'''
Seqlogo requires me to have input sequences of all the same size. Since all cleavage windows have 100nts downstream of the cleavage site, I first find the length of the shortest sequences in the input--let's call it m. If a given sequence has length n>m, I slice the first n-m characters from the beginning of the string so that it has length m. I then recreate the input file format [fasta format] with the new sequences
'''

infile = open("subtilis_partially_dependent.txt", "r")
outfile = open("shortened_subtilis_partially_dependent.txt", "w")

fastaheaders = [] #contains all lines that serve as fasta headers, e.g. '>atpI'
sequences = [] #contains all sequences
shortenedsequences = [] #contains all shortened sequences

for line in infile: 
	if ">" in line: #i.e. if it's a fasta header
		fastaheaders.append(line[:-1]) #the -1 gets rid of the \n at the end of the line
	else: 
		sequences.append(line[:-1])

shortestlength = len(min(sequences, key=len)) #the min() just gets the sequence of smallest length in the array sequences

for i in range(len(sequences)): 
	sequence = sequences[i]
	if len(sequence)>shortestlength: 
		lengthtocut = len(sequence)-shortestlength #so if the shortest length is 100nts and this sequence is 102nts, we want to cut 2nt
		shortenedsequence = sequence[lengthtocut:] #we're cutting those 2nt from the beginning of the sequence
		
		shortenedsequences.append(shortenedsequence)
	else: 
		shortenedsequences.append(sequence)

for i in range(len(fastaheaders)): 
	outfile.write(fastaheaders[i]+"\n")
	outfile.write(shortenedsequences[i]+"\n")


outfile.close()
infile.close()
