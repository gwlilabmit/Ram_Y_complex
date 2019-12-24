#this is a really dumb script. It literally takes in a fasta format file with multiple entries
#e.g. >gene_1 \n sequence_1 \n >gene_2 \n sequence_2 ...
#and it exports each >gene \n sequence into a separate file named based on the salient header

infile = open("random_mfe.txt", "r") #choose your input file

for line in infile: 
	if ">" in line: 
		name = line[1:-1] #get rid of > and \n before and after then name
	else: 
		outfile = open(name+"_random.txt", "w") # e.g. ">cggR" is written into cggr.txt
		outfile.write(">"+name+"\n"+line) #the output file is just >cggr \n cggr_sequence 
		outfile.close()

infile.close()
