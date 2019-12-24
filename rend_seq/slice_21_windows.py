#Literally slices specific sequences in 21_cleavage_windows_final.txt based on indices you provide

infile = open("21_cleavage_windows_final.txt", "r") #input file with full cleavage sequences
outfile = open("21_downstream_hairpins.txt", "w") #file to write out sliced sequences

indices = {"atpI":[104, 127], "cggR":[102,128], "cwlO":[97,121], "dacA":[26,63], "ddl":[48,69], "menA":[11,31], "metN":[110,142],\
"pyrG":[98,123], "rpoB":[103,143], "rpsL":[15,63], "ybfG":[70,109], "yonR":[113,147]} #indices that I'd like to have sliced out

for line in infile:
	if ">" in line: #if the line is a fasta header--we assume fasta-formatted txt files
		name = line[1:-1] #>atpi\n becomes atpi
	else: 
		if name in indices.keys(): #if that sequence needs to be sliced
			roi = indices[name] #this gives you an array with the start/end indices of the region of interest
			roistart = roi[0]-1 #first index is the start, subtract -1 to account for 0 indexing 
			roiend = roi[1] #don't subtract -1 because string slicing automatically subtracts -1 from the ending index
			outfile.write(">"+name+"\n") #e.g. write >atpi\n
			outfile.write(line[roistart:roiend]+"\n")
			
 
outfile.close()
infile.close() 
