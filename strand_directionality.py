from base_complement import * 
from string_reversal import * 

#open genome file and extract window

bedfile = "cleavage_windows.bed12"  #bedfile, change name as needed
bed = open(bedfile, 'r')

old_cleavage = "cleavage_sites.txt" #old getfasta output, change name as needed
old_sites = open(old_cleavage, 'r')

cleavage_site_seqs = "cleavage_seqs.txt" #no need to change, intermediate file
stripped = open(cleavage_site_seqs, 'w')

cleavage_site_names = "cleavage_names.txt" #no need to change, intermediate file
names = open(cleavage_site_names, 'w')

new_cleavage = "21_cleavage_newer.txt" #new output, change name as needed
new_sites = open(new_cleavage, 'w')

for line in old_sites: 
	if ">" in line: 
		names.write(line)
	else: 
		stripped.write(line)

stripped.close()
names.close()
stripped = open(cleavage_site_seqs, 'r')
names = open(cleavage_site_names, 'r')

for line1, line2, line3 in zip(bed, stripped, names): 
	new_sites.write(line3)
	if '-' in line1:
		comp = reversal(seqcomplement(line2, 'DNA'))
		new_sites.write(comp+'\n')
	else: 
		new_sites.write(line2) 


bed.close()
old_sites.close()
stripped.close()
names.close()
new_sites.close()

