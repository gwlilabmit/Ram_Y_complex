from base_complement import * 
from string_reversal import * 

#open genome file and extract window

bedfile = "metq_correct_evolutionary_longer.bed12"  #edit path as needed
bed = open(bedfile, 'r')

old_cleavage = "metq_evolutionary_longer_old.txt"
old_sites = open(old_cleavage, 'r')

cleavage_site_seqs = "metq_cleavage_seqs_shortened.txt"
stripped = open(cleavage_site_seqs, 'w')

cleavage_site_names = "metq_cleavage_names_shortened.txt"
names = open(cleavage_site_names, 'w')

new_cleavage = "metq_evolutionary_longer.txt"
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

