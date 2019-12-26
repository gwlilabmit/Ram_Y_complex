from base_complement import * 

#open genome file and extract window

bedfile = "glna_cleavage.bed12"  #edit path as needed
bed = open(bedfile, 'r')

old_cleavage = "glna_windows.txt"
old_sites = open(old_cleavage, 'r')

new_cleavage = "glna_windows_newer.txt"
new_sites = open(new_cleavage, 'w')

for line1 in bed:
	for line2 in old_sites: 
		if '>' in line2:
			new_sites.write(line2)
		else: 
			if '+' in line1: 
				new_sites.write(line2)
			else: 
				comp = seqcomplement(line2, 'DNA')
				new_sites.write(comp) 

bed.close()
old_sites.close()
new_sites.close()

