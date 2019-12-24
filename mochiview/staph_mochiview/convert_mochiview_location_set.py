'''
Converts location sets from the NCBI format [e.g. https://www.ncbi.nlm.nih.gov/genome/proteins/154?genome_assembly_id=299272] into the Mochiview format
'''

infile = open("NC_007795_cds.txt", "r")
outfile = open("NC_007795_mochi_cds.txt", "w")

#write header line to outfile
outfile.write('SEQ_NAME\tSTART\tEND\tSTRAND\tFEATURE_NAME\tGENE_NAME\tALIASES\tDESCRIPTION\tTXN_START\tTXN_END\tCDS_START\tCDS_END\tEXON_COUNT\tEXON_STARTS\tEXON_ENDS\n')


for line in infile: 
	if "#" not in line: #the header row starts with a #
		splitline = line.split('\t')
		#we want the genome name, start, stop and strand
		seqname = splitline[1]
		start = splitline[2]
		end = splitline[3]
		strand = splitline[4]
		name = splitline[6]

		outfile.write(seqname+"\t"+start+"\t"+end+"\t"+strand+"\t"+name+"\t"+name+"\t"+' '+"\t"+' '+'\t') #write seq_name->description columns
		outfile.write(start+"\t"+end+"\t"+start+"\t"+end+"\t"+'1'+"\t"+start+"\t"+end+'\n') #write the rest of the columns--exon_count is set to 1
	

outfile.close()
infile.close()

