#sequence = raw_input("Enter the sequence you'd like the complement of: \n")
#seqtype = raw_input("Is this DNA or RNA? \n")

sequence = ''
seqtype = 'DNA'

def seqcomplement(sequence, seqtype): 
	complement = ''
	for i in range(len(sequence)):
		char = sequence[i]
		if char == 'A' or char == 'a':
			if 'R' in seqtype or 'r' in seqtype: 
				complement += 'U'
			else: 
				complement += 'T'

		elif char == 'T' or char == 't': 
			if 'R' in seqtype or 'r' in seqtype: 
				answer = raw_input("Just letting you know, this sequence has T's despite you saying it's RNA. Do you want me to treat it as DNA? [Answer yes/no]\n")
				if 'y' in answer or 'Y' in answer:
					seqtype = 'DNA'
					for j in range(len(complement)): 
						if complement[j] == 'U':
							complementlist = list(complement)
							complementlist[j] = 'T'
							complement = "".join(complementlist)								
				complement += 'A'
			else:
				complement += 'A'
		elif char == 'U' or char == 'u': 
			if 'D' in seqtype or 'd' in seqtype: 
				answer = raw_input("Just letting you know, this sequence has U's despite you saying it's DNA. Do you want me to treat it as RNA? [Answer yes/no]\n")
				if 'y' in answer or 'Y' in answer: 
					seqtype = 'RNA'
					for j in range(len(complement)): 
						if complement[j]== 'T':
							complementlist = list(complement)
							complementlist[j] = 'U' 
							complement = "".join(complementlist)
				complement += 'A'

		elif char == 'C' or char == 'c':
			complement += 'G'

		elif char == 'G' or char == 'g': 
			complement += 'C'

		elif char == '\n' or char == ' ': 
			continue 
		else: 
			print "There's some non-A/U/T/C/G base here. Check index", i, "for the problem. There may be more issues later on in the string but this is the first issue I'm detecting."
			complement = ""
			break 
		i+=1

	return complement

seqcomplement(sequence, seqtype)
