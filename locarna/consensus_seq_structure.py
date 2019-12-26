#sequence and structure produced by locaRNA
#getting rid of the spaces so I can produce a structure ready for visualization

consensus_seq = 'G-------------------------------------------------------CAAUUGCCGGAGC-CAUUGGUUCAGGAGGCCUUGGAAACCUUGCGUACGU-----------------------------AGAAGGAUA----UCAGUCUAAUAACGCUGAUGUU-------------------ACCUUCG--------------------------------UUGCUACUGUGUUUAUUUUAAUCAUCGUGUUUAUUAUUCAAAUCAUCGGUGAUCUAAUAACAAAUAUUAUAGACAAACGAU------------------------------------------------------AAGGGAGGAUUUACAUUGAAAAAGCUAUUUU'

consensus_struct = '........................................................((((.((((((((......)).)))..)))))))..............................................(((((..(....(((((........)))))..................).....))))).............................................(((((((.(((((((...................))))))))))))))..........................................................................................................'

seq = ""
struct = ""

for i in range(len(consensus_seq)): 
	if consensus_seq[i]!='-': 
		seq+= consensus_seq[i]
		struct+= consensus_struct[i]

print seq
print struct


