'''
Gets all unique kmers for a given value of k for all 21 sites and calculates their frequency as well as p values
I iterate through k = 3->10 and identify the most enriched kmers, then calculate their p values
'''

from scipy.misc import comb


#############################
#completely global variables
#############################

pvalthresh = 0.001 #threshold for p value significance

def ncr(n, r): 
	return comb(n, r, exact=0) #approximate the combination--this lets me deal with larger numbers 

def kmerdist(k): 

	#################
	#input files
	#################

	cleavagesites = open("21_cleavage_final.txt", 'r') #windows for all 21 sites in subtilis

	###############################
	#global vars and dictionaries
	###############################

	kmers = set() #a set of all unique kmers, produced by get_kmers() so I don't redo it on every iter of enrichment()--that's why it's a set, not a list

	subtiseqs = [] #an array with all subtilis sequence windows

	subtiliscounts = {} #key: kmer, value: number of times it occurs in all subtilis sequences

	subtiprobs = {} #key: kmer, value: probability of seeing it subtiliscounts[kmer] times in the 21 sites

	pvalues = {} #key: kmer, value: p value of subtiprobs[kmer]

	numkmers = 0 #number of kmers total in all sequences 

	###probabilities from nt_probabilities.py--I took the 168 genome and counted the A's, G's etc. 
	proba = 0.281827333959
	probg = 0.217077212624
	probc = 0.218066868678
	probt = 0.28302858474

	##################################################
	#get all sequences for background sets
	##################################################

	#subtilis 

	for line in cleavagesites: 
		if ">" not in line: #i.e. if the line isn't a fasta header
			subtiseqs.append(line[:-1]) #the [:-1] gets rid of the \n at the end

	##########################################################################
	#get number of kmers in total for all sites and identify all unique kmers
	##########################################################################

	#take a rolling window of size k through each sequence in subtiseqs--if it's unique add it to kmers

	for sequence in subtiseqs: 

		#for a sequence of length L, there are L-k kmers--add this number to numkmers
			
		numkmers += len(sequence)-k+1

		for i in range(len(sequence)-k+1):
				
			#identify each kmer and add to the set kmers		
	 
			kmer = sequence[i:i+k]
			kmers.add(kmer)

	#we don't need kmers to be a set anymore, convert to a list so we can manipulate it later
	kmers = list(kmers)

	#################################################################
	#calculate frequency of each kmer in subtilis+background windows
	#################################################################

	for kmer in kmers: 
		subtiliscounts[kmer] = 0 #add it as a key to subtiliscounts
		
		for sequence in subtiseqs: 
			frequency = sequence.count(kmer)
			subtiliscounts[kmer]+= frequency

		
	###############################################################################################
	#calculate the probability of seeing each kmer subtiliscounts[kmer] times and add to subtiprobs
	################################################################################################

	#We take this to be p(n times) = (N choose n)*p^n(1-p)^{N-n}, where p = p(kmer), n=freq and N = numkmers
	#say the kmer is AGG--in that case p = at*gc*gc 

	#additionally, we take the p value to be sum_n^N (N choose n)p^n(1-p)^{N-n}

	for kmer in subtiliscounts.keys():

		#############
		#probability
		#############
		
		p = 1
		freq = subtiliscounts[kmer]
 
		for i in range(len(kmer)): #calculate probabilities for each kmer
			if kmer[i] == "A": 
				p *= proba 
			elif kmer[i]=="T": 
				p *= probt
			elif kmer[i] == "G": 
				p *= probg
			else: 
				p *= probc

		comb = ncr(numkmers, freq)
		pfreq = p**freq
		npfreq = (1-p)**(numkmers-freq)	
		prob = comb*pfreq*npfreq
	
		subtiprobs[kmer] = prob

		##########
		#p value
		##########
		
		pvalue = 0

		for newfreq in range(freq, numkmers+1): #+1 because python's range function takes -1 from the second argument
			newpfreq = p**newfreq
			newnpfreq = (1-p)**(numkmers-newfreq)
			if newpfreq == 0.0 or newnpfreq == 0:
				newprob = 0.0
			else:  
				newcomb = ncr(numkmers, newfreq)
				newprob = newcomb*newpfreq*newnpfreq
			
			pvalue += newprob

		
		pvalues[kmer] = pvalue*numkmers #multiplying by numkmers for bonferroni correction
		

	cleavagesites.close()

	print sorted(subtiprobs.items(), key=lambda x: x[1]) #uncomment to get a list of all kmers and their associated p values

	return pvalues, subtiliscounts


##########################################################################################
#iterate over all kmers of size 3-10, find the lowest p value kmers and their frequencies
##########################################################################################

kmerpvals = {} #key: kmer, value: p value [for all values of k]
kmerfreqs = {} #key: kmer, value: frequency [for all values of k]

pvalsandfreqs = {} #key: kmer, value: [p>0.001, frequency]

#kmerdist returns (pvalues, subtiliscounts), so kmerdist(k)[0] = pvalues, kmerdist(k)[1] = subtiliscounts

for k in range(3, 11): #using 11 so we capture all kmers of size 3->10
	kpvals = kmerdist(k)[0]
	kmerpvals.update(kpvals)

	kfreq = kmerdist(k)[1]
	kmerfreqs.update(kfreq)


###################################
#get all kmers with p<=pvalthresh
###################################

#from https://stackoverflow.com/questions/18807079/selecting-elements-of-a-python-dictionary-greater-than-a-certain-value/18807120
significantkmers = dict((k, v) for k, v in kmerpvals.items() if v <= pvalthresh)

for key in significantkmers.keys(): 
	pval = significantkmers[key]
	freq = kmerfreqs[key]

	pvalsandfreqs[key] = [pval, freq]

#from https://thispointer.com/python-how-to-sort-a-dictionary-by-key-or-value/ -- prints kmercounts in ascending order of probabilities [values]
#print sorted(pvalsandfreqs.items() , key=lambda x: x[1])
	

